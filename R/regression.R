# ============================================================================
#              EOG Regression: Ocular Artifact Removal by Regression
# ============================================================================
#
# A fast, deterministic alternative to ICA for removing eye blinks and eye
# movements. Works on continuous recordings AND on epoched (trial-cut) data,
# via one shared model object:
#
#   subtract_evoked()       hide each trial's evoked response before fitting
#   fit_eog_regression()    learn one weight per EEG channel per EOG channel
#   apply_eog_regression()  subtract weight x EOG from each EEG channel
#   new_eog_regression()    constructor / validator for the model object
#
# Ported from MNE-Python's mne/preprocessing/_regress.py (EOGRegression and
# regress_artifact), extended here with the epoch/Gratton path that source
# only documents as a usage recipe (fit on evoked-subtracted epochs, apply
# to the originals - see subtract_evoked()'s docs). Design points, all
# checked against that source for the continuous case:
#
#  - Ordinary least squares, all EOG channels fitted jointly, one fit per
#    target channel. The EOG and each target are mean-removed first, so no
#    intercept is stored and a channel's overall level (DC offset) is left
#    alone when the weights are applied. Continuous data is demeaned once,
#    over the whole recording; epoched data is demeaned per trial (a
#    trial's own resting level says nothing about the blink, only its
#    within-trial fluctuation does) - see .eog_solve_weights() and the
#    epoch branches of fit_eog_regression()/apply_eog_regression().
#  - Every sample is used and annotations are NOT consulted directly (MNE's
#    regression has no reject_by_annotation either, unlike ICA fitting) -
#    for epochs, annotation-based rejection already happened earlier, inside
#    epoch_eeg(), so by the time data reaches this file that step is done.
#  - The EEG must already be re-referenced: the weights depend on the
#    reference (MNE raises an error otherwise, and so does this file).
#  - Fit and apply are separate steps, so a model can be inspected, saved
#    with saveRDS() and reused on another recording - or on the other data
#    shape entirely (a model fit on epochs works on continuous data and vice
#    versa, since apply matches channels BY NAME, the same way get_sources()
#    / apply_ica() do; MNE itself requires identical channel order).
#  - Never touches ICA or ica$exclude: a separate choice from ica1.R /
#    ica_detect.R for the same artifact.
#
# Phase 1: continuous 'eeg' objects. Phase 2: epoched 'eeg_epochs' objects,
# via subtract_evoked() (Gratton et al. 1983 - subtract the response, fit on
# what is left; Croft & Barry 2000 - fit directly on averaged blinks, apply
# to continuous data, still possible here by fitting with no subtraction)
# and apply_eog_regression()'s reapply_baseline. Planned next: a weights
# topography plot and a bipolar VEOG/HEOG helper.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# Status: Phase 1 (continuous data) and Phase 2 (epoched data) built.
# Tested: see tests/testthat/test-regression.R
# ============================================================================

# ----------------------------------------------------------------------------
# new_eog_regression() - constructor / validator for the model object
# ----------------------------------------------------------------------------
#' Create an EOG Regression Model
#'
#' Builds and validates the small model object that
#' \code{\link{fit_eog_regression}} returns and
#' \code{\link{apply_eog_regression}} consumes. You rarely need to call it
#' yourself. Use it to bring in weights computed elsewhere (for example the
#' \code{betas} that MNE-Python's \code{regress_artifact()} returns, saved to
#' CSV and read back with \code{read.csv()}) so they can be used with
#' \code{apply_eog_regression()}.
#'
#' @param coef Numeric matrix of regression weights: one row per target
#'   channel, one column per EOG channel (MNE's \code{coef_}). Its dimnames
#'   are replaced by \code{ch_names} and \code{ch_names_artifact}.
#' @param ch_names Character vector of target channel names, in the row order
#'   of \code{coef}.
#' @param ch_names_artifact Character vector of EOG channel names, in the
#'   column order of \code{coef}.
#' @param reference Character string: the reference scheme the weights were
#'   fit in (\code{eeg$reference} at fit time), or \code{NA} if unknown.
#'   \code{\link{apply_eog_regression}} warns when the data's reference
#'   differs from it.
#' @param fit_on \code{"continuous"} (default) or \code{"epochs"}: what kind
#'   of data the weights were fit on. Informational only.
#' @param n_samples Number of samples the fit used, or \code{NA} if unknown.
#'   For an epoch fit this is the flattened total (time points per trial x
#'   number of trials).
#'
#' @return An object of class \code{eeg_eog_regression}: a list with
#'  \describe{
#'    \item{coef_}{Numeric matrix, target channels x EOG channels, with
#'      dimnames. Row i, column j is how much of EOG channel j is subtracted
#'      from target channel i.}
#'    \item{ch_names}{Target channel names.}
#'    \item{ch_names_artifact}{EOG channel names.}
#'    \item{reference}{Reference scheme at fit time (\code{NA} if unknown).}
#'    \item{fit_on}{\code{"continuous"} or \code{"epochs"}.}
#'    \item{n_samples_}{Number of samples used in the fit.}
#'  }
#'
#' @examples
#' \dontrun{
#'   # weights from elsewhere: 3 target channels x 2 EOG channels
#'   betas <- matrix(c(0.60, 0.35, 0.10,
#'                     0.20, 0.05, 0.00), nrow = 3)
#'   model <- new_eog_regression(betas,
#'                               ch_names = c("Fp1", "Fz", "Cz"),
#'                               ch_names_artifact = c("VEOG", "HEOG"),
#'                               reference = "Common Average")
#'   eeg_clean <- apply_eog_regression(model, eeg)
#' }
#'
#' @seealso \code{\link{fit_eog_regression}}, \code{\link{apply_eog_regression}}
#'
#' @export
new_eog_regression <- function(coef,
                               ch_names,
                               ch_names_artifact,
                               reference = NA_character_,
                               fit_on = c("continuous", "epochs"),
                               n_samples = NA_real_) {

  fit_on <- match.arg(fit_on)

  # ========== VALIDATE coef ==========

  if (!is.matrix(coef) || !is.numeric(coef)) {
    stop("ERROR: 'coef' must be a numeric matrix (target channels x EOG ",
         "channels).", call. = FALSE)
  }
  if (!all(is.finite(coef))) {
    stop("ERROR: 'coef' contains NA, NaN or Inf values.", call. = FALSE)
  }

  # ========== VALIDATE channel names ==========

  valid_names <- function(x) {
    is.character(x) && length(x) > 0 && !anyNA(x) && !anyDuplicated(x)
  }
  if (!valid_names(ch_names) || !valid_names(ch_names_artifact)) {
    stop("ERROR: 'ch_names' and 'ch_names_artifact' must be character ",
         "vectors of unique, non-missing names.", call. = FALSE)
  }
  if (nrow(coef) != length(ch_names) || ncol(coef) != length(ch_names_artifact)) {
    stop("ERROR: 'coef' has ", nrow(coef), " row(s) and ", ncol(coef),
         " column(s), but 'ch_names' has ", length(ch_names),
         " name(s) and 'ch_names_artifact' has ", length(ch_names_artifact),
         ".", call. = FALSE)
  }
  overlap <- intersect(ch_names, ch_names_artifact)
  if (length(overlap) > 0) {
    stop("ERROR: channel(s) ", paste(overlap, collapse = ", "),
         " appear in both 'ch_names' and 'ch_names_artifact'; a channel ",
         "cannot be regressed on itself.", call. = FALSE)
  }

  # ========== VALIDATE the descriptive fields ==========

  if (!is.character(reference) || length(reference) != 1) {
    stop("ERROR: 'reference' must be a single character string (or NA).",
         call. = FALSE)
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1) {
    stop("ERROR: 'n_samples' must be a single number (or NA).", call. = FALSE)
  }

  # ========== BUILD ==========

  storage.mode(coef) <- "double"
  dimnames(coef) <- list(ch_names, ch_names_artifact)

  structure(
    list(
      coef_ = coef,
      ch_names = ch_names,
      ch_names_artifact = ch_names_artifact,
      reference = reference,
      fit_on = fit_on,
      n_samples_ = n_samples
    ),
    class = "eeg_eog_regression"
  )
}

# ----------------------------------------------------------------------------
# print.eeg_eog_regression() - readable summary of a model
# ----------------------------------------------------------------------------
#' Print Method for EOG Regression Models
#'
#' Displays what a \code{eeg_eog_regression} model was fit on (data type,
#' reference, EOG channels, target channels) and, for each EOG channel, the
#' target channel it is subtracted from most strongly. The full weight table
#' is \code{model$coef_}.
#'
#' @param x An object of class \code{eeg_eog_regression}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x} (standard R print method convention).
#'
#' @examples
#' \dontrun{
#'   model <- fit_eog_regression(eeg)
#'   print(model)  # Calls this method automatically
#' }
#'
#' @export
print.eeg_eog_regression <- function(x, ...) {

  n_tgt <- length(x$ch_names)
  n_art <- length(x$ch_names_artifact)

  # ========== HEADER ==========
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("EOG Regression Model\n")
  cat(strrep("=", 70), "\n\n")

  # ========== FIT INFORMATION ==========
  cat("FIT:\n")
  cat("  Fitted on:        ", x$fit_on, " data",
      if (!is.na(x$n_samples_)) paste0(" (", x$n_samples_, " samples)") else "",
      "\n", sep = "")
  cat("  Reference at fit: ",
      if (is.na(x$reference)) "unknown" else x$reference, "\n", sep = "")

  # ========== CHANNELS ==========
  cat("\nCHANNELS:\n")
  cat("  EOG regressors:   ", n_art, " (",
      paste(x$ch_names_artifact, collapse = ", "), ")\n", sep = "")

  tgt_display <- paste(x$ch_names[seq_len(min(5, n_tgt))], collapse = ", ")
  if (n_tgt > 5) {
    tgt_display <- paste0(tgt_display, " ... (+", n_tgt - 5, " more)")
  }
  cat("  Target channels:  ", n_tgt, " (", tgt_display, ")\n", sep = "")

  # ========== STRONGEST WEIGHT PER EOG CHANNEL ==========
  cat("\nSTRONGEST WEIGHT PER EOG CHANNEL:\n")
  for (j in seq_len(n_art)) {
    i <- which.max(abs(x$coef_[, j]))
    cat(sprintf("  %-16s -> %-10s %+.3f\n",
                x$ch_names_artifact[j], x$ch_names[i], x$coef_[i, j]))
  }

  # ========== FOOTER ==========
  cat("\n", strrep("=", 70), "\n\n", sep = "")

  invisible(x)
}

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

#' Has No EEG Reference Been Applied Yet? (internal)
#'
#' The package has two "nothing re-referenced yet" labels: \code{new_eeg()}
#' defaults \code{eeg$reference} to \code{"original"}, and
#' \code{read_bdf_native()} sets \code{"Biosemi CMS/DRL"}.
#' \code{eeg_rereference()} replaces either with a real scheme
#' (\code{"Common Average"}, \code{"M1+M2"}, ...). A missing or \code{NA}
#' reference counts as "not applied" as well.
#'
#' @param reference \code{eeg$reference} (works the same for an \code{eeg} or
#'   an \code{eeg_epochs} object - both carry this field the same way).
#' @return \code{TRUE} if no real reference has been applied, else \code{FALSE}.
#' @keywords internal
.eog_reference_missing <- function(reference) {
  if (is.null(reference) || length(reference) != 1 || is.na(reference)) {
    return(TRUE)
  }
  tolower(trimws(reference)) == "original" ||
    grepl("cms", reference, ignore.case = TRUE)
}

#' Stop If the EEG Has No Reference Applied (internal)
#'
#' MNE refuses to run EOG regression on EEG without a reference
#' (\code{_needs_eeg_average_ref_proj}); the weights are only valid in the
#' reference frame they were fit in.
#'
#' @param eeg An object of class 'eeg' or 'eeg_epochs'.
#' @return Invisibly \code{TRUE}; otherwise stops with a message that says how
#'   to fix it.
#' @keywords internal
.eog_check_reference <- function(eeg) {
  if (.eog_reference_missing(eeg$reference)) {
    shown <- if (is.null(eeg$reference) || length(eeg$reference) != 1) {
      "NULL"
    } else {
      as.character(eeg$reference)
    }
    stop("ERROR: no EEG reference has been applied yet (eeg$reference is '",
         shown, "'). EOG regression weights depend on the reference, so ",
         "re-reference first, e.g. eeg_rereference(eeg, ref = \"average\", ",
         "exclude = <EOG and Status channels>).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Turn Channel Picks Into Channel Names (internal)
#'
#' @param picks Character vector of channel names, or numeric vector of
#'   channel indices.
#' @param eeg An object of class 'eeg' or 'eeg_epochs'.
#' @param arg Name of the calling argument, used in error messages.
#' @return Character vector of channel names, in the order given.
#' @keywords internal
.eog_resolve_picks <- function(picks, eeg, arg) {

  if (is.character(picks)) {
    idx <- match(picks, eeg$channels)
    if (anyNA(idx)) {
      stop("ERROR: '", arg, "' contains channel name(s) not found in ",
           "eeg$channels: ", paste(picks[is.na(idx)], collapse = ", "),
           call. = FALSE)
    }
  } else if (is.numeric(picks)) {
    idx <- suppressWarnings(as.integer(picks))
    if (anyNA(idx) || any(picks != idx) ||
        any(idx < 1L | idx > length(eeg$channels))) {
      stop("ERROR: '", arg, "' must hold valid channel indices (whole ",
           "numbers from 1 to ", length(eeg$channels), ").", call. = FALSE)
    }
  } else {
    stop("ERROR: '", arg, "' must be a character vector of channel names ",
         "or a numeric vector of channel indices.", call. = FALSE)
  }

  if (length(idx) == 0) {
    stop("ERROR: '", arg, "' is empty.", call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    stop("ERROR: '", arg, "' lists the same channel more than once.",
         call. = FALSE)
  }

  eeg$channels[idx]
}

# ----------------------------------------------------------------------------
# .eog_solve_weights() - shared last step of fit_eog_regression()
# ----------------------------------------------------------------------------
#' Solve the Normal Equations for EOG Regression Weights (internal)
#'
#' The shared last step of \code{\link{fit_eog_regression}} for both
#' continuous and epoched data: given the EOG channels' Gram matrix and their
#' cross-product with each (already demeaned) target channel, solves for the
#' least-squares weights. The two input shapes build \code{G}/\code{B}
#' differently (continuous: demeaned once over the whole recording; epochs:
#' demeaned per trial, accumulated trial by trial) but from here on the
#' maths - and the collinearity error/warning - are identical either way.
#'
#' @param G Numeric matrix, k x k (EOG channels): \code{R \%*\% t(R)} for
#'   whatever demeaned EOG matrix \code{R} the caller built (or the sum of
#'   that quantity across trials, for epochs).
#' @param B Numeric matrix, k x p (EOG channels x target channels), column j
#'   is \code{R \%*\% demeaned_target_j} (summed across trials, for epochs).
#' @return Numeric matrix, p x k (target channels x EOG channels) - ready to
#'   store as \code{coef_} (see \code{\link{new_eog_regression}}).
#' @keywords internal
.eog_solve_weights <- function(G, B) {
  coef_t <- tryCatch(
    solve(G, B),                            # k x p
    error = function(e) {
      stop("ERROR: the EOG channels are collinear (one is an exact or ",
           "near-exact linear combination of the others), so the weights ",
           "cannot be estimated. Use fewer EOG channels in ",
           "'picks_artifact', or replace redundant ones with a bipolar ",
           "derivation such as upper minus lower. (", conditionMessage(e),
           ")", call. = FALSE)
    }
  )

  rc <- rcond(G)
  if (rc < 1e-10) {
    warning("The EOG channels are nearly collinear (reciprocal condition ",
            "number ", signif(rc, 2), "), so the weights may be unstable. ",
            "Consider using fewer EOG channels.", call. = FALSE)
  }

  t(coef_t)                                 # p x k
}

# ----------------------------------------------------------------------------
# .reapply_baseline() - put corrected channels' baseline back at zero
# ----------------------------------------------------------------------------
#' Re-Baseline Specific Channels of an Epoch Array (internal)
#'
#' After \code{\link{apply_eog_regression}} subtracts weight x EOG from each
#' target channel of an \code{eeg_epochs} object, a trial's baseline-window
#' average can drift away from zero (the correction is not itself
#' baseline-anchored). This re-applies the SAME baseline window and method
#' \code{epoch_eeg()} already used to build the epochs
#' (\code{epochs$baseline}, \code{epochs$baseline_method}) - but only to the
#' given channel rows; every other channel is left exactly as
#' \code{apply_eog_regression()} produced it.
#'
#' @param data Numeric array, channels x times x trials (\code{epochs$data}
#'   after the correction has already been subtracted).
#' @param idx Integer vector - which rows of \code{data} to re-baseline (the
#'   model's target channels).
#' @param times Numeric vector, length \code{dim(data)[2]} (\code{epochs$times}).
#' @param baseline Numeric \code{c(start, end)}, in seconds
#'   (\code{epochs$baseline}).
#' @param baseline_method \code{"mean"} or \code{"median"}
#'   (\code{epochs$baseline_method}).
#' @return \code{data}, with rows \code{idx} re-baselined; every other row
#'   unchanged.
#' @keywords internal
.reapply_baseline <- function(data, idx, times, baseline, baseline_method) {

  bl_min <- which.min(abs(times - baseline[1]))
  bl_max <- which.min(abs(times - baseline[2]))
  bl_idx <- bl_min:bl_max

  fun <- if (identical(baseline_method, "median")) median else mean
  p <- length(idx)

  for (n in seq_len(dim(data)[3])) {
    trial  <- matrix(data[idx, , n], nrow = p)                  # p x n_times
    bl_val <- apply(trial[, bl_idx, drop = FALSE], 1, fun, na.rm = TRUE)
    data[idx, , n] <- trial - bl_val
  }

  data
}

# ----------------------------------------------------------------------------
# subtract_evoked() - hide each trial's evoked response before fitting
# ----------------------------------------------------------------------------
#' Subtract Each Trial's Own Average Response (Gratton-Style Pre-Processing)
#'
#' Removes the repeatable, time-locked brain response from every trial of an
#' \code{eeg_epochs} object, leaving mostly what varies trial to trial - a
#' stray eye movement chief among it, since a blink happens at a random
#' moment rather than locked to the event. Feeding this to
#' \code{\link{fit_eog_regression}} instead of the raw epochs keeps the real
#' evoked response from being mistaken for eye-artifact leakage while
#' fitting, following Gratton, Coles & Donchin (1983).
#'
#' @param epochs An \code{eeg_epochs} object (see \code{\link{epoch_eeg}}),
#'   fit with \code{preload = TRUE} (the default), so \code{epochs$data} is
#'   loaded.
#' @param by \code{"event_type"} (default) or \code{"all"}.
#'   \code{"event_type"} subtracts each condition's own average from its own
#'   trials (see Details); \code{"all"} subtracts one grand average from
#'   every trial, matching what a single-condition analysis reduces to.
#'
#' @return A new object of class \code{eeg_epochs}, the same shape as
#'   \code{epochs}, with \code{data} replaced by the evoked-subtracted
#'   version and a note appended to \code{preprocessing_history}. Every other
#'   field (\code{events}, \code{channels}, \code{bads}, ...) is copied over
#'   unchanged. This is scratch data for fitting only - pass the
#'   \strong{original} \code{epochs} (not this function's output) to
#'   \code{\link{apply_eog_regression}}, so the real evoked response stays in
#'   the cleaned result.
#'
#' @details
#' Which average is used matters when conditions have genuinely different
#' responses: subtracting one grand average leaves a leftover trace of each
#' condition's own response behind (the gap between the true response and
#' the blended average), which the regression can mistake for part of the
#' eye artifact. Subtracting each condition's own average removes that
#' leftover almost completely. When blinks are large this difference barely
#' matters; when blinks are small or infrequent, per-condition subtraction
#' recovers noticeably more accurate weights.
#'
#' A condition with only a single trial has that trial's own data as its
#' "average", so subtracting it leaves that trial at exactly zero - a
#' warning is issued when this happens, since a zeroed-out trial contributes
#' nothing useful to the fit.
#'
#' @examples
#' \dontrun{
#'   epochs   <- epoch_eeg(eeg, events = "all", tmin = -0.2, tmax = 0.8)
#'   learn_on <- subtract_evoked(epochs)                # per condition
#'   model    <- fit_eog_regression(learn_on)
#'   epochs   <- apply_eog_regression(model, epochs)    # the ORIGINAL epochs
#' }
#'
#' @seealso \code{\link{fit_eog_regression}}, \code{\link{apply_eog_regression}},
#'   \code{\link{epoch_eeg}}
#'
#' @export
subtract_evoked <- function(epochs, by = c("event_type", "all")) {

  by <- match.arg(by)

  # ========== VALIDATE inputs ==========

  if (!inherits(epochs, "eeg_epochs")) {
    stop("ERROR: 'epochs' must be an object of class 'eeg_epochs' (see ",
         "epoch_eeg()).", call. = FALSE)
  }
  if (is.null(epochs$data)) {
    stop("ERROR: epochs$data is not loaded - re-run epoch_eeg() with ",
         "preload = TRUE.", call. = FALSE)
  }

  n_trials <- dim(epochs$data)[3]
  if (is.null(n_trials) || n_trials < 2) {
    stop("ERROR: at least 2 trials are required to subtract an average, ",
         "got ", if (is.null(n_trials)) 0 else n_trials, ".", call. = FALSE)
  }

  # ========== GROUP TRIALS ==========

  groups <- if (by == "all") {
    list(all = seq_len(n_trials))
  } else {
    split(seq_len(n_trials), epochs$events$type)
  }

  singleton <- names(groups)[lengths(groups) == 1]
  if (length(singleton) > 0) {
    warning("Condition(s) with only 1 trial (", paste(singleton, collapse = ", "),
            ") have their own trial subtracted as the 'average', leaving ",
            "them exactly zero.", call. = FALSE)
  }

  # ========== SUBTRACT EACH GROUP'S OWN AVERAGE ==========

  out_data <- epochs$data
  for (idx in groups) {
    grp_mean <- rowMeans(epochs$data[, , idx, drop = FALSE], dims = 2, na.rm = TRUE)
    for (i in idx) {
      out_data[, , i] <- out_data[, , i] - grp_mean
    }
  }

  # ========== BUILD THE RESULT ==========

  out <- epochs
  out$data <- out_data
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Evoked response subtracted (by = \"", by, "\"): ",
                length(groups), " group(s), ", n_trials, " trial(s) total"))
  )

  out
}

# ----------------------------------------------------------------------------
# fit_eog_regression() - learn the weights
# ----------------------------------------------------------------------------
#' Fit an EOG Regression Model (Ocular Artifact Removal)
#'
#' Works out how strongly each EEG channel follows the EOG channel(s): one
#' weight per EEG channel per EOG channel, fitted together when there are
#' several (say vertical and horizontal). The result is a small model that
#' \code{\link{apply_eog_regression}} uses to subtract that share of the EOG
#' from the EEG. This is a fast, deterministic alternative to ICA for blinks
#' and eye movements, and a port of MNE-Python's \code{EOGRegression.fit()},
#' extended here to also fit on epoched data.
#'
#' @param eeg A continuous \code{eeg} object (see \code{new_eeg()}) or an
#'   \code{eeg_epochs} object (see \code{\link{epoch_eeg}}), either way
#'   already re-referenced (see \code{\link{eeg_rereference}}). For epochs,
#'   fitting on \code{\link{subtract_evoked}}'s output rather than the raw
#'   epochs is usually the better choice - see its docs.
#' @param picks \code{NULL} (default), a character vector of channel names, or
#'   a numeric vector of channel indices: the channels to compute weights for
#'   (the targets). If \code{NULL}, every channel with
#'   \code{channel_types == "eeg"} is used, except channels listed in
#'   \code{eeg$bads} and the EOG channels themselves - so the status channel
#'   and other external channels (ECG, EMG, ...) are left out. An explicit
#'   \code{picks} is taken as-is (it can include bad channels) but must not
#'   include any EOG channel.
#' @param picks_artifact \code{NULL} (default), a character vector of channel
#'   names, or a numeric vector of channel indices: the EOG channels used as
#'   predictors. If \code{NULL}, they are found by name, using the same lookup
#'   as \code{\link{find_bads_eog}} (external channels whose name contains
#'   \code{"EOG"}, e.g. \code{"EOG_L (EXG1)"} or \code{"VEOG"}), and any of
#'   those marked bad in \code{eeg$bads} are skipped. An explicit
#'   \code{picks_artifact} is taken as-is.
#'
#' @return An object of class \code{eeg_eog_regression} (see
#'   \code{\link{new_eog_regression}}). Its \code{coef_} matrix has one row
#'   per target channel and one column per EOG channel. It is a plain list, so
#'   it can be inspected, saved with \code{saveRDS()} and reused - on either
#'   data shape, since \code{\link{apply_eog_regression}} matches channels by
#'   name.
#'
#' @details
#' \strong{How the weights are found.} For each target channel, ordinary
#' least squares of the channel on all EOG channels together, after taking
#' each signal's average level out, so only the eye-related fluctuations
#' count. No intercept is stored; a channel's overall level is untouched when
#' the weights are applied. For continuous data the average level is taken
#' over the whole recording; for epoched data it is taken \strong{per trial}
#' (Gratton/Croft-Barry convention: a trial's own resting level says nothing
#' about the blink, only the fluctuation within that trial does).
#'
#' \strong{Reference.} The weights depend on the reference, so the EEG must
#' already be re-referenced; an object whose \code{eeg$reference} is still
#' \code{"original"} or \code{"Biosemi CMS/DRL"} is refused. When calling
#' \code{eeg_rereference()}, pass the EOG and status channels to
#' \code{exclude} so they stay out of the average.
#'
#' \strong{Every sample counts.} As in MNE, all samples are used; for
#' continuous data \code{eeg$annotations} is not consulted (for epoched data,
#' \code{\link{epoch_eeg}}'s own \code{reject_by_annotation} already handled
#' that earlier). Data must be complete: an NA, NaN or Inf in a used channel
#' is an error.
#'
#' \strong{Caveats.} This needs real EOG electrodes. EOG electrodes also pick
#' up some frontal brain activity, so regression can over-correct frontal
#' channels. If two EOG channels are (nearly) copies of each other the weights
#' cannot be estimated reliably: an error is raised for exact collinearity and
#' a warning for near-collinearity. Use it after filtering, and pick either
#' this or ICA for the eyes rather than reaching for both by habit; it never
#' touches an \code{eeg_ica} object or \code{ica$exclude}.
#'
#' @examples
#' \dontrun{
#'   # continuous data, after filtering, bad-channel handling and re-referencing
#'   eeg <- eeg_rereference(eeg, ref = "average",
#'                          exclude = c("EOG_L (EXG1)", "EOG_R (EXG2)", "Status"))
#'
#'   model <- fit_eog_regression(eeg)        # EOG channels found by name
#'   print(model)
#'   model$coef_                             # the full weight table
#'
#'   # or name the EOG channels yourself, and keep the model for later
#'   model <- fit_eog_regression(eeg, picks_artifact = c("VEOG", "HEOG"))
#'   saveRDS(model, "eog_model.rds")
#'
#'   eeg_clean <- apply_eog_regression(model, eeg)
#'
#'   # epoched data: fit on the evoked-subtracted residual, apply to the epochs
#'   epochs   <- epoch_eeg(eeg, events = "all", tmin = -0.2, tmax = 0.8)
#'   model    <- fit_eog_regression(subtract_evoked(epochs))
#'   epochs   <- apply_eog_regression(model, epochs)
#' }
#'
#' @seealso \code{\link{apply_eog_regression}}, \code{\link{new_eog_regression}},
#'   \code{\link{subtract_evoked}}, \code{\link{eeg_rereference}},
#'   \code{\link{find_bads_eog}}, \code{\link{apply_ica}}
#'
#' @export
fit_eog_regression <- function(eeg, picks = NULL, picks_artifact = NULL) {

  # ========== VALIDATE inputs ==========

  is_epochs <- inherits(eeg, "eeg_epochs")
  if (!is_epochs && !inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' or 'eeg_epochs' ",
         "(see new_eeg()/epoch_eeg()).", call. = FALSE)
  }
  .eog_check_reference(eeg)

  if (is_epochs && is.null(eeg$data)) {
    stop("ERROR: epochs$data is not loaded - re-run epoch_eeg() with ",
         "preload = TRUE.", call. = FALSE)
  }

  X <- eeg$data
  if (is_epochs) {
    if (!is.array(X) || length(dim(X)) != 3) {
      stop("ERROR: epochs$data must be a 3D array (channels x times x ",
           "trials).", call. = FALSE)
    }
    n_times  <- dim(X)[2]
    n_trials <- dim(X)[3]
    if (n_times * n_trials < 2) {
      stop("ERROR: at least 2 samples are required to fit, got ",
           n_times * n_trials, ".", call. = FALSE)
    }
  } else {
    if (!is.matrix(X) || !is.numeric(X)) {
      stop("ERROR: eeg$data must be a numeric matrix (channels x time).",
           call. = FALSE)
    }
    if (ncol(X) < 2) {
      stop("ERROR: at least 2 samples are required to fit, got ", ncol(X),
           ".", call. = FALSE)
    }
  }

  # ========== RESOLVE THE EOG (PREDICTOR) CHANNELS ==========
  # (shape-independent: only reads eeg$channels/channel_types/bads)

  if (is.null(picks_artifact)) {
    # same lookup find_bads_eog() uses (R/ica_detect.R)
    art <- tryCatch(.resolve_reference_channels(eeg, NULL, "EOG"),
                    error = function(e) NULL)
    if (is.null(art)) {
      ext <- eeg$channels[eeg$channel_types == "external"]
      stop("ERROR: no EOG channel found among the external channels (",
           if (length(ext) > 0) paste(ext, collapse = ", ") else "none",
           "). Pass 'picks_artifact' explicitly, e.g. ",
           "picks_artifact = c(\"VEOG\", \"HEOG\").", call. = FALSE)
    }
    art <- art[!(art %in% eeg$bads)]
    if (length(art) == 0) {
      stop("ERROR: every EOG channel found is marked bad in eeg$bads. Pass ",
           "'picks_artifact' explicitly to use one anyway.", call. = FALSE)
    }
  } else {
    art <- .eog_resolve_picks(picks_artifact, eeg, "picks_artifact")
  }

  # ========== RESOLVE THE TARGET CHANNELS ==========
  # (shape-independent: only reads eeg$channels/channel_types/bads)

  if (is.null(picks)) {
    if (is.null(eeg$channel_types)) {
      stop("ERROR: eeg$channel_types is missing, so the target channels ",
           "cannot be chosen automatically. Pass 'picks' explicitly.",
           call. = FALSE)
    }
    is_target <- eeg$channel_types == "eeg" &
      !(eeg$channels %in% eeg$bads) &
      !(eeg$channels %in% art)
    tgt <- eeg$channels[is_target]
  } else {
    tgt <- .eog_resolve_picks(picks, eeg, "picks")
    overlap <- intersect(tgt, art)
    if (length(overlap) > 0) {
      stop("ERROR: 'picks' includes channel(s) that are also EOG ",
           "regressors (", paste(overlap, collapse = ", "), "); a channel ",
           "cannot be regressed on itself.", call. = FALSE)
    }
  }
  if (length(tgt) == 0) {
    stop("ERROR: no target channels left to fit (need at least one good ",
         "EEG channel).", call. = FALSE)
  }

  art_idx <- match(art, eeg$channels)
  tgt_idx <- match(tgt, eeg$channels)
  k <- length(art_idx)
  p <- length(tgt_idx)

  if (is_epochs) {

    # ---- EOG: check finiteness over the whole block up front ----

    R_full <- X[art_idx, , , drop = FALSE]
    if (!all(is.finite(R_full))) {
      bad <- art[!apply(is.finite(R_full), 1, all)]
      stop("ERROR: non-finite values (NA/NaN/Inf) found in EOG channel(s): ",
           paste(bad, collapse = ", "), ". Regression needs complete data.",
           call. = FALSE)
    }

    # ---- ACCUMULATE G AND B ONE TRIAL AT A TIME ----
    #
    # Each trial's EOG and target channels are demeaned using THAT TRIAL's
    # own mean, then that trial's contribution is added to the running
    # totals - so the full data is never flattened into one big matrix, and
    # each trial's demeaned EOG is computed once and reused across every
    # target channel rather than recomputed per channel.

    G <- matrix(0, k, k)
    B <- matrix(0, k, p)
    not_finite <- character(0)

    for (n in seq_len(n_trials)) {
      eog_trial <- matrix(X[art_idx, , n], nrow = k)          # k x n_times
      eog_trial <- eog_trial - rowMeans(eog_trial)
      G <- G + tcrossprod(eog_trial)

      for (j in seq_len(p)) {
        y <- X[tgt_idx[j], , n]
        if (!all(is.finite(y))) {
          not_finite <- union(not_finite, tgt[j])
          next
        }
        B[, j] <- B[, j] + eog_trial %*% (y - mean(y))
      }
    }
    if (length(not_finite) > 0) {
      stop("ERROR: non-finite values (NA/NaN/Inf) found in target ",
           "channel(s): ", paste(not_finite, collapse = ", "), ". ",
           "Regression needs complete data; interpolate or crop the ",
           "affected stretch first (see annotate_nan()).", call. = FALSE)
    }

    n_samples_used <- n_times * n_trials
    fit_on <- "epochs"

  } else {

    # ---- EOG: check, then remove each channel's mean over the recording ----

    R <- X[art_idx, , drop = FALSE]
    if (!all(is.finite(R))) {
      stop("ERROR: non-finite values (NA/NaN/Inf) found in EOG channel(s): ",
           paste(art[!apply(is.finite(R), 1, all)], collapse = ", "),
           ". Regression needs complete data.", call. = FALSE)
    }
    R <- R - rowMeans(R)                    # k x n, every row now sums to ~0

    # ---- targets: cross-product with the EOG, one channel at a time ----
    #
    # B[, j] = R %*% (y_j - mean(y_j)), the right-hand side of the normal
    # equations for target j. Done channel by channel (as MNE does) so the
    # full data matrix is never copied.

    B <- matrix(0, nrow = k, ncol = p)
    not_finite <- character(0)
    for (j in seq_len(p)) {
      y <- X[tgt_idx[j], ]
      if (!all(is.finite(y))) {
        not_finite <- c(not_finite, tgt[j])
        next
      }
      B[, j] <- R %*% (y - mean(y))
    }
    if (length(not_finite) > 0) {
      stop("ERROR: non-finite values (NA/NaN/Inf) found in target ",
           "channel(s): ", paste(not_finite, collapse = ", "), ". ",
           "Regression needs complete data; interpolate or crop the ",
           "affected stretch first (see annotate_nan()).", call. = FALSE)
    }

    G <- tcrossprod(R)                      # k x k
    n_samples_used <- ncol(X)
    fit_on <- "continuous"
  }

  # ========== SOLVE THE NORMAL EQUATIONS (shared) AND BUILD THE MODEL ==========

  new_eog_regression(
    coef = .eog_solve_weights(G, B),        # p x k: targets x EOG
    ch_names = tgt,
    ch_names_artifact = art,
    reference = eeg$reference,
    fit_on = fit_on,
    n_samples = n_samples_used
  )
}

# ----------------------------------------------------------------------------
# apply_eog_regression() - subtract weight x EOG from each target channel
# ----------------------------------------------------------------------------
#' Remove EOG Artifacts With a Fitted Regression Model
#'
#' Subtracts, from each target channel, its weight times the EOG signal (with
#' the EOG's average level taken out first), and returns the cleaned object.
#' A channel's overall level is left alone, only the eye-related fluctuation
#' is removed. The EOG channels themselves and any channel not in the model
#' (status, bad channels left out of the fit, ...) pass through unchanged. A
#' port of MNE-Python's \code{EOGRegression.apply()}, extended here to also
#' apply to epoched data.
#'
#' @param model A fitted \code{eeg_eog_regression} object, from
#'   \code{\link{fit_eog_regression}} or \code{\link{new_eog_regression}}.
#' @param eeg A continuous \code{eeg} object or an \code{eeg_epochs} object to
#'   clean. It does not have to be the recording the model was fit on, or
#'   even the same data shape (a model fit on epochs can clean continuous
#'   data and vice versa), but it must contain every channel named in the
#'   model, and it must already be re-referenced.
#' @param reapply_baseline Logical, default \code{TRUE}. \code{eeg_epochs}
#'   input only (ignored for continuous data): subtracting the eye artifact
#'   can shift a trial's baseline-window average away from zero, since the
#'   correction is not itself baseline-anchored. When \code{TRUE}, the
#'   target channels are re-baselined afterward using the SAME window and
#'   method \code{\link{epoch_eeg}} already recorded on \code{eeg}
#'   (\code{eeg$baseline}, \code{eeg$baseline_method}); a no-op if
#'   \code{eeg$baseline} is \code{NULL} (no baseline was set).
#'
#' @return A new object (see \code{\link{new_eeg}}/\code{\link{epoch_eeg}}),
#'   the same class as \code{eeg}, with the model's target channels cleaned
#'   and a note appended to \code{preprocessing_history}. Since R does not
#'   change arguments in place, reassign the result (\code{eeg <-
#'   apply_eog_regression(model, eeg)}); the input object itself is left
#'   untouched.
#'
#' @details
#' Channels are matched to the model by name, so the channel order in
#' \code{eeg} does not matter (MNE requires identical order). If a channel
#' the model needs is missing, an error lists it. A warning is issued when
#' \code{eeg$reference} differs from the reference the model was fit in,
#' since the weights only strictly apply in that frame.
#'
#' For continuous data, each channel's average level is taken over the whole
#' recording. For epoched data, it is taken \strong{per trial}, matching how
#' \code{\link{fit_eog_regression}} treats epochs.
#'
#' Where the EOG is NaN, the cleaned value is NaN too (the value is unknown);
#' every other sample is cleaned normally.
#'
#' Applying the same model twice removes the EOG twice, so apply it once per
#' recording. This function never touches ICA; see
#' \code{\link{apply_ica}} for that route.
#'
#' @examples
#' \dontrun{
#'   model     <- fit_eog_regression(eeg)
#'   eeg_clean <- apply_eog_regression(model, eeg)
#'
#'   # reuse a saved model on another recording (same channel names)
#'   model <- readRDS("eog_model.rds")
#'   eeg2_clean <- apply_eog_regression(model, eeg2)
#'
#'   # epoched data: fit on the residual, apply to the original epochs
#'   model  <- fit_eog_regression(subtract_evoked(epochs))
#'   epochs <- apply_eog_regression(model, epochs)     # reapply_baseline = TRUE
#' }
#'
#' @seealso \code{\link{fit_eog_regression}}, \code{\link{new_eog_regression}},
#'   \code{\link{subtract_evoked}}, \code{\link{apply_ica}}
#'
#' @export
apply_eog_regression <- function(model, eeg, reapply_baseline = TRUE) {

  # ========== VALIDATE inputs ==========

  if (!inherits(model, "eeg_eog_regression")) {
    stop("ERROR: 'model' must be an object of class 'eeg_eog_regression' ",
         "(see fit_eog_regression()).", call. = FALSE)
  }
  is_epochs <- inherits(eeg, "eeg_epochs")
  if (!is_epochs && !inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' or 'eeg_epochs' ",
         "(see new_eeg()/epoch_eeg()).", call. = FALSE)
  }
  .eog_check_reference(eeg)

  if (is_epochs && is.null(eeg$data)) {
    stop("ERROR: epochs$data is not loaded - re-run epoch_eeg() with ",
         "preload = TRUE.", call. = FALSE)
  }

  # ========== MATCH CHANNELS BY NAME (shape-independent) ==========

  tgt_idx <- match(model$ch_names, eeg$channels)
  art_idx <- match(model$ch_names_artifact, eeg$channels)
  if (anyNA(tgt_idx)) {
    stop("ERROR: 'eeg' is missing target channel(s) that 'model' was fit ",
         "on: ", paste(model$ch_names[is.na(tgt_idx)], collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(art_idx)) {
    stop("ERROR: 'eeg' is missing EOG channel(s) that 'model' was fit on: ",
         paste(model$ch_names_artifact[is.na(art_idx)], collapse = ", "),
         call. = FALSE)
  }

  if (!is.na(model$reference) && !identical(model$reference, eeg$reference)) {
    warning("The model was fit on data referenced to '", model$reference,
            "' but 'eeg' is referenced to '", eeg$reference, "'. The ",
            "weights depend on the reference, so the correction may be off.",
            call. = FALSE)
  }

  cf <- unname(model$coef_)                 # p x k

  if (is_epochs) {

    # ========== SUBTRACT weight x (per-trial demeaned EOG), TRIAL BY TRIAL ==========

    data <- eeg$data                        # channels x times x trials
    n_trials <- dim(data)[3]
    k <- length(art_idx)
    p <- length(tgt_idx)

    for (n in seq_len(n_trials)) {
      eog_trial <- matrix(data[art_idx, , n], nrow = k)
      eog_trial <- eog_trial - rowMeans(eog_trial, na.rm = TRUE)
      tgt_trial <- matrix(data[tgt_idx, , n], nrow = p)
      data[tgt_idx, , n] <- tgt_trial - cf %*% eog_trial
    }

    # ========== OPTIONALLY PUT THE BASELINE WINDOW BACK AT ZERO ==========

    if (isTRUE(reapply_baseline) && !is.null(eeg$baseline)) {
      data <- .reapply_baseline(data, tgt_idx, eeg$times, eeg$baseline,
                                eeg$baseline_method)
    }

  } else {

    # ========== SUBTRACT weight x (mean-removed EOG), IN BLOCKS ==========
    #
    # Done in blocks of samples so the temporary matrices stay small on long
    # recordings. na.rm = TRUE keeps a NaN in the EOG local to its own
    # samples instead of poisoning the mean (and with it the whole channel).

    data <- eeg$data
    R <- data[art_idx, , drop = FALSE]
    R <- R - rowMeans(R, na.rm = TRUE)

    n <- ncol(data)
    block <- 32768L
    for (start in seq.int(1L, n, by = block)) {
      cols <- start:min(n, start + block - 1L)
      data[tgt_idx, cols] <- data[tgt_idx, cols, drop = FALSE] -
        cf %*% R[, cols, drop = FALSE]
    }
  }

  # ========== BUILD THE RESULT ==========

  out <- eeg
  out$data <- data
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("EOG regression applied", if (is_epochs) " (epochs)" else "",
                ": ", length(tgt_idx), " channel(s) regressed on ",
                paste(model$ch_names_artifact, collapse = ", "),
                if (!is.na(model$reference)) {
                  paste0(" (reference: ", model$reference, ")")
                } else {
                  ""
                }))
  )

  out
}
