# ============================================================================
#              EOG Regression: Ocular Artifact Removal by Regression
# ============================================================================
#
# A fast, deterministic alternative to ICA for removing eye blinks and eye
# movements. Two functions that each do one job, plus a small model object:
#
#   fit_eog_regression()    learn one weight per EEG channel per EOG channel
#   apply_eog_regression()  subtract weight x EOG from each EEG channel
#   new_eog_regression()    constructor / validator for the model object
#
# Ported from MNE-Python's mne/preprocessing/_regress.py (EOGRegression and
# regress_artifact). Design points, all checked against that source:
#
#  - Ordinary least squares, all EOG channels fitted jointly, one fit per
#    target channel. The EOG and each target are mean-removed first, so no
#    intercept is stored and a channel's overall level (DC offset) is left
#    alone when the weights are applied.
#  - Every sample is used and annotations are NOT consulted (MNE's regression
#    has no reject_by_annotation either, unlike ICA fitting).
#  - The EEG must already be re-referenced: the weights depend on the
#    reference (MNE raises an error otherwise, and so does this file).
#  - Fit and apply are separate steps, so a model can be inspected, saved with
#    saveRDS() and reused on another recording. apply matches channels BY NAME
#    (MNE requires identical order), the same way get_sources() / apply_ica()
#    do.
#  - Never touches ICA or ica$exclude: a separate choice from ica1.R /
#    ica_detect.R for the same artifact.
#
# Phase 1 (this file): continuous 'eeg' objects. Planned next: epoched data
# (Gratton et al. 1983 and Croft & Barry 2000, via a subtract_evoked() step),
# a weights topography plot, and a bipolar VEOG/HEOG helper.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# Status: Phase 1 (continuous data) built.
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
#' @param reference \code{eeg$reference}.
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
#' @param eeg An object of class 'eeg'.
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
#' @param eeg An object of class 'eeg'.
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
# fit_eog_regression() - learn the weights
# ----------------------------------------------------------------------------
#' Fit an EOG Regression Model (Ocular Artifact Removal)
#'
#' Works out how strongly each EEG channel follows the EOG channel(s): one
#' weight per EEG channel per EOG channel, fitted together when there are
#' several (say vertical and horizontal). The result is a small model that
#' \code{\link{apply_eog_regression}} uses to subtract that share of the EOG
#' from the EEG. This is a fast, deterministic alternative to ICA for blinks
#' and eye movements, and a port of MNE-Python's \code{EOGRegression.fit()}.
#'
#' @param eeg A continuous \code{eeg} object (see \code{new_eeg()}) that has
#'   already been re-referenced (see \code{\link{eeg_rereference}}).
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
#'   it can be inspected, saved with \code{saveRDS()} and reused.
#'
#' @details
#' \strong{How the weights are found.} For each target channel, ordinary
#' least squares of the channel on all EOG channels together, after taking
#' each signal's average level out, so only the eye-related fluctuations
#' count. No intercept is stored; a channel's overall level is untouched when
#' the weights are applied.
#'
#' \strong{Reference.} The weights depend on the reference, so the EEG must
#' already be re-referenced; an object whose \code{eeg$reference} is still
#' \code{"original"} or \code{"Biosemi CMS/DRL"} is refused. When calling
#' \code{eeg_rereference()}, pass the EOG and status channels to
#' \code{exclude} so they stay out of the average.
#'
#' \strong{Every sample counts.} As in MNE, all samples are used and
#' \code{eeg$annotations} is not consulted. Data must be complete: an NA, NaN
#' or Inf in a used channel is an error.
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
#'   # after filtering, bad-channel handling and re-referencing
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
#' }
#'
#' @seealso \code{\link{apply_eog_regression}}, \code{\link{new_eog_regression}},
#'   \code{\link{eeg_rereference}}, \code{\link{find_bads_eog}},
#'   \code{\link{apply_ica}}
#'
#' @export
fit_eog_regression <- function(eeg, picks = NULL, picks_artifact = NULL) {

  # ========== VALIDATE inputs ==========

  if (inherits(eeg, "eeg_epochs")) {
    stop("ERROR: epoched data ('eeg_epochs') is not supported yet - ",
         "fit_eog_regression() currently works on continuous 'eeg' objects ",
         "only.", call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }
  .eog_check_reference(eeg)

  X <- eeg$data
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("ERROR: eeg$data must be a numeric matrix (channels x time).",
         call. = FALSE)
  }
  if (ncol(X) < 2) {
    stop("ERROR: at least 2 samples are required to fit, got ", ncol(X), ".",
         call. = FALSE)
  }

  # ========== RESOLVE THE EOG (PREDICTOR) CHANNELS ==========

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

  # ========== EOG: CHECK, THEN REMOVE EACH CHANNEL'S MEAN ==========

  R <- X[art_idx, , drop = FALSE]
  if (!all(is.finite(R))) {
    stop("ERROR: non-finite values (NA/NaN/Inf) found in EOG channel(s): ",
         paste(art[!apply(is.finite(R), 1, all)], collapse = ", "),
         ". Regression needs complete data.", call. = FALSE)
  }
  R <- R - rowMeans(R)                      # k x n, every row now sums to ~0

  # ========== TARGETS: CROSS-PRODUCT WITH THE EOG, ONE CHANNEL AT A TIME ==========
  #
  # B[, j] = R %*% (y_j - mean(y_j)), the right-hand side of the normal
  # equations for target j. Done channel by channel (as MNE does) so the full
  # data matrix is never copied.

  B <- matrix(0, nrow = length(art_idx), ncol = length(tgt_idx))
  not_finite <- character(0)
  for (j in seq_along(tgt_idx)) {
    y <- X[tgt_idx[j], ]
    if (!all(is.finite(y))) {
      not_finite <- c(not_finite, tgt[j])
      next
    }
    B[, j] <- R %*% (y - mean(y))
  }
  if (length(not_finite) > 0) {
    stop("ERROR: non-finite values (NA/NaN/Inf) found in target channel(s): ",
         paste(not_finite, collapse = ", "), ". Regression needs complete ",
         "data; interpolate or crop the affected stretch first ",
         "(see annotate_nan()).", call. = FALSE)
  }

  # ========== SOLVE THE NORMAL EQUATIONS: (R R') coef' = R Y' ==========

  G <- tcrossprod(R)                        # k x k
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

  # ========== BUILD THE MODEL ==========

  new_eog_regression(
    coef = t(coef_t),                       # p x k: targets x EOG
    ch_names = tgt,
    ch_names_artifact = art,
    reference = eeg$reference,
    fit_on = "continuous",
    n_samples = ncol(X)
  )
}

# ----------------------------------------------------------------------------
# apply_eog_regression() - subtract weight x EOG from each target channel
# ----------------------------------------------------------------------------
#' Remove EOG Artifacts With a Fitted Regression Model
#'
#' Subtracts, from each target channel, its weight times the EOG signal (with
#' the EOG's average level taken out first), and returns the cleaned
#' \code{eeg} object. A channel's overall level is left alone, only the
#' eye-related fluctuation is removed. The EOG channels themselves and any
#' channel not in the model (status, bad channels left out of the fit, ...)
#' pass through unchanged. A port of MNE-Python's \code{EOGRegression.apply()}.
#'
#' @param model A fitted \code{eeg_eog_regression} object, from
#'   \code{\link{fit_eog_regression}} or \code{\link{new_eog_regression}}.
#' @param eeg A continuous \code{eeg} object to clean. It does not have to be
#'   the recording the model was fit on, but it must contain every channel
#'   named in the model, and it must already be re-referenced.
#'
#' @return A new \code{eeg} object (see \code{\link{new_eeg}}) with the
#'   model's target channels cleaned and a note appended to
#'   \code{preprocessing_history}. Since R does not change arguments in
#'   place, reassign the result (\code{eeg <- apply_eog_regression(model,
#'   eeg)}); the input object itself is left untouched.
#'
#' @details
#' Channels are matched to the model by name, so the channel order in
#' \code{eeg} does not matter (MNE requires identical order). If a channel
#' the model needs is missing, an error lists it. A warning is issued when
#' \code{eeg$reference} differs from the reference the model was fit in,
#' since the weights only strictly apply in that frame.
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
#' }
#'
#' @seealso \code{\link{fit_eog_regression}}, \code{\link{new_eog_regression}},
#'   \code{\link{apply_ica}}
#'
#' @export
apply_eog_regression <- function(model, eeg) {

  # ========== VALIDATE inputs ==========

  if (!inherits(model, "eeg_eog_regression")) {
    stop("ERROR: 'model' must be an object of class 'eeg_eog_regression' ",
         "(see fit_eog_regression()).", call. = FALSE)
  }
  if (inherits(eeg, "eeg_epochs")) {
    stop("ERROR: epoched data ('eeg_epochs') is not supported yet - ",
         "apply_eog_regression() currently works on continuous 'eeg' ",
         "objects only.", call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }
  .eog_check_reference(eeg)

  # ========== MATCH CHANNELS BY NAME ==========

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

  # ========== SUBTRACT weight x (mean-removed EOG) ==========
  #
  # Done in blocks of samples so the temporary matrices stay small on long
  # recordings. na.rm = TRUE keeps a NaN in the EOG local to its own samples
  # instead of poisoning the mean (and with it the whole channel).

  data <- eeg$data
  R <- data[art_idx, , drop = FALSE]
  R <- R - rowMeans(R, na.rm = TRUE)
  cf <- unname(model$coef_)

  n <- ncol(data)
  block <- 32768L
  for (start in seq.int(1L, n, by = block)) {
    cols <- start:min(n, start + block - 1L)
    data[tgt_idx, cols] <- data[tgt_idx, cols, drop = FALSE] -
      cf %*% R[, cols, drop = FALSE]
  }

  # ========== BUILD THE RESULT ==========

  out <- eeg
  out$data <- data
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("EOG regression applied: ", length(tgt_idx),
                " channel(s) regressed on ",
                paste(model$ch_names_artifact, collapse = ", "),
                if (!is.na(model$reference)) {
                  paste0(" (reference: ", model$reference, ")")
                } else {
                  ""
                }))
  )

  out
}
