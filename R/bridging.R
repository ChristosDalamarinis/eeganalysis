# ============================================================================
#                  Electrode Bridging: Detection and Repair
# ============================================================================
#
# Gel caps can short two neighbouring electrodes together through a puddle of
# conductive gel ("bridging"). The two channels then record almost the same
# signal, which is invisible when each channel is looked at on its own. This
# module finds such pairs and, if wanted, rebuilds them.
#
# find_bridged_electrodes() uses the "electrical distance" of Tenke & Kayser
# (2001): the variance of the difference between two channels, in microvolts
# squared, computed for every pair of EEG channels in short windows of the
# 0.5-30 Hz filtered data. A bridged pair has an almost flat difference, so its
# distance is close to zero; the cutoff is found from the data (Greischar et
# al., 2004). interpolate_bridged_electrodes() rebuilds the bridged electrodes
# by spherical-spline interpolation, with a virtual electrode halfway between
# each bridged group as an extra donor so that the brain signal the pair still
# carries is not thrown away. The repair reuses make_interpolation_matrix()
# from R/interpolate.R.
#
# In short, detection filters the EEG channels to 0.5-30 Hz and cuts them into
# consecutive 2-second windows (a window that overlaps a BAD annotation is
# skipped). It computes the distance of every pair in every window, takes the
# cutoff to be the valley between the near-zero pile that bridging produces and
# the rest of the small distances, and calls a pair bridged when it is below
# that cutoff in more than half of the windows (by default). The repair merges
# pairs that share an electrode into groups, gives each group a virtual
# electrode (the group's average signal, at its average position pushed back
# onto the scalp sphere) and rebuilds all bridged electrodes together with one
# interpolation matrix.
#
# Author: Christos Dalamarinis
# Date: Oct - 2026
# ============================================================================

#' Find Bridged Electrodes
#'
#' Looks for pairs of EEG electrodes that have been shorted together by a
#' puddle of conductive gel ("bridging"), a common problem with gel caps such
#' as BioSemi's. Two bridged electrodes record almost the same signal, so the
#' difference between them is nearly flat. The function measures that with the
#' "electrical distance" (the variance of the difference between two channels,
#' in microvolts squared) for every pair of channels in short windows of the
#' filtered data, and calls a pair bridged when its distance is below a
#' data-driven cutoff in most windows - by default, in more than half of them
#' (Tenke & Kayser, 2001; Greischar et al., 2004).
#'
#' @param eeg_obj An object of class \code{'eeg'} (continuous data). Only
#'   channels classified as \code{"eeg"} in \code{eeg_obj$channel_types} and
#'   not listed in \code{eeg_obj$bads} are examined. Epoched data are not
#'   supported: run this before epoching.
#'
#' @param lm_cutoff Numeric, positive. Only electrical distances below this
#'   value (in microvolts squared) are searched for the pile of near-zero
#'   values that bridging produces. A value of 16 is conservative (it is based
#'   on the distributions reported by Greischar et al., 2004); EEGLAB uses 5.
#'   Default \code{16}.
#'
#' @param epoch_threshold Numeric from 0 up to (not including) 1. A pair is
#'   called bridged when its distance is below the data-driven cutoff in MORE
#'   than this proportion of the windows. It is also used as a first check: if
#'   the average number of distances below \code{lm_cutoff} per window is
#'   smaller than this, nothing is searched. Default \code{0.5}.
#'
#' @param l_freq,h_freq Numeric, the band (in Hz) the data are filtered to
#'   before the distances are computed. Defaults \code{0.5} and \code{30}.
#'   \code{h_freq} must be below half the sampling rate.
#'
#' @param epoch_duration Numeric, positive, length in seconds of the windows
#'   the recording is cut into. Default \code{2}.
#'
#' @param verbose Logical. If \code{TRUE} (default), print a one-line summary.
#'
#' @return An object of class \code{eeg_bridges}, a list with:
#'   \describe{
#'     \item{bridged}{Data frame, one row per bridged pair: \code{channel_1}
#'       and \code{channel_2} (in the order of \code{eeg_obj$channels}),
#'       \code{fraction_below} (proportion of windows with a distance below
#'       the cutoff) and \code{median_ed} (median distance over the windows, in
#'       microvolts squared). Zero rows if nothing was found.}
#'     \item{ed}{Numeric matrix of electrical distances in microvolts squared:
#'       one row per channel pair (see \code{ed_pairs}), one column per window
#'       used.}
#'     \item{ed_pairs}{Data frame with \code{channel_1} and \code{channel_2},
#'       saying which pair each row of \code{ed} belongs to.}
#'     \item{local_minimum}{The data-driven cutoff in microvolts squared (the
#'       valley between the near-zero pile of bridged pairs and the rest), or
#'       \code{NA} if the first check found too few small distances to search.}
#'     \item{n_windows, n_windows_dropped}{Windows used, and windows left out
#'       because they overlap a \code{"BAD"} annotation.}
#'     \item{channels}{Names of the channels examined.}
#'     \item{lm_cutoff, epoch_threshold, l_freq, h_freq, epoch_duration}{The
#'       settings used.}
#'   }
#'   \code{eeg_obj} itself is not changed, and nothing is written to
#'   \code{eeg_obj$bads}: bridged electrodes still carry brain signal.
#'
#' @details
#' \strong{Steps.} (1) The examined channels are filtered to
#' \code{l_freq}-\code{h_freq} with \code{\link{eeg_bandpass}}. (2) The
#' recording is cut into consecutive windows of \code{epoch_duration} seconds;
#' a partial window at the end is dropped, and so is any window that overlaps a
#' row of \code{eeg_obj$annotations} whose description starts with "BAD" (any
#' case). (3) For every pair of channels and every window, the electrical
#' distance is the variance of the difference of the two channels. (4) If the
#' average number of distances below \code{lm_cutoff} per window is smaller
#' than \code{epoch_threshold}, no bridging is suspected and the search stops.
#' Otherwise a smooth density curve is fitted to the distances below
#' \code{lm_cutoff} and its valley (\code{local_minimum}) is used as the
#' cutoff. (5) A pair is bridged if its distance is below that cutoff in more
#' than \code{epoch_threshold} of the windows.
#'
#' \strong{Unit.} The data are assumed to be in microvolts (as read by
#' \code{\link{read_bdf_native}}), so the distances are in microvolts squared.
#'
#' \strong{Why the reference does not matter.} A re-reference subtracts the
#' same signal from both channels of a pair, so it cancels in their
#' difference: the result is the same before and after re-referencing.
#'
#' \strong{Chains.} Electrodes can be bridged in a chain (A-B and B-C); then all
#' pairs among the connected electrodes are usually reported.
#' \code{\link{interpolate_bridged_electrodes}} merges them into one group.
#'
#' \strong{Order of steps.} Run this early, on the raw recording, before
#' re-referencing, ICA or CSD: two bridged channels carry one signal twice,
#' which also removes a degree of freedom from ICA.
#'
#' @examples
#' \dontrun{
#'   bridges <- find_bridged_electrodes(eeg)
#'   bridges                 # print method: channels, cutoffs, pairs
#'   bridges$bridged         # the pairs as a data frame
#'
#'   # Rebuild them (needs a montage)
#'   eeg <- set_montage(eeg, create_montage())
#'   eeg <- interpolate_bridged_electrodes(eeg, bridges)
#' }
#'
#' @seealso \code{\link{interpolate_bridged_electrodes}},
#'   \code{\link{find_bad_channels}}, \code{\link{eeg_bandpass}}
#'
#' @export
find_bridged_electrodes <- function(eeg_obj,
                                    lm_cutoff       = 16,
                                    epoch_threshold = 0.5,
                                    l_freq          = 0.5,
                                    h_freq          = 30,
                                    epoch_duration  = 2,
                                    verbose         = TRUE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg' (continuous data). ",
         "Epoched data are not supported - run this before epoching.", call. = FALSE)
  }

  if (!is.matrix(eeg_obj$data)) {
    stop("ERROR: eeg_obj$data must be a numeric matrix (channels x timepoints).",
         call. = FALSE)
  }

  is_num1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)

  if (!is_num1(lm_cutoff) || lm_cutoff <= 0) {
    stop("ERROR: 'lm_cutoff' must be a single positive number.", call. = FALSE)
  }

  if (!is_num1(epoch_threshold) || epoch_threshold < 0 || epoch_threshold >= 1) {
    stop("ERROR: 'epoch_threshold' must be a single number from 0 up to (but not including) 1.",
         call. = FALSE)
  }

  if (!is_num1(l_freq) || !is_num1(h_freq) || l_freq <= 0 || l_freq >= h_freq) {
    stop("ERROR: 'l_freq' and 'h_freq' must be single numbers with 0 < l_freq < h_freq.",
         call. = FALSE)
  }

  sr <- eeg_obj$sampling_rate
  if (h_freq >= sr / 2) {
    stop("ERROR: 'h_freq' must be below the Nyquist frequency (", sr / 2, " Hz).",
         call. = FALSE)
  }

  if (!is_num1(epoch_duration) || epoch_duration <= 0) {
    stop("ERROR: 'epoch_duration' must be a single positive number (seconds).",
         call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("ERROR: 'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  # ========== PICK CHANNELS ==========

  picks <- which(eeg_obj$channel_types == "eeg" & !(eeg_obj$channels %in% eeg_obj$bads))
  ch    <- eeg_obj$channels[picks]

  if (length(ch) < 2) {
    stop("ERROR: At least 2 good EEG channels are needed to look for bridges (found ",
         length(ch), ").", call. = FALSE)
  }

  if (!all(is.finite(eeg_obj$data[picks, , drop = FALSE]))) {
    stop("ERROR: EEG data contains NA or non-finite values (the filter would spread them). ",
         "Remove or interpolate them first.", call. = FALSE)
  }

  # ========== FILTER ==========

  filtered <- eeg_bandpass(eeg_obj, l_freq = l_freq, h_freq = h_freq,
                           channels = ch, verbose = FALSE)$data[picks, , drop = FALSE]

  # ========== WINDOWS ==========

  win_len <- round(epoch_duration * sr)
  n_all   <- floor(ncol(filtered) / win_len)

  if (n_all < 1) {
    stop("ERROR: The recording is too short for even one window of ", epoch_duration,
         " s (it lasts ", round(ncol(filtered) / sr, 2), " s).", call. = FALSE)
  }

  starts <- (seq_len(n_all) - 1) * win_len + 1
  keep   <- .bridge_keep_windows(eeg_obj$annotations, starts, win_len, sr)

  if (!any(keep)) {
    stop("ERROR: Every window overlaps a BAD annotation - nothing left to analyse.",
         call. = FALSE)
  }

  # ========== ELECTRICAL DISTANCES ==========

  ed  <- .bridge_ed(filtered, starts[keep], win_len)
  idx <- .bridge_pairs_index(length(ch))
  ed_pairs <- data.frame(channel_1 = ch[idx[, "i"]], channel_2 = ch[idx[, "j"]],
                         stringsAsFactors = FALSE)

  # ========== CUTOFF AND PAIRS ==========

  local_minimum <- .bridge_cutoff(ed, lm_cutoff, epoch_threshold)

  if (is.na(local_minimum)) {
    hit            <- integer(0)
    fraction_below <- rep(NA_real_, nrow(ed))
  } else {
    fraction_below <- rowMeans(ed < local_minimum)
    hit            <- which(fraction_below > epoch_threshold)      # strictly greater
  }

  bridged <- data.frame(
    channel_1        = ed_pairs$channel_1[hit],
    channel_2        = ed_pairs$channel_2[hit],
    fraction_below   = fraction_below[hit],
    median_ed        = apply(ed[hit, , drop = FALSE], 1, median),
    stringsAsFactors = FALSE)

  result <- structure(
    list(bridged           = bridged,
         ed                = ed,
         ed_pairs          = ed_pairs,
         local_minimum     = local_minimum,
         n_windows         = sum(keep),
         n_windows_dropped = sum(!keep),
         channels          = ch,
         lm_cutoff         = lm_cutoff,
         epoch_threshold   = epoch_threshold,
         l_freq            = l_freq,
         h_freq            = h_freq,
         epoch_duration    = epoch_duration),
    class = "eeg_bridges")

  # ========== SUMMARY ==========

  if (verbose) {
    message("find_bridged_electrodes(): ", length(ch), " EEG channel(s), ", nrow(ed),
            " pair(s), ", sum(keep), " window(s) of ", epoch_duration, " s",
            if (sum(!keep) > 0) paste0(" (", sum(!keep), " dropped because of BAD annotations)") else "",
            ". ",
            if (is.na(local_minimum)) {
              "Too few small distances to suspect bridging. "
            } else {
              paste0("Data-driven cutoff: ", format(round(local_minimum, 3), nsmall = 3), " uV^2. ")
            },
            if (nrow(bridged) == 0) {
              "No bridged electrodes found."
            } else {
              paste0(nrow(bridged), " bridged pair(s): ",
                     paste(bridged$channel_1, bridged$channel_2, sep = " - ", collapse = ", "), ".")
            })
  }

  result
}


#' Print Method for Bridging Results
#'
#' Displays what \code{\link{find_bridged_electrodes}} examined (channels,
#' windows, filter), the cutoffs it used and the bridged pairs it found. The
#' full table of electrical distances is \code{x$ed}.
#'
#' @param x An object of class \code{eeg_bridges}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x} (standard R print method convention).
#'
#' @examples
#' \dontrun{
#'   bridges <- find_bridged_electrodes(eeg)
#'   print(bridges)  # Calls this method automatically
#' }
#'
#' @export
print.eeg_bridges <- function(x, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Electrode Bridging Check\n")
  cat(strrep("=", 70), "\n\n")

  cat("DATA:\n")
  cat("  EEG channels examined: ", length(x$channels), " (", nrow(x$ed), " pairs)\n", sep = "")
  cat("  Windows used:          ", x$n_windows, " of ", x$epoch_duration, " s",
      if (x$n_windows_dropped > 0) {
        paste0(" (", x$n_windows_dropped, " dropped: BAD annotation)")
      } else "", "\n", sep = "")
  cat("  Filter:                ", x$l_freq, " - ", x$h_freq, " Hz\n", sep = "")

  cat("\nCUTOFFS:\n")
  cat("  lm_cutoff:             ", x$lm_cutoff, " uV^2\n", sep = "")
  cat("  Data-driven cutoff:    ",
      if (is.na(x$local_minimum)) {
        "not needed (too few small distances)"
      } else {
        paste0(format(round(x$local_minimum, 3), nsmall = 3), " uV^2")
      }, "\n", sep = "")
  cat("  epoch_threshold:       ", x$epoch_threshold, "\n", sep = "")

  cat("\nBRIDGED PAIRS:\n")
  if (nrow(x$bridged) == 0) {
    cat("  none found\n")
  } else {
    b <- x$bridged
    for (i in seq_len(nrow(b))) {
      cat(sprintf("  %-8s - %-8s below the cutoff in %3.0f%% of windows (median distance %.3f uV^2)\n",
                  b$channel_1[i], b$channel_2[i], 100 * b$fraction_below[i], b$median_ed[i]))
    }
  }
  cat("\n")

  invisible(x)
}


#' Interpolate Bridged Electrodes
#'
#' Rebuilds electrodes that are shorted together by gel (see
#' \code{\link{find_bridged_electrodes}}). Because bridged electrodes still
#' contain brain signal - it is only smeared between them - they are not simply
#' thrown away: for each group of bridged electrodes a "virtual" electrode is
#' placed halfway between them, carrying their average signal, and the bridged
#' electrodes are then rebuilt by spherical-spline interpolation from all good
#' channels plus these virtual electrodes. Pairs that share an electrode are
#' first merged into groups, and all groups are rebuilt together in one step.
#'
#' @param eeg_obj An object of class \code{'eeg'}, with a montage attached (see
#'   \code{\link{set_montage}}).
#'
#' @param bridged Either the result of \code{\link{find_bridged_electrodes}}, or
#'   a table (data frame or matrix) with two columns of channel names, one
#'   bridged pair per row. Pairs that share a channel are merged into one group
#'   (A-B and B-C make the group A, B, C).
#'
#' @param bad_limit Integer, at least 1. The largest number of electrodes
#'   allowed in one bridged group. A larger group is an error, because
#'   interpolating a large area from the electrodes around it is inaccurate.
#'   Default \code{4}.
#'
#' @param origin Numeric vector of length 3 (x, y, z, in the same units as the
#'   montage positions) giving the centre of the sphere the electrodes lie on.
#'   Default \code{c(0, 0, 0)}, correct for any montage built with
#'   \code{\link{create_montage}} - see \code{\link{interpolate_bads}}.
#'
#' @return \code{eeg_obj} with the data of the bridged channels replaced and a
#'   step appended to \code{$preprocessing_history}. Everything else, including
#'   \code{$bads}, is unchanged (the repaired channels are not marked bad). If
#'   \code{bridged} contains no pairs, \code{eeg_obj} is returned unchanged,
#'   with a message.
#'
#' @details
#' \strong{Donors.} The channels the bridged electrodes are rebuilt from are the
#' EEG channels that have a montage position, are not part of a bridged group
#' and are not listed in \code{eeg_obj$bads}, plus one virtual electrode per
#' group. All groups are rebuilt together with one interpolation matrix, so
#' every virtual electrode helps every bridged electrode. Channels already
#' marked bad are deliberately not used as donors, because their signal cannot
#' be trusted - the same rule as in \code{\link{interpolate_bads}}.
#'
#' \strong{Virtual electrode.} Its position is the average of the group's
#' positions, pushed back onto the scalp sphere at the group's average radius;
#' its signal is the average of the group's channels, sample by sample.
#'
#' \strong{Requirements.} Every bridged channel must be an EEG channel with a
#' position in the montage and must not be marked bad (use
#' \code{\link{interpolate_bads}} for those). The repaired channels are an
#' estimate from their surroundings, not recorded data.
#'
#' @examples
#' \dontrun{
#'   eeg <- set_montage(eeg, create_montage())
#'   bridges <- find_bridged_electrodes(eeg)
#'   eeg <- interpolate_bridged_electrodes(eeg, bridges)
#'
#'   # Or name the pairs yourself, one pair per row
#'   eeg <- interpolate_bridged_electrodes(eeg, rbind(c("Cz", "CPz"), c("P3", "P5")))
#' }
#'
#' @seealso \code{\link{find_bridged_electrodes}}, \code{\link{interpolate_bads}},
#'   \code{\link{create_montage}}, \code{\link{set_montage}}
#'
#' @export
interpolate_bridged_electrodes <- function(eeg_obj,
                                           bridged,
                                           bad_limit = 4,
                                           origin = c(0, 0, 0)) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }

  if (!is.matrix(eeg_obj$data)) {
    stop("ERROR: eeg_obj$data must be a numeric matrix (channels x timepoints).",
         call. = FALSE)
  }

  if (!is.numeric(bad_limit) || length(bad_limit) != 1 || !is.finite(bad_limit) ||
      bad_limit < 1 || bad_limit != round(bad_limit)) {
    stop("ERROR: 'bad_limit' must be a single whole number of at least 1.", call. = FALSE)
  }

  if (!is.numeric(origin) || length(origin) != 3 || !all(is.finite(origin))) {
    stop("ERROR: 'origin' must be a numeric vector of length 3 (x, y, z).", call. = FALSE)
  }

  # The pairs, with the zero-row trap handled first (as.matrix() of a zero-row
  # data frame is a logical matrix).
  pairs <- NULL
  if (inherits(bridged, "eeg_bridges")) {
    pairs <- if (nrow(bridged$bridged) == 0) {
      matrix(character(0), 0, 2)
    } else {
      as.matrix(bridged$bridged[, c("channel_1", "channel_2")])
    }
  } else if ((is.matrix(bridged) || is.data.frame(bridged)) && ncol(bridged) == 2) {
    pairs <- if (nrow(bridged) == 0) matrix(character(0), 0, 2) else as.matrix(bridged)
  }

  if (is.null(pairs) || !is.character(pairs)) {
    stop("ERROR: 'bridged' must be the result of find_bridged_electrodes() or a ",
         "two-column table of channel names.", call. = FALSE)
  }

  if (nrow(pairs) > 0 && any(pairs[, 1] == pairs[, 2])) {
    stop("ERROR: Each row of 'bridged' must name two different channels.", call. = FALSE)
  }

  if (nrow(pairs) == 0) {
    message("interpolate_bridged_electrodes(): no bridged pairs given - returning eeg_obj unchanged.")
    return(eeg_obj)
  }

  named <- unique(as.vector(pairs))

  eeg_ch <- eeg_obj$channels[eeg_obj$channel_types == "eeg"]
  not_eeg <- named[!(named %in% eeg_ch)]
  if (length(not_eeg) > 0) {
    stop("ERROR: These channels in 'bridged' are not EEG channels of eeg_obj: ",
         paste(not_eeg, collapse = ", "), ".", call. = FALSE)
  }

  if (is.null(eeg_obj$montage) || !inherits(eeg_obj$montage, "montage")) {
    stop("ERROR: No montage attached - repairing bridged electrodes needs channel ",
         "scalp positions. Attach one with set_montage().", call. = FALSE)
  }

  no_pos <- named[!(named %in% eeg_obj$montage$channels)]
  if (length(no_pos) > 0) {
    stop("ERROR: Bridged channel(s) with no position in the montage: ",
         paste(no_pos, collapse = ", "), ".", call. = FALSE)
  }

  already_bad <- named[named %in% eeg_obj$bads]
  if (length(already_bad) > 0) {
    stop("ERROR: These bridged channels are already marked bad: ",
         paste(already_bad, collapse = ", "),
         ". Interpolate them with interpolate_bads() instead.", call. = FALSE)
  }

  # ========== GROUPS ==========

  # In channel order, so the history entry and the donor order do not depend on
  # how the pairs were listed.
  groups <- .bridge_groups(pairs)
  groups <- lapply(groups, function(g) g[order(match(g, eeg_obj$channels))])
  groups <- groups[order(vapply(groups, function(g) match(g[1], eeg_obj$channels), numeric(1)))]
  members <- unlist(groups)

  too_big <- which(vapply(groups, length, integer(1)) > bad_limit)
  if (length(too_big) > 0) {
    g <- groups[[too_big[1]]]
    stop("ERROR: The channels ", paste(g, collapse = ", "), " are bridged together and form a group of ",
         length(g), " electrodes (limit: ", bad_limit, ") - interpolation would be inaccurate. ",
         "Raise 'bad_limit' to override.", call. = FALSE)
  }

  # ========== DONORS AND POSITIONS ==========

  donors <- setdiff(eeg_ch[eeg_ch %in% eeg_obj$montage$channels & !(eeg_ch %in% eeg_obj$bads)],
                    members)

  if (length(donors) == 0) {
    stop("ERROR: No good EEG channels with a montage position are available to interpolate from.",
         call. = FALSE)
  }

  mont_pos <- eeg_obj$montage$positions
  get_pos  <- function(chs) {
    p <- as.matrix(mont_pos[match(chs, mont_pos$channel), c("x", "y", "z")])
    sweep(p, 2, as.numeric(origin), "-")
  }
  pos_donors  <- get_pos(donors)
  pos_members <- get_pos(members)

  radii <- sqrt(rowSums(rbind(pos_donors, pos_members)^2))
  if (max(abs(radii / mean(radii) - 1)) > 0.1) {
    warning("interpolate_bridged_electrodes(): channel positions are not close to ",
            "spherical around 'origin' - results may be inaccurate.",
            call. = FALSE, immediate. = TRUE)
  }

  pos_virtual <- t(vapply(groups, function(g) .bridge_centroid(get_pos(g)), numeric(3)))

  no_place <- which(is.na(pos_virtual[, 1]))
  if (length(no_place) > 0) {
    stop("ERROR: Cannot place a virtual electrode for ",
         paste(groups[[no_place[1]]], collapse = ", "),
         ": they lie on opposite sides of the head.", call. = FALSE)
  }

  # ========== REBUILD ==========

  donor_idx  <- match(donors, eeg_obj$channels)
  member_idx <- match(members, eeg_obj$channels)
  virt_data  <- t(vapply(groups, function(g) {
    colMeans(eeg_obj$data[match(g, eeg_obj$channels), , drop = FALSE])
  }, numeric(ncol(eeg_obj$data))))

  # One matrix for all groups; the virtual electrodes are the last columns.
  W        <- make_interpolation_matrix(rbind(pos_donors, pos_virtual), pos_members)
  n_donors <- length(donors)

  eeg_obj$data[member_idx, ] <-
    W[, seq_len(n_donors), drop = FALSE] %*% eeg_obj$data[donor_idx, , drop = FALSE] +
    W[, n_donors + seq_along(groups), drop = FALSE] %*% virt_data

  group_text <- vapply(groups, function(g) paste(g, collapse = "+"), character(1))
  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history, list(paste0(
    "interpolate_bridged_electrodes(): rebuilt ", length(members), " bridged channel(s) in ",
    length(groups), " group(s) (", paste(group_text, collapse = ", "), ") from ",
    n_donors, " good channel(s) plus ", length(groups), " virtual electrode(s).")))

  eeg_obj
}


#' Index of All Channel Pairs (internal)
#'
#' Internal helper listing every pair of channels \code{i < j} of \code{n}
#' channels, ordered by \code{i} and then \code{j} (row by row: all pairs of
#' channel 1 first, then all pairs of channel 2, and so on - not the
#' column-by-column order of \code{upper.tri()}).
#'
#' @param n Integer, number of channels, at least 2.
#'
#' @return Integer matrix with \code{n * (n - 1) / 2} rows and columns
#'   \code{i} and \code{j}.
#'
#' @seealso \code{\link{find_bridged_electrodes}}
#' @keywords internal
.bridge_pairs_index <- function(n) {
  i <- rep(seq_len(n - 1), times = (n - 1):1)
  j <- unlist(lapply(seq_len(n - 1), function(a) (a + 1):n))
  cbind(i = i, j = j)
}


#' Electrical Distances for All Pairs and Windows (internal)
#'
#' Internal helper computing, for every pair of channels and every window, the
#' variance of the difference between the two channels (the "electrical
#' distance"). It uses var(a - b) = var(a) + var(b) - 2 cov(a, b), so one
#' matrix product per window gives all pairs at once.
#'
#' @param data Numeric matrix, channels x time, already filtered.
#' @param win_start Integer vector, 1-based first sample of each window.
#' @param win_len Integer, window length in samples.
#'
#' @return Numeric matrix, one row per pair (in the order of
#'   \code{.bridge_pairs_index}), one column per window, in squared data units
#'   (microvolts squared). Values are clamped at 0 (rounding can make a distance
#'   a hair negative).
#'
#' @seealso \code{\link{find_bridged_electrodes}}
#' @keywords internal
.bridge_ed <- function(data, win_start, win_len) {
  n_ch <- nrow(data)
  idx  <- .bridge_pairs_index(n_ch)
  ed   <- matrix(NA_real_, nrow(idx), length(win_start))

  for (w in seq_along(win_start)) {
    xw <- data[, win_start[w] + seq_len(win_len) - 1L, drop = FALSE]
    xc <- xw - rowMeans(xw)
    cv <- tcrossprod(xc) / win_len          # var(a - b) = var(a) + var(b) - 2 cov(a, b)
    d  <- diag(cv)
    ed[, w] <- (outer(d, d, "+") - 2 * cv)[idx]
  }

  pmax(ed, 0)
}


#' Windows Not Overlapping a BAD Annotation (internal)
#'
#' Internal helper telling which windows to keep: a window is dropped if it
#' overlaps a row of the annotations table whose description starts with "BAD"
#' (any case). The window spans \code{[start, start + win_len / sampling_rate)}
#' seconds, the annotation \code{[onset, onset + duration)}.
#'
#' @param annotations Data frame with \code{onset}, \code{duration} and
#'   \code{description} (see \code{eeg_obj$annotations}), or \code{NULL}.
#' @param win_start Integer vector, 1-based first sample of each window.
#' @param win_len Integer, window length in samples.
#' @param sampling_rate Numeric, samples per second.
#'
#' @return Logical vector, one entry per window, \code{TRUE} = keep.
#'
#' @seealso \code{\link{find_bridged_electrodes}}
#' @keywords internal
.bridge_keep_windows <- function(annotations, win_start, win_len, sampling_rate) {
  keep <- rep(TRUE, length(win_start))

  if (is.null(annotations) || nrow(annotations) == 0) return(keep)
  bad <- annotations[startsWith(toupper(annotations$description), "BAD"), , drop = FALSE]
  if (nrow(bad) == 0) return(keep)

  t_start <- (win_start - 1) / sampling_rate
  t_end   <- t_start + win_len / sampling_rate      # the window spans [t_start, t_end)
  for (k in seq_len(nrow(bad))) {
    keep <- keep & !(bad$onset[k] < t_end & (bad$onset[k] + bad$duration[k]) > t_start)
  }
  keep
}


#' Data-Driven Cutoff for Bridging (internal)
#'
#' Internal helper finding the electrical-distance cutoff below which a pair
#' counts as bridged. If the average number of distances below
#' \code{lm_cutoff} per window is smaller than \code{epoch_threshold}, nothing
#' is suspicious and \code{NA} is returned. Otherwise a Gaussian kernel
#' density (bandwidth sd * n^(-1/5)) is fitted to the distances below
#' \code{lm_cutoff} and its minimum over (0, \code{lm_cutoff}) is returned:
#' the valley between the near-zero pile of bridged pairs and the rest.
#'
#' @param ed Numeric matrix of electrical distances (pairs x windows).
#' @param lm_cutoff Numeric, only distances below this are used.
#' @param epoch_threshold Numeric, the first-check threshold.
#'
#' @return A single number, or \code{NA_real_}.
#'
#' @seealso \code{\link{find_bridged_electrodes}}
#' @keywords internal
.bridge_cutoff <- function(ed, lm_cutoff, epoch_threshold) {
  below <- ed[ed < lm_cutoff]

  if (length(below) / ncol(ed) < epoch_threshold) return(NA_real_)

  if (length(below) < 2 || sd(below) == 0) {
    stop("ERROR: Too few small electrical distances to estimate the cutoff (found ",
         length(below), " below lm_cutoff, need at least 2 that differ). ",
         "Use a longer recording.", call. = FALSE)
  }

  bw  <- sd(below) * length(below)^(-1 / 5)         # SciPy's default (Scott's rule)
  kde <- function(x) vapply(x, function(v) mean(stats::dnorm(v, below, bw)), numeric(1))

  stats::optimize(kde, interval = c(0, lm_cutoff), tol = 1e-8)$minimum
}


#' Group Bridged Pairs (internal)
#'
#' Internal helper merging pairs that share a channel into groups (A-B and B-C
#' make the group A, B, C).
#'
#' @param pairs Character matrix with two columns, one pair per row.
#'
#' @return A list of character vectors, one per group (members and groups in
#'   the order they were first met).
#'
#' @seealso \code{\link{interpolate_bridged_electrodes}}
#' @keywords internal
.bridge_groups <- function(pairs) {
  groups <- list()
  for (k in seq_len(nrow(pairs))) {
    a   <- pairs[k, 1]
    b   <- pairs[k, 2]
    hit <- which(vapply(groups, function(g) any(c(a, b) %in% g), logical(1)))
    merged <- unique(c(a, b, unlist(groups[hit])))
    if (length(hit) > 0) groups <- groups[-hit]     # groups[-integer(0)] would be empty
    groups <- c(groups, list(merged))
  }
  groups
}


#' Position of the Virtual Electrode of a Bridged Group (internal)
#'
#' Internal helper placing the virtual electrode of a group: the mean of the
#' member positions, pushed back onto the sphere at the group's mean radius.
#' (The plain mean of points on a sphere lies inside the sphere, which would
#' put the virtual electrode inside the head; moving it outwards along its own
#' direction, until its distance from the centre is the members' average
#' distance, keeps it on the scalp.)
#'
#' @param pos Numeric matrix, \code{n x 3}, positions centred on the sphere's
#'   origin.
#'
#' @return Numeric vector of length 3, or three \code{NA} if the electrodes lie
#'   on opposite sides of the head (their mean is at the centre).
#'
#' @seealso \code{\link{interpolate_bridged_electrodes}}
#' @keywords internal
.bridge_centroid <- function(pos) {
  mid <- colMeans(pos)
  if (sqrt(sum(mid^2)) < 1e-9 * max(sqrt(rowSums(pos^2)))) return(rep(NA_real_, 3))
  mid / sqrt(sum(mid^2)) * mean(sqrt(rowSums(pos^2)))
}
