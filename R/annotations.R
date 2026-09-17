#' ============================================================================
#'                       Time-Based Data Annotations
#' ============================================================================
#'
#' This module detects bad TIME RANGES in EEG data automatically and writes
#' them to eeg_obj$annotations (see R/eeg_class.R) - the time-domain sibling
#' of eeg_obj$bads (see R/bad_channels.R). Where `bads` disqualifies a whole
#' channel, an annotation disqualifies only a stretch of time, so downstream
#' steps (epoching, ICA fitting) have something to check against before
#' touching a given moment - without throwing away the whole channel over it.
#'
#' Four independent detectors:
#'   - annotate_amplitude() - consecutive-sample jumps that are too small
#'     (flat/disconnected) or too large (spikes), sustained for longer than
#'     a minimum duration. A channel bad for most of the recording is added
#'     to eeg_obj$bads instead of annotated (mirrors find_bad_channels()) -
#'     that's a channel problem, not a moment problem.
#'   - annotate_muscle()    - jaw-clench/EMG bursts, via a 110-140 Hz
#'     envelope-and-z-score method.
#'   - annotate_nan()       - amplifier dropouts (NA runs), per channel.
#'   - annotate_break()     - dead time between experimental blocks: gaps
#'     between existing annotations (or, with use_events = TRUE, between
#'     eeg_obj$events) that are long enough to count as a break.
#'
#' Unlike eeg_obj$bads, there is no "fixer" function for annotations - a bad
#' moment in time can't be rebuilt from its neighbors the way a bad channel
#' can be interpolated from its spatial neighbors (see R/interpolate.R). All
#' a detector here can do is mark a stretch of time for downstream steps to
#' skip.
#'
#' eeg_obj$annotations is a data frame with columns:
#'   onset       - start time in seconds from the beginning of the recording
#'   duration    - length in seconds
#'   description - reason, e.g. "BAD_muscle", "BAD_flat", "BAD_peak",
#'                 "BAD_NAN", "BAD_break"
#'   channel     - channel name, if the row is channel-specific (only
#'                 annotate_nan() produces these); NA otherwise
#'
#' R/eeg_class.R does not define this field yet, so every exported function
#' below initializes eeg_obj$annotations to a zero-row frame the first time
#' it is touched (.ensure_annotations()) - this file works standalone. When
#' eeg_class.R is next updated, giving new_eeg() the same empty-data-frame
#' default it already gives `bads` would make that initialization redundant
#' (harmless either way - it only fires when the field is NULL).
#'
#' Author: Christos Dalamarinis
#' Date: Sep - 2026
#' Status: Ready
#' ============================================================================
#
#
# ============================================================================
#                     SHARED INTERNAL HELPERS
# ============================================================================

#' Empty Annotations Data Frame (internal)
#'
#' Zero-row \code{onset, duration, description, channel} data frame shared by
#' every detector in this file, so none of them have to special-case a
#' \code{NULL}/empty result.
#'
#' @return A zero-row data frame with columns \code{onset, duration}
#'   (numeric) and \code{description, channel} (character).
#' @keywords internal
.empty_annotations <- function() {
  data.frame(onset = numeric(0), duration = numeric(0),
             description = character(0), channel = character(0),
             stringsAsFactors = FALSE)
}

#' Ensure eeg_obj$annotations Exists (internal)
#'
#' \code{eeg_obj$annotations} is not yet a field \code{\link{new_eeg}}
#' initializes (see file header), so every exported function in this file
#' calls this first to default it to a zero-row frame if missing.
#'
#' @param eeg_obj An object of class 'eeg'.
#' @return \code{eeg_obj}, with \code{$annotations} guaranteed non-NULL.
#' @keywords internal
.ensure_annotations <- function(eeg_obj) {
  if (is.null(eeg_obj$annotations)) {
    eeg_obj$annotations <- .empty_annotations()
  }
  eeg_obj
}

#' Append New Annotation Rows, Kept Sorted by Onset (internal)
#'
#' Combines \code{new_rows} onto \code{eeg_obj$annotations} and re-sorts the
#' result by \code{onset}, so the table always reads in chronological order.
#'
#' @param eeg_obj An object of class 'eeg' (already passed through
#'   \code{\link{.ensure_annotations}}).
#' @param new_rows A data frame with the same columns as
#'   \code{eeg_obj$annotations} (see \code{\link{.empty_annotations}}).
#' @return \code{eeg_obj}, with \code{$annotations} updated.
#' @keywords internal
.append_annotations <- function(eeg_obj, new_rows) {
  combined <- rbind(eeg_obj$annotations, new_rows)
  combined <- combined[order(combined$onset), , drop = FALSE]
  rownames(combined) <- NULL
  eeg_obj$annotations <- combined
  eeg_obj
}

#' Resolve a `channels` Argument to Row Indices (internal)
#'
#' Shared channel-picking logic - \code{NULL} falls back to
#' \code{default_idx}, a character vector is matched against
#' \code{eeg_obj$channels}, a numeric vector is used directly. Mirrors the
#' inline logic already used by \code{eeg_bandpass()} (see R/filter1.R).
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param channels \code{NULL}, a character vector of channel names, or an
#'   integer vector of channel indices.
#' @param default_idx Integer vector to use when \code{channels} is
#'   \code{NULL}.
#' @return Integer vector of row indices into \code{eeg_obj$data}.
#' @keywords internal
.resolve_channel_idx <- function(eeg_obj, channels, default_idx) {

  n_channels <- length(eeg_obj$channels)

  if (is.null(channels)) {
    return(default_idx)
  }
  if (is.character(channels)) {
    idx <- match(channels, eeg_obj$channels)
    if (any(is.na(idx))) {
      stop("Channel(s) not found: ",
           paste(channels[is.na(idx)], collapse = ", "), call. = FALSE)
    }
    return(idx)
  }
  if (is.numeric(channels)) {
    idx <- as.integer(channels)
    if (any(idx < 1L) || any(idx > n_channels)) {
      stop("Channel indices out of range [1, ", n_channels, "].", call. = FALSE)
    }
    return(idx)
  }
  stop("'channels' must be NULL, a character vector, or an integer vector.",
       call. = FALSE)
}

#' Group a Boolean Mask into Onset/Offset Sample Pairs (internal)
#'
#' Finds the rising and falling edges of \code{TRUE} runs in \code{mask}.
#' \code{onsets} are 1-based and inclusive (the first \code{TRUE} sample of
#' a run); \code{offsets} are 1-based and exclusive, i.e.
#' \code{mask[onset:(offset - 1)]} is the run.
#'
#' @param mask Logical vector.
#' @return A list with integer vectors \code{onsets} and \code{offsets}, the
#'   same length (zero-length if \code{mask} has no \code{TRUE} runs).
#' @keywords internal
.mask_to_onsets_offsets <- function(mask) {

  n <- length(mask)
  if (n == 0) return(list(onsets = integer(0), offsets = integer(0)))

  m <- as.integer(mask)
  d <- diff(m)

  onsets <- which(d > 0) + 1L
  if (mask[1]) onsets <- c(1L, onsets)

  offsets <- which(d < 0) + 1L
  if (mask[n]) offsets <- c(offsets, n + 1L)

  list(onsets = onsets, offsets = offsets)
}

#' Convert a Boolean Mask into Annotation Rows (internal)
#'
#' Runs \code{\link{.mask_to_onsets_offsets}} on \code{mask} and converts
#' each run to seconds, assuming sample 1 falls at t = 0 (matching
#' \code{\link{new_eeg}}'s auto-generated \code{times}).
#'
#' @param mask Logical vector, one element per sample.
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param description Character scalar. Value for the \code{description}
#'   column of every row produced.
#' @param channel Character scalar. Value for the \code{channel} column of
#'   every row produced. Default \code{NA_character_} (not channel-specific).
#' @return A data frame with columns \code{onset, duration, description,
#'   channel} (see \code{\link{.empty_annotations}}), zero rows if
#'   \code{mask} has no \code{TRUE} runs.
#' @keywords internal
.mask_to_annotation_rows <- function(mask, sfreq, description,
                                      channel = NA_character_) {

  oo <- .mask_to_onsets_offsets(mask)
  if (length(oo$onsets) == 0) return(.empty_annotations())

  data.frame(
    onset       = (oo$onsets - 1L) / sfreq,
    duration    = (oo$offsets - oo$onsets) / sfreq,
    description = description,
    channel     = channel,
    stringsAsFactors = FALSE
  )
}

#' Flip Short Runs of a Given Value (internal)
#'
#' Run-length-encodes \code{mask} and flips any run equal to \code{run_value}
#' shorter than \code{min_samples} samples to \code{flip_to}. Used two ways
#' in this file: \code{\link{annotate_muscle}} uses it to fold short "good"
#' gaps into the surrounding flagged stretch (\code{run_value = FALSE,
#' flip_to = TRUE}); \code{\link{annotate_amplitude}} uses it to discard
#' flagged runs that don't last long enough to count (\code{run_value =
#' TRUE, flip_to = FALSE}).
#'
#' @param mask Logical vector.
#' @param min_samples Numeric. Runs of \code{run_value} shorter than this
#'   many samples are flipped.
#' @param run_value Logical. Which run value is subject to flipping.
#' @param flip_to Logical. What a too-short run of \code{run_value} becomes.
#' @return Logical vector, same length as \code{mask}.
#' @keywords internal
.flip_short_runs <- function(mask, min_samples, run_value, flip_to) {

  if (length(mask) == 0 || min_samples <= 0) return(mask)

  r <- rle(mask)
  hit <- (r$values == run_value) & (r$lengths < min_samples)
  if (any(hit)) r$values[hit] <- flip_to

  inverse.rle(r)
}


# ============================================================================
#                     annotate_amplitude()
# ============================================================================

#' Annotate Segments with Flat or Excessive Peak-to-Peak Amplitude
#'
#' Scans consecutive-sample differences per channel and flags stretches
#' where the signal barely moves (\code{flat}) or jumps too hard
#' (\code{peak}) for at least \code{min_duration} seconds. A channel
#' flagged for \code{bad_percent} of the recording or more is added to
#' \code{eeg_obj$bads} instead (mirrors \code{\link{find_bad_channels}}) -
#' the problem is the channel, not a moment in it. Everything below that
#' threshold is written to \code{eeg_obj$annotations} as \code{"BAD_flat"} /
#' \code{"BAD_peak"} rows, unioned across every channel that tripped it.
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param peak Numeric or \code{NULL}. Flag stretches where
#'   \code{abs(diff(x))} is at or above this value (same units as
#'   \code{eeg_obj$data}, typically microV). \code{NULL} (default) disables
#'   the peak check. At least one of \code{peak}/\code{flat} must be set.
#' @param flat Numeric or \code{NULL}. Flag stretches where
#'   \code{abs(diff(x))} is at or below this value. \code{NULL} (default)
#'   disables the flat check.
#' @param bad_percent Numeric in \verb{[0, 100]}. A channel flagged for this
#'   percentage of the recording (or more) is added to \code{eeg_obj$bads}
#'   instead of annotated. Default: 5.
#' @param min_duration Numeric. Minimum duration in seconds a stretch of
#'   consecutive flagged samples must last to count. Default: 0.005 (5 ms).
#' @param channels Character or integer vector, or \code{NULL}. Channels to
#'   check. \code{NULL} (default) uses every channel classified \code{"eeg"}
#'   that is not already in \code{eeg_obj$bads}.
#' @param return_details Logical. If \code{FALSE} (default), returns just
#'   the updated \code{eeg_obj}. If \code{TRUE}, returns
#'   \code{list(eeg_obj, annotations)}, where \code{annotations} is a data
#'   frame of just the rows this call added.
#'
#' @return If \code{return_details = FALSE} (default), the input
#'   \code{eeg_obj} with new rows appended to \code{$annotations}, any
#'   channel bad for most of the recording added to \code{$bads}, and a step
#'   appended to \code{$preprocessing_history}. If \code{return_details =
#'   TRUE}, a list with elements \code{eeg_obj} (as above) and
#'   \code{annotations} (data frame with columns \code{onset, duration,
#'   description, channel} - just the rows added by this call).
#'
#' @details
#' A channel's flagged-sample count is bumped by 1 before converting to a
#' percentage, correcting for \code{diff()} being one sample shorter than
#' the channel - a per-channel, not per-run, correction.
#'
#' This only detects edges where consecutive samples change abruptly; a slow
#' drift or a plateau reached gradually will not trigger \code{peak} or
#' \code{flat} the way the names alone might suggest.
#'
#' @examples
#' \dontrun{
#'   eeg <- annotate_amplitude(eeg, peak = 500, flat = 0.5)
#'   eeg$annotations
#'   eeg$bads
#' }
#'
#' @seealso \code{\link{find_bad_channels}}, \code{\link{annotate_muscle}},
#'   \code{\link{annotate_nan}}, \code{\link{annotate_break}}
#'
#' @export
annotate_amplitude <- function(eeg_obj,
                                peak           = NULL,
                                flat           = NULL,
                                bad_percent    = 5,
                                min_duration   = 0.005,
                                channels       = NULL,
                                return_details = FALSE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (is.null(peak) && is.null(flat)) {
    stop("ERROR: at least one of 'peak' or 'flat' must be supplied.",
         call. = FALSE)
  }
  if (!is.null(peak) && (!is.numeric(peak) || length(peak) != 1 || peak < 0)) {
    stop("ERROR: 'peak' must be a single non-negative numeric value.",
         call. = FALSE)
  }
  if (!is.null(flat) && (!is.numeric(flat) || length(flat) != 1 || flat < 0)) {
    stop("ERROR: 'flat' must be a single non-negative numeric value.",
         call. = FALSE)
  }
  if (!is.numeric(bad_percent) || length(bad_percent) != 1 ||
      bad_percent < 0 || bad_percent > 100) {
    stop("ERROR: 'bad_percent' must be a single numeric value in [0, 100].",
         call. = FALSE)
  }

  eeg_obj <- .ensure_annotations(eeg_obj)

  sfreq        <- eeg_obj$sampling_rate
  n_times      <- ncol(eeg_obj$data)
  rec_duration <- n_times / sfreq

  if (!is.numeric(min_duration) || length(min_duration) != 1 ||
      min_duration < 0 || min_duration >= rec_duration) {
    stop("ERROR: 'min_duration' must be a single non-negative numeric value ",
         "shorter than the recording (", rec_duration, " s).", call. = FALSE)
  }
  min_duration_samples <- round(min_duration * sfreq)

  # ========== RESOLVE CHANNELS ==========

  default_idx <- which(eeg_obj$channel_types == "eeg" &
                          !(eeg_obj$channels %in% eeg_obj$bads))
  ch_idx <- .resolve_channel_idx(eeg_obj, channels, default_idx)

  if (length(ch_idx) == 0) {
    stop("ERROR: No channels available to check - either none are ",
         "classified as EEG, or all of them are already in eeg_obj$bads.",
         call. = FALSE)
  }

  # ========== PER-CHANNEL SCAN ==========

  any_flat  <- rep(FALSE, n_times - 1L)
  any_peak  <- rep(FALSE, n_times - 1L)
  newly_bad <- character(0)

  for (i in ch_idx) {

    ch <- eeg_obj$channels[i]
    d  <- abs(diff(eeg_obj$data[i, ]))

    if (!is.null(flat)) {
      flat_mask <- .flip_short_runs(d <= flat, min_duration_samples,
                                     run_value = TRUE, flip_to = FALSE)
      cnt <- sum(flat_mask)
      if (cnt > 0) cnt <- cnt + 1L
      pct <- cnt / n_times * 100
      if (pct >= bad_percent) {
        newly_bad <- c(newly_bad, ch)
      } else if (pct > 0) {
        any_flat <- any_flat | flat_mask
      }
    }

    if (!is.null(peak)) {
      peak_mask <- .flip_short_runs(d >= peak, min_duration_samples,
                                     run_value = TRUE, flip_to = FALSE)
      cnt <- sum(peak_mask)
      if (cnt > 0) cnt <- cnt + 1L
      pct <- cnt / n_times * 100
      if (pct >= bad_percent) {
        newly_bad <- c(newly_bad, ch)
      } else if (pct > 0) {
        any_peak <- any_peak | peak_mask
      }
    }
  }

  # ========== BUILD ANNOTATIONS + WRITE bads ==========

  new_rows <- rbind(
    .mask_to_annotation_rows(any_flat, sfreq, "BAD_flat"),
    .mask_to_annotation_rows(any_peak, sfreq, "BAD_peak")
  )

  newly_bad    <- unique(newly_bad)
  eeg_obj$bads <- union(eeg_obj$bads, newly_bad)

  eeg_obj <- .append_annotations(eeg_obj, new_rows)

  history_entry <- paste0(
    "annotate_amplitude(): ", nrow(new_rows), " new annotation(s) added",
    if (length(newly_bad) > 0) {
      paste0("; ", length(newly_bad), " channel(s) newly marked bad (",
             paste(newly_bad, collapse = ", "), ")")
    } else {
      ""
    }, ".")

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  if (isTRUE(return_details)) {
    return(list(eeg_obj = eeg_obj, annotations = new_rows))
  }

  eeg_obj
}


# ============================================================================
#                     annotate_muscle()
# ============================================================================

#' Annotate Muscle (EMG) Artifacts via High-Frequency Envelope Z-Score
#'
#' Band-pass filters the picked channels to \code{filter_freq} (jaw-clench
#' and other EMG activity shows up here, well above EEG's own frequency
#' range), takes the Hilbert envelope of each, z-scores every channel's
#' envelope across time, sums the z-scores across channels (divided by
#' \code{sqrt(n_channels)}), and low-pass filters that combined trace at
#' 4 Hz to suppress spurious transient peaks. Samples where the smoothed
#' trace exceeds \code{threshold} are flagged; "good" gaps between flagged
#' stretches shorter than \code{min_length_good} are absorbed into the
#' surrounding flagged stretch rather than left as a sliver of good data.
#' Flagged stretches are written to \code{eeg_obj$annotations} as
#' \code{"BAD_muscle"} rows.
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param threshold Numeric. Z-score threshold above which a sample counts
#'   as muscle activity. Default: 4.
#' @param filter_freq Numeric vector of length 2, \verb{c(low, high)} in Hz.
#'   Band tested for EMG envelope activity. Default: \code{c(110, 140)} -
#'   requires \code{eeg_obj$sampling_rate} well above
#'   \code{2 * filter_freq[2]} to leave room for the filter's transition
#'   band; lower sampling rates error out (see Details).
#' @param min_length_good Numeric. Shortest allowed run of "good" data (in
#'   seconds) between two flagged stretches; shorter runs are folded into
#'   the surrounding flagged stretch. Default: 0.1.
#' @param channels Character or integer vector, or \code{NULL}. Channels to
#'   use. \code{NULL} (default) uses every channel classified \code{"eeg"}
#'   that is not already in \code{eeg_obj$bads}.
#' @param return_details Logical. If \code{FALSE} (default), returns just
#'   the updated \code{eeg_obj}. If \code{TRUE}, returns
#'   \code{list(eeg_obj, annotations, scores)}, where \code{annotations} is
#'   a data frame of just the rows this call added and \code{scores} is the
#'   smoothed combined z-score trace (length \code{ncol(eeg_obj$data)}),
#'   useful for plotting.
#'
#' @return If \code{return_details = FALSE} (default), the input
#'   \code{eeg_obj} with new \code{"BAD_muscle"} rows appended to
#'   \code{$annotations} and a step appended to
#'   \code{$preprocessing_history}. If \code{return_details = TRUE}, a list
#'   with elements \code{eeg_obj} (as above), \code{annotations} (just the
#'   rows added by this call), and \code{scores} (the z-score trace).
#'
#' @details
#' This function bandpass- and lowpass-filters data via the same FFT-based
#' overlap-add convolution \code{\link{eeg_bandpass}} uses (see
#' R/filter1.R), which has no concept of skipping over already-bad
#' stretches of a continuous signal - it filters straight through. Because
#' both steps are FFT-based (global support), a single \code{NA} anywhere
#' in the picked channels would silently turn the entire filtered trace to
#' \code{NA} rather than just the affected stretch - so this function
#' errors instead if it finds one. Run
#' \code{\link{annotate_nan}} first to locate such gaps; there is no
#' automatic way to filter around them yet.
#'
#' @examples
#' \dontrun{
#'   eeg <- annotate_muscle(eeg)
#'   eeg$annotations
#'
#'   result <- annotate_muscle(eeg, return_details = TRUE)
#'   plot(eeg$times, result$scores, type = "l")
#'   abline(h = 4, col = "red", lty = 2)
#' }
#'
#' @seealso \code{\link{annotate_amplitude}}, \code{\link{annotate_nan}},
#'   \code{\link{annotate_break}}
#'
#' @export
annotate_muscle <- function(eeg_obj,
                             threshold       = 4,
                             filter_freq     = c(110, 140),
                             min_length_good = 0.1,
                             channels        = NULL,
                             return_details  = FALSE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("ERROR: 'threshold' must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(filter_freq) || length(filter_freq) != 2 ||
      filter_freq[1] <= 0 || filter_freq[2] <= filter_freq[1]) {
    stop("ERROR: 'filter_freq' must be a numeric vector c(low, high) with ",
         "0 < low < high.", call. = FALSE)
  }
  if (!is.numeric(min_length_good) || length(min_length_good) != 1 ||
      min_length_good < 0) {
    stop("ERROR: 'min_length_good' must be a single non-negative numeric ",
         "value.", call. = FALSE)
  }

  eeg_obj <- .ensure_annotations(eeg_obj)

  sfreq   <- eeg_obj$sampling_rate
  nyquist <- sfreq / 2

  if (filter_freq[2] >= nyquist) {
    stop("ERROR: 'filter_freq' upper bound (", filter_freq[2],
         " Hz) must be below the Nyquist frequency (", nyquist,
         " Hz) - eeg_obj$sampling_rate (", sfreq, " Hz) is too low for ",
         "this filter_freq.", call. = FALSE)
  }

  # ========== RESOLVE CHANNELS ==========

  default_idx <- which(eeg_obj$channel_types == "eeg" &
                          !(eeg_obj$channels %in% eeg_obj$bads))
  ch_idx <- .resolve_channel_idx(eeg_obj, channels, default_idx)

  if (length(ch_idx) == 0) {
    stop("ERROR: No channels available to check - either none are ",
         "classified as EEG, or all of them are already in eeg_obj$bads.",
         call. = FALSE)
  }
  if (anyNA(eeg_obj$data[ch_idx, , drop = FALSE])) {
    stop("ERROR: picked channel(s) contain NA values - FFT-based filtering ",
         "cannot run past them (see Details in ?annotate_muscle). Use ",
         "annotate_nan() to locate the gap(s) first.", call. = FALSE)
  }

  # ========== BANDPASS + HILBERT ENVELOPE PER CHANNEL ==========

  n_picked <- length(ch_idx)
  n_times  <- ncol(eeg_obj$data)
  envelope <- matrix(0, nrow = n_picked, ncol = n_times)

  for (i in seq_len(n_picked)) {
    x_filt <- .fir_filter_vector(eeg_obj$data[ch_idx[i], ], sfreq,
                                  l_freq = filter_freq[1],
                                  h_freq = filter_freq[2])
    envelope[i, ] <- .hilbert_envelope(x_filt)
  }

  # ========== Z-SCORE, COMBINE, SMOOTH ==========

  z        <- .zscore_rows(envelope)
  combined <- colSums(z) / sqrt(n_picked)
  smoothed <- .fir_filter_vector(combined, sfreq, h_freq = 4)

  art_mask <- smoothed > threshold

  min_samps <- min_length_good * sfreq
  art_mask  <- .flip_short_runs(art_mask, min_samps,
                                 run_value = FALSE, flip_to = TRUE)

  # ========== BUILD ANNOTATIONS ==========

  new_rows <- .mask_to_annotation_rows(art_mask, sfreq, "BAD_muscle")
  eeg_obj  <- .append_annotations(eeg_obj, new_rows)

  history_entry <- paste0(
    "annotate_muscle(): ", nrow(new_rows), " new BAD_muscle annotation(s) ",
    "added, covering ", round(sum(new_rows$duration), 2), " s (",
    n_picked, " channel(s), filter_freq = ", filter_freq[1], "-",
    filter_freq[2], " Hz).")

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  if (isTRUE(return_details)) {
    return(list(eeg_obj = eeg_obj, annotations = new_rows, scores = smoothed))
  }

  eeg_obj
}

#' Design-and-Apply a Bandpass or Lowpass FIR Filter to One Vector (internal)
#'
#' A single-vector counterpart to \code{\link{eeg_bandpass}} (see
#' R/filter1.R): builds the same freq/gain kernel specification
#' \code{eeg_bandpass()} does (bandpass when both \code{l_freq}/\code{h_freq}
#' are given, lowpass when only \code{h_freq} is given), then applies it with
#' the same \code{.firwin_kernel()} / \code{.overlap_add_filter()} primitives
#' - without constructing a full \code{eeg} object, which is all
#' \code{\link{annotate_muscle}} needs for its two internal filtering steps
#' (110-140 Hz bandpass on each channel, 4 Hz lowpass on the combined
#' z-score trace). Highpass-only is intentionally not implemented - no
#' caller in this file needs it.
#'
#' @param x Numeric vector. The signal to filter.
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param l_freq Numeric or \code{NULL}. Lower passband edge in Hz.
#' @param h_freq Numeric or \code{NULL}. Upper passband edge in Hz.
#' @param fir_window,phase,pad Passed through to \code{.firwin_kernel()} /
#'   \code{.overlap_add_filter()}. Defaults match \code{eeg_bandpass()}.
#'
#' @return Numeric vector, same length as \code{x}.
#' @keywords internal
.fir_filter_vector <- function(x, sfreq, l_freq = NULL, h_freq = NULL,
                                fir_window = "hamming", phase = "zero",
                                pad = "reflect_limited") {

  nyquist <- sfreq / 2.0

  if (is.null(l_freq) && is.null(h_freq)) {
    stop(".fir_filter_vector(): at least one of 'l_freq' or 'h_freq' must ",
         "be supplied.", call. = FALSE)
  }
  if (!is.null(l_freq) && is.null(h_freq)) {
    stop(".fir_filter_vector(): highpass-only is not implemented - this ",
         "file only needs bandpass (both freqs) or lowpass (h_freq only).",
         call. = FALSE)
  }

  l_trans_bw <- if (is.null(l_freq)) NULL else min(max(0.25 * l_freq, 2.0), l_freq)
  h_trans_bw <- min(max(0.25 * h_freq, 2.0), nyquist - h_freq)

  N <- .auto_filter_length(l_trans_bw, h_trans_bw, sfreq, window = fir_window)

  if (!is.null(l_freq)) {
    # Bandpass
    f_s1 <- (l_freq - l_trans_bw) / nyquist
    f_p1 <- l_freq / nyquist
    f_p2 <- h_freq / nyquist
    f_s2 <- (h_freq + h_trans_bw) / nyquist
    freq <- c(f_s1, f_p1, f_p2, f_s2)
    gain <- c(   0,    1,    1,    0)
    if (f_s1 != 0) { freq <- c(0, freq); gain <- c(0, gain) }
    if (f_s2 != 1) { freq <- c(freq, 1); gain <- c(gain, 0) }
  } else {
    # Lowpass only
    f_p  <- h_freq / nyquist
    f_s  <- (h_freq + h_trans_bw) / nyquist
    freq <- c(0, f_p, f_s)
    gain <- c(1,   1,   0)
    if (f_s != 1) { freq <- c(freq, 1); gain <- c(gain, 0) }
  }

  freq <- pmax(pmin(freq, 1.0), 0.0)
  h <- .firwin_kernel(N, freq, gain, window = fir_window)
  .overlap_add_filter(x, h, phase = phase, pad = pad)
}

#' Hilbert Envelope via FFT (internal)
#'
#' Computes the analytic signal of \code{x} via the standard FFT recipe
#' (double the positive frequencies, zero the negative ones, keep DC and
#' Nyquist unscaled) and returns its magnitude, i.e. the envelope.
#'
#' @param x Numeric vector.
#' @return Numeric vector, same length as \code{x} - the envelope.
#' @keywords internal
.hilbert_envelope <- function(x) {

  n  <- length(x)
  Xf <- fft(x)
  h  <- numeric(n)

  if (n %% 2 == 0) {
    h[1] <- 1
    h[n / 2 + 1] <- 1
    if (n > 2) h[2:(n / 2)] <- 2
  } else {
    h[1] <- 1
    if (n > 1) h[2:((n + 1) / 2)] <- 2
  }

  analytic <- fft(Xf * h, inverse = TRUE) / n
  Mod(analytic)
}

#' Per-Row Z-Score, Population SD (internal)
#'
#' Z-scores each row of \code{x} independently using the population
#' standard deviation (divide by \code{n}, not \code{n - 1}) -
#' \code{\link{annotate_muscle}} relies on this exact convention.
#'
#' @param x Numeric matrix, rows to z-score independently (channels x time).
#' @return Numeric matrix, same shape as \code{x}.
#' @keywords internal
.zscore_rows <- function(x) {
  mu       <- rowMeans(x)
  centered <- sweep(x, 1, mu, FUN = "-")
  popsd    <- sqrt(rowMeans(centered^2))
  sweep(centered, 1, popsd, FUN = "/")
}


# ============================================================================
#                     annotate_nan()
# ============================================================================

#' Annotate Segments with NA (Amplifier Dropouts)
#'
#' Scans every picked channel independently for runs of \code{NA} and
#' writes each run to \code{eeg_obj$annotations} as a channel-specific
#' \code{"BAD_NAN"} row - this is the one detector in this file that fills
#' in the \code{channel} column, since a dropout is a fact about one
#' amplifier line, not the whole recording, and shouldn't block channels
#' that stayed connected.
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param channels Character or integer vector, or \code{NULL}. Channels to
#'   scan. \code{NULL} (default) scans every channel in
#'   \code{eeg_obj$channels}, regardless of type or \code{bads} status - a
#'   dropout is a hardware fact, independent of how a channel was later
#'   classified.
#' @param return_details Logical. If \code{FALSE} (default), returns just
#'   the updated \code{eeg_obj}. If \code{TRUE}, returns
#'   \code{list(eeg_obj, annotations)}, where \code{annotations} is a data
#'   frame of just the rows this call added.
#'
#' @return If \code{return_details = FALSE} (default), the input
#'   \code{eeg_obj} with new \code{"BAD_NAN"} rows appended to
#'   \code{$annotations} and a step appended to
#'   \code{$preprocessing_history}. If \code{return_details = TRUE}, a list
#'   with elements \code{eeg_obj} (as above) and \code{annotations} (just
#'   the rows added by this call).
#'
#' @examples
#' \dontrun{
#'   eeg <- annotate_nan(eeg)
#'   eeg$annotations
#' }
#'
#' @seealso \code{\link{annotate_amplitude}}, \code{\link{annotate_muscle}},
#'   \code{\link{annotate_break}}
#'
#' @export
annotate_nan <- function(eeg_obj,
                          channels       = NULL,
                          return_details = FALSE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }

  eeg_obj <- .ensure_annotations(eeg_obj)
  sfreq   <- eeg_obj$sampling_rate

  # ========== RESOLVE CHANNELS ==========

  ch_idx <- .resolve_channel_idx(eeg_obj, channels,
                                  default_idx = seq_along(eeg_obj$channels))

  # ========== SCAN EACH CHANNEL ==========

  new_rows <- .empty_annotations()

  for (i in ch_idx) {
    ch   <- eeg_obj$channels[i]
    rows <- .mask_to_annotation_rows(is.na(eeg_obj$data[i, ]), sfreq,
                                      "BAD_NAN", channel = ch)
    new_rows <- rbind(new_rows, rows)
  }

  eeg_obj <- .append_annotations(eeg_obj, new_rows)

  history_entry <- paste0(
    "annotate_nan(): ", nrow(new_rows), " new BAD_NAN annotation(s) added",
    if (nrow(new_rows) > 0) {
      paste0(" across ", length(unique(new_rows$channel)), " channel(s).")
    } else {
      "."
    })

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  if (isTRUE(return_details)) {
    return(list(eeg_obj = eeg_obj, annotations = new_rows))
  }

  eeg_obj
}


# ============================================================================
#                     annotate_break()
# ============================================================================

#' Annotate Dead Time Between Experimental Blocks
#'
#' Looks for gaps at least \code{min_break_duration} seconds long between
#' known-occupied stretches of the recording, and writes each gap to
#' \code{eeg_obj$annotations} as a \code{"BAD_break"} row - trimmed on both
#' ends so the annotation doesn't start or stop right on top of a real
#' event. Unlike the other three detectors here, this one never looks at
#' \code{eeg_obj$data} - only at time already accounted for.
#'
#' "Occupied" stretches come from one of two places:
#' \itemize{
#'   \item \code{use_events = FALSE} (default): every existing row in
#'     \code{eeg_obj$annotations} whose \code{description} does not start
#'     with one of \code{ignore} (case-insensitive) - typically whatever
#'     \code{\link{annotate_amplitude}}/\code{\link{annotate_muscle}}/
#'     \code{\link{annotate_nan}} already wrote, though by default those all
#'     start with \code{"BAD_"} and so are themselves excluded by
#'     \code{ignore = "bad"}. Pass \code{ignore = character(0)} to count
#'     every annotation instead.
#'   \item \code{use_events = TRUE}: every trigger in \code{eeg_obj$events}
#'     (see R/extract_bdf_events.R), each treated as an instantaneous
#'     marker rather than a stretch. Natural for finding dead time between
#'     experimental blocks delimited by real trigger codes.
#' }
#' Overlapping occupied stretches are merged before gaps are computed, so
#' this is safe to call with a busy, overlapping \code{annotations} table.
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param use_events Logical. See above. Default: \code{FALSE}.
#' @param min_break_duration Numeric. Minimum gap, in seconds, between two
#'   occupied stretches (or events) to call it a break. This is the minimum
#'   size of the gap itself, not of the resulting annotation - see
#'   \code{t_start_after_previous}/\code{t_stop_before_next}. Default: 15.
#' @param t_start_after_previous,t_stop_before_next Numeric. How far the
#'   break annotation is trimmed back from the stretches bounding it on each
#'   side, so it doesn't start or end right on top of a real event. Both
#'   default to 5. \code{min_break_duration - t_start_after_previous -
#'   t_stop_before_next} must be positive (otherwise every gap found would
#'   produce a non-positive-duration annotation).
#' @param ignore Character vector of description prefixes to exclude,
#'   matched case-insensitively, when \code{use_events = FALSE}. Default:
#'   \code{"bad"}. This package has no concept of "edge" annotations from
#'   concatenating separate recordings together, so there is nothing else
#'   ignored by default. Has no effect when \code{use_events = TRUE}. Pass
#'   \code{character(0)} to keep everything.
#' @param return_details Logical. If \code{FALSE} (default), returns just
#'   the updated \code{eeg_obj}. If \code{TRUE}, returns
#'   \code{list(eeg_obj, annotations)}, where \code{annotations} is a data
#'   frame of just the rows this call added.
#'
#' @return If \code{return_details = FALSE} (default), the input
#'   \code{eeg_obj} with new \code{"BAD_break"} rows appended to
#'   \code{$annotations} and a step appended to
#'   \code{$preprocessing_history}. If \code{return_details = TRUE}, a list
#'   with elements \code{eeg_obj} (as above) and \code{annotations} (just
#'   the rows added by this call).
#'
#' @examples
#' \dontrun{
#'   eeg <- annotate_amplitude(eeg, peak = 500)
#'   eeg <- annotate_break(eeg)
#'
#'   # Or, based on real trigger events instead:
#'   eeg <- annotate_break(eeg, use_events = TRUE)
#' }
#'
#' @seealso \code{\link{annotate_amplitude}}, \code{\link{annotate_muscle}},
#'   \code{\link{annotate_nan}}
#'
#' @export
annotate_break <- function(eeg_obj,
                            use_events             = FALSE,
                            min_break_duration     = 15,
                            t_start_after_previous = 5,
                            t_stop_before_next     = 5,
                            ignore                 = "bad",
                            return_details         = FALSE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }

  annot_dur <- min_break_duration - t_start_after_previous - t_stop_before_next
  if (annot_dur <= 0) {
    stop("ERROR: min_break_duration - t_start_after_previous - ",
         "t_stop_before_next must be positive, got ", annot_dur, ".",
         call. = FALSE)
  }

  eeg_obj <- .ensure_annotations(eeg_obj)

  # ========== BUILD OCCUPIED INTERVALS ==========

  if (isTRUE(use_events)) {

    if (is.null(eeg_obj$events) || nrow(eeg_obj$events) == 0) {
      stop("ERROR: annotate_break(use_events = TRUE) needs a non-empty ",
           "eeg_obj$events.", call. = FALSE)
    }
    starts    <- sort(eeg_obj$events$onset_time)
    intervals <- cbind(starts, starts, deparse.level = 0)

  } else {

    ann <- eeg_obj$annotations
    if (length(ignore) > 0 && nrow(ann) > 0) {
      ignore_lc <- tolower(ignore)
      descr_lc  <- tolower(ann$description)
      excluded  <- vapply(descr_lc, function(d) any(startsWith(d, ignore_lc)),
                          logical(1))
      ann <- ann[!excluded, , drop = FALSE]
    }
    if (nrow(ann) == 0) {
      stop("ERROR: annotate_break() found no usable rows in ",
           "eeg_obj$annotations (empty, or every row excluded by 'ignore'). ",
           "Run one of annotate_amplitude()/annotate_muscle()/annotate_nan() ",
           "first, or pass use_events = TRUE to use eeg_obj$events instead.",
           call. = FALSE)
    }
    ord       <- order(ann$onset)
    starts    <- ann$onset[ord]
    intervals <- cbind(starts, starts + ann$duration[ord], deparse.level = 0)
  }

  # ========== MERGE OVERLAPPING INTERVALS ==========

  merged <- .merge_intervals(intervals)

  # ========== FIND GAPS >= min_break_duration ==========

  rec_end        <- eeg_obj$times[length(eeg_obj$times)]
  break_onset    <- numeric(0)
  break_duration <- numeric(0)

  if (merged[1, 1] > 0 && merged[1, 1] >= min_break_duration) {
    onset  <- 0
    offset <- merged[1, 1] - t_stop_before_next
    break_onset    <- c(break_onset, onset)
    break_duration <- c(break_duration, offset - onset)
  }

  if (nrow(merged) > 1) {
    for (i in 2:nrow(merged)) {
      this_start <- merged[i, 1]
      prev_stop  <- merged[i - 1, 2]
      if (this_start - prev_stop < min_break_duration) next
      onset  <- prev_stop + t_start_after_previous
      offset <- this_start - t_stop_before_next
      break_onset    <- c(break_onset, onset)
      break_duration <- c(break_duration, offset - onset)
    }
  }

  last_stop <- merged[nrow(merged), 2]
  if (rec_end > last_stop && (rec_end - last_stop) >= min_break_duration) {
    onset  <- last_stop + t_start_after_previous
    offset <- rec_end
    break_onset    <- c(break_onset, onset)
    break_duration <- c(break_duration, offset - onset)
  }

  # ========== BUILD ANNOTATIONS ==========

  new_rows <- if (length(break_onset) == 0) {
    .empty_annotations()
  } else {
    data.frame(onset = break_onset, duration = break_duration,
               description = "BAD_break", channel = NA_character_,
               stringsAsFactors = FALSE)
  }

  eeg_obj <- .append_annotations(eeg_obj, new_rows)

  history_entry <- paste0(
    "annotate_break(): ", nrow(new_rows), " new BAD_break annotation(s) ",
    "added, covering ", round(sum(new_rows$duration), 2), " s (source: ",
    if (isTRUE(use_events)) "events" else "annotations", ").")

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  if (isTRUE(return_details)) {
    return(list(eeg_obj = eeg_obj, annotations = new_rows))
  }

  eeg_obj
}

#' Merge Overlapping (or Touching) Intervals (internal)
#'
#' Sweeps \code{intervals} (already sorted by start) left to right, merging
#' any interval that overlaps or touches the previous merged interval.
#'
#' @param intervals A 2-column numeric matrix, one row per interval,
#'   \code{[, 1]} = start, \code{[, 2]} = stop, sorted by start ascending.
#' @return A 2-column numeric matrix of merged, non-overlapping intervals.
#' @keywords internal
.merge_intervals <- function(intervals) {

  merged <- intervals[1, , drop = FALSE]

  for (i in seq_len(nrow(intervals))) {
    last_stop <- merged[nrow(merged), 2]
    start     <- intervals[i, 1]
    stop_at   <- intervals[i, 2]

    if (stop_at < last_stop) {
      next
    } else if (start <= last_stop) {
      merged[nrow(merged), 2] <- stop_at
    } else {
      merged <- rbind(merged, c(start, stop_at))
    }
  }

  merged
}
