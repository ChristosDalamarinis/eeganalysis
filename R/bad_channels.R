#' ============================================================================
#'                       Bad EEG Channel Detection
#' ============================================================================
#'
#' This module detects bad EEG channels automatically and writes them to
#' eeg_obj$bads (see R/eeg_class.R), instead of only reporting them the way
#' eeg_summary() does. find_bad_channels() runs four independent checks:
#'   - Flat / disconnected, excessive amplitude, outlier noise level (the
#'     same three eeg_summary() already reports, print-only)
#'   - Low correlation with physically-nearby channels (uses the montage
#'     attached via set_montage(); skipped if none is attached)
#'   - Global Local Outlier Factor, mirroring MNE-Python's
#'     find_bad_channels_lof() (no montage needed - compares raw channel
#'     signals directly, with no scalp position involved at all)
#'
#' Flagged channels are added to eeg_obj$bads additively - existing entries
#' are never removed. The checks run independently against the same
#' starting eeg_obj$bads snapshot and are combined only at the end, so the
#' order they run in has no effect on the result.
#'
#' Author: Christos Dalamarinis
#' Date: Sep - 2026
#' Status: Ready
#' ============================================================================
#'
#' Automatically Detect and Flag Bad EEG Channels
#'
#' Runs four independent checks against an \code{eeg} object - flat /
#' disconnected, excessive amplitude, outlier noise level, low correlation
#' with physically nearby channels, and a global Local Outlier Factor (LOF)
#' check - and adds every channel any check flags to \code{eeg_obj$bads}.
#' Unlike \code{\link{eeg_summary}}, which only prints its flags, this
#' function writes them.
#'
#' @param eeg_obj An object of class 'eeg'.
#'
#' @param flag_flat_threshold Numeric. Channels with std below this value
#'   (in uV) are flagged as potentially flat / disconnected. Default: 0.5.
#'
#' @param flag_amplitude_threshold Numeric. Channels whose absolute
#'   amplitude, measured from the channel's own level (see
#'   \code{amplitude_window}), ever exceeds this value (in uV) are flagged as
#'   potentially noisy or saturated. Default: 500.
#'
#' @param amplitude_window Numeric. Length in seconds of the stretches in
#'   which each channel's own level (its average) is measured and subtracted
#'   before the amplitude check looks at it. Raw BioSemi data keeps a
#'   constant offset on every electrode (often tens of mV) that slowly
#'   drifts, and without this step every channel would look far above
#'   \code{flag_amplitude_threshold}. The level is measured again in every
#'   stretch because it creeps; a last stretch shorter than the window joins
#'   the stretch before it. \code{Inf} measures one level for the whole
#'   recording. Must be at least 2 samples long. Default: 2.
#'
#' @param flag_outlier_sd_multiplier Numeric. A channel whose std exceeds
#'   this multiple of the median std across all EEG channels is flagged
#'   as an outlier. Default: 3.
#'
#' @param k Integer. Number of physically-nearest neighbor channels (by
#'   Euclidean distance on the attached montage's x/y/z positions) used to
#'   build each channel's local reference signal for the neighbor-
#'   correlation check. Only used when \code{eeg_obj$montage} is attached.
#'   Default: 4.
#'
#' @param correlation_threshold Numeric. A channel whose correlation with
#'   its neighbor reference signal (see \code{k}) falls below this value is
#'   flagged. Default: 0.45 - a starting point inspired by the PREP
#'   pipeline's channel-correlation criterion, not a fixed rule; adjust per
#'   dataset.
#'
#' @param robust Logical. If \code{TRUE} (default), the neighbor reference
#'   signal is the per-timepoint median of a channel's \code{k} neighbors
#'   (robust to one bad neighbor). If \code{FALSE}, the per-timepoint mean.
#'
#' @param lof_n_neighbors Integer. Number of comparison channels used by the
#'   global Local Outlier Factor check, mirroring MNE-Python's
#'   \code{find_bad_channels_lof(n_neighbors = ...)}. Default: 20
#'   (automatically capped if fewer EEG channels are available).
#'
#' @param lof_threshold Numeric. A channel whose LOF score is at or above
#'   this value is flagged, mirroring MNE-Python's
#'   \code{find_bad_channels_lof(threshold = ...)}. Default: 1.5.
#'
#' @param return_details Logical. If \code{FALSE} (default), returns just
#'   the updated \code{eeg} object. If \code{TRUE}, returns
#'   \code{list(eeg_obj, flags)}, where \code{flags} is a data frame of
#'   every channel/reason/value combination flagged this call.
#'
#' @return If \code{return_details = FALSE} (default), the input
#'   \code{eeg_obj} with newly-detected channels added to \code{$bads}
#'   (existing entries are never removed) and a step appended to
#'   \code{$preprocessing_history}. If \code{return_details = TRUE}, a list
#'   with elements \code{eeg_obj} (as above) and \code{flags} (data frame
#'   with columns \code{channel, reason, value}).
#'
#' @details
#' The flat/excessive-amplitude/outlier-noise checks always run and use the
#' same default thresholds \code{\link{eeg_summary}} already reports - this
#' function writes them to \code{bads} instead of only printing them. The
#' global LOF check also always runs.
#'
#' The amplitude check measures every channel from its own level, so it also
#' works on raw data with large DC offsets (see \code{amplitude_window}). The
#' neighbor-correlation and LOF checks use the data as it is: large offsets
#' or a strong slow drift can disturb them, so for raw data a high-pass
#' filter (for example 1 Hz, see \code{\link{eeg_bandpass}}) before this
#' function gives the most reliable result.
#'
#' The neighbor-correlation check only runs if \code{eeg_obj$montage} is
#' attached (see \code{\link{set_montage}}); if not, it is skipped with a
#' message and the other three checks still run.
#'
#' The global LOF check mirrors MNE-Python's
#' \code{mne.preprocessing.find_bad_channels_lof()}: every eligible
#' channel's raw signal is compared directly against every other's
#' (Euclidean distance, via \code{\link[dbscan]{lof}}), with no scalp
#' position involved at all - deliberately independent of the neighbor-
#' correlation check above, which is the opposite: purely spatial, no
#' data-similarity comparison across the whole cap. See Breunig et al.
#' (2000), "LOF: Identifying Density-Based Local Outliers".
#'
#' Bridged electrodes (two electrodes shorted together by gel) are not
#' caught by any of these checks: they record the same signal twice, so
#' they look unusually alike rather than bad. Use
#' \code{\link{find_bridged_electrodes}} for that.
#'
#' @examples
#' \dontrun{
#'   eeg <- find_bad_channels(eeg)
#'   eeg$bads
#'
#'   # With a montage attached, the neighbor-correlation check also runs
#'   eeg <- set_montage(eeg, create_montage())
#'   eeg <- find_bad_channels(eeg)
#'
#'   # Inspect what was flagged and why
#'   result <- find_bad_channels(eeg, return_details = TRUE)
#'   result$flags
#' }
#'
#' @seealso \code{\link{create_montage}}, \code{\link{set_montage}},
#'   \code{\link{eeg_summary}}, \code{\link{eeg_rereference}},
#'   \code{\link{find_bridged_electrodes}}
#'
#' @export
find_bad_channels <- function(eeg_obj,
                               flag_flat_threshold         = 0.5,
                               flag_amplitude_threshold    = 500,
                               amplitude_window            = 2,
                               flag_outlier_sd_multiplier  = 3,
                               k                           = 4,
                               correlation_threshold       = 0.45,
                               robust                      = TRUE,
                               lof_n_neighbors             = 20,
                               lof_threshold               = 1.5,
                               return_details              = FALSE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!is.matrix(eeg_obj$data)) {
    stop("ERROR: eeg_obj$data must be a numeric matrix (channels x timepoints).",
         call. = FALSE)
  }

  eeg_idx <- which(eeg_obj$channel_types == "eeg" &
                      !(eeg_obj$channels %in% eeg_obj$bads))

  if (length(eeg_idx) == 0) {
    stop("ERROR: No EEG channels available to check - either none are ",
         "classified as EEG, or all of them are already in eeg_obj$bads.",
         call. = FALSE)
  }

  eeg_channels <- eeg_obj$channels[eeg_idx]

  # ========== CHECK 1: FLAT / AMPLITUDE / OUTLIER ==========

  quality_flags <- flag_channel_quality(
    eeg_obj,
    flag_flat_threshold        = flag_flat_threshold,
    flag_amplitude_threshold   = flag_amplitude_threshold,
    amplitude_window           = amplitude_window,
    flag_outlier_sd_multiplier = flag_outlier_sd_multiplier)

  # ========== CHECK 2+3: NEIGHBOR CORRELATION (MONTAGE-GATED) ==========

  has_montage <- !is.null(eeg_obj$montage) && inherits(eeg_obj$montage, "montage")

  if (!has_montage) {

    message("find_bad_channels(): no montage attached - skipping the ",
            "neighbor-correlation check. Attach one with set_montage() to ",
            "include it.")
    corr_flags <- .empty_flags()

  } else {

    montage_eligible <- intersect(eeg_channels, eeg_obj$montage$channels)
    missing_montage   <- setdiff(eeg_channels, eeg_obj$montage$channels)

    if (length(missing_montage) > 0) {
      warning("find_bad_channels(): ", length(missing_montage),
              " EEG channel(s) have no position in the attached montage ",
              "and are excluded from the neighbor-correlation check: ",
              paste(missing_montage, collapse = ", "),
              call. = FALSE, immediate. = TRUE)
    }

    if (length(montage_eligible) < 2) {

      warning("find_bad_channels(): fewer than 2 EEG channels have a ",
              "montage position - skipping the neighbor-correlation check.",
              call. = FALSE, immediate. = TRUE)
      corr_flags <- .empty_flags()

    } else {

      positions_sub <- eeg_obj$montage$positions[
        eeg_obj$montage$positions$channel %in% montage_eligible, , drop = FALSE]
      neighbor_map <- find_channel_neighbors(positions_sub, k = k)

      data_idx <- match(montage_eligible, eeg_obj$channels)
      corr_flags <- flag_low_neighbor_correlation(
        data                  = eeg_obj$data[data_idx, , drop = FALSE],
        channels              = montage_eligible,
        neighbor_map          = neighbor_map,
        correlation_threshold = correlation_threshold,
        robust                = robust)
    }
  }

  # ========== CHECK 4: GLOBAL LOF ==========

  if (length(eeg_idx) < 2) {
    warning("find_bad_channels(): fewer than 2 eligible EEG channels - ",
            "skipping the Local Outlier Factor check.",
            call. = FALSE, immediate. = TRUE)
    lof_flags <- .empty_flags()
  } else {
    lof_flags <- flag_lof_outlier_channels(
      data        = eeg_obj$data[eeg_idx, , drop = FALSE],
      channels    = eeg_channels,
      n_neighbors = lof_n_neighbors,
      threshold   = lof_threshold)
  }

  # ========== COMBINE AND WRITE TO bads ==========

  all_flags     <- rbind(quality_flags, corr_flags, lof_flags)
  newly_flagged <- unique(all_flags$channel)

  eeg_obj$bads <- union(eeg_obj$bads, newly_flagged)

  history_entry <- paste0(
    "find_bad_channels(): ", length(newly_flagged), " channel(s) newly flagged",
    if (length(newly_flagged) > 0) {
      paste0(" (", paste(newly_flagged, collapse = ", "), ")")
    } else {
      ""
    },
    " - ", length(eeg_obj$bads), " total in bads.")

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  if (isTRUE(return_details)) {
    return(list(eeg_obj = eeg_obj, flags = all_flags))
  }

  eeg_obj
}

#' Flag Flat, Excessive-Amplitude, and Outlier-Noise Channels
#'
#' Internal helper implementing the same three per-channel quality checks
#' \code{\link{eeg_summary}} reports (flat/disconnected, excessive
#' amplitude, outlier noise level), but returning the flags instead of
#' printing them, so \code{\link{find_bad_channels}} can write them to
#' \code{eeg_obj$bads}. Designed to also be called from \code{eeg_summary()}
#' itself in place of its own copy of this logic (same reuse rationale as
#' \code{classify_channels()} in \code{R/eeg_class.R}).
#'
#' @param eeg_obj An object of class 'eeg'.
#' @param flag_flat_threshold Numeric, std (uV) below which a channel is
#'   flagged as flat. Default: 0.5.
#' @param flag_amplitude_threshold Numeric, peak absolute amplitude (uV)
#'   above which a channel is flagged. The amplitude is measured from the
#'   channel's own level (see \code{amplitude_window}), not from zero.
#'   Default: 500.
#' @param amplitude_window Numeric, length in seconds of the stretches in
#'   which each channel's own level (its average) is measured and subtracted
#'   before the amplitude check. Needed for raw BioSemi data, whose
#'   electrodes sit on constant offsets of tens of mV that also drift
#'   slowly. The last stretch, if shorter than the window, joins the one
#'   before it; \code{Inf} uses one level for the whole recording. Must be a
#'   single number above 0 and at least 2 samples long. Default: 2.
#' @param flag_outlier_sd_multiplier Numeric, multiple of the median EEG-
#'   channel std above which a channel is flagged as an outlier. Default: 3.
#'
#' @return A data frame with columns \code{channel, reason, value}
#'   (character), one row per flag raised (a channel may collect more than
#'   one). Zero rows (never \code{NULL}) if nothing is flagged.
#'
#' @keywords internal
flag_channel_quality <- function(eeg_obj,
                                  flag_flat_threshold        = 0.5,
                                  flag_amplitude_threshold   = 500,
                                  amplitude_window           = 2,
                                  flag_outlier_sd_multiplier = 3) {

  # The window for measuring each channel's own level (checked first, so a
  # wrong value is reported even when no channel is eligible)
  if (!is.numeric(amplitude_window) || length(amplitude_window) != 1L ||
      is.na(amplitude_window) || amplitude_window <= 0) {
    stop("ERROR: 'amplitude_window' must be a single number above 0 (the ",
         "length in seconds of the stretches in which each channel's own ",
         "level is measured), or Inf.", call. = FALSE)
  }
  window_samples <- if (is.infinite(amplitude_window)) {
    ncol(eeg_obj$data)
  } else {
    round(amplitude_window * eeg_obj$sampling_rate)
  }
  if (window_samples < 2) {
    stop("ERROR: 'amplitude_window' (", amplitude_window, " s) is shorter ",
         "than 2 samples at ", eeg_obj$sampling_rate, " Hz.", call. = FALSE)
  }

  eeg_idx <- which(eeg_obj$channel_types == "eeg" &
                      !(eeg_obj$channels %in% eeg_obj$bads))

  if (length(eeg_idx) == 0) {
    return(.empty_flags())
  }

  channels <- eeg_obj$channels[eeg_idx]
  data     <- eeg_obj$data[eeg_idx, , drop = FALSE]

  # The amplitude is measured from each channel's own level, re-measured in
  # every stretch (see .centre_by_stretch()). The std below does not need
  # this: a constant level does not change it.
  centred <- .centre_by_stretch(data, window_samples)

  # Round before thresholding - matches eeg_summary()'s existing order
  std_uv <- round(apply(data, 1, sd), 2)
  min_uv <- round(apply(centred, 1, min), 2)
  max_uv <- round(apply(centred, 1, max), 2)

  level_text <- if (is.infinite(amplitude_window)) {
    "its own level"
  } else {
    paste0("its own level, re-measured every ", format(amplitude_window), " s")
  }

  median_std <- median(std_uv, na.rm = TRUE)

  flags <- .empty_flags()

  for (i in seq_along(channels)) {

    ch   <- channels[i]
    std  <- std_uv[i]
    amax <- max(abs(min_uv[i]), abs(max_uv[i]))

    if (std < flag_flat_threshold) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Flat / possibly disconnected",
        value   = paste0("std = ", std, " uV  (threshold: < ",
                          flag_flat_threshold, " uV)"),
        stringsAsFactors = FALSE))
    }

    if (amax > flag_amplitude_threshold) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Excessive amplitude",
        value   = paste0("peak |amplitude| from ", level_text, " = ",
                          round(amax, 2), " uV  (threshold: > ",
                          flag_amplitude_threshold, " uV)"),
        stringsAsFactors = FALSE))
    }

    if (std > flag_outlier_sd_multiplier * median_std) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Outlier noise level",
        value   = paste0("std = ", std, " uV  (", flag_outlier_sd_multiplier,
                          " x median std = ",
                          round(flag_outlier_sd_multiplier * median_std, 2),
                          " uV)"),
        stringsAsFactors = FALSE))
    }
  }

  flags
}

#' Subtract Each Channel's Own Level, Stretch by Stretch (internal)
#'
#' Cuts the recording into consecutive stretches of \code{window_samples}
#' samples and, in every stretch, subtracts each channel's own average from
#' that channel. What is left is how much the channel wiggles around where it
#' stands. The level is measured again in every stretch because it creeps
#' (raw BioSemi data keeps a constant offset on every electrode, often tens of
#' mV, that also drifts slowly). As many whole stretches are used as fit into
#' the recording (at least one), so a last piece shorter than one stretch
#' joins the stretch before it, and a recording shorter than one stretch is a
#' single stretch.
#'
#' @param data Numeric matrix, channels x timepoints.
#' @param window_samples Number of samples in one stretch (at least 1).
#' @return A matrix of the same size and dimnames as \code{data}. In every
#'   stretch each row has mean 0.
#' @keywords internal
.centre_by_stretch <- function(data, window_samples) {

  n <- ncol(data)
  if (n == 0L) {
    return(data)
  }

  n_stretch <- max(1, n %/% window_samples)
  out       <- data

  for (s in seq_len(n_stretch)) {
    first <- (s - 1) * window_samples + 1
    last  <- if (s == n_stretch) n else s * window_samples
    idx   <- first:last
    block <- data[, idx, drop = FALSE]
    out[, idx] <- block - rowMeans(block)
  }

  out
}

#' Find Each Channel's Nearest Physical Neighbors
#'
#' Internal helper computing, for every channel in \code{positions}, its
#' \code{k} nearest other channels by Euclidean distance on the montage's
#' Cartesian x/y/z coordinates (millimeters). Pure geometry - no signal data
#' is touched - so this is also meant to be reusable by other channel-
#' geometry needs (e.g. a future bad-channel interpolation function), not
#' just \code{\link{find_bad_channels}}.
#'
#' @param positions A data frame with at least columns \code{channel, x, y,
#'   z} - typically \code{montage$positions} or a subset of it (see
#'   \code{\link{create_montage}}). Extra columns (e.g. \code{theta, phi,
#'   radius}) are ignored - \code{radius} in particular is a fixed constant
#'   for every channel in this package's electrode database, so it carries
#'   no spatial information.
#' @param k Integer, number of nearest neighbors to find per channel.
#'   Default: 4. Automatically capped (with a warning) if fewer than
#'   \code{k} other channels are available.
#'
#' @return A data frame in long format, one row per (channel, neighbor)
#'   pair, with columns \code{channel, neighbor} (character), \code{rank}
#'   (integer, 1 = nearest), \code{distance} (numeric, mm). Zero rows if
#'   fewer than 2 channels are supplied.
#'
#' @keywords internal
find_channel_neighbors <- function(positions, k = 4) {

  n <- nrow(positions)

  if (n < 2) {
    return(data.frame(channel = character(0), neighbor = character(0),
                       rank = integer(0), distance = numeric(0),
                       stringsAsFactors = FALSE))
  }

  k_eff <- min(k, n - 1)
  if (k_eff < k) {
    warning("find_channel_neighbors(): only ", n - 1, " other channel(s) ",
            "available - using k = ", k_eff, " instead of the requested ",
            k, ".", call. = FALSE, immediate. = TRUE)
  }

  dist_mat <- as.matrix(dist(positions[, c("x", "y", "z")], method = "euclidean"))
  diag(dist_mat) <- Inf
  rownames(dist_mat) <- positions$channel
  colnames(dist_mat) <- positions$channel

  neighbor_rows <- lapply(positions$channel, function(ch) {
    nearest <- sort(dist_mat[ch, ])[seq_len(k_eff)]
    data.frame(
      channel  = ch,
      neighbor = names(nearest),
      rank     = seq_len(k_eff),
      distance = as.numeric(nearest),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, neighbor_rows)
}

#' Middle Value of Every Column, All at Once (internal)
#'
#' Does the same job as \code{apply(m, 2, median)} without calling
#' \code{median()} once per column: every column is sorted in a single step,
#' then the middle row is taken (an odd number of rows) or the mean of the two
#' middle rows (an even number). For the neighbor-correlation check this is
#' the difference between about 30 seconds and a fraction of a second for 20
#' seconds of 64-channel data, because the middle value is needed at every
#' sample of every channel. A column that holds a missing value (\code{NA} or
#' \code{NaN}) gives \code{NA}, like \code{median()} with its default
#' \code{na.rm = FALSE}.
#'
#' @param m Numeric matrix, rows = the values to take the middle of (for
#'   example the neighbor channels), columns = the samples.
#' @return A numeric vector with one value per column of \code{m}
#'   (\code{numeric(0)} for a matrix without columns).
#' @keywords internal
.column_medians <- function(m) {

  k <- nrow(m)
  n <- ncol(m)

  if (n == 0L) {
    return(numeric(0))
  }
  if (k == 0L) {
    return(rep(NA_real_, n))
  }

  sorted <- matrix(m[order(col(m), m)], nrow = k)

  middle <- if (k %% 2L == 1L) {
    sorted[(k + 1L) / 2L, ]
  } else {
    (sorted[k / 2L, ] + sorted[k / 2L + 1L, ]) / 2
  }

  middle[colSums(is.na(m)) > 0L] <- NA_real_
  middle
}

#' Flag Channels That Disagree With Their Physical Neighbors
#'
#' Internal helper: for each channel, compares its signal against the
#' (robust) average of the neighbor channels found by
#' \code{\link{find_channel_neighbors}}, and flags channels whose
#' correlation with that local reference falls below
#' \code{correlation_threshold}. Loosely inspired by the PREP pipeline's
#' channel-correlation criterion, restricted to physically nearby channels
#' rather than the whole cap - deliberately simple, not an attempt to
#' reimplement Local Outlier Factor (see \code{\link{flag_lof_outlier_channels}}
#' for that, a separate and intentionally unrestricted check).
#'
#' @param data Numeric matrix, channels x timepoints, in the same row order
#'   as \code{channels}.
#' @param channels Character vector of channel names, one per row of
#'   \code{data}.
#' @param neighbor_map A data frame as returned by
#'   \code{\link{find_channel_neighbors}} (columns \code{channel, neighbor,
#'   rank, distance}).
#' @param correlation_threshold Numeric, correlation below which a channel
#'   is flagged. Default: 0.45.
#' @param robust Logical. If \code{TRUE} (default), the neighbor reference
#'   signal is the per-timepoint median of the neighbor channels; if
#'   \code{FALSE}, the per-timepoint mean.
#'
#' @return A data frame with columns \code{channel, reason, value}
#'   (character), one row per flag raised. Zero rows (never \code{NULL}) if
#'   nothing is flagged. A channel with no neighbors available in
#'   \code{channels}/\code{data} is silently skipped (not flagged, not an
#'   error).
#'
#' @keywords internal
flag_low_neighbor_correlation <- function(data, channels, neighbor_map,
                                           correlation_threshold = 0.45,
                                           robust = TRUE) {

  flags <- .empty_flags()

  for (ch in channels) {

    neighbor_names <- neighbor_map$neighbor[neighbor_map$channel == ch]
    neighbor_names <- intersect(neighbor_names, channels)

    if (length(neighbor_names) == 0) next

    ch_row             <- data[match(ch, channels), ]
    neighbor_data_rows <- data[match(neighbor_names, channels), , drop = FALSE]

    reference <- if (isTRUE(robust)) {
      .column_medians(neighbor_data_rows)
    } else {
      colMeans(neighbor_data_rows)
    }

    # A zero-variance channel or neighbor-reference makes cor() return NA
    # and emit "the standard deviation is zero" - expected here (flat
    # channels are exactly what this package flags elsewhere), and already
    # handled deliberately below, so the raw warning is suppressed.
    correlation <- suppressWarnings(cor(ch_row, reference))

    if (is.na(correlation) || correlation < correlation_threshold) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Low correlation with neighboring channels",
        value   = paste0(
          "correlation = ", if (is.na(correlation)) "NA" else round(correlation, 3),
          "  (threshold: < ", correlation_threshold, "; ",
          length(neighbor_names), " neighbor(s): ",
          paste(neighbor_names, collapse = ", "), ")"),
        stringsAsFactors = FALSE))
    }
  }

  flags
}

#' Flag Globally Outlying Channels via Local Outlier Factor
#'
#' Internal helper mirroring MNE-Python's
#' \code{mne.preprocessing.find_bad_channels_lof()}: computes a Local
#' Outlier Factor (LOF) score for every channel by comparing its raw signal
#' directly against every other channel's (Euclidean distance across the
#' full time series), with no scalp position involved at all. Uses
#' \code{\link[dbscan]{lof}}. See Breunig et al. (2000), "LOF: Identifying
#' Density-Based Local Outliers".
#'
#' @param data Numeric matrix, channels x timepoints, in the same row order
#'   as \code{channels}. Rows are channels, columns are timepoints - do not
#'   transpose (this must match MNE's own \code{raw.get_data()}
#'   orientation, or the check silently scores "outlier timepoint" instead
#'   of "outlier channel", with no error to catch it).
#' @param channels Character vector of channel names, one per row of
#'   \code{data}.
#' @param n_neighbors Integer, number of comparison channels, mirroring
#'   MNE's \code{n_neighbors} (default 20; automatically capped if fewer
#'   channels are available). Passed to \code{\link[dbscan]{lof}} as
#'   \code{minPts = n_neighbors + 1}, since \code{dbscan::lof()}'s
#'   \code{minPts} counts the point itself while MNE/sklearn's
#'   \code{n_neighbors} does not.
#' @param threshold Numeric, a channel is flagged if its LOF score is at or
#'   above this value. Default: 1.5, matching MNE's default. Both
#'   \code{dbscan::lof()} and MNE's score use the same convention (~1.0 =
#'   inlier, higher = more outlier-like), so no rescaling is needed.
#'
#' @return A data frame with columns \code{channel, reason, value}
#'   (character), one row per flagged channel. Zero rows (never \code{NULL})
#'   if nothing is flagged.
#'
#' @seealso \code{\link{find_bad_channels}}
#'
#' @importFrom dbscan lof
#' @keywords internal
flag_lof_outlier_channels <- function(data, channels,
                                       n_neighbors = 20,
                                       threshold   = 1.5) {

  n <- nrow(data)

  if (n < 2) {
    stop("flag_lof_outlier_channels(): need at least 2 channels, got ",
         n, ".", call. = FALSE)
  }

  n_neighbors_eff <- min(n_neighbors, n - 1)
  if (n_neighbors_eff < n_neighbors) {
    warning("flag_lof_outlier_channels(): only ", n - 1, " other channel(s) ",
            "available - using n_neighbors = ", n_neighbors_eff,
            " instead of the requested ", n_neighbors, ".",
            call. = FALSE, immediate. = TRUE)
  }

  lof_scores <- dbscan::lof(data, minPts = n_neighbors_eff + 1)

  flags <- .empty_flags()

  for (i in seq_along(channels)) {
    if (lof_scores[i] >= threshold) {
      flags <- rbind(flags, data.frame(
        channel = channels[i],
        reason  = "Local Outlier Factor (global)",
        value   = paste0("LOF score = ", round(lof_scores[i], 3),
                          "  (threshold: >= ", threshold,
                          "; n_neighbors = ", n_neighbors_eff, ")"),
        stringsAsFactors = FALSE))
    }
  }

  flags
}

#' Empty Flags Data Frame (internal)
#'
#' Zero-row \code{channel, reason, value} data frame shared by every check
#' in this file, so \code{\link{find_bad_channels}} never has to
#' special-case a \code{NULL} result.
#'
#' @return A zero-row data frame with columns \code{channel, reason, value}
#'   (all character).
#' @keywords internal
.empty_flags <- function() {
  data.frame(channel = character(0), reason = character(0),
             value = character(0), stringsAsFactors = FALSE)
}
