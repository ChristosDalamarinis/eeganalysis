#' ============================================================================
#'                  ICA Phase 2: Automatic Bad-Component Detection
#' ============================================================================
#'
#' This module is Phase 2 of the three-phase ICA pipeline documented in
#' notes/mne-ica-pipeline-reference.md: Phase 1 (fit - R/ica1.R) produces
#' independent components; this file scores each component against a
#' reference (an EOG/ECG channel, or the component's own shape) and flags
#' outliers as candidate artifacts; Phase 3 (apply - R/ica1.R) removes
#' whichever components were actually chosen for exclusion.
#'
#' "Detect, don't decide": every function here only ever writes to
#' ica$labels_ (see new_ica(), R/ica1.R) - never to ica$exclude. Deciding
#' what actually gets removed stays a separate, explicit set_exclude() call,
#' mirroring MNE's own contract (ICA.find_bads_*() populates ica.labels_ but
#' leaves ica.exclude untouched unless the caller extends it themselves).
#'
#' Ported from the vendored MNE source (python/ica/bads.py, python/ica/eog.py,
#' python/ica/ecg.py, python/ica/ctps_.py, python/ica/ica.py lines 1385-2172
#' and 1934-2078), not just the prose summary in
#' notes/mne-ica-pipeline-reference.md, so constants and merge logic match
#' exactly.
#'
#' Scope (see the "R status" table in notes/mne-ica-pipeline-reference.md):
#' find_bads_eog(), find_bads_ecg() (both method = "correlation" and
#' method = "ctps"), find_bads_muscle(). Deferred: corrmap() (cross-subject
#' template matching).
#'
#' Author: Christos Dalamarinis
#' Date: Sep - 2026
#' Status: find_bads_eog(), find_bads_ecg() (correlation + ctps),
#'         find_bads_muscle() built. corrmap() deferred.
#' Tested: Tested
#' ============================================================================
#
#
# ============================================================================
#                     SHARED INTERNAL HELPERS
# ============================================================================

# ----------------------------------------------------------------------------
# .find_outliers() - iterative adaptive z-score outlier rule
# ----------------------------------------------------------------------------
#' Find Outliers via Iterated Z-Scoring (internal)
#'
#' Ported from \code{python/ica/bads.py:_find_outliers} (the shared decision
#' rule behind MNE's EOG/ECG-correlation/MEG-reference detectors). Z-scores
#' the elements not yet flagged, flags anything with \code{|z| > threshold}
#' (or one-sided, via \code{tail}), then repeats - excluding already-flagged
#' elements from the z-score computation each time - since one huge outlier
#' can otherwise deflate the mean/sd enough to hide a second, smaller-but-real
#' one. Stops early once an iteration flags nothing new.
#'
#' Uses \code{\link{.population_sd}} (ddof = 0), matching
#' \code{scipy.stats.zscore}'s default exactly - R's own \code{sd()}/
#' \code{scale()} use ddof = 1 and would silently diverge from MNE's numbers.
#'
#' @param scores Numeric vector of per-component scores.
#' @param threshold Numeric. Flag \code{|z| > threshold} (or one-sided, per
#'   \code{tail}). Default: \code{3.0}.
#' @param max_iter Integer. Maximum number of re-scoring passes. Default:
#'   \code{2}.
#' @param tail \code{0}, \code{1}, or \code{-1}. \code{0} (default) flags
#'   outliers on both ends; \code{1} only large-positive; \code{-1} only
#'   large-negative.
#' @return Integer vector of indices into \code{scores} flagged as outliers.
#' @keywords internal
.find_outliers <- function(scores, threshold = 3.0, max_iter = 2, tail = 0) {

  n   <- length(scores)
  bad <- logical(n)

  for (i in seq_len(max_iter)) {
    idx <- which(!bad)
    if (length(idx) < 2) break  # z-score of <2 points is undefined

    x <- scores[idx]
    z <- (x - mean(x)) / .population_sd(x)

    this_z <- if (tail == 0) {
      abs(z)
    } else if (tail == 1) {
      z
    } else if (tail == -1) {
      -z
    } else {
      stop("ERROR: 'tail' must be 0, 1, or -1, got ", tail, ".", call. = FALSE)
    }

    local_bad <- this_z > threshold
    local_bad[is.na(local_bad)] <- FALSE  # sd=0 among the remainder -> NaN
    # z; numpy's `nan > threshold` is False, so match that instead of R's NA
    if (!any(local_bad)) break
    bad[idx[local_bad]] <- TRUE
  }

  which(bad)
}

# ----------------------------------------------------------------------------
# .resolve_reference_channels() - locate EOG/ECG channel(s) by name
# ----------------------------------------------------------------------------
#' Resolve Reference Channel Name(s) for Artifact Detection (internal)
#'
#' MNE auto-detects EOG/ECG channels via a dedicated channel \emph{type}
#' (\code{eog}/\code{ecg}) it stores in its own info structure
#' (\code{_get_eog_channel_index()}, \code{python/ica/eog.py:186}). This
#' package has no such sub-type - \code{eeg_obj$channel_types} only
#' distinguishes \code{"eeg"}/\code{"external"}/\code{"status"} (see
#' \code{classify_channels()}, R/eeg_class.R) - the specific physiological
#' role instead lives in the channel \emph{name} itself, via
#' \code{\link{apply_external_labels}} (e.g. \code{"EOG_L (EXG1)"}). So
#' auto-detection here is a case-insensitive substring match of
#' \code{pattern} among \code{channel_types == "external"} channel names.
#'
#' @param eeg An object of class 'eeg'.
#' @param ch_name \code{NULL}, or a character vector of explicit channel
#'   name(s) to use instead of auto-detection. Every name must exist in
#'   \code{eeg$channels}.
#' @param pattern Character scalar. Case-insensitive substring to search for
#'   among external channel names when \code{ch_name} is \code{NULL} (e.g.
#'   \code{"EOG"}, \code{"ECG"}).
#' @return Character vector of resolved channel name(s), length >= 1.
#' @keywords internal
.resolve_reference_channels <- function(eeg, ch_name, pattern) {

  if (!is.null(ch_name)) {
    idx <- match(ch_name, eeg$channels)
    if (anyNA(idx)) {
      stop("ERROR: 'ch_name' contains channel name(s) not found in ",
           "eeg$channels: ", paste(ch_name[is.na(idx)], collapse = ", "),
           call. = FALSE)
    }
    return(eeg$channels[idx])
  }

  ext_idx <- which(eeg$channel_types == "external")
  matches <- ext_idx[grepl(pattern, eeg$channels[ext_idx], ignore.case = TRUE)]

  if (length(matches) == 0) {
    stop("ERROR: no channel name matching '", pattern, "' found among ",
         "eeg$channels' external channels (",
         paste(eeg$channels[ext_idx], collapse = ", "),
         "). Pass 'ch_name' explicitly.", call. = FALSE)
  }

  eeg$channels[matches]
}

# ----------------------------------------------------------------------------
# .score_sources() - correlate every IC against a reference signal
# ----------------------------------------------------------------------------
#' Score Every Independent Component Against a Reference Signal (internal)
#'
#' Ported from \code{ica.py:score_sources} (1385) + \code{_band_pass_filter},
#' restricted to Pearson correlation against a fixed target (MNE's more
#' general \code{score_func}/no-target-skewness path is not needed by any
#' detector in this file). Band-pass filters every IC time-course and the
#' target signal to \code{[l_freq, h_freq]} via the existing
#' \code{\link{.fir_filter_vector}} (R/annotations.R - already used by
#' \code{\link{annotate_muscle}} for the same kind of per-vector bandpass),
#' then correlates each filtered IC against the filtered target.
#'
#' @param ica A fitted \code{eeg_ica} object.
#' @param eeg An \code{eeg} object matching \code{ica$ch_names} (see
#'   \code{\link{get_sources}}).
#' @param target_data Numeric vector, length \code{ncol(eeg$data)} - the
#'   reference signal (e.g. one row of \code{eeg$data}).
#' @param l_freq,h_freq Numeric. Band-pass edges in Hz.
#' @return Named numeric vector, length \code{ica$n_components_}, one Pearson
#'   correlation per component, named \code{"IC1"}, \code{"IC2"}, ...
#' @keywords internal
.score_sources <- function(ica, eeg, target_data, l_freq, h_freq) {

  sources <- get_sources(ica, eeg)
  sfreq   <- eeg$sampling_rate

  target_filt <- .fir_filter_vector(target_data, sfreq, l_freq = l_freq, h_freq = h_freq)

  scores <- vapply(seq_len(nrow(sources)), function(i) {
    ic_filt <- .fir_filter_vector(sources[i, ], sfreq, l_freq = l_freq, h_freq = h_freq)
    cor(ic_filt, target_filt)
  }, numeric(1))

  names(scores) <- rownames(sources)
  scores
}

# ----------------------------------------------------------------------------
# .find_bads_ch() - shared EOG/ECG(correlation) detection + merge logic
# ----------------------------------------------------------------------------
#' Detect Components Correlated With One or More Reference Channels (internal)
#'
#' Ported from \code{ica.py:_find_bads_ch} (1503-1576) - the shared machinery
#' behind both \code{\link{find_bads_eog}} and
#' \code{\link{find_bads_ecg}(method = "correlation")}. For each channel in
#' \code{chs}: scores every component (\code{\link{.score_sources}}), flags
#' outliers via \code{\link{.find_outliers}} (\code{measure = "zscore"}) or a
#' direct threshold on \code{abs(score)} (\code{measure = "correlation"}),
#' and records the per-channel result into
#' \code{ica$labels_[["<prefix>/<i>/<ch>"]]}. Results across every channel in
#' \code{chs} are then merged: sorted by \strong{descending
#' \code{abs(score)}}, deduplicated keeping the first (highest-scoring)
#' occurrence of a repeated component - exactly reproducing the Python merge
#' loop at lines 1561-1574 - and written to \code{ica$labels_[[prefix]]}.
#' Never touches \code{ica$exclude} (see file header).
#'
#' @param ica A fitted \code{eeg_ica} object.
#' @param eeg An \code{eeg} object matching \code{ica$ch_names}.
#' @param chs Character vector of one or more reference channel names (must
#'   exist in \code{eeg$channels}).
#' @param threshold Numeric. Passed to \code{\link{.find_outliers}} (as
#'   \code{measure = "zscore"}) or compared directly against
#'   \code{abs(score)} (\code{measure = "correlation"}).
#' @param l_freq,h_freq Numeric. Band-pass edges in Hz, passed to
#'   \code{\link{.score_sources}}.
#' @param prefix Character scalar, e.g. \code{"eog"} or \code{"ecg"} - the
#'   \code{labels_} key the merged result is stored under.
#' @param measure \code{"zscore"} or \code{"correlation"}.
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with
#'   \code{labels_} updated.
#' @keywords internal
.find_bads_ch <- function(ica, eeg, chs, threshold, l_freq, h_freq, prefix, measure) {

  all_scores <- vector("list", length(chs))
  all_idx    <- vector("list", length(chs))

  for (i in seq_along(chs)) {
    ch          <- chs[i]
    target_data <- eeg$data[match(ch, eeg$channels), ]
    scores      <- .score_sources(ica, eeg, target_data, l_freq, h_freq)

    this_idx <- if (identical(measure, "zscore")) {
      .find_outliers(scores, threshold = threshold)
    } else {
      unname(which(abs(scores) > threshold))
    }

    ica$labels_[[paste0(prefix, "/", i, "/", ch)]] <- this_idx
    all_scores[[i]] <- scores
    all_idx[[i]]    <- this_idx
  }

  flat_idx <- unlist(all_idx)

  if (length(flat_idx) > 0) {
    flat_scores <- unlist(lapply(seq_along(chs), function(i) all_scores[[i]][all_idx[[i]]]))
    merged <- flat_idx[order(-abs(flat_scores))]
    merged <- merged[!duplicated(merged)]
  } else {
    merged <- integer(0)
  }

  ica$labels_[[prefix]] <- merged
  ica
}

# ----------------------------------------------------------------------------
# .qrs_detector() - locate heartbeats in an ECG channel
# ----------------------------------------------------------------------------
#' Detect QRS (Heartbeat) Peaks in an ECG Channel (internal)
#'
#' Ported from \code{python/ica/ecg.py:qrs_detector()} (18-154). Slides a
#' half-second window over \code{abs(ecg)}; whenever a window's first sample
#' exceeds an adaptive threshold, records the window's peak position, keeping
#' only candidates whose window RMS is below \code{mean(rms) + levels *
#' population_sd(rms)} (see \code{\link{.population_sd}}, R/ica1.R, for exact
#' ddof = 0 parity) and whose in-window threshold-crossing count is below
#' \code{n_thresh}. Assumes \code{ecg} is already appropriately band-pass
#' filtered by the caller (see \code{\link{find_bads_ecg}}'s
#' \code{method = "ctps"} branch) - unlike the vendored source, this function
#' takes no \code{l_freq}/\code{h_freq} of its own, since in the one call
#' path that reaches it, filtering always happens upstream and its own
#' filtering step is always disabled.
#'
#' When \code{thresh_value = "auto"} (default), tries 16 threshold
#' multipliers (\code{seq(0.30, 1.05, by = 0.05)}) and keeps whichever run's
#' implied heart rate lands closest to the median of any in-range (40-160
#' bpm) result, or closest to 80 bpm if none qualify - ported exactly,
#' including that fallback.
#'
#' @param ecg Numeric vector. The (already filtered) ECG channel.
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param thresh_value \code{"auto"} (default), or a single numeric threshold
#'   multiplier.
#' @param levels Numeric. RMS-based noise rejection strictness. Default:
#'   \code{2.5}.
#' @param n_thresh Integer. Max in-window threshold-crossing count to accept.
#'   Default: \code{3}.
#' @param tstart Numeric. Seconds to skip at the start before detecting.
#'   Default: \code{0}.
#' @return Integer vector of detected heartbeat sample positions (1-indexed,
#'   into \code{ecg}), empty if none found.
#' @keywords internal
.qrs_detector <- function(ecg, sfreq, thresh_value = "auto", levels = 2.5,
                           n_thresh = 3, tstart = 0) {

  win_size <- round((60 * sfreq) / 120)  # = round(sfreq / 2)
  init     <- round(sfreq)               # 1 second, in samples

  n_samples_start <- round(sfreq * tstart)
  ecg_abs  <- abs(ecg)[(n_samples_start + 1):length(ecg)]
  n_points <- length(ecg_abs)

  if (n_points < 3 * init) {
    stop("ERROR: .qrs_detector() needs at least 3 seconds of data after ",
         "'tstart' to seed its initial threshold, got ",
         round(n_points / sfreq, 2), " s.", call. = FALSE)
  }

  maxpt <- c(
    max(ecg_abs[1:init]),
    max(ecg_abs[(init + 1):(2 * init)]),
    max(ecg_abs[(2 * init + 1):(3 * init)])
  )
  init_max <- mean(maxpt)

  thresh_runs <- if (identical(thresh_value, "auto")) {
    seq(0.30, 1.05, by = 0.05)
  } else {
    thresh_value
  }

  clean_events <- list()

  for (tv in thresh_runs) {
    thresh1  <- init_max * tv
    numcross <- integer(0)
    time_pos <- integer(0)
    rms      <- numeric(0)

    ii <- 0L  # kept 0-indexed throughout, mirroring the source exactly;
    # only ever offset by +1 when actually subsetting an R vector below.
    while (ii < (n_points - win_size)) {
      window <- ecg_abs[(ii + 1):(ii + win_size)]
      if (window[1] > thresh1) {
        max_time <- which.max(window) - 1L  # 0-indexed, matches np.argmax
        time_pos <- c(time_pos, ii + max_time)
        numcross <- c(numcross, sum(diff(as.integer(window > thresh1))))
        rms      <- c(rms, sqrt(mean(window^2)))
        ii <- ii + win_size
      } else {
        ii <- ii + 1L
      }
    }

    if (length(rms) == 0) {
      # No candidates at all for this threshold - a dummy entry that can
      # never pass the RMS/crossing checks below, so this run contributes
      # no events (avoids mean()/sd() of an empty vector, matches the
      # vendored source's equivalent empty-safe fallback).
      rms      <- 0
      time_pos <- 0L
      numcross <- n_thresh
    }

    rms_thresh <- mean(rms) + .population_sd(rms) * levels
    keep <- which(rms < rms_thresh & numcross < n_thresh)
    ce   <- time_pos[keep] + n_samples_start  # still 0-indexed

    if (length(ce) > 0) clean_events[[length(clean_events) + 1]] <- ce
  }

  if (length(clean_events) > 0) {
    rates <- vapply(clean_events, function(cev) {
      60 * length(cev) / (length(ecg) / sfreq)
    }, numeric(1))

    in_range   <- which(rates <= 160 & rates >= 40)
    ideal_rate <- if (length(in_range) > 0) stats::median(rates[in_range]) else 80

    events_0indexed <- clean_events[[which.min(abs(rates - ideal_rate))]]
  } else {
    events_0indexed <- integer(0)
  }

  events_0indexed + 1L  # 0-indexed -> 1-indexed R sample positions
}

# ----------------------------------------------------------------------------
# .hilbert_phase() - normalized instantaneous phase via FFT
# ----------------------------------------------------------------------------
#' Hilbert Instantaneous Phase via FFT (internal)
#'
#' Sibling of \code{\link{.hilbert_envelope}} (R/annotations.R) - identical
#' analytic-signal FFT construction, but returns the normalized
#' \emph{phase} (\code{[0, 1)}, one full cycle) instead of the magnitude.
#' Matches \code{python/ica/ctps_.py:_compute_normalized_phase()}.
#'
#' @param x Numeric vector.
#' @return Numeric vector, same length as \code{x}, values in \code{[0, 1)}.
#' @keywords internal
.hilbert_phase <- function(x) {

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
  (Arg(analytic) + pi) / (2 * pi)
}

# ----------------------------------------------------------------------------
# .kuiper_test() - Kuiper's test of uniformity, per time-column
# ----------------------------------------------------------------------------
#' Kuiper's Statistic for Each Column of a Phase Matrix (internal)
#'
#' Ported from \code{python/ica/ctps_.py:kuiper()} (82-118). For each column
#' of \code{phase_matrix} (one within-epoch time point), sorts the epoch
#' values ascending and compares them against the ideal uniform-CDF
#' reference points, returning the Kuiper statistic \code{d1 + d2}.
#'
#' @param phase_matrix Numeric matrix, \code{n_epochs x n_times}: one row
#'   per heartbeat epoch, values in \code{[0, 1)} (see
#'   \code{\link{.hilbert_phase}}).
#' @return Numeric vector, length \code{ncol(phase_matrix)}.
#' @keywords internal
.kuiper_test <- function(phase_matrix) {

  n_trials <- nrow(phase_matrix)
  n_times  <- ncol(phase_matrix)

  sorted <- apply(phase_matrix, 2, sort)
  dim(sorted) <- c(n_trials, n_times)  # apply() simplifies to a vector when n_times == 1

  j1 <- seq_len(n_trials) / n_trials        # 1/n, 2/n, ..., n/n
  j2 <- (seq_len(n_trials) - 1) / n_trials  # 0/n, 1/n, ..., (n-1)/n

  d1 <- apply(j1 - sorted, 2, max)
  d2 <- apply(sorted - j2, 2, max)

  d1 + d2
}

# ----------------------------------------------------------------------------
# .prob_kuiper() - Kuiper statistic -> normalized [0,1] significance
# ----------------------------------------------------------------------------
#' Normalized Significance of a Kuiper Statistic (internal)
#'
#' Ported from \code{python/ica/ctps_.py:_prob_kuiper()} (121-167) - the
#' asymptotic significance formula (Stephens 1970; Kuiper 1962), evaluated as
#' a 100-term series per \code{d} value. Uses a numerically stable
#' "subtract the max exponent before exponentiating" reduction (the same
#' role \code{scipy.special.logsumexp} plays in the source), written out
#' directly here since this is its only call site - not worth a
#' general-purpose \code{logsumexp()} utility for one use.
#'
#' @param d Numeric vector. Kuiper statistic(s) (see \code{\link{.kuiper_test}}).
#' @param n_eff Integer. Effective number of trials (epochs) the statistic
#'   was computed from.
#' @return Numeric vector, same length as \code{d}, values in \code{[0, 1]} -
#'   \code{0} forced whenever the statistic is too small to be
#'   distinguishable from a uniform distribution (\code{k_lambda < 0.4}).
#' @keywords internal
.prob_kuiper <- function(d, n_eff) {

  n_points <- 100

  en       <- sqrt(n_eff)
  k_lambda <- (en + 0.155 + 0.24 / en) * d
  l2       <- k_lambda^2

  j2 <- seq_len(n_points)^2                                   # length 100
  a    <- outer(j2, l2, function(jj, ll) -2 * jj * ll)         # 100 x length(d)
  fact <- outer(j2, l2, function(jj, ll) 4 * jj * ll - 1)
  b    <- 2 * fact

  max_a  <- apply(a, 2, max)
  scaled <- sweep(a, 2, max_a, FUN = "-")
  lse    <- max_a + log(colSums(b * exp(scaled)))

  pk_norm <- -lse / (2 * n_eff)

  pk_norm[k_lambda < 0.4] <- 0  # no difference from uniform
  pk_norm[pk_norm > 1]    <- 1  # round-off guard

  pk_norm
}

# ----------------------------------------------------------------------------
# .get_ctps_threshold() - automatic CTPS significance cutoff
# ----------------------------------------------------------------------------
#' Automatic Threshold for CTPS Detection (internal)
#'
#' Ported from \code{ica.py:_get_ctps_threshold()} (1578-1600). Searches 99
#' candidate Kuiper-index values for whichever one's implied significance is
#' closest to \code{10^-pk_threshold} - what
#' \code{find_bads_ecg(threshold = "auto", method = "ctps")} resolves to.
#'
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param pk_threshold Numeric. Target significance exponent. Default:
#'   \code{20} (i.e. target \code{1e-20}), matching MNE's own default -
#'   deliberately strict, since CTPS is a trial-by-trial test and false
#'   positives compound across many components.
#' @return A single numeric value in \code{(0, 1)}.
#' @keywords internal
.get_ctps_threshold <- function(sfreq, pk_threshold = 20) {

  Vs <- seq_len(99) / 100  # 0.01, 0.02, ..., 0.99
  C  <- sqrt(sfreq) + 0.155 + 0.24 / sqrt(sfreq)
  Pks <- 2 * (4 * (Vs * C)^2 - 1) * exp(-2 * (Vs * C)^2)

  Vs[which.min(abs(Pks - 10^(-pk_threshold)))]
}


# ============================================================================
#                    find_bads_eog() - PUBLIC: EYE-MOVEMENT DETECTION
# ============================================================================

#' Detect Eye-Movement (EOG) Related Components
#'
#' Correlates every independent component against an EOG channel (or every
#' matching channel, if more than one) and flags outliers as candidate eye-
#' movement artifacts. Ported from \code{ica.py:find_bads_eog} (2081-2172).
#' \strong{Only writes to \code{ica$labels_}} - never to \code{ica$exclude}
#' (see file header); call \code{\link{set_exclude}} yourself to actually
#' remove anything.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{\link{fit_ica}}).
#' @param eeg The \code{eeg} object \code{ica} was fit on (or one with
#'   matching channels - see \code{\link{get_sources}}).
#' @param ch_name \code{NULL} (default), or a character vector of explicit
#'   EOG channel name(s). If \code{NULL}, every \code{channel_types ==
#'   "external"} channel whose name contains \code{"EOG"} (case-insensitive)
#'   is used (see \code{\link{.resolve_reference_channels}}) - errors if none
#'   is found.
#' @param threshold Numeric, or \code{"auto"} (default). \code{"auto"}
#'   resolves to \code{3.0} if \code{measure = "zscore"}, or \code{0.9} if
#'   \code{measure = "correlation"}.
#' @param l_freq,h_freq Numeric. Band-pass edges in Hz applied to both the IC
#'   time-courses and the EOG channel before correlating. Default:
#'   \code{c(1, 10)} - the eye-movement band.
#' @param measure \code{"zscore"} (default) - iterative adaptive z-scoring of
#'   the correlations (\code{\link{.find_outliers}}) - or
#'   \code{"correlation"} - a direct threshold on \code{abs(correlation)}.
#'
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with
#'   \code{labels_[["eog"]]} (merged, sorted by strength) and
#'   \code{labels_[["eog/<i>/<ch_name>"]]} (per-channel) populated. Since R
#'   does not mutate arguments in place, the caller must reassign the result
#'   (\code{ica <- find_bads_eog(ica, eeg)}).
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- find_bads_eog(ica, eeg)
#'   ica$labels_$eog
#'   ica <- set_exclude(ica, ica$labels_$eog)
#' }
#'
#' @seealso \code{\link{find_bads_ecg}}, \code{\link{find_bads_muscle}},
#'   \code{\link{set_exclude}}
#'
#' @export
find_bads_eog <- function(ica, eeg, ch_name = NULL, threshold = "auto",
                           l_freq = 1, h_freq = 10, measure = "zscore") {

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).",
         call. = FALSE)
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.", call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }
  if (!identical(measure, "zscore") && !identical(measure, "correlation")) {
    stop("ERROR: 'measure' must be 'zscore' or 'correlation', got '",
         measure, "'.", call. = FALSE)
  }

  if (identical(threshold, "auto")) {
    threshold <- if (identical(measure, "zscore")) 3.0 else 0.9
  }

  eog_chs <- .resolve_reference_channels(eeg, ch_name, pattern = "EOG")

  .find_bads_ch(ica, eeg, eog_chs, threshold = threshold, l_freq = l_freq,
                h_freq = h_freq, prefix = "eog", measure = measure)
}


# ============================================================================
#                    find_bads_ecg() - PUBLIC: HEARTBEAT DETECTION
# ============================================================================

#' Detect Heartbeat (ECG) Related Components
#'
#' Correlates every independent component against an ECG channel (or every
#' matching channel, if more than one) and flags outliers as candidate
#' heartbeat artifacts. Ported from \code{ica.py:find_bads_ecg}'s
#' \code{method = "correlation"} branch (1603-1760). \strong{Only writes to
#' \code{ica$labels_}} - never to \code{ica$exclude} (see file header); call
#' \code{\link{set_exclude}} yourself to actually remove anything.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{\link{fit_ica}}).
#' @param eeg The \code{eeg} object \code{ica} was fit on (or one with
#'   matching channels - see \code{\link{get_sources}}).
#' @param ch_name \code{NULL} (default), or a character vector of explicit
#'   ECG channel name(s). If \code{NULL}, every \code{channel_types ==
#'   "external"} channel whose name contains \code{"ECG"} (case-insensitive)
#'   is used (see \code{\link{.resolve_reference_channels}}) - errors if none
#'   is found.
#' @param threshold Numeric, or \code{"auto"} (default). \code{"auto"}
#'   resolves to \code{3.0} if \code{measure = "zscore"}, or \code{0.9} if
#'   \code{measure = "correlation"}.
#' @param l_freq,h_freq Numeric. Band-pass edges in Hz applied to both the IC
#'   time-courses and the ECG channel before correlating. Default:
#'   \code{c(8, 16)} - the heartbeat band.
#' @param method \code{"correlation"} (default) or \code{"ctps"} - MNE's
#'   actual default, and a meaningfully more powerful phase-locking method:
#'   instead of asking whether a component roughly moves with the ECG
#'   channel, it asks whether the component's instantaneous phase lands in
#'   nearly the same spot, every single heartbeat, reliably. Ported from
#'   \code{ica.py}'s ctps branch (1707-1741) plus \code{python/ica/ecg.py}'s
#'   \code{qrs_detector()} and \code{python/ica/ctps_.py} in full - see
#'   \code{\link{.qrs_detector}}, \code{\link{.hilbert_phase}},
#'   \code{\link{.kuiper_test}}, \code{\link{.prob_kuiper}},
#'   \code{\link{.get_ctps_threshold}}. \code{threshold = "auto"} resolves
#'   via \code{\link{.get_ctps_threshold}} instead of the correlation path's
#'   \code{3.0}/\code{0.9}. \code{l_freq}/\code{h_freq} only filter the ECG
#'   channel for heartbeat detection here - unlike the correlation path, the
#'   IC epochs analyzed are never filtered (matches the vendored source's own
#'   assumption that "the sources are already appropriately filtered" before
#'   \code{ctps()} ever sees them). \code{measure} has no effect when
#'   \code{method = "ctps"} (still validated, just unused - matches MNE).
#' @param measure \code{"zscore"} (default) - iterative adaptive z-scoring of
#'   the correlations (\code{\link{.find_outliers}}) - or
#'   \code{"correlation"} - a direct threshold on \code{abs(correlation)}.
#'   Only applies to \code{method = "correlation"}.
#'
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with
#'   \code{labels_[["ecg"]]} (sorted by strength) populated. For
#'   \code{method = "correlation"}, per-channel detail is also stored at
#'   \code{labels_[["ecg/<i>/<ch_name>"]]} (one or more reference channels);
#'   for \code{method = "ctps"} (always exactly one ECG channel), at
#'   \code{labels_[["ecg/<ch_name>"]]} instead (matches MNE's own naming in
#'   each path). Since R does not mutate arguments in place, the caller must
#'   reassign the result (\code{ica <- find_bads_ecg(ica, eeg)}).
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- find_bads_ecg(ica, eeg, method = "ctps")
#'   ica$labels_$ecg
#'   ica <- set_exclude(ica, ica$labels_$ecg)
#' }
#'
#' @seealso \code{\link{find_bads_eog}}, \code{\link{find_bads_muscle}},
#'   \code{\link{set_exclude}}
#'
#' @export
find_bads_ecg <- function(ica, eeg, ch_name = NULL, threshold = "auto",
                           l_freq = 8, h_freq = 16, method = "correlation",
                           measure = "zscore") {

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).",
         call. = FALSE)
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.", call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }
  if (!identical(measure, "zscore") && !identical(measure, "correlation")) {
    stop("ERROR: 'measure' must be 'zscore' or 'correlation', got '",
         measure, "'.", call. = FALSE)
  }
  if (identical(method, "ctps")) {

    ecg_chs <- .resolve_reference_channels(eeg, ch_name, pattern = "ECG")
    if (length(ecg_chs) > 1) {
      warning("More than one ECG-like channel found for method = 'ctps' ",
              "(which only supports one) - using '", ecg_chs[1], "' only. ",
              "Pass 'ch_name' explicitly to choose a different one.",
              call. = FALSE)
    }
    ecg_ch <- ecg_chs[1]
    sfreq  <- eeg$sampling_rate

    ecg_filt <- .fir_filter_vector(eeg$data[match(ecg_ch, eeg$channels), ],
                                    sfreq, l_freq = l_freq, h_freq = h_freq)
    beats <- .qrs_detector(ecg_filt, sfreq)

    if (identical(threshold, "auto")) {
      threshold <- .get_ctps_threshold(sfreq)
    }

    tmin_samp <- round(-0.5 * sfreq)
    tmax_samp <- round( 0.5 * sfreq)
    n_total   <- ncol(eeg$data)

    starts <- beats + tmin_samp
    stops  <- beats + tmax_samp
    keep   <- which(starts >= 1 & stops <= n_total)

    if (length(keep) == 0) {
      stop("ERROR: no complete heartbeat epochs found (either no heartbeats ",
           "were detected in '", ecg_ch, "', or every detected one was too ",
           "close to the start/end of the recording). Consider changing ",
           "'l_freq'/'h_freq', or check the ECG channel.", call. = FALSE)
    }
    starts   <- starts[keep]
    stops    <- stops[keep]
    n_epochs <- length(starts)
    n_times  <- tmax_samp - tmin_samp + 1

    sources <- get_sources(ica, eeg)
    n_comp  <- ica$n_components_

    scores <- vapply(seq_len(n_comp), function(k) {
      epoch_mat <- t(vapply(seq_len(n_epochs), function(e) {
        sources[k, starts[e]:stops[e]]
      }, numeric(n_times)))
      phase_mat <- t(apply(epoch_mat, 1, .hilbert_phase))
      d  <- .kuiper_test(phase_mat)
      pk <- .prob_kuiper(d, n_epochs)
      max(pk)
    }, numeric(1))

    ecg_idx <- which(scores >= threshold)
    ecg_idx <- ecg_idx[order(-scores[ecg_idx])]

    ica$labels_[["ecg"]] <- ecg_idx
    ica$labels_[[paste0("ecg/", ecg_ch)]] <- ecg_idx

    return(ica)
  }
  if (!identical(method, "correlation")) {
    stop("ERROR: 'method' must be 'correlation' or 'ctps', got '",
         method, "'.", call. = FALSE)
  }

  if (identical(threshold, "auto")) {
    threshold <- if (identical(measure, "zscore")) 3.0 else 0.9
  }

  ecg_chs <- .resolve_reference_channels(eeg, ch_name, pattern = "ECG")

  .find_bads_ch(ica, eeg, ecg_chs, threshold = threshold, l_freq = l_freq,
                h_freq = h_freq, prefix = "ecg", measure = measure)
}


# ============================================================================
#                    find_bads_muscle() - PUBLIC: MUSCLE DETECTION
# ============================================================================

#' Detect Muscle (EMG) Related Components
#'
#' Flags components whose own shape looks like muscle activity - no
#' reference channel needed. Ported from \code{ica.py:find_bads_muscle}
#' (1934-2078) with its exact constants, combining up to three criteria
#' (each squashed through a logistic curve into \code{[0, 1]}, then
#' multiplied together):
#' \enumerate{
#'   \item \strong{Spectral slope} (always computed): the component's
#'     log-log power-spectrum slope from \code{l_freq} to \code{h_freq}.
#'     Muscle activity is broadband and roughly flat/rising at high
#'     frequency (slope near \code{+0.15}); genuine neural/artifact
#'     components typically fall off. Logistic centered at \code{-0.5}.
#'   \item \strong{Peripherality} (needs electrode positions): how far the
#'     component's topography is weighted toward the edge of the montage
#'     (jaw/neck/temple), vs. the vertex. Logistic centered at 65\% of the
#'     maximum electrode distance from the fitted channels' own centroid.
#'   \item \strong{Spatial smoothness} (needs electrode positions): how
#'     focal/spiky (vs. smooth) the topography is. Logistic, inverted,
#'     centered at a smoothness value of 300.
#' }
#' If fewer than 3 of the fitted channels (\code{ica$ch_names}) have a
#' montage position, only criterion 1 is used and a warning is issued -
#' \code{channel_types == "external"} channels never have montage positions
#' in this package by design (see \code{\link{create_montage}}), so an
#' MNE-strict "every channel needs a position" check would trigger this
#' fallback far more often than useful; instead, criteria 2-3 use whichever
#' subset of \code{ica$ch_names} does have one. The combined threshold is
#' raised to the power of however many criteria were actually used, so the
#' default \code{0.5} behaves consistently whether 1 or 3 are active.
#'
#' \strong{Only writes to \code{ica$labels_}} - never to \code{ica$exclude}
#' (see file header); call \code{\link{set_exclude}} yourself to actually
#' remove anything.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{\link{fit_ica}}).
#' @param eeg The \code{eeg} object \code{ica} was fit on (or one with
#'   matching channels - see \code{\link{get_sources}}).
#' @param threshold Numeric. Combined-score cutoff, raised to the power of
#'   the number of criteria used. Default: \code{0.5}.
#' @param l_freq,h_freq Numeric. Frequency band in Hz for the spectral-slope
#'   criterion. Default: \code{c(7, 45)}.
#' @param montage An object of class \code{'montage'}. If \code{NULL}
#'   (default), \code{eeg$montage} is used (same convention as
#'   \code{\link{plot_ica_topography}}).
#'
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with
#'   \code{labels_[["muscle"]]} populated. Since R does not mutate arguments
#'   in place, the caller must reassign the result
#'   (\code{ica <- find_bads_muscle(ica, eeg)}).
#'
#' @details
#' Two deliberate, disclosed divergences from MNE, in the spirit of the
#' existing \code{.compute_pre_whitener()} \code{noise_cov} precedent ("match
#' exactly or knowingly diverge and document why"):
#' \enumerate{
#'   \item The spectral-slope criterion uses a single whole-segment
#'     periodogram (\code{\link{.fft_one_sided}}, R/fourier.R), not MNE's
#'     Welch-averaged \code{compute_psd()} - reasonable since only a robust
#'     log-log \emph{slope} is needed, not a precisely calibrated spectrum
#'     (a constant scale factor shifts \code{log10(power)} but never changes
#'     the fitted slope).
#'   \item Criteria 2-3 use raw 3D \code{(x, y, z)} montage positions, not
#'     MNE's 2D azimuthal-equidistant topomap projection - captures the same
#'     peripherality/smoothness concept without porting MNE's 2D projection
#'     algorithm.
#' }
#'
#' @examples
#' \dontrun{
#'   eeg <- set_montage(eeg, create_montage())
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- find_bads_muscle(ica, eeg)
#'   ica$labels_$muscle
#' }
#'
#' @seealso \code{\link{find_bads_eog}}, \code{\link{find_bads_ecg}},
#'   \code{\link{set_exclude}}
#'
#' @export
find_bads_muscle <- function(ica, eeg, threshold = 0.5, l_freq = 7, h_freq = 45,
                              montage = NULL) {

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).",
         call. = FALSE)
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.", call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  sfreq   <- eeg$sampling_rate
  sources <- get_sources(ica, eeg)
  n_comp  <- ica$n_components_

  # ========== CRITERION 1: LOG-LOG PSD SLOPE (always computed) ==========

  n_fft <- ncol(sources)
  slopes <- vapply(seq_len(n_comp), function(k) {
    psd  <- .fft_one_sided(sources[k, ], sfreq, n_fft)
    band <- which(psd$freqs >= l_freq & psd$freqs <= h_freq)
    if (length(band) < 2) {
      stop("ERROR: fewer than 2 frequency bins fall within [", l_freq, ", ",
           h_freq, "] Hz - check l_freq/h_freq against eeg$sampling_rate.",
           call. = FALSE)
    }
    fit <- lm(log10(psd$power[band]) ~ log10(psd$freqs[band]))
    unname(stats::coef(fit)[2])
  }, numeric(1))

  # Typical muscle slope ~ +0.15, non-muscle negative: logistic shift -0.5,
  # slope 0.25, so -0.5 -> 0.5 and 0 -> 1 (matches ica.py:2013-2015 exactly).
  slope_score <- 1 / (1 + exp(-((slopes + 0.5) / 0.25)))

  # ========== CRITERIA 2-3: NEED ELECTRODE POSITIONS ==========

  if (is.null(montage)) montage <- eeg$montage

  pos_idx <- match(ica$ch_names, montage$positions$channel)
  has_pos <- !is.na(pos_idx)

  if (sum(has_pos) < 3) {
    warning("No (or too few, <3) sensor positions found for the fitted ",
            "channels - scores for bad muscle components are only based ",
            "on the 'slope' criterion.", call. = FALSE)
    ica$labels_[["muscle"]] <- which(slope_score > threshold)
    return(ica)
  }

  components <- vapply(seq_len(n_comp), function(k) {
    get_component_topography(ica, k)
  }, numeric(length(ica$ch_names)))

  pos <- as.matrix(montage$positions[pos_idx[has_pos], c("x", "y", "z")])
  pos <- sweep(pos, 2, colMeans(pos), FUN = "-")  # center on the fitted subset's own mean

  comp_sub  <- components[has_pos, , drop = FALSE]
  comp_norm <- sweep(abs(comp_sub), 2, apply(abs(comp_sub), 2, max), FUN = "/")

  # Metric #2: distance from the centroid, weighted by each component's
  # (normalized) loading at that electrode (matches ica.py:2033-2050).
  dists <- sqrt(rowSums(pos^2))
  dists <- dists / max(dists)
  focus_dists <- as.vector(dists %*% comp_norm)
  focus_score <- 1 / (1 + exp(-((focus_dists - 0.65) / 0.1)))

  # Metric #3: spatial smoothness (matches ica.py:2052-2063).
  geo_dist <- as.matrix(dist(pos))
  geo_dist <- 1 - (geo_dist / max(geo_dist))

  smoothnesses <- vapply(seq_len(n_comp), function(k) {
    amp_dist <- as.matrix(dist(comp_sub[, k]))
    amp_dist <- amp_dist / max(amp_dist)
    sum(geo_dist * amp_dist)
  }, numeric(1))

  smoothness_score <- 1 - 1 / (1 + exp(-((smoothnesses - 300) / 100)))

  # ========== COMBINE (matches ica.py:2065-2077) ==========

  scores <- slope_score * focus_score * smoothness_score
  ica$labels_[["muscle"]] <- which(scores > threshold^3)

  ica
}
