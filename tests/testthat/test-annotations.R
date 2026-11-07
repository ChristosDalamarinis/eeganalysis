# ============================================================================
#                       Test File for annotations.R
# ============================================================================
#
# Functions tested:
#   1. annotate_amplitude()       - Exported: flat/peak detection
#   2. annotate_muscle()          - Exported: EMG envelope z-score detection
#   3. annotate_nan()             - Exported: per-channel NA runs
#   4. annotate_break()           - Exported: dead time between blocks
#   5. Shared internal helpers    - mask/interval/z-score/Hilbert utilities
#
# Test suites:
#   1. Shared internal helpers
#   2. annotate_amplitude()
#   3. annotate_muscle()
#   4. annotate_nan()
#   5. annotate_break()
#
# Author: Christos Dalamarinis
# Date: Sep 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

.chans <- c("Cz", "Pz", "Oz", "Fz")

# Plain multi-channel eeg object, IID noise on every channel - no artifacts.
.make_plain_eeg <- function(n_tp = 2000, sfreq = 256, sd = 10, seed = 1) {
  set.seed(seed)
  data <- matrix(rnorm(length(.chans) * n_tp, sd = sd), nrow = length(.chans))
  new_eeg(data = data, channels = .chans, sampling_rate = sfreq)
}


# ============================================================================
# TEST SUITE 1: Shared internal helpers
# ============================================================================

test_that(".mask_to_onsets_offsets finds runs and round-trips", {
  m <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  oo <- eeganalysis:::.mask_to_onsets_offsets(m)
  expect_equal(oo$onsets, c(3L, 8L))
  expect_equal(oo$offsets, c(6L, 9L))

  rebuilt <- rep(FALSE, length(m))
  for (i in seq_along(oo$onsets)) {
    rebuilt[oo$onsets[i]:(oo$offsets[i] - 1)] <- TRUE
  }
  expect_equal(rebuilt, m)
})

test_that(".mask_to_onsets_offsets handles all-TRUE, all-FALSE, and empty masks", {
  oo_true  <- eeganalysis:::.mask_to_onsets_offsets(rep(TRUE, 5))
  expect_equal(oo_true$onsets, 1L)
  expect_equal(oo_true$offsets, 6L)

  oo_false <- eeganalysis:::.mask_to_onsets_offsets(rep(FALSE, 5))
  expect_length(oo_false$onsets, 0)
  expect_length(oo_false$offsets, 0)

  oo_empty <- eeganalysis:::.mask_to_onsets_offsets(logical(0))
  expect_length(oo_empty$onsets, 0)
})

test_that(".mask_to_annotation_rows converts sample runs to seconds", {
  m <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  rows <- eeganalysis:::.mask_to_annotation_rows(m, sfreq = 10, description = "TEST")
  expect_equal(nrow(rows), 2)
  expect_equal(rows$onset, c(0.2, 0.7))
  expect_equal(rows$duration, c(0.3, 0.1))
  expect_true(all(rows$description == "TEST"))
  expect_true(all(is.na(rows$channel)))

  empty <- eeganalysis:::.mask_to_annotation_rows(rep(FALSE, 5), sfreq = 10,
                                                   description = "TEST")
  expect_equal(nrow(empty), 0)
  expect_named(empty, c("onset", "duration", "description", "channel"))
})

test_that(".flip_short_runs flips only runs shorter than min_samples", {
  # short TRUE runs -> FALSE (annotate_amplitude's "reject short" direction)
  mask <- c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  flipped <- eeganalysis:::.flip_short_runs(mask, min_samples = 3,
                                             run_value = TRUE, flip_to = FALSE)
  expect_equal(flipped, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))

  # short FALSE runs -> TRUE (annotate_muscle's "merge good gaps" direction)
  mask2 <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)
  flipped2 <- eeganalysis:::.flip_short_runs(mask2, min_samples = 2,
                                              run_value = FALSE, flip_to = TRUE)
  expect_equal(flipped2, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))

  # min_samples <= 0 is a no-op
  expect_equal(eeganalysis:::.flip_short_runs(mask, 0, TRUE, FALSE), mask)
})

test_that(".hilbert_envelope recovers constant amplitude of a pure sinusoid", {
  sfreq <- 1000
  t <- seq(0, 1 - 1 / sfreq, by = 1 / sfreq)
  amp <- 3.7
  x <- amp * sin(2 * pi * 50 * t)
  env <- eeganalysis:::.hilbert_envelope(x)

  # ignore edges, where the FFT-based analytic signal has boundary artifacts
  mid <- env[100:900]
  expect_equal(mean(mid), amp, tolerance = 0.01)
  expect_lt(sd(mid), 0.05)
})

test_that(".zscore_rows z-scores each row with population sd", {
  set.seed(1)
  mat <- rbind(rnorm(1000, mean = 5, sd = 2), rnorm(1000, mean = -3, sd = 10))
  z <- eeganalysis:::.zscore_rows(mat)

  popsd <- function(v) sqrt(mean((v - mean(v))^2))
  expect_equal(rowMeans(z), c(0, 0), tolerance = 1e-8)
  expect_equal(apply(z, 1, popsd), c(1, 1), tolerance = 1e-8)
})

test_that(".merge_intervals merges overlapping/touching intervals only", {
  iv <- rbind(c(0, 5), c(4, 8), c(20, 25), c(24, 24), c(30, 31))
  merged <- eeganalysis:::.merge_intervals(iv)
  expect_equal(unname(merged),
               matrix(c(0, 20, 30, 8, 25, 31), ncol = 2))
  expect_null(dimnames(merged)[[1]])
})

test_that(".resolve_channel_idx handles NULL/character/numeric/invalid input", {
  eeg <- .make_plain_eeg()
  expect_equal(eeganalysis:::.resolve_channel_idx(eeg, NULL, c(2L, 3L)), c(2L, 3L))
  expect_equal(eeganalysis:::.resolve_channel_idx(eeg, c("Pz", "Oz")), c(2L, 3L))
  expect_equal(eeganalysis:::.resolve_channel_idx(eeg, c(1L, 4L)), c(1L, 4L))
  expect_error(eeganalysis:::.resolve_channel_idx(eeg, "NotAChannel"), "not found")
  expect_error(eeganalysis:::.resolve_channel_idx(eeg, 99L), "out of range")
  expect_error(eeganalysis:::.resolve_channel_idx(eeg, TRUE), "character vector, or an integer")
})


# ============================================================================
# TEST SUITE 2: annotate_amplitude()
# ============================================================================

test_that("annotate_amplitude validates its inputs", {
  eeg <- .make_plain_eeg()
  expect_error(annotate_amplitude(list(), peak = 500), "class 'eeg'")
  expect_error(annotate_amplitude(eeg), "peak' or 'flat'")
  expect_error(annotate_amplitude(eeg, peak = -1), "non-negative")
  expect_error(annotate_amplitude(eeg, flat = -1), "non-negative")
  expect_error(annotate_amplitude(eeg, peak = 500, bad_percent = 200), "\\[0, 100\\]")
  expect_error(annotate_amplitude(eeg, peak = 500, min_duration = -1), "min_duration")
})

test_that("annotate_amplitude flags a short flat stretch as BAD_flat, not bads", {
  set.seed(1)
  n_tp <- 2000
  data <- matrix(rnorm(length(.chans) * n_tp, sd = 10), nrow = length(.chans))
  data[2, 900:939] <- 5  # ~40-sample flat run on Pz, well under bad_percent
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = 256)

  result <- annotate_amplitude(eeg, flat = 0.5, min_duration = 0.02,
                                return_details = TRUE)

  expect_equal(nrow(result$annotations), 1)
  expect_equal(result$annotations$description, "BAD_flat")
  expect_equal(result$annotations$onset, 899 / 256, tolerance = 1e-6)
  expect_length(result$eeg_obj$bads, 0)
  expect_true(any(grepl("annotate_amplitude",
                        unlist(result$eeg_obj$preprocessing_history))))
})

test_that("annotate_amplitude routes a mostly-flat channel to bads instead", {
  set.seed(2)
  n_tp <- 500
  data <- matrix(rnorm(length(.chans) * n_tp, sd = 10), nrow = length(.chans))
  data[2, 100:199] <- 5  # 100/500 samples flat -> ~20%, over default bad_percent
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = 256)

  result <- annotate_amplitude(eeg, flat = 0.5, min_duration = 0.02,
                                return_details = TRUE)

  expect_equal(result$eeg_obj$bads, "Pz")
  expect_equal(nrow(result$annotations), 0)
})

test_that("annotate_amplitude flags a sharp spike stretch as BAD_peak", {
  set.seed(3)
  n_tp <- 2000
  data <- matrix(rnorm(length(.chans) * n_tp, sd = 10), nrow = length(.chans))
  data[1, 500:539] <- rep(c(1000, -1000), 20)  # sawtooth -> huge diffs
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = 256)

  result <- annotate_amplitude(eeg, peak = 500, return_details = TRUE)

  expect_equal(nrow(result$annotations), 1)
  expect_equal(result$annotations$description, "BAD_peak")
  expect_true(result$annotations$onset >= 499 / 256 - 0.01)
  expect_true(result$annotations$onset <= 500 / 256 + 0.01)
})

test_that("annotate_amplitude initializes eeg_obj$annotations when absent", {
  eeg <- .make_plain_eeg()
  eeg$annotations <- NULL  # simulate an object predating the annotations field
  result <- annotate_amplitude(eeg, peak = 1e6)  # threshold nothing can hit
  expect_s3_class(result$annotations, "data.frame")
  expect_named(result$annotations, c("onset", "duration", "description", "channel"))
  expect_equal(nrow(result$annotations), 0)
})

test_that("new_eeg initializes annotations to an empty frame by default", {
  eeg <- .make_plain_eeg()
  expect_s3_class(eeg$annotations, "data.frame")
  expect_named(eeg$annotations, c("onset", "duration", "description", "channel"))
  expect_equal(nrow(eeg$annotations), 0)
})


# ============================================================================
# TEST SUITE 3: annotate_muscle()
# ============================================================================

test_that("annotate_muscle validates its inputs", {
  eeg <- .make_plain_eeg(sfreq = 1024)
  expect_error(annotate_muscle(list()), "class 'eeg'")
  expect_error(annotate_muscle(eeg, filter_freq = c(140, 110)), "filter_freq")
  expect_error(annotate_muscle(eeg, filter_freq = c(-10, 20)), "filter_freq")
  expect_error(annotate_muscle(eeg, min_length_good = -1), "min_length_good")

  eeg_slow <- .make_plain_eeg(sfreq = 256)
  expect_error(annotate_muscle(eeg_slow), "Nyquist")
})

test_that("annotate_muscle errors on NA in the picked channels rather than corrupting output", {
  eeg <- .make_plain_eeg(sfreq = 1024)
  eeg$data[1, 10] <- NA
  expect_error(annotate_muscle(eeg), "NA values")
})

test_that("annotate_muscle detects an injected 110-140 Hz burst and stays quiet without one", {
  sfreq <- 1024
  dur_s <- 30
  n_tp  <- sfreq * dur_s
  t     <- (0:(n_tp - 1)) / sfreq
  burst_idx <- which(t >= 14 & t < 15)

  set.seed(42)
  data <- t(vapply(.chans, function(ch) {
    x <- rnorm(n_tp, sd = 5)
    x[burst_idx] <- x[burst_idx] + 60 * sin(2 * pi * 125 * t[burst_idx])
    x
  }, numeric(n_tp)))
  eeg_burst <- new_eeg(data = data, channels = .chans, sampling_rate = sfreq)

  result <- annotate_muscle(eeg_burst, return_details = TRUE)

  expect_gte(nrow(result$annotations), 1)
  expect_true(all(result$annotations$description == "BAD_muscle"))
  span_start <- min(result$annotations$onset)
  span_end   <- max(result$annotations$onset + result$annotations$duration)
  expect_true(span_start >= 13.5 && span_start <= 14.5)
  expect_true(span_end   >= 14.5 && span_end   <= 15.5)
  expect_length(result$scores, n_tp)

  # No burst at all: should stay well under threshold, no annotations.
  set.seed(43)
  data_clean <- t(vapply(.chans, function(ch) rnorm(n_tp, sd = 5), numeric(n_tp)))
  eeg_clean  <- new_eeg(data = data_clean, channels = .chans, sampling_rate = sfreq)
  result_clean <- annotate_muscle(eeg_clean, return_details = TRUE)

  expect_equal(nrow(result_clean$annotations), 0)
  expect_lt(max(result_clean$scores), 4)
})

test_that("annotate_muscle's score is smooth: the 110-140 Hz strength is low-passed at 4 Hz", {
  sfreq <- 1024
  n_tp  <- sfreq * 30

  set.seed(44)
  data <- t(vapply(.chans, function(ch) rnorm(n_tp, sd = 5), numeric(n_tp)))
  eeg  <- new_eeg(data = data, channels = .chans, sampling_rate = sfreq)

  scores <- annotate_muscle(eeg, return_details = TRUE)$scores

  # A 4 Hz low-pass leaves almost nothing above 8 Hz (its stop band starts at
  # 6 Hz). The strength of noise in a 30 Hz wide band has a lot between 8 and
  # 30 Hz, so a weaker smoothing would leave a large share of the variance there.
  fast <- filter_data(scores, sfreq, l_freq = 8, h_freq = 100)
  expect_lt(var(fast) / var(scores), 1e-3)
})

# --- annotate_muscle() and notes that are already in eeg$annotations ---------
#
# Fixture: 30 s at 1024 Hz, 4 channels of noise (sd 5), plus a 125 Hz burst on
# every channel: burst A at 8-9 s (amplitude 200, strong) and burst B at 22-23 s
# (amplitude 40, weak). With nothing marked, A makes the "usual loudness" so
# large that B stays under the threshold (top score about 1.9, threshold 4) and
# only A is flagged. With A already marked as bad, A is left out of the
# z-scoring and B is found (top score about 11). Every note below starts and
# ends on a whole sample (1024 Hz: multiples of 1/1024 s), so the stretch that
# is left out can be counted by hand: 1 s = 1024 samples.

.make_muscle_eeg <- function(seed = 61, amp_a = 200, amp_b = 40) {
  sfreq <- 1024
  n_tp  <- sfreq * 30
  t     <- (0:(n_tp - 1)) / sfreq
  set.seed(seed)
  data <- t(vapply(.chans, function(ch) {
    x <- rnorm(n_tp, sd = 5)
    a <- which(t >= 8 & t < 9)
    b <- which(t >= 22 & t < 23)
    x[a] <- x[a] + amp_a * sin(2 * pi * 125 * t[a])
    x[b] <- x[b] + amp_b * sin(2 * pi * 125 * t[b])
    x
  }, numeric(n_tp)))
  list(eeg = new_eeg(data = data, channels = .chans, sampling_rate = sfreq), t = t)
}

.note <- function(onset, duration, description = "BAD_manual", channel = NA_character_) {
  data.frame(onset = onset, duration = duration, description = description,
             channel = channel, stringsAsFactors = FALSE)
}

test_that("annotate_muscle with no earlier annotation gives the same trace as before (pinned numbers)", {
  # The numbers were taken from the code before notes were taken into account.
  f <- .make_muscle_eeg()
  r <- annotate_muscle(f$eeg, return_details = TRUE)
  expect_equal(max(r$scores[f$t >= 8 & f$t < 9]),   11.5612517067, tolerance = 1e-8)
  expect_equal(max(r$scores[f$t >= 22 & f$t < 23]),  1.8704928273, tolerance = 1e-8)
  expect_equal(sum(r$scores), -2.5475917103, tolerance = 1e-8)
  expect_equal(sd(r$scores),   1.9845882269, tolerance = 1e-8)
  expect_false(anyNA(r$scores))
  expect_equal(nrow(r$annotations), 1)          # only the strong burst A
  expect_gt(r$annotations$onset, 7.9)
  expect_lt(r$annotations$onset + r$annotations$duration, 9.1)
})

test_that("annotate_muscle does not flag a stretch marked BAD_ by hand again, and then finds the weaker burst", {
  f   <- .make_muscle_eeg()
  eeg <- f$eeg
  eeg$annotations <- .note(8, 1)                 # burst A, marked by hand
  r <- annotate_muscle(eeg, return_details = TRUE)

  # no new row touches 8-9 s; exactly one row covers burst B (22-23 s)
  expect_equal(nrow(r$annotations), 1)
  expect_true(all(r$annotations$description == "BAD_muscle"))
  expect_gt(r$annotations$onset, 21.5)
  expect_lt(r$annotations$onset, 22.1)
  expect_gt(r$annotations$onset + r$annotations$duration, 22.9)
  expect_lt(r$annotations$onset + r$annotations$duration, 23.5)
  expect_gt(max(r$scores[f$t >= 22 & f$t < 23]), 4)    # was about 1.9 when A counted
  # the table now holds the hand-made note and the new row, in time order
  expect_equal(r$eeg_obj$annotations$description, c("BAD_manual", "BAD_muscle"))
})

test_that("annotate_muscle's scores are NA exactly inside the marked stretch and numbers everywhere else", {
  f   <- .make_muscle_eeg()
  eeg <- f$eeg
  eeg$annotations <- .note(8, 1)
  scores <- annotate_muscle(eeg, return_details = TRUE)$scores

  expect_type(scores, "double")
  expect_length(scores, ncol(eeg$data))
  expect_equal(which(is.na(scores)), which(f$t >= 8 & f$t < 9))   # samples 8193 to 9216
  expect_equal(sum(is.na(scores)), 1024)                           # 1 s at 1024 Hz
})

test_that("annotate_muscle's scores equal a step-by-step recomputation that leaves the marked stretch out", {
  f   <- .make_muscle_eeg()
  eeg <- f$eeg
  eeg$annotations <- .note(8, 1)
  scores <- annotate_muscle(eeg, return_details = TRUE)$scores

  # 1. band-pass + envelope on EVERY sample; 2. drop the marked samples;
  # 3. z-score each channel over the kept samples; 4. add up, divide by
  # sqrt(4 channels); 5. smooth the glued kept samples at 4 Hz; 6. put back.
  env <- t(vapply(seq_len(nrow(eeg$data)), function(i) {
    .hilbert_envelope(filter_data(eeg$data[i, ], 1024, l_freq = 110, h_freq = 140))
  }, numeric(ncol(eeg$data))))
  keep <- !(f$t >= 8 & f$t < 9)
  z    <- .zscore_rows(env[, keep])
  expected <- rep(NA_real_, ncol(eeg$data))
  expected[keep] <- filter_data(colSums(z) / sqrt(4), 1024, h_freq = 4)

  expect_equal(scores, expected, tolerance = 1e-12)
})

test_that("annotate_muscle leaves out only notes whose description starts with 'bad' (any letter case)", {
  f <- .make_muscle_eeg()
  plain <- annotate_muscle(f$eeg, return_details = TRUE)$scores

  for (d in c("BAD_manual", "bad_x", "Bad Segment", "BAD_peak", "badminton")) {
    eeg <- f$eeg
    eeg$annotations <- .note(8, 1, d)
    scores <- annotate_muscle(eeg, return_details = TRUE)$scores
    expect_equal(sum(is.na(scores)), 1024, info = d)
  }
  # other descriptions are just labels: nothing is left out, trace unchanged
  for (d in c("Block1", "xBAD", "Not_bad", "Baseline", "good")) {
    eeg <- f$eeg
    eeg$annotations <- .note(8, 1, d)
    scores <- annotate_muscle(eeg, return_details = TRUE)$scores
    expect_false(anyNA(scores), info = d)
    expect_equal(scores, plain, info = d)
  }
})

test_that("annotate_muscle ignores the channel column of an earlier note", {
  f <- .make_muscle_eeg()
  eeg_all <- f$eeg
  eeg_all$annotations <- .note(8, 1)
  eeg_cz <- f$eeg
  eeg_cz$annotations <- .note(8, 1, "BAD_manual", "Cz")      # marked for one channel only
  s_all <- annotate_muscle(eeg_all, return_details = TRUE)$scores
  s_cz  <- annotate_muscle(eeg_cz,  return_details = TRUE)$scores
  expect_equal(which(is.na(s_cz)), which(f$t >= 8 & f$t < 9))
  expect_equal(s_cz, s_all)
})

test_that("annotate_muscle called on its own result does not flag a marked burst again", {
  f <- .make_muscle_eeg(amp_b = 0)               # burst A only
  first <- annotate_muscle(f$eeg)
  expect_equal(nrow(first$annotations), 1)
  r1 <- first$annotations[1, ]
  second <- annotate_muscle(first, return_details = TRUE)$annotations

  # no new row overlaps the first one ...
  expect_false(any(second$onset < r1$onset + r1$duration &
                   r1$onset < second$onset + second$duration))
  # ... the only new rows are two short stretches that touch its two edges: the
  # fringe of the burst, which the first call left unflagged and which counts
  # as loud once burst A is out of the z-scoring
  expect_equal(nrow(second), 2)
  expect_true(all(second$duration < 0.05))
  expect_equal(second$onset[1] + second$duration[1], r1$onset, tolerance = 1e-9)
  expect_equal(second$onset[2], r1$onset + r1$duration, tolerance = 1e-9)
})

test_that("annotate_muscle on its own result finds the weaker burst that the first call could not see", {
  g  <- .make_muscle_eeg()
  g1 <- annotate_muscle(g$eeg)                   # flags A only, B is hidden by A
  expect_equal(nrow(g1$annotations), 1)
  a_row <- g1$annotations[1, ]
  g2 <- annotate_muscle(g1, return_details = TRUE)$annotations

  # burst B (22-23 s) is now covered by one new row ...
  expect_true(any(g2$onset > 21.5 & g2$onset < 22.1 & g2$onset + g2$duration > 22.9))
  # ... and nothing overlaps burst A's row again
  expect_false(any(g2$onset < a_row$onset + a_row$duration &
                   a_row$onset < g2$onset + g2$duration))
})

test_that("annotate_muscle errors when fewer than 2 samples are left outside earlier BAD notes", {
  f <- .make_muscle_eeg()
  eeg <- f$eeg
  eeg$annotations <- .note(0, 30, "BAD_all")                       # every sample
  expect_error(annotate_muscle(eeg), "fewer than 2 samples")
  eeg$annotations <- .note(0, 30719 / 1024, "BAD_all")             # all but the last one
  expect_error(annotate_muscle(eeg), "fewer than 2 samples")
})

test_that("annotate_muscle's history line names the seconds left out only when some were", {
  f <- .make_muscle_eeg()
  eeg <- f$eeg
  eeg$annotations <- .note(8, 1)
  h <- annotate_muscle(eeg)$preprocessing_history
  expect_match(h[[length(h)]], "1 s inside earlier BAD annotation(s) left out", fixed = TRUE)

  h0 <- annotate_muscle(f$eeg)$preprocessing_history
  expect_false(grepl("left out", h0[[length(h0)]], fixed = TRUE))
})

test_that("annotate_muscle counts a short left-out stretch between two flagged stretches as good data", {
  # One continuous burst at 8-10 s with a 0.05 s (51 samples) stretch marked in
  # its middle. min_length_good = 0.1 s (102 samples): the marked stretch is a
  # shorter "good" run, so it is folded into the flagged stretch -> ONE row.
  # min_length_good = 0.01 s (10 samples): it stays a gap -> TWO rows.
  sfreq <- 1024; n_tp <- sfreq * 30; t <- (0:(n_tp - 1)) / sfreq
  set.seed(62)
  data <- t(vapply(.chans, function(ch) {
    x <- rnorm(n_tp, sd = 5)
    a <- which(t >= 8 & t < 10)
    x[a] <- x[a] + 100 * sin(2 * pi * 125 * t[a])
    x
  }, numeric(n_tp)))
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = sfreq)
  eeg$annotations <- .note(9, 51 / 1024, "BAD_gap")

  one <- annotate_muscle(eeg, min_length_good = 0.1)$annotations
  two <- annotate_muscle(eeg, min_length_good = 0.01)$annotations
  expect_equal(sum(one$description == "BAD_muscle"), 1)
  expect_equal(sum(two$description == "BAD_muscle"), 2)
  # in the two-row case the gap between the new rows is exactly the marked stretch
  m <- two[two$description == "BAD_muscle", ]
  expect_equal(m$onset[2], 9 + 51 / 1024, tolerance = 1e-9)
  expect_equal(m$onset[1] + m$duration[1], 9, tolerance = 1e-9)
})


# ============================================================================
# TEST SUITE 4: annotate_nan()
# ============================================================================

test_that("annotate_nan finds per-channel NA runs and leaves clean data alone", {
  set.seed(4)
  n_tp <- 500
  data <- matrix(rnorm(length(.chans) * n_tp), nrow = length(.chans))
  data[1, 100:149] <- NA  # Cz: 50-sample dropout
  data[3, 300:302] <- NA  # Oz: 3-sample dropout
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = 256)

  result <- annotate_nan(eeg, return_details = TRUE)

  expect_equal(nrow(result$annotations), 2)
  expect_true(all(result$annotations$description == "BAD_NAN"))
  expect_setequal(result$annotations$channel, c("Cz", "Oz"))

  cz_row <- result$annotations[result$annotations$channel == "Cz", ]
  expect_equal(cz_row$onset, 99 / 256, tolerance = 1e-6)
  expect_equal(cz_row$duration, 50 / 256, tolerance = 1e-6)

  clean_eeg <- .make_plain_eeg()
  clean_result <- annotate_nan(clean_eeg, return_details = TRUE)
  expect_equal(nrow(clean_result$annotations), 0)
})

test_that("annotate_nan respects an explicit channels argument", {
  set.seed(5)
  n_tp <- 300
  data <- matrix(rnorm(length(.chans) * n_tp), nrow = length(.chans))
  data[1, 10:19] <- NA
  data[2, 10:19] <- NA
  eeg <- new_eeg(data = data, channels = .chans, sampling_rate = 256)

  result <- annotate_nan(eeg, channels = "Pz", return_details = TRUE)
  expect_equal(result$annotations$channel, "Pz")
})


# ============================================================================
# TEST SUITE 5: annotate_break()
# ============================================================================

.make_break_eeg <- function(dur_s = 200, sfreq = 256) {
  n_tp <- dur_s * sfreq
  data <- matrix(rnorm(length(.chans) * n_tp), nrow = length(.chans))
  new_eeg(data = data, channels = .chans, sampling_rate = sfreq)
}

test_that("annotate_break validates duration parameters and requires source rows", {
  eeg <- .make_break_eeg()
  expect_error(
    annotate_break(eeg, min_break_duration = 5,
                    t_start_after_previous = 3, t_stop_before_next = 3),
    "must be positive")
  expect_error(annotate_break(eeg), "no usable rows")
  expect_error(annotate_break(eeg, use_events = TRUE), "eeg_obj\\$events")
})

test_that("annotate_break finds the gap between two occupied stretches", {
  eeg <- .make_break_eeg(dur_s = 200)
  eeg$annotations <- data.frame(
    onset = c(5, 100), duration = c(2, 3),
    description = c("stimulus", "stimulus"), channel = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- annotate_break(eeg, return_details = TRUE)

  expect_equal(nrow(result$annotations), 2)
  expect_equal(result$annotations$onset, c(12, 108))
  expect_equal(result$annotations$duration, c(83, 91.99609375), tolerance = 1e-4)
  expect_true(all(result$annotations$description == "BAD_break"))
})

test_that("annotate_break's default ignore excludes BAD_-prefixed rows", {
  eeg <- .make_break_eeg(dur_s = 200)
  eeg$annotations <- data.frame(
    onset = c(5, 100), duration = c(2, 3),
    description = c("BAD_flat", "stimulus"), channel = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- annotate_break(eeg, return_details = TRUE)

  # BAD_flat ignored -> the only occupied stretch is [100, 103], so the
  # "before first" gap now starts at 0, not at the BAD_flat-bounded 12.
  expect_equal(result$annotations$onset, c(0, 108))

  # Passing an empty ignore list brings BAD_flat back into consideration.
  result_all <- annotate_break(eeg, ignore = character(0), return_details = TRUE)
  expect_equal(result_all$annotations$onset, c(12, 108))
})

test_that("annotate_break can use eeg_obj$events instead of annotations", {
  eeg <- .make_break_eeg(dur_s = 200)
  eeg$events <- data.frame(
    onset = c(1, 2), onset_time = c(5, 150),
    type = c(1, 1), description = c("stim", "stim"),
    stringsAsFactors = FALSE
  )

  result <- annotate_break(eeg, use_events = TRUE, return_details = TRUE)

  expect_equal(nrow(result$annotations), 2)
  expect_equal(result$annotations$onset, c(10, 155))
  expect_equal(result$annotations$duration, c(135, 44.99609375), tolerance = 1e-4)
})

test_that("annotate_break merges overlapping occupied stretches before finding gaps", {
  eeg <- .make_break_eeg(dur_s = 100)
  eeg$annotations <- data.frame(
    onset = c(1, 3), duration = c(4, 4),         # [1,5] and [3,7] overlap -> [1,7]
    description = "stimulus", channel = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- annotate_break(eeg, min_break_duration = 15,
                            t_start_after_previous = 5, t_stop_before_next = 5,
                            return_details = TRUE)

  # Only one merged occupied stretch [1,7] -> only the trailing gap qualifies.
  expect_equal(nrow(result$annotations), 1)
  expect_equal(result$annotations$onset, 12)
})

test_that("annotate_break rows are appended and sorted into eeg_obj$annotations", {
  eeg <- .make_break_eeg(dur_s = 200)
  eeg$annotations <- data.frame(
    onset = c(5, 100), duration = c(2, 3),
    description = c("stimulus", "stimulus"), channel = NA_character_,
    stringsAsFactors = FALSE
  )

  eeg2 <- annotate_break(eeg)
  expect_equal(nrow(eeg2$annotations), 4)
  expect_equal(eeg2$annotations$onset, sort(eeg2$annotations$onset))
})
