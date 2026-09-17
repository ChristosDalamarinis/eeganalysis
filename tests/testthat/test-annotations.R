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
