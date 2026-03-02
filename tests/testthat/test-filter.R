# ============================================================================
#                       Test File for filter.R
# ============================================================================
#
# This test file provides comprehensive testing for the filter.R script, which
# implements three EEG filtering functions for the eeganalysis package.
#
# Functions tested:
#   1. eeg_bandpass()  - Butterworth / FIR bandpass, highpass, and lowpass
#                        filtering with channel selection, zero-phase option,
#                        mirror-padding, and an optional in-line notch
#   2. eeg_notch()     - Band-stop filtering for line-noise removal, with
#                        optional harmonic suppression
#   3. eeg_filter()    - Unified wrapper that calls eeg_bandpass() then
#                        eeg_notch() in the correct order
#
# Test structure:
#   - A shared make_eeg() helper creates minimal but realistic eeg objects
#   - Signal-quality tests synthesise pure sine-wave signals and verify
#     that the filtered output has appropriate power relative to the input
#     (stopband power < 10 % of input power; passband power > 90 %)
#   - Each test block opens with a "WHAT THIS TESTS" comment that explains
#     exactly which line(s) of source code the test exercises
#
# Author: Christos Dalamarinis
# Date: March 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                         SHARED HELPER FUNCTIONS
# ============================================================================

#' Build a minimal mock eeg object.
#'
#' @param n_channels   Number of EEG channels (rows in data matrix).
#' @param n_timepoints Number of time samples (columns in data matrix).
#' @param sr           Sampling rate in Hz.
#' @param data_values  Optional numeric matrix; random noise used if NULL.
#' @return An S3 object of class "eeg".
make_eeg <- function(n_channels   = 3,
                     n_timepoints = 2560,
                     sr           = 512,
                     data_values  = NULL) {
  
  if (is.null(data_values)) {
    set.seed(42)
    data_values <- matrix(rnorm(n_channels * n_timepoints),
                          nrow = n_channels, ncol = n_timepoints)
  }
  channels <- paste0("Ch", seq_len(n_channels))
  times    <- seq(0, (n_timepoints - 1) / sr, length.out = n_timepoints)
  
  structure(
    list(
      data                  = data_values,
      channels              = channels,
      sampling_rate         = sr,
      times                 = times,
      events                = data.frame(onset        = integer(0),
                                         onset_time   = numeric(0),
                                         type         = character(0),
                                         description  = character(0),
                                         stringsAsFactors = FALSE),
      metadata              = list(),
      reference             = "original",
      preprocessing_history = list()
    ),
    class = "eeg"
  )
}

#' Build a pure sine-wave signal matrix.
#'
#' @param freq_hz    Frequency of the sine in Hz.
#' @param sr         Sampling rate in Hz.
#' @param n_samp     Number of samples.
#' @param n_channels Number of identical channel rows.
#' @return A numeric matrix (n_channels x n_samp).
sine_matrix <- function(freq_hz, sr = 512, n_samp = 2560, n_channels = 2) {
  t   <- seq(0, (n_samp - 1) / sr, length.out = n_samp)
  sig <- sin(2 * pi * freq_hz * t)
  matrix(rep(sig, n_channels), nrow = n_channels, byrow = TRUE)
}

#' Signal power (mean squared amplitude).
signal_power <- function(x) mean(x^2)


# ============================================================================
#           TEST SUITE 1: eeg_bandpass() - Input Validation
# ============================================================================
#
# These tests exercise the input-validation block that runs before any
# filtering. Each test checks one specific guard condition.

# ----------------------------------------------------------------------------
# Test 1.1: Rejects non-eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The very first guard in eeg_bandpass():
#   if (!inherits(eeg_obj, "eeg")) stop(...)
# Passing a plain list must raise an error containing "class 'eeg'".
test_that("eeg_bandpass() rejects a non-eeg object", {
  not_eeg <- list(data = matrix(1:10, 2, 5), channels = c("A", "B"))
  expect_error(
    eeg_bandpass(not_eeg, low_freq = 1, high_freq = 40, verbose = FALSE),
    regexp = "class 'eeg'"
  )
})

# ----------------------------------------------------------------------------
# Test 1.2: Rejects when both low_freq and high_freq are NULL
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (is.null(low_freq) && is.null(high_freq)) stop(...)
# At least one cutoff must be specified; omitting both must error.
test_that("eeg_bandpass() errors when both low_freq and high_freq are NULL", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = NULL, high_freq = NULL, verbose = FALSE),
    regexp = "At least one"
  )
})

# ----------------------------------------------------------------------------
# Test 1.3: Rejects low_freq at or above Nyquist
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (low_freq >= nyquist) stop(...)
# Nyquist for 512 Hz is 256 Hz; low_freq = 256 must error.
test_that("eeg_bandpass() errors when low_freq equals Nyquist", {
  eeg <- make_eeg(sr = 512)   # Nyquist = 256 Hz
  expect_error(
    eeg_bandpass(eeg, low_freq = 256, high_freq = NULL, verbose = FALSE),
    regexp = "Nyquist"
  )
})

# ----------------------------------------------------------------------------
# Test 1.4: Rejects high_freq at or above Nyquist
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (high_freq >= nyquist) stop(...)
test_that("eeg_bandpass() errors when high_freq equals Nyquist", {
  eeg <- make_eeg(sr = 512)
  expect_error(
    eeg_bandpass(eeg, low_freq = NULL, high_freq = 256, verbose = FALSE),
    regexp = "Nyquist"
  )
})

# ----------------------------------------------------------------------------
# Test 1.5: Rejects low_freq >= high_freq
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.null(low_freq) && !is.null(high_freq) && low_freq >= high_freq)
#     stop(...)
# A lowpass cutoff lower than or equal to the highpass cutoff is invalid.
test_that("eeg_bandpass() errors when low_freq >= high_freq", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 40, high_freq = 1, verbose = FALSE),
    regexp = "low_freq.*high_freq|high_freq.*low_freq"
  )
  expect_error(
    eeg_bandpass(eeg, low_freq = 20, high_freq = 20, verbose = FALSE),
    regexp = "low_freq.*high_freq|high_freq.*low_freq"
  )
})

# ----------------------------------------------------------------------------
# Test 1.6: Rejects non-positive low_freq
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.numeric(low_freq) || length(low_freq) != 1 || low_freq <= 0)
#     stop(...)
test_that("eeg_bandpass() errors when low_freq is zero or negative", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 0, high_freq = 40, verbose = FALSE),
    regexp = "positive"
  )
  expect_error(
    eeg_bandpass(eeg, low_freq = -1, high_freq = 40, verbose = FALSE),
    regexp = "positive"
  )
})

# ----------------------------------------------------------------------------
# Test 1.7: Rejects non-integer filter_order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.numeric(filter_order) || filter_order < 1 || filter_order %% 1 != 0)
#     stop(...)
test_that("eeg_bandpass() errors when filter_order is not a positive integer", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 filter_order = 0, verbose = FALSE),
    regexp = "positive integer"
  )
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 filter_order = 2.5, verbose = FALSE),
    regexp = "positive integer"
  )
})

# ----------------------------------------------------------------------------
# Test 1.8: Rejects non-logical zero_phase
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.logical(zero_phase) || length(zero_phase) != 1) stop(...)
test_that("eeg_bandpass() errors when zero_phase is not logical", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 zero_phase = "yes", verbose = FALSE),
    regexp = "zero_phase"
  )
})

# ----------------------------------------------------------------------------
# Test 1.9: Rejects non-logical verbose
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.logical(verbose) || length(verbose) != 1) stop(...)
test_that("eeg_bandpass() errors when verbose is not logical", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 verbose = 1),
    regexp = "verbose"
  )
})

# ----------------------------------------------------------------------------
# Test 1.10: Rejects data containing NA values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (any(is.na(eeg_obj$data))) stop(...)
test_that("eeg_bandpass() errors on NA values in data", {
  eeg           <- make_eeg()
  eeg$data[1, 5] <- NA
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE),
    regexp = "NA"
  )
})

# ----------------------------------------------------------------------------
# Test 1.11: Rejects data containing Inf values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (any(!is.finite(eeg_obj$data))) stop(...)
test_that("eeg_bandpass() errors on Inf values in data", {
  eeg           <- make_eeg()
  eeg$data[2, 10] <- Inf
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE),
    regexp = "non-finite|Inf"
  )
})

# ----------------------------------------------------------------------------
# Test 1.12: Rejects invalid method argument
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   method <- match.arg(method)
# Passing an unrecognised string must trigger match.arg's error.
test_that("eeg_bandpass() errors on unrecognised method string", {
  eeg <- make_eeg()
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 method = "chebyshev", verbose = FALSE),
    regexp = "arg"
  )
})

# ----------------------------------------------------------------------------
# Test 1.13: Rejects character channel names not present in the eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   bad <- channels[is.na(ch_idx)]; stop("Channel(s) not found in eeg_obj: ...")
test_that("eeg_bandpass() errors on unknown channel names", {
  eeg <- make_eeg(n_channels = 3)   # channels: Ch1, Ch2, Ch3
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 channels = c("Ch1", "BadChan"), verbose = FALSE),
    regexp = "BadChan"
  )
})

# ----------------------------------------------------------------------------
# Test 1.14: Rejects out-of-range integer channel indices
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (any(ch_idx < 1) || any(ch_idx > n_channels)) stop(...)
test_that("eeg_bandpass() errors on out-of-range integer channel index", {
  eeg <- make_eeg(n_channels = 3)
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 channels = c(1, 99), verbose = FALSE),
    regexp = "out of range"
  )
})

# ----------------------------------------------------------------------------
# Test 1.15: Warns when padding >= n_timepoints
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (padding > 0 && padding >= n_timepoints) warning(...)
# The function must warn but continue (not error) by disabling padding.
test_that("eeg_bandpass() warns when padding >= n_timepoints", {
  eeg <- make_eeg(n_timepoints = 100)
  expect_warning(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 padding = 200, verbose = FALSE),
    regexp = "padding"
  )
})

# ----------------------------------------------------------------------------
# Test 1.16: Rejects notch_freq at or above Nyquist (inline notch)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The notch_freq validation block inside eeg_bandpass():
#   if (notch_freq >= nyquist) stop(...)
test_that("eeg_bandpass() errors when inline notch_freq equals Nyquist", {
  eeg <- make_eeg(sr = 512)
  expect_error(
    eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                 notch_freq = 256, verbose = FALSE),
    regexp = "Nyquist|notch"
  )
})


# ============================================================================
#       TEST SUITE 2: eeg_bandpass() - Output Structure and Metadata
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Returns an object of class "eeg"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The final return value is built with new_eeg() and must
# carry the class attribute "eeg" for all downstream functions to work.
test_that("eeg_bandpass() returns an object of class 'eeg'", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expect_s3_class(result, "eeg")
})

# ----------------------------------------------------------------------------
# Test 2.2: Returned data matrix has same dimensions as input
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The filtered data is written back into a new_eeg() call
# using the original data matrix shape. Dimensions must not change.
test_that("eeg_bandpass() preserves data matrix dimensions", {
  eeg    <- make_eeg(n_channels = 4, n_timepoints = 1024)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expect_equal(dim(result$data), dim(eeg$data))
})

# ----------------------------------------------------------------------------
# Test 2.3: Returned object preserves all required eeg fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: new_eeg() is called with all original fields forwarded.
# The eight required fields must be present in the output.
test_that("eeg_bandpass() result contains all eight eeg fields", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expected_fields <- c("data", "channels", "sampling_rate", "times",
                       "events", "metadata", "reference", "preprocessing_history")
  expect_true(all(expected_fields %in% names(result)))
})

# ----------------------------------------------------------------------------
# Test 2.4: Sampling rate, channels, events, metadata, and reference are
#           passed through unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The new_eeg() call at the end of eeg_bandpass() forwards
# all non-data fields verbatim. None of them should be mutated.
test_that("eeg_bandpass() forwards non-data fields unchanged", {
  eeg              <- make_eeg()
  eeg$metadata     <- list(subject = "S01", session = 2L)
  eeg$reference    <- "average"
  result           <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                                   verbose = FALSE)
  expect_equal(result$sampling_rate, eeg$sampling_rate)
  expect_equal(result$channels,      eeg$channels)
  expect_equal(result$metadata,      eeg$metadata)
  expect_equal(result$reference,     eeg$reference)
  expect_equal(result$events,        eeg$events)
})

# ----------------------------------------------------------------------------
# Test 2.5: Appends exactly one entry to preprocessing_history
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: At the end of eeg_bandpass():
#   new_history <- c(eeg_obj$preprocessing_history, list(history_entry))
# Starting from an empty history, output must have exactly one entry.
test_that("eeg_bandpass() appends exactly one preprocessing history entry", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expect_equal(length(result$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
})

# ----------------------------------------------------------------------------
# Test 2.6: Appends to a pre-existing preprocessing history (does not overwrite)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: c(eeg_obj$preprocessing_history, list(history_entry))
# A pre-existing history entry must survive the call unchanged.
test_that("eeg_bandpass() appends to pre-existing preprocessing history", {
  eeg <- make_eeg()
  eeg$preprocessing_history <- list("Prior step: re-referenced to average")
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expect_equal(length(result$preprocessing_history), 2)
  expect_equal(result$preprocessing_history[[1]],
               "Prior step: re-referenced to average")
})

# ----------------------------------------------------------------------------
# Test 2.7: Bandpass history entry contains both cutoff frequencies
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The history_entry string built for a bandpass call:
#   paste0("Bandpass filtered: ", low_freq, " - ", high_freq, " Hz | ...")
# Both frequency values and the word "Bandpass" must appear.
test_that("eeg_bandpass() bandpass history entry contains cutoff frequencies", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("1",  entry))
  expect_true(grepl("40", entry))
  expect_true(grepl("(?i)bandpass|band.pass", entry, perl = TRUE))
})

# ----------------------------------------------------------------------------
# Test 2.8: Highpass-only history entry says "Highpass"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   filter_desc <- if (!is.null(low_freq) && is.null(high_freq))
#                    paste0("Highpass filtered: ", low_freq, " Hz")
test_that("eeg_bandpass() highpass-only history entry is labelled correctly", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = NULL, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("(?i)highpass|high.pass", entry, perl = TRUE))
  expect_true(grepl("1", entry))
})

# ----------------------------------------------------------------------------
# Test 2.9: Lowpass-only history entry says "Lowpass"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   filter_desc <- paste0("Lowpass filtered: ", high_freq, " Hz")
test_that("eeg_bandpass() lowpass-only history entry is labelled correctly", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = NULL, high_freq = 40, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("(?i)lowpass|low.pass", entry, perl = TRUE))
  expect_true(grepl("40", entry))
})

# ----------------------------------------------------------------------------
# Test 2.10: History entry records the method, order, and zero-phase flag
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The history string includes:
#   " | method: butter, order: 4, zero-phase: yes, channels: ..."
test_that("eeg_bandpass() history entry records method, order, and zero-phase", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         method = "butter", filter_order = 4,
                         zero_phase = TRUE, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("butter", entry))
  expect_true(grepl("4",      entry))
  expect_true(grepl("yes",    entry))
})


# ============================================================================
#     TEST SUITE 3: eeg_bandpass() - Channel Selection Behaviour
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: NULL channels filters all channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (is.null(channels)) ch_idx <- seq_len(n_channels)
# Every row of data must be modified when channels = NULL.
test_that("eeg_bandpass() filters all channels when channels = NULL", {
  set.seed(1)
  eeg    <- make_eeg(n_channels = 4)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         channels = NULL, verbose = FALSE)
  # Every channel should differ from input after filtering noise
  for (i in seq_len(4)) {
    expect_false(identical(result$data[i, ], eeg$data[i, ]))
  }
})

# ----------------------------------------------------------------------------
# Test 3.2: Integer channel indices — only selected channels are modified
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   ch_idx <- as.integer(channels)   # numeric branch
#   for (i in ch_idx) { ... }        # loop only processes ch_idx
# Channels NOT in ch_idx must be byte-for-byte identical to the input.
test_that("eeg_bandpass() with integer indices only modifies selected channels", {
  eeg    <- make_eeg(n_channels = 4)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         channels = c(1L, 3L), verbose = FALSE)
  # Channels 2 and 4 must be untouched
  expect_equal(result$data[2, ], eeg$data[2, ])
  expect_equal(result$data[4, ], eeg$data[4, ])
  # Channels 1 and 3 must be different (noise input is broadband)
  expect_false(identical(result$data[1, ], eeg$data[1, ]))
  expect_false(identical(result$data[3, ], eeg$data[3, ]))
})

# ----------------------------------------------------------------------------
# Test 3.3: Character channel names — selected channels are modified
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   ch_idx <- match(channels, chan_names)   # character branch
test_that("eeg_bandpass() with character names only modifies named channels", {
  eeg    <- make_eeg(n_channels = 3)  # Ch1, Ch2, Ch3
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         channels = "Ch2", verbose = FALSE)
  expect_equal(result$data[1, ], eeg$data[1, ])
  expect_equal(result$data[3, ], eeg$data[3, ])
  expect_false(identical(result$data[2, ], eeg$data[2, ]))
})


# ============================================================================
#     TEST SUITE 4: eeg_bandpass() - Signal Quality (Filtering Efficacy)
# ============================================================================
#
# These tests verify that filtering actually attenuates the correct frequencies.
# A pure sine wave is created at a target frequency, filtered, and the power
# ratio (filtered / original) is inspected.

# ----------------------------------------------------------------------------
# Test 4.1: Highpass filtering attenuates sub-cutoff signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The .apply_filter() / signal::filtfilt() call for the
# highpass filter object. A 0.1 Hz sine (below a 1 Hz cutoff) must lose most
# of its power after filtering.
test_that("eeg_bandpass() highpass filter strongly attenuates sub-cutoff sine", {
  sr     <- 512; n <- 5120    # 10 seconds — long enough for stable low-freq filter
  sine   <- sine_matrix(freq_hz = 0.1, sr = sr, n_samp = n, n_channels = 2)
  eeg    <- make_eeg(n_channels = 2, n_timepoints = n, sr = sr,
                     data_values = sine)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = NULL, verbose = FALSE)
  
  # Power ratio < 5 % = strongly attenuated
  power_in  <- signal_power(eeg$data[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_lt(power_out / power_in, 0.05)
})

# ----------------------------------------------------------------------------
# Test 4.2: Highpass filtering passes above-cutoff signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The same highpass filter must leave a 10 Hz sine mostly
# intact (passband). Power ratio must exceed 90 %.
test_that("eeg_bandpass() highpass filter preserves above-cutoff sine", {
  sr     <- 512; n <- 2560
  sine   <- sine_matrix(freq_hz = 10, sr = sr, n_samp = n, n_channels = 2)
  eeg    <- make_eeg(n_channels = 2, n_timepoints = n, sr = sr,
                     data_values = sine)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = NULL, verbose = FALSE)
  
  power_in  <- signal_power(eeg$data[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_gt(power_out / power_in, 0.90)
})

# ----------------------------------------------------------------------------
# Test 4.3: Lowpass filtering attenuates above-cutoff signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The lowpass filter object applied via .apply_filter().
# A 80 Hz sine (above a 20 Hz cutoff) must be strongly attenuated.
test_that("eeg_bandpass() lowpass filter strongly attenuates above-cutoff sine", {
  sr     <- 512; n <- 2560
  sine   <- sine_matrix(freq_hz = 80, sr = sr, n_samp = n, n_channels = 2)
  eeg    <- make_eeg(n_channels = 2, n_timepoints = n, sr = sr,
                     data_values = sine)
  result <- eeg_bandpass(eeg, low_freq = NULL, high_freq = 20, verbose = FALSE)
  
  power_in  <- signal_power(eeg$data[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_lt(power_out / power_in, 0.05)
})

# ----------------------------------------------------------------------------
# Test 4.4: Lowpass filtering preserves below-cutoff signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The lowpass filter must pass a 5 Hz sine (below 20 Hz
# cutoff) with > 90 % power retention.
test_that("eeg_bandpass() lowpass filter preserves below-cutoff sine", {
  sr     <- 512; n <- 2560
  sine   <- sine_matrix(freq_hz = 5, sr = sr, n_samp = n, n_channels = 2)
  eeg    <- make_eeg(n_channels = 2, n_timepoints = n, sr = sr,
                     data_values = sine)
  result <- eeg_bandpass(eeg, low_freq = NULL, high_freq = 20, verbose = FALSE)
  
  power_in  <- signal_power(eeg$data[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_gt(power_out / power_in, 0.90)
})

# ----------------------------------------------------------------------------
# Test 4.5: Bandpass filtering attenuates both out-of-band ends
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Combined application of highpass + lowpass in sequence.
# A 0.2 Hz (below passband 1-40 Hz) and a 100 Hz (above passband) tone must
# both be strongly attenuated.
test_that("eeg_bandpass() bandpass attenuates low-freq and high-freq sines", {
  sr <- 512; n <- 5120
  
  # Low-freq tone: 0.2 Hz
  low_sine  <- sine_matrix(0.2,  sr, n, 1)
  eeg_low   <- make_eeg(1, n, sr, data_values = low_sine)
  res_low   <- eeg_bandpass(eeg_low, 1, 40, verbose = FALSE)
  expect_lt(signal_power(res_low$data[1, ]) / signal_power(low_sine[1, ]), 0.05)
  
  # High-freq tone: 100 Hz
  high_sine <- sine_matrix(100, sr, n, 1)
  eeg_high  <- make_eeg(1, n, sr, data_values = high_sine)
  res_high  <- eeg_bandpass(eeg_high, 1, 40, verbose = FALSE)
  expect_lt(signal_power(res_high$data[1, ]) / signal_power(high_sine[1, ]), 0.05)
})

# ----------------------------------------------------------------------------
# Test 4.6: Bandpass filtering preserves in-band signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Same bandpass (1-40 Hz) must preserve a 10 Hz sine.
test_that("eeg_bandpass() bandpass preserves in-band sine", {
  sr   <- 512; n <- 2560
  sine <- sine_matrix(10, sr, n, 2)
  eeg  <- make_eeg(2, n, sr, data_values = sine)
  res  <- eeg_bandpass(eeg, 1, 40, verbose = FALSE)
  expect_gt(signal_power(res$data[1, ]) / signal_power(sine[1, ]), 0.90)
})

# ----------------------------------------------------------------------------
# Test 4.7: FIR method also attenuates out-of-band signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The FIR branch in the "design filters" block:
#   } else {  hp_filt <- signal::fir1(fir_len - 1, W_low, type = "high") }
# Verifies the FIR code path is reachable and functional.
test_that("eeg_bandpass() FIR method attenuates sub-cutoff sine", {
  sr     <- 512; n <- 5120
  sine   <- sine_matrix(0.2, sr, n, 2)
  eeg    <- make_eeg(2, n, sr, data_values = sine)
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = NULL,
                         method = "fir", verbose = FALSE)
  power_in  <- signal_power(sine[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_lt(power_out / power_in, 0.10)
})


# ============================================================================
#     TEST SUITE 5: eeg_bandpass() - Padding and Zero-Phase Options
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Explicit padding = 0 disables padding
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (padding > 0) { left_pad <- ... }   (the mirror-pad block is skipped)
# With padding = 0, the result must still be a valid filtered eeg object.
test_that("eeg_bandpass() works correctly with padding = 0", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         padding = 0, verbose = FALSE)
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
})

# ----------------------------------------------------------------------------
# Test 5.2: Auto-padding (NULL) produces correct output dimensions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (is.null(padding)) padding <- 3L * as.integer(filter_order)
# Mirror-padded samples are added and then removed; data shape must be
# identical to the input.
test_that("eeg_bandpass() auto-padding (NULL) does not change data dimensions", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         padding = NULL, verbose = FALSE)
  expect_equal(dim(result$data), dim(eeg$data))
})

# ----------------------------------------------------------------------------
# Test 5.3: zero_phase = FALSE (causal filtering) still returns valid output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (zero_phase) { signal::filtfilt(...) } else { signal::filter(...) }
# The causal code path must produce a numeric output of the right shape.
test_that("eeg_bandpass() zero_phase = FALSE returns valid filtered output", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         zero_phase = FALSE, verbose = FALSE)
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
  expect_true(all(is.finite(result$data)))
})

# ----------------------------------------------------------------------------
# Test 5.4: Inline notch via notch_freq argument in eeg_bandpass()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The optional notch block inside eeg_bandpass():
#   if (!is.null(notch_freq)) { notch_filt <- butter(..., type="stop"); ... }
# A 50 Hz sine passed through the inline notch should be attenuated.
test_that("eeg_bandpass() inline notch_freq attenuates the target frequency", {
  sr   <- 512; n <- 4096
  sine <- sine_matrix(50, sr, n, 2)
  eeg  <- make_eeg(2, n, sr, data_values = sine)
  res  <- eeg_bandpass(eeg, low_freq = 1, high_freq = 100,
                       notch_freq = 50, notch_bandwidth = 4,
                       verbose = FALSE)
  power_in  <- signal_power(sine[1, ])
  power_out <- signal_power(res$data[1, ])
  expect_lt(power_out / power_in, 0.10)
})

# ----------------------------------------------------------------------------
# Test 5.5: Inline notch history entry mentions the notch frequency
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The history_entry construction when notch_freq is supplied:
#   if (!is.null(notch_freq)) paste0(", notch: ", notch_freq, " Hz")
test_that("eeg_bandpass() history entry includes inline notch frequency", {
  eeg    <- make_eeg()
  result <- eeg_bandpass(eeg, low_freq = 1, high_freq = 40,
                         notch_freq = 50, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("50", entry))
  expect_true(grepl("notch", entry, ignore.case = TRUE))
})


# ============================================================================
#         TEST SUITE 6: eeg_notch() - Input Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: Rejects non-eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!inherits(eeg_obj, "eeg")) stop("'eeg_obj' must be an object of class 'eeg'.")
test_that("eeg_notch() rejects a non-eeg object", {
  not_eeg <- list(data = matrix(1:12, 2, 6))
  expect_error(
    eeg_notch(not_eeg, freq = 50, verbose = FALSE),
    regexp = "class 'eeg'"
  )
})

# ----------------------------------------------------------------------------
# Test 6.2: Rejects non-positive freq
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.numeric(freq) || length(freq) != 1 || freq <= 0) stop(...)
test_that("eeg_notch() errors on non-positive freq", {
  eeg <- make_eeg()
  expect_error(eeg_notch(eeg, freq = 0,  verbose = FALSE), regexp = "positive")
  expect_error(eeg_notch(eeg, freq = -5, verbose = FALSE), regexp = "positive")
})

# ----------------------------------------------------------------------------
# Test 6.3: Rejects freq at or above Nyquist
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (freq >= nyquist) stop("'freq' ... must be less than the Nyquist frequency")
test_that("eeg_notch() errors when freq >= Nyquist", {
  eeg <- make_eeg(sr = 256)   # Nyquist = 128 Hz
  expect_error(
    eeg_notch(eeg, freq = 128, verbose = FALSE),
    regexp = "Nyquist"
  )
})

# ----------------------------------------------------------------------------
# Test 6.4: Rejects non-positive bandwidth
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.numeric(bandwidth) || bandwidth <= 0) stop(...)
test_that("eeg_notch() errors on non-positive bandwidth", {
  eeg <- make_eeg()
  expect_error(eeg_notch(eeg, freq = 50, bandwidth = 0,  verbose = FALSE),
               regexp = "positive")
  expect_error(eeg_notch(eeg, freq = 50, bandwidth = -1, verbose = FALSE),
               regexp = "positive")
})

# ----------------------------------------------------------------------------
# Test 6.5: Rejects bandwidth that pushes lower notch edge to zero or below
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   notch_low <- freq - bandwidth / 2
#   if (notch_low <= 0) stop("Notch lower edge ... must be above 0 Hz.")
test_that("eeg_notch() errors when notch lower edge <= 0", {
  eeg <- make_eeg()
  # freq = 2, bandwidth = 6 → lower edge = 2 - 3 = -1 Hz
  expect_error(
    eeg_notch(eeg, freq = 2, bandwidth = 6, verbose = FALSE),
    regexp = "lower edge|Reduce|bandwidth"
  )
})

# ----------------------------------------------------------------------------
# Test 6.6: Rejects bandwidth that pushes upper notch edge to Nyquist or above
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   notch_high <- freq + bandwidth / 2
#   if (notch_high >= nyquist) stop("Notch upper edge ...")
test_that("eeg_notch() errors when notch upper edge >= Nyquist", {
  eeg <- make_eeg(sr = 256)   # Nyquist = 128
  # freq = 126, bandwidth = 6 → upper edge = 126 + 3 = 129 >= 128
  expect_error(
    eeg_notch(eeg, freq = 126, bandwidth = 6, verbose = FALSE),
    regexp = "upper edge|Nyquist|bandwidth"
  )
})

# ----------------------------------------------------------------------------
# Test 6.7: Rejects non-integer filter_order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.numeric(filter_order) || filter_order < 1 || filter_order %% 1 != 0)
#     stop(...)
test_that("eeg_notch() errors when filter_order is not a positive integer", {
  eeg <- make_eeg()
  expect_error(eeg_notch(eeg, freq = 50, filter_order = 0,   verbose = FALSE),
               regexp = "positive integer")
  expect_error(eeg_notch(eeg, freq = 50, filter_order = 1.5, verbose = FALSE),
               regexp = "positive integer")
})

# ----------------------------------------------------------------------------
# Test 6.8: Rejects non-logical harmonics
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.logical(harmonics) || length(harmonics) != 1) stop(...)
test_that("eeg_notch() errors when harmonics is not logical", {
  eeg <- make_eeg()
  expect_error(eeg_notch(eeg, freq = 50, harmonics = "yes", verbose = FALSE),
               regexp = "harmonics")
})

# ----------------------------------------------------------------------------
# Test 6.9: Rejects data with NA values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (any(is.na(eeg_obj$data))) stop(...)
test_that("eeg_notch() errors on NA in data", {
  eeg            <- make_eeg()
  eeg$data[1, 1] <- NA
  expect_error(eeg_notch(eeg, freq = 50, verbose = FALSE), regexp = "NA")
})

# ----------------------------------------------------------------------------
# Test 6.10: Rejects data with Inf values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (any(!is.finite(eeg_obj$data))) stop(...)
test_that("eeg_notch() errors on Inf in data", {
  eeg            <- make_eeg()
  eeg$data[2, 3] <- -Inf
  expect_error(eeg_notch(eeg, freq = 50, verbose = FALSE), regexp = "non-finite|Inf")
})

# ----------------------------------------------------------------------------
# Test 6.11: Warns when padding >= n_timepoints
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (padding > 0 && padding >= n_timepoints) warning(...)
test_that("eeg_notch() warns when padding >= n_timepoints", {
  eeg <- make_eeg(n_timepoints = 50)
  expect_warning(
    eeg_notch(eeg, freq = 50, padding = 100, verbose = FALSE),
    regexp = "padding"
  )
})


# ============================================================================
#     TEST SUITE 7: eeg_notch() - Output Structure and Metadata
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: Returns an object of class "eeg"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The returned object is built with new_eeg() and must carry
# the "eeg" class tag for dispatch by downstream functions.
test_that("eeg_notch() returns an object of class 'eeg'", {
  eeg    <- make_eeg()
  result <- eeg_notch(eeg, freq = 50, verbose = FALSE)
  expect_s3_class(result, "eeg")
})

# ----------------------------------------------------------------------------
# Test 7.2: Data dimensions are unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The filtered_data matrix is initialised as a copy of
# eeg_obj$data. Padding is added and removed; dimensions must be invariant.
test_that("eeg_notch() preserves data matrix dimensions", {
  eeg    <- make_eeg(n_channels = 5, n_timepoints = 1024)
  result <- eeg_notch(eeg, freq = 50, verbose = FALSE)
  expect_equal(dim(result$data), dim(eeg$data))
})

# ----------------------------------------------------------------------------
# Test 7.3: Non-data fields are forwarded unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: channels, sampling_rate, times, events, metadata, reference
# are all passed directly to new_eeg() without modification.
test_that("eeg_notch() forwards all non-data fields unchanged", {
  eeg           <- make_eeg()
  eeg$metadata  <- list(subject = "S02")
  eeg$reference <- "linked_mastoids"
  result        <- eeg_notch(eeg, freq = 50, verbose = FALSE)
  expect_equal(result$sampling_rate, eeg$sampling_rate)
  expect_equal(result$channels,      eeg$channels)
  expect_equal(result$metadata,      eeg$metadata)
  expect_equal(result$reference,     eeg$reference)
})

# ----------------------------------------------------------------------------
# Test 7.4: Appends exactly one entry to preprocessing history
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   new_history <- c(eeg_obj$preprocessing_history, list(history_entry))
test_that("eeg_notch() appends exactly one preprocessing history entry", {
  eeg    <- make_eeg()
  result <- eeg_notch(eeg, freq = 50, verbose = FALSE)
  expect_equal(length(result$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
})

# ----------------------------------------------------------------------------
# Test 7.5: History entry contains the notch frequency
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   history_entry <- paste0("Notch filtered: ", freq, " Hz ...")
test_that("eeg_notch() history entry contains the notch frequency", {
  eeg    <- make_eeg()
  result <- eeg_notch(eeg, freq = 50, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("50",    entry))
  expect_true(grepl("(?i)notch", entry, perl = TRUE))
})

# ----------------------------------------------------------------------------
# Test 7.6: History entry contains bandwidth, order, and zero-phase flag
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The full history_entry string includes:
#   " | bandwidth: 2 Hz, order: 4, zero-phase: yes, ..."
test_that("eeg_notch() history entry records bandwidth, order, and zero-phase", {
  eeg    <- make_eeg()
  result <- eeg_notch(eeg, freq = 50, bandwidth = 2,
                      filter_order = 4, zero_phase = TRUE, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  expect_true(grepl("2",   entry))   # bandwidth
  expect_true(grepl("4",   entry))   # order
  expect_true(grepl("yes", entry))   # zero-phase
})


# ============================================================================
#     TEST SUITE 8: eeg_notch() - Harmonics Logic
# ============================================================================
#
# These tests verify the harmonic frequency expansion block:
#   if (harmonics) { harmonic_mult <- 2; while ((freq * harmonic_mult) < nyquist) ... }

# ----------------------------------------------------------------------------
# Test 8.1: harmonics = FALSE (default) creates only one notch filter object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Without harmonics, notch_freqs = c(freq) and only one
# filter is designed. The output data matrix must still be the correct shape.
test_that("eeg_notch() with harmonics = FALSE applies only the fundamental notch", {
  eeg    <- make_eeg()
  result <- eeg_notch(eeg, freq = 50, harmonics = FALSE, verbose = FALSE)
  # Dimensions and class
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
  # History must not mention additional harmonic frequencies
  entry <- result$preprocessing_history[[1]]
  expect_false(grepl("100|150", entry))
})

# ----------------------------------------------------------------------------
# Test 8.2: harmonics = TRUE attenuates the harmonic frequencies
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With harmonics = TRUE and a 50 Hz notch on a 512 Hz
# recording (Nyquist = 256), harmonics at 100, 150, 200 Hz are added to
# notch_freqs and a separate filter is applied for each. A pure 100 Hz sine
# must therefore be significantly attenuated.
test_that("eeg_notch() with harmonics = TRUE attenuates harmonic tones", {
  sr   <- 512; n <- 4096
  # 100 Hz = first harmonic of 50 Hz notch
  sine_100 <- sine_matrix(100, sr, n, 2)
  eeg      <- make_eeg(2, n, sr, data_values = sine_100)
  result   <- eeg_notch(eeg, freq = 50, harmonics = TRUE, verbose = FALSE)
  
  power_in  <- signal_power(sine_100[1, ])
  power_out <- signal_power(result$data[1, ])
  expect_lt(power_out / power_in, 0.10)
})

# ----------------------------------------------------------------------------
# Test 8.3: harmonics = TRUE history entry lists harmonic frequencies
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: In the history_entry construction:
#   if (harmonics) paste0(", harmonics: ", paste(notch_freqs[-1], collapse=", "), " Hz")
# The entry must include at least one harmonic frequency value.
test_that("eeg_notch() harmonics history entry lists harmonic frequencies", {
  eeg    <- make_eeg(sr = 512)
  result <- eeg_notch(eeg, freq = 50, harmonics = TRUE, verbose = FALSE)
  entry  <- result$preprocessing_history[[1]]
  # At least the first harmonic (100 Hz) must appear
  expect_true(grepl("100", entry))
  expect_true(grepl("(?i)harmonic", entry, perl = TRUE))
})


# ============================================================================
#     TEST SUITE 9: eeg_notch() - Signal Quality (Filtering Efficacy)
# ============================================================================

# ----------------------------------------------------------------------------
# Test 9.1: 50 Hz sine is strongly attenuated by a 50 Hz notch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The .make_notch_filter() and .apply_filter() pipeline:
#   signal::butter(filter_order, c(W_low, W_high), type = "stop")
#   signal::filtfilt(filt_obj, signal_vec)
test_that("eeg_notch() strongly attenuates a 50 Hz sine", {
  sr   <- 512; n <- 4096
  sine <- sine_matrix(50, sr, n, 2)
  eeg  <- make_eeg(2, n, sr, data_values = sine)
  res  <- eeg_notch(eeg, freq = 50, bandwidth = 4, verbose = FALSE)
  
  power_in  <- signal_power(sine[1, ])
  power_out <- signal_power(res$data[1, ])
  expect_lt(power_out / power_in, 0.05)
})

# ----------------------------------------------------------------------------
# Test 9.2: An out-of-band sine (10 Hz) is preserved by the 50 Hz notch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Only the 50 Hz band is removed; frequencies well outside
# the notch must retain > 90 % of their power.
test_that("eeg_notch() preserves out-of-band signal (10 Hz)", {
  sr   <- 512; n <- 2560
  sine <- sine_matrix(10, sr, n, 2)
  eeg  <- make_eeg(2, n, sr, data_values = sine)
  res  <- eeg_notch(eeg, freq = 50, bandwidth = 2, verbose = FALSE)
  
  power_in  <- signal_power(sine[1, ])
  power_out <- signal_power(res$data[1, ])
  expect_gt(power_out / power_in, 0.90)
})

# ----------------------------------------------------------------------------
# Test 9.3: Channel selectivity — unselected channels are unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   filtered_data <- eeg_obj$data   # copy; unselected channels untouched
#   for (i in ch_idx) { ... }
test_that("eeg_notch() does not modify unselected channels", {
  eeg    <- make_eeg(n_channels = 4)
  result <- eeg_notch(eeg, freq = 50, channels = c(1L, 3L), verbose = FALSE)
  expect_equal(result$data[2, ], eeg$data[2, ])
  expect_equal(result$data[4, ], eeg$data[4, ])
})


# ============================================================================
#     TEST SUITE 10: eeg_filter() - Input Validation
# ============================================================================
#
# eeg_filter() is a thin wrapper. Its own validation is limited to two checks:
# class membership and the all-NULL frequency guard. Parameter-level validation
# is delegated to eeg_bandpass() and eeg_notch().

# ----------------------------------------------------------------------------
# Test 10.1: Rejects non-eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!inherits(eeg_obj, "eeg")) stop("'eeg_obj' must be an object of class 'eeg'.")
test_that("eeg_filter() rejects a non-eeg object", {
  not_eeg <- list(data = matrix(1, 2, 10))
  expect_error(
    eeg_filter(not_eeg, low_freq = 1, high_freq = 40, verbose = FALSE),
    regexp = "class 'eeg'"
  )
})

# ----------------------------------------------------------------------------
# Test 10.2: Rejects when all three frequency arguments are NULL
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (is.null(low_freq) && is.null(high_freq) && is.null(notch_freq)) stop(...)
# This is the only validation unique to eeg_filter().
test_that("eeg_filter() errors when low_freq, high_freq, and notch_freq are all NULL", {
  eeg <- make_eeg()
  expect_error(
    eeg_filter(eeg, low_freq = NULL, high_freq = NULL,
               notch_freq = NULL, verbose = FALSE),
    regexp = "At least one"
  )
})


# ============================================================================
#     TEST SUITE 11: eeg_filter() - Pipeline Logic
# ============================================================================
#
# These tests verify that eeg_filter() correctly routes arguments to
# eeg_bandpass() and eeg_notch() and that the output reflects both stages.

# ----------------------------------------------------------------------------
# Test 11.1: Bandpass only (notch_freq = NULL) runs only the bandpass stage
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.null(low_freq) || !is.null(high_freq)) eeg_obj <- eeg_bandpass(...)
#   if (!is.null(notch_freq)) eeg_obj <- eeg_notch(...)   ← skipped
# With notch_freq = NULL only one history entry should be added.
test_that("eeg_filter() bandpass-only call adds exactly one history entry", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq   = 1,
                       high_freq  = 40,
                       notch_freq = NULL,
                       verbose    = FALSE)
  expect_equal(length(result$preprocessing_history), 1)
})

# ----------------------------------------------------------------------------
# Test 11.2: Bandpass + notch call adds exactly two history entries
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Both if-branches execute:
#   eeg_obj <- eeg_bandpass(...)   → history entry 1
#   eeg_obj <- eeg_notch(...)      → history entry 2
# The wrapper must not lose or overwrite either entry.
test_that("eeg_filter() bandpass + notch call adds exactly two history entries", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq   = 1,
                       high_freq  = 40,
                       notch_freq = 50,
                       verbose    = FALSE)
  expect_equal(length(result$preprocessing_history), 2)
})

# ----------------------------------------------------------------------------
# Test 11.3: First history entry is from bandpass, second from notch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The ordering is bandpass → notch. The first history entry
# must mention the bandpass frequencies; the second must mention "Notch".
test_that("eeg_filter() history order is bandpass then notch", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq   = 1,
                       high_freq  = 40,
                       notch_freq = 50,
                       verbose    = FALSE)
  first  <- result$preprocessing_history[[1]]
  second <- result$preprocessing_history[[2]]
  expect_true(grepl("(?i)bandpass|highpass|lowpass|filtered", first,  perl = TRUE))
  expect_true(grepl("(?i)notch",                              second, perl = TRUE))
})

# ----------------------------------------------------------------------------
# Test 11.4: Notch only (both low_freq and high_freq are NULL) skips bandpass
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   if (!is.null(low_freq) || !is.null(high_freq)) ...   ← skipped
#   if (!is.null(notch_freq)) eeg_obj <- eeg_notch(...)  ← runs
# Only one history entry should be added, and it must be from eeg_notch().
test_that("eeg_filter() notch-only call skips bandpass and runs notch", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq   = NULL,
                       high_freq  = NULL,
                       notch_freq = 50,
                       verbose    = FALSE)
  expect_equal(length(result$preprocessing_history), 1)
  entry <- result$preprocessing_history[[1]]
  expect_true(grepl("(?i)notch", entry, perl = TRUE))
})

# ----------------------------------------------------------------------------
# Test 11.5: Returns an object of class "eeg"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The return value of eeg_filter() is the eeg object returned
# by the last stage (eeg_bandpass or eeg_notch), which is always an "eeg".
test_that("eeg_filter() returns an object of class 'eeg'", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg, low_freq = 1, high_freq = 40, verbose = FALSE)
  expect_s3_class(result, "eeg")
})

# ----------------------------------------------------------------------------
# Test 11.6: Data dimensions are unchanged after the full pipeline
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Both stages use new_eeg() which preserves the data matrix
# shape. The final object must have the same shape as the input.
test_that("eeg_filter() preserves data matrix dimensions across both stages", {
  eeg    <- make_eeg(n_channels = 5, n_timepoints = 1024)
  result <- eeg_filter(eeg,
                       low_freq   = 1,
                       high_freq  = 40,
                       notch_freq = 50,
                       verbose    = FALSE)
  expect_equal(dim(result$data), dim(eeg$data))
})

# ----------------------------------------------------------------------------
# Test 11.7: notch_harmonics = TRUE is forwarded to eeg_notch()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   eeg_obj <- eeg_notch(eeg_obj, ..., harmonics = notch_harmonics, ...)
# When notch_harmonics = TRUE the notch history entry must mention harmonics.
test_that("eeg_filter() forwards notch_harmonics = TRUE to eeg_notch()", {
  eeg    <- make_eeg(sr = 512)
  result <- eeg_filter(eeg,
                       low_freq        = 1,
                       high_freq       = 40,
                       notch_freq      = 50,
                       notch_harmonics = TRUE,
                       verbose         = FALSE)
  notch_entry <- result$preprocessing_history[[2]]
  expect_true(grepl("(?i)harmonic|100", notch_entry, perl = TRUE))
})

# ----------------------------------------------------------------------------
# Test 11.8: channels argument is forwarded to both stages
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   eeg_bandpass(..., channels = channels, ...)
#   eeg_notch(...,   channels = channels, ...)
# Unselected channels must be identical to the input in the final output.
test_that("eeg_filter() channels argument leaves unselected channels unchanged", {
  eeg    <- make_eeg(n_channels = 4)
  result <- eeg_filter(eeg,
                       low_freq   = 1,
                       high_freq  = 40,
                       notch_freq = 50,
                       channels   = c(1L, 3L),
                       verbose    = FALSE)
  expect_equal(result$data[2, ], eeg$data[2, ])
  expect_equal(result$data[4, ], eeg$data[4, ])
})

# ----------------------------------------------------------------------------
# Test 11.9: notch_bandwidth is forwarded to eeg_notch()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS:
#   eeg_obj <- eeg_notch(eeg_obj, ..., bandwidth = notch_bandwidth, ...)
# A wider bandwidth should remove a broader swath around 50 Hz. We verify
# the history entry records the custom bandwidth value.
test_that("eeg_filter() forwards custom notch_bandwidth to eeg_notch()", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq        = 1,
                       high_freq       = 40,
                       notch_freq      = 50,
                       notch_bandwidth = 4,
                       verbose         = FALSE)
  notch_entry <- result$preprocessing_history[[2]]
  expect_true(grepl("4", notch_entry))
})

# ----------------------------------------------------------------------------
# Test 11.10: verbose = FALSE suppresses all console output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: verbose = FALSE is forwarded to both eeg_bandpass() and
# eeg_notch(), which each guard their cat() calls with the verbose flag.
# No output should be produced.
test_that("eeg_filter() with verbose = FALSE produces no console output", {
  eeg <- make_eeg()
  expect_silent(
    eeg_filter(eeg,
               low_freq   = 1,
               high_freq  = 40,
               notch_freq = 50,
               verbose    = FALSE)
  )
})


# ============================================================================
#     TEST SUITE 12: eeg_filter() - Signal Quality (End-to-End)
# ============================================================================

# ----------------------------------------------------------------------------
# Test 12.1: Full pipeline attenuates both out-of-band and notch frequency
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The combined eeg_bandpass() + eeg_notch() pipeline run by
# eeg_filter(). A signal composed of three sines (sub-cutoff + in-band + notch
# tone) is filtered, and only the in-band component should survive.
test_that("eeg_filter() full pipeline attenuates out-of-band and notch tones", {
  sr  <- 512; n <- 5120
  t   <- seq(0, (n - 1) / sr, length.out = n)
  
  # Three-component signal: 0.2 Hz (sub-HP) + 10 Hz (passband) + 50 Hz (notch)
  composite <- sin(2 * pi * 0.2  * t) +
    sin(2 * pi * 10   * t) +
    sin(2 * pi * 50   * t)
  data_mat  <- matrix(rep(composite, 2), nrow = 2, byrow = TRUE)
  eeg       <- make_eeg(2, n, sr, data_values = data_mat)
  
  result    <- eeg_filter(eeg,
                          low_freq   = 1,
                          high_freq  = 40,
                          notch_freq = 50,
                          verbose    = FALSE)
  
  filtered  <- result$data[1, ]
  
  # Power of filtered output should be much lower than sum-of-three input
  # (two of three components removed) but definitely non-zero (10 Hz retained)
  power_in  <- signal_power(composite)
  power_out <- signal_power(filtered)
  
  # Output power should be substantially less than input (two sines removed)
  expect_lt(power_out, power_in * 0.70)
  # But not near zero — the 10 Hz component must still be present
  expect_gt(power_out, 0.05)
})

# ----------------------------------------------------------------------------
# Test 12.2: Highpass-only call in eeg_filter() (high_freq = NULL) works
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The branch
#   if (!is.null(low_freq) || !is.null(high_freq)) eeg_obj <- eeg_bandpass(...)
# is triggered by low_freq alone. The result must be a valid eeg.
test_that("eeg_filter() with only low_freq (highpass-only) returns valid eeg", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq  = 1,
                       high_freq = NULL,
                       verbose   = FALSE)
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
  expect_equal(length(result$preprocessing_history), 1)
})

# ----------------------------------------------------------------------------
# Test 12.3: Lowpass-only call in eeg_filter() (low_freq = NULL) works
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The same bandpass branch triggered by high_freq alone.
test_that("eeg_filter() with only high_freq (lowpass-only) returns valid eeg", {
  eeg    <- make_eeg()
  result <- eeg_filter(eeg,
                       low_freq  = NULL,
                       high_freq = 40,
                       verbose   = FALSE)
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
  expect_equal(length(result$preprocessing_history), 1)
})