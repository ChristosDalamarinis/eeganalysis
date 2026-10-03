# ============================================================================
#                       Test File for downsample.R
# ============================================================================
#
# This test file provides comprehensive testing for the downsample.R script,
# which implements EEG downsampling with anti-aliasing filters.
#
# Function tested:
#   downsample() - Reduces sampling rate of EEG data with anti-aliasing
#
# Dependencies:
#   - testthat
#   - eeganalysis
#   - signal (for the "decimate" method)
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# HELPER FUNCTION: Create test EEG objects
# ============================================================================
# A small utility to quickly build EEG objects for each test, so we don't
# repeat the same boilerplate over and over. All arguments have sensible
# defaults so individual tests only need to override what they care about.

make_eeg <- function(n_channels   = 3,
                     n_timepoints = 1024,
                     sampling_rate = 512,
                     with_events   = FALSE,
                     event_onsets  = NULL,
                     data_values   = NULL) {
  
  if (is.null(data_values)) {
    set.seed(42)
    data_values <- matrix(rnorm(n_channels * n_timepoints),
                          nrow = n_channels, ncol = n_timepoints)
  }
  
  channels <- paste0("Ch", seq_len(n_channels))
  times    <- seq(0, (n_timepoints - 1) / sampling_rate, length.out = n_timepoints)
  
  events <- if (with_events && !is.null(event_onsets)) {
    data.frame(
      onset       = as.integer(event_onsets),
      onset_time  = times[event_onsets],
      type        = rep("S1", length(event_onsets)),
      description = rep("stimulus", length(event_onsets))
    )
  } else {
    data.frame(
      onset       = integer(0),
      onset_time  = numeric(0),
      type        = character(0),
      description = character(0)
    )
  }
  
  new_eeg(
    data              = data_values,
    channels          = channels,
    sampling_rate     = sampling_rate,
    times             = times,
    events            = events
  )
}


# ============================================================================
# TEST SUITE 1: Input Validation - eeg_obj argument
# ============================================================================

test_that("downsample() errors when eeg_obj is not an eeg object", {
  # WHAT THIS TESTS:
  # The function must reject objects that are not of class 'eeg'. Passing a
  # plain list, a matrix, or a data.frame should throw an informative error.
  # This guards against silent mis-use.
  
  expect_error(downsample(list(data = matrix(1:100, 10, 10)), target_rate = 128),
               regexp = NULL)   # any error is acceptable
  
  expect_error(downsample(matrix(1:100, 10, 10), target_rate = 128),
               regexp = NULL)
  
  expect_error(downsample("not_an_eeg", target_rate = 128),
               regexp = NULL)
})


# ============================================================================
# TEST SUITE 2: Input Validation - target_rate argument
# ============================================================================

test_that("downsample() errors when target_rate is not a single positive numeric", {
  # WHAT THIS TESTS:
  # target_rate must be a single, finite, positive number.  The function
  # should reject: negative values, zero, character strings, vectors, NULL.
  
  eeg <- make_eeg()
  
  expect_error(downsample(eeg, target_rate = -100,   verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 0,      verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = "256",  verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = c(128, 256), verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = NULL,   verbose = FALSE), regexp = NULL)
})

test_that("downsample() errors when target_rate >= current sampling_rate", {
  # WHAT THIS TESTS:
  # Downsampling by definition must reduce the sampling rate.  Requesting a
  # target_rate that equals or exceeds the original rate (512 Hz here) is
  # logically impossible and must raise an error.
  
  eeg <- make_eeg(sampling_rate = 512)
  
  expect_error(downsample(eeg, target_rate = 512, verbose = FALSE),
               regexp = "must be less than current sampling rate")
  
  expect_error(downsample(eeg, target_rate = 1024, verbose = FALSE),
               regexp = "must be less than current sampling rate")
})

test_that("downsample() errors when downsampling factor is too small (< 1.5)", {
  # WHAT THIS TESTS:
  # The calculated factor = sampling_rate / target_rate must be at least 1.5
  # to constitute meaningful downsampling.  E.g. going from 512 to 400 Hz
  # gives factor = 1.28 which is below the threshold and must error.
  
  eeg <- make_eeg(sampling_rate = 512)
  
  expect_error(downsample(eeg, target_rate = 400, verbose = FALSE),
               regexp = "Downsampling factor too small")
})


# ============================================================================
# TEST SUITE 3: Input Validation - filter_order argument
# ============================================================================

test_that("downsample() errors when filter_order is invalid", {
  # WHAT THIS TESTS:
  # filter_order must be a single positive integer.  Non-integers, zero,
  # negative values, and vectors must all be rejected.
  
  eeg <- make_eeg()
  
  expect_error(downsample(eeg, target_rate = 256, filter_order = 0,    verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, filter_order = -4,   verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, filter_order = 3.5,  verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, filter_order = c(4, 8), verbose = FALSE), regexp = NULL)
})


# ============================================================================
# TEST SUITE 4: Input Validation - cutoff_ratio argument
# ============================================================================

test_that("downsample() errors when cutoff_ratio is outside (0, 1]", {
  # WHAT THIS TESTS:
  # cutoff_ratio defines the low-pass filter cutoff as a fraction of the new
  # Nyquist frequency, so it must lie strictly between 0 and 1 (inclusive of 1).
  # Zero, negatives, and values > 1 are physically meaningless and must error.
  
  eeg <- make_eeg()
  
  expect_error(downsample(eeg, target_rate = 256, cutoff_ratio = 0,    verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, cutoff_ratio = -0.5, verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, cutoff_ratio = 1.1,  verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, cutoff_ratio = 2,    verbose = FALSE), regexp = NULL)
})


# ============================================================================
# TEST SUITE 5: Input Validation - preserve_bands argument
# ============================================================================

test_that("downsample() errors when preserve_bands is malformed", {
  # WHAT THIS TESTS:
  # preserve_bands must be a numeric vector of exactly 2 elements where
  # element[1] < element[2] and both are positive. The function validates
  # this strictly before any processing begins.
  
  eeg <- make_eeg()
  
  # Single value instead of pair
  expect_error(downsample(eeg, target_rate = 256, preserve_bands = 50, verbose = FALSE),
               regexp = "numeric vector of length 2")
  
  # Three values
  expect_error(downsample(eeg, target_rate = 256, preserve_bands = c(1, 30, 50), verbose = FALSE),
               regexp = NULL)
  
  # min >= max
  expect_error(downsample(eeg, target_rate = 256, preserve_bands = c(50, 30), verbose = FALSE),
               regexp = NULL)
  
  # Negative values
  expect_error(downsample(eeg, target_rate = 256, preserve_bands = c(-1, 30), verbose = FALSE),
               regexp = NULL)
})

test_that("downsample() errors when preserve_bands frequencies exceed filter cutoff", {
  # WHAT THIS TESTS:
  # If the user asks to preserve frequencies up to 100 Hz but the target
  # rate of 256 Hz only supports ~115 Hz Nyquist (cutoff ~103 Hz at ratio 0.9),
  # preservation of up to 100 Hz is fine — but requesting preservation up to
  # 120 Hz with a target of 256 Hz (cutoff ~115 Hz) must error because 120 Hz
  # cannot survive the anti-aliasing filter.
  
  eeg <- make_eeg(sampling_rate = 512)
  
  # 256 Hz target, cutoff_ratio 0.9 -> cutoff = (256/2)*0.9 = 115.2 Hz
  # Requesting preservation up to 120 Hz must fail.
  expect_error(
    downsample(eeg, target_rate = 256, preserve_bands = c(1, 120),
               cutoff_ratio = 0.9, verbose = FALSE),
    regexp = "Cannot preserve frequencies"
  )
})

test_that("downsample() succeeds when preserve_bands frequencies are within cutoff", {
  # WHAT THIS TESTS:
  # Conversely, when the user requests preservation within valid bounds
  # (e.g., up to 50 Hz with a 256 Hz target giving ~115 Hz cutoff), the
  # function must NOT error on this argument.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  
  expect_no_error(
    downsample(eeg, target_rate = 256, preserve_bands = c(0.5, 50),
               verbose = FALSE)
  )
})


# ============================================================================
# TEST SUITE 6: Input Validation - logical flags
# ============================================================================

test_that("downsample() errors when remove_outbound_events is not a single logical", {
  # WHAT THIS TESTS:
  # remove_outbound_events must be exactly TRUE or FALSE.
  
  eeg <- make_eeg()
  
  expect_error(downsample(eeg, target_rate = 256,
                          remove_outbound_events = "yes", verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256,
                          remove_outbound_events = 1, verbose = FALSE), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256,
                          remove_outbound_events = c(TRUE, FALSE), verbose = FALSE), regexp = NULL)
})

test_that("downsample() errors when verbose is not a single logical", {
  # WHAT THIS TESTS:
  # verbose must be exactly TRUE or FALSE.
  
  eeg <- make_eeg()
  
  expect_error(downsample(eeg, target_rate = 256, verbose = "true"), regexp = NULL)
  expect_error(downsample(eeg, target_rate = 256, verbose = 1),      regexp = NULL)
})


# ============================================================================
# TEST SUITE 7: Input Validation - data quality (NA / Inf values)
# ============================================================================

test_that("downsample() errors when EEG data contains NA values", {
  # WHAT THIS TESTS:
  # NA values in EEG data would silently corrupt the anti-aliasing filter
  # output.  The function must detect and reject data with NAs before any
  # filtering occurs.
  
  data_with_na <- matrix(rnorm(3 * 1024), nrow = 3, ncol = 1024)
  data_with_na[2, 500] <- NA
  
  eeg <- make_eeg(data_values = data_with_na)
  
  expect_error(downsample(eeg, target_rate = 256, verbose = FALSE),
               regexp = "NA values")
})

test_that("downsample() errors when EEG data contains Inf values", {
  # WHAT THIS TESTS:
  # Similarly, infinite values (e.g., from division by zero during preprocessing)
  # are non-finite and cannot be filtered.  The function must catch these too.
  
  data_with_inf <- matrix(rnorm(3 * 1024), nrow = 3, ncol = 1024)
  data_with_inf[1, 300] <- Inf
  
  eeg <- make_eeg(data_values = data_with_inf)
  
  expect_error(downsample(eeg, target_rate = 256, verbose = FALSE),
               regexp = "non-finite")
})


# ============================================================================
# TEST SUITE 8: Core output - sampling rate and data dimensions
# ============================================================================

test_that("downsample() produces correct new sampling rate with exact integer factor", {
  # WHAT THIS TESTS:
  # Going from 512 Hz to 256 Hz gives an exact integer factor of 2.
  # The resulting eeg object must report sampling_rate = 256 Hz.
  # This is the most fundamental correctness check.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(result$sampling_rate, 256)
})

test_that("downsample() produces correct data matrix dimensions", {
  # WHAT THIS TESTS:
  # With 1024 samples at 512 Hz and a target of 256 Hz (factor = 2),
  # we expect ceil(1024 / 2) = 512 samples per channel after decimation.
  # The number of channels must remain unchanged.
  
  eeg    <- make_eeg(n_channels = 4, sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  # Channels preserved
  expect_equal(nrow(result$data), 4)
  
  # Time points halved (factor 2)
  expected_samples <- length(seq(1, 1024, by = 2))
  expect_equal(ncol(result$data), expected_samples)
})

test_that("downsample() preserves channel names exactly", {
  # WHAT THIS TESTS:
  # The channels vector in the output must be identical to the input.
  # Downsampling only changes the time dimension; channel metadata is untouched.
  
  eeg    <- make_eeg(n_channels = 3, sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(result$channels, eeg$channels)
})


# ============================================================================
# TEST SUITE 9: Core output - time vector
# ============================================================================

test_that("downsample() returns the correct time vector length", {
  # WHAT THIS TESTS:
  # The new times vector must have exactly as many entries as there are
  # columns in the downsampled data matrix.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(length(result$times), ncol(result$data))
})

test_that("downsample() time vector values match original time points at decimated indices", {
  # WHAT THIS TESTS:
  # The new time vector is built by selecting every nth entry from the
  # original time vector.  So new_times[k] must equal old_times[1 + (k-1)*factor].
  # Here we check the first and last values as sentinels.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  factor <- 512 / 256  # = 2
  expected_times <- eeg$times[seq(1, 1024, by = factor)]
  
  expect_equal(result$times, expected_times)
})

test_that("downsample() time spacing matches the new sampling rate", {
  # WHAT THIS TESTS:
  # The mean interval between consecutive time points in the output must
  # equal 1 / new_sampling_rate (within floating-point tolerance).
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  mean_dt <- mean(diff(result$times))
  expected_dt <- 1 / 256
  
  expect_equal(mean_dt, expected_dt, tolerance = 1e-10)
})


# ============================================================================
# TEST SUITE 10: Event onset adjustment - all four strategies
# ============================================================================

# HOW THE KEPT SAMPLES LINE UP (needed to read the tests below):
# downsample() keeps original samples 1, 1 + f, 1 + 2f, ... (f = the factor),
# so new sample k is original sample 1 + (k - 1) * f.  With f = 2 the kept
# original samples are 1, 3, 5, ... so an event on original sample 101 sits
# exactly on new sample 51, while an event on original sample 100 sits half-way
# between new samples 50 (original 99) and 51 (original 101).
# An event's new position is therefore (onset - 1) / f + 1, and the strategy
# decides how that position is turned into a whole sample.

test_that("event_strategy = 'round' correctly rescales onset indices", {
  # WHAT THIS TESTS:
  # 'round' moves each event to the closest kept sample.  With f = 4 the kept
  # original samples are 1, 5, 9, ... 101, 105, ...
  #   original 102: neighbours 101 (1 away) and 105 (3 away) -> 101 = new sample 26
  #   original 104: neighbours 101 (3 away) and 105 (1 away) -> 105 = new sample 27
  # Neither event is a tie, so the answer does not depend on rounding rules.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(102, 104))

  result <- downsample(eeg, target_rate = 128, event_strategy = "round", verbose = FALSE)

  expect_equal(result$events$onset, c(26, 27))
})

test_that("event_strategy = 'floor' always rounds onset down", {
  # WHAT THIS TESTS:
  # 'floor' picks the closest kept sample at or BEFORE the event.  With f = 2:
  #   original 100 (between kept 99 and 101)  -> 99  = new sample 50
  #   original 101 (on a kept sample)         -> 101 = new sample 51
  #   original 200 (between kept 199 and 201) -> 199 = new sample 100
  #   original 201 (on a kept sample)         -> 201 = new sample 101

  onsets <- c(100, 101, 200, 201)
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = onsets)

  result <- downsample(eeg, target_rate = 256, event_strategy = "floor", verbose = FALSE)

  expect_equal(result$events$onset, c(50, 51, 100, 101))

  # The same thing in words: the kept sample is never after the event.
  kept_original_sample <- 1 + (result$events$onset - 1) * 2
  expect_true(all(kept_original_sample <= onsets))
})

test_that("event_strategy = 'ceiling' always rounds onset up", {
  # WHAT THIS TESTS:
  # 'ceiling' picks the closest kept sample at or AFTER the event.  With f = 2:
  #   original 100 (between kept 99 and 101)  -> 101 = new sample 51
  #   original 101 (on a kept sample)         -> 101 = new sample 51
  #   original 200 (between kept 199 and 201) -> 201 = new sample 101
  #   original 201 (on a kept sample)         -> 201 = new sample 101

  onsets <- c(100, 101, 200, 201)
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = onsets)

  result <- downsample(eeg, target_rate = 256, event_strategy = "ceiling", verbose = FALSE)

  expect_equal(result$events$onset, c(51, 51, 101, 101))

  # The same thing in words: the kept sample is never before the event.
  kept_original_sample <- 1 + (result$events$onset - 1) * 2
  expect_true(all(kept_original_sample >= onsets))
})

test_that("event_strategy = 'nearest' produces same result as 'round'", {
  # WHAT THIS TESTS:
  # The 'nearest' strategy is implemented identically to 'round' in the
  # source code (both round the same new position).  This test confirms the
  # two produce identical output.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(100, 300, 600))

  result_nearest <- downsample(eeg, target_rate = 256, event_strategy = "nearest", verbose = FALSE)
  result_round   <- downsample(eeg, target_rate = 256, event_strategy = "round",   verbose = FALSE)

  expect_equal(result_nearest$events$onset, result_round$events$onset)
})

test_that("an event sitting on a kept sample points at that sample's data, for every strategy", {
  # WHAT THIS TESTS:
  # The check that matters most: after downsampling, an event marker must name
  # the column of the data where the event really happened.  We put a spike
  # (value 1000) into the data at original samples 101, 201 and 401 (all kept
  # samples when f = 2), mark an event on each spike, downsample without
  # filtering (method = "simple", so the spike is not smeared), and look at
  # which column now holds each spike.  The event markers must name those columns.

  n_tp     <- 1024
  spikes   <- c(101, 201, 401)
  data_mat <- matrix(0, nrow = 3, ncol = n_tp)
  data_mat[, spikes] <- 1000

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = n_tp,
                  with_events = TRUE, event_onsets = spikes, data_values = data_mat)

  for (strategy in c("round", "nearest", "floor", "ceiling")) {
    result <- suppressWarnings(
      downsample(eeg, target_rate = 256, method = "simple",
                 event_strategy = strategy, verbose = FALSE)
    )
    expect_equal(result$events$onset, which(result$data[1, ] == 1000), info = strategy)
  }
})

test_that("event times stay within half a new-sample step of the true times", {
  # WHAT THIS TESTS:
  # An event can only sit on a kept sample, so with 'round' it can never move
  # by more than half a new-sample step (f / 2 original samples), whatever the
  # factor and whatever the onset (odd, even, anywhere in between).  We check
  # factors 2, 4 and 8 with many onsets, using the original onset_time as the truth.

  onsets <- seq(11, 1000, by = 7)

  for (target in c(256, 128, 64)) {
    f   <- 512 / target
    eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                    with_events = TRUE, event_onsets = onsets)

    result <- downsample(eeg, target_rate = target, event_strategy = "round", verbose = FALSE)

    expect_equal(nrow(result$events), length(onsets), info = paste("factor", f))
    shift <- abs(result$events$onset_time - eeg$events$onset_time)
    expect_true(all(shift <= (f / 2) / 512 + 1e-9), info = paste("factor", f))
  }
})


# ============================================================================
# TEST SUITE 11: Event onset_time update
# ============================================================================

test_that("downsample() updates onset_time to match new time vector", {
  # WHAT THIS TESTS:
  # After rescaling onset sample indices, the function updates the onset_time
  # column by looking up each new onset index in the new time vector:
  #   new_events$onset_time <- new_times[new_events$onset]
  # We verify this directly: for each event, onset_time must equal
  # result$times[result$events$onset].
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(100, 200, 400))
  
  result <- downsample(eeg, target_rate = 256, event_strategy = "round", verbose = FALSE)
  
  expected_times <- result$times[result$events$onset]
  expect_equal(result$events$onset_time, expected_times)
})


# ============================================================================
# TEST SUITE 12: Out-of-bounds event handling
# ============================================================================

test_that("remove_outbound_events = TRUE removes events outside new time range", {
  # WHAT THIS TESTS:
  # With f = 4 and 1024 samples the kept original samples are 1, 5, ... 1021,
  # so the downsampled recording ends at new sample 256 (original 1021).  An
  # event on original sample 1024 lies 3 original samples past that end: the
  # sample closest to it would be new sample 257, which does not exist.  With
  # remove_outbound_events = TRUE that event must be dropped and the early
  # event must stay.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(100, 1024))

  result <- suppressWarnings(
    downsample(eeg, target_rate = 128, event_strategy = "round",
               remove_outbound_events = TRUE, verbose = FALSE)
  )

  expect_equal(nrow(result$events), 1)
  expect_equal(result$events$onset, 26)   # original 100 -> kept original 101 = new sample 26
})

test_that("remove_outbound_events = FALSE clamps out-of-bounds events instead of removing", {
  # WHAT THIS TESTS:
  # When remove_outbound_events = FALSE, events that fall outside the new
  # range are clamped to [1, n_new_samples] rather than removed.  The number
  # of events in the output must equal the number in the input (no events lost).
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(100, 300, 700))
  
  result_remove <- suppressWarnings(
    downsample(eeg, target_rate = 256, remove_outbound_events = TRUE,  verbose = FALSE)
  )
  result_clamp  <- suppressWarnings(
    downsample(eeg, target_rate = 256, remove_outbound_events = FALSE, verbose = FALSE)
  )
  
  # When events are in-bounds, both strategies should retain all events
  expect_equal(nrow(result_clamp$events), nrow(eeg$events))
  expect_true(all(result_clamp$events$onset >= 1))
  expect_true(all(result_clamp$events$onset <= ncol(result_clamp$data)))
})

test_that("downsample() produces a warning when remove_outbound_events = TRUE removes events", {
  # WHAT THIS TESTS:
  # Dropping an event must never be silent: the user has to be told.  We use
  # the same past-the-end event as above (original sample 1024 with f = 4) and
  # check that the warning says why.  Events that are well inside the
  # recording must not trigger it.

  eeg_edge <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                       with_events = TRUE, event_onsets = c(100, 1024))
  eeg_ok   <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                       with_events = TRUE, event_onsets = c(100, 900))

  expect_warning(
    downsample(eeg_edge, target_rate = 128,
               remove_outbound_events = TRUE, verbose = FALSE),
    regexp = "outside the downsampled time range"
  )

  expect_no_warning(
    downsample(eeg_ok, target_rate = 128,
               remove_outbound_events = TRUE, verbose = FALSE)
  )
})

test_that("remove_outbound_events = FALSE clamps a past-the-end event onto the last sample", {
  # WHAT THIS TESTS:
  # Same past-the-end event as above.  With remove_outbound_events = FALSE the
  # event is kept and pulled back onto the last new sample (256), with a
  # warning that says so.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE, event_onsets = c(100, 1024))

  expect_warning(
    result <- downsample(eeg, target_rate = 128, event_strategy = "round",
                         remove_outbound_events = FALSE, verbose = FALSE),
    regexp = "Clamping"
  )

  expect_equal(nrow(result$events), 2)
  expect_equal(result$events$onset[2], ncol(result$data))
  expect_equal(result$events$onset_time[2], tail(result$times, 1))
})


# ============================================================================
# TEST SUITE 13: EEG with no events
# ============================================================================

test_that("downsample() handles EEG objects with no events gracefully", {
  # WHAT THIS TESTS:
  # An EEG object may have an empty events data frame (0 rows).  The function
  # must not crash when iterating over events, and must return an eeg object
  # with an empty events data frame.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024, with_events = FALSE)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_s3_class(result, "eeg")
  expect_equal(nrow(result$events), 0)
})


# ============================================================================
# TEST SUITE 14: Method = "simple"
# ============================================================================

test_that("downsample() with method = 'simple' produces a warning about aliasing", {
  # WHAT THIS TESTS:
  # The "simple" method skips anti-aliasing filtering, which is scientifically
  # inadvisable.  The function must warn the user about potential spectral
  # aliasing so they can make an informed decision.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  
  expect_warning(
    downsample(eeg, target_rate = 256, method = "simple", verbose = FALSE),
    regexp = "aliasing"
  )
})

test_that("downsample() with method = 'simple' returns correct dimensions", {
  # WHAT THIS TESTS:
  # Even though 'simple' skips filtering, the decimation step (select every
  # nth sample) must still produce correctly-sized output.  With factor 2
  # and 1024 samples, we expect 512 output samples and the correct new rate.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  
  result <- suppressWarnings(
    downsample(eeg, target_rate = 256, method = "simple", verbose = FALSE)
  )
  
  expect_equal(result$sampling_rate, 256)
  expected_samples <- length(seq(1, 1024, by = 2))
  expect_equal(ncol(result$data), expected_samples)
})

test_that("downsample() with method = 'simple' selects correct data samples", {
  # WHAT THIS TESTS:
  # Simple downsampling is purely selecting every nth column of the data matrix.
  # We can verify this exactly: result$data must equal eeg$data[, seq(1, n, by=factor)].
  # This test uses a small, deterministic data matrix to make verification easy.
  
  n_tp <- 100
  set.seed(1)
  data_mat <- matrix(seq_len(3 * n_tp), nrow = 3, ncol = n_tp)
  eeg <- new_eeg(
    data          = data_mat,
    channels      = c("A", "B", "C"),
    sampling_rate = 200,
    times         = seq(0, (n_tp - 1) / 200, length.out = n_tp)
  )
  
  result <- suppressWarnings(
    downsample(eeg, target_rate = 100, method = "simple", verbose = FALSE)
  )
  
  expected_data <- data_mat[, seq(1, n_tp, by = 2), drop = FALSE]
  expect_equal(result$data, expected_data)
})


# ============================================================================
# TEST SUITE 15: Method = "decimate" (anti-aliasing filter applied)
# ============================================================================

test_that("downsample() with method = 'decimate' returns an eeg object", {
  # WHAT THIS TESTS:
  # The 'decimate' method invokes the signal package's Butterworth filter and
  # filtfilt.  If all dependencies are available, the function must complete
  # successfully and return a proper eeg object.
  
  skip_if_not_installed("signal")
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  
  expect_s3_class(result, "eeg")
})

test_that("downsample() with method = 'decimate' returns correct sampling rate", {
  # WHAT THIS TESTS:
  # After anti-aliasing filtering and decimation, the sampling_rate field of
  # the returned object must reflect the new rate (256 Hz).
  
  skip_if_not_installed("signal")
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  
  expect_equal(result$sampling_rate, 256)
})

test_that("downsample() with method = 'decimate' returns correct number of samples", {
  # WHAT THIS TESTS:
  # After filtering, the decimation step still selects every nth sample.
  # The number of output samples must match seq(1, n, by=factor).
  
  skip_if_not_installed("signal")
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  
  expected_n <- length(seq(1, 1024, by = 2))
  expect_equal(ncol(result$data), expected_n)
})

test_that("downsample() 'decimate' preserves low-frequency signals better than 'simple'", {
  # WHAT THIS TESTS:
  # Anti-aliasing filtering in 'decimate' should not distort frequencies well
  # below the cutoff.  We create a pure sine wave at 5 Hz (well below the
  # 115 Hz cutoff when going 512->256 Hz) and verify the RMS energy is nearly
  # identical between filtered and unfiltered methods for low-frequency content.
  
  skip_if_not_installed("signal")
  
  sr  <- 512
  n   <- 2048
  t   <- seq(0, (n - 1) / sr, length.out = n)
  sine_5hz <- sin(2 * pi * 5 * t)
  data_mat <- matrix(sine_5hz, nrow = 1, ncol = n)
  
  eeg <- new_eeg(data = data_mat, channels = "Cz", sampling_rate = sr, times = t)
  
  result_decimate <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  result_simple   <- suppressWarnings(
    downsample(eeg, target_rate = 256, method = "simple", verbose = FALSE)
  )
  
  rms_dec <- sqrt(mean(result_decimate$data^2))
  rms_sim <- sqrt(mean(result_simple$data^2))
  
  # Both should have RMS close to 1/sqrt(2) ≈ 0.707 (RMS of a unit sine wave)
  expect_equal(rms_dec, rms_sim, tolerance = 0.05)
})


# ============================================================================
# TEST SUITE 16: Non-integer downsampling factor warning
# ============================================================================

test_that("downsample() warns when downsampling factor is not an integer", {
  # WHAT THIS TESTS:
  # Going from 512 Hz to 200 Hz gives factor = 2.56, which is not an integer.
  # The function must emit a warning and proceed by rounding the factor to 3,
  # producing an actual rate of 512/3 ≈ 170.67 Hz instead of 200 Hz.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1536)
  
  expect_warning(
    downsample(eeg, target_rate = 200, verbose = FALSE),
    regexp = "not an integer"
  )
})

test_that("downsample() rounds non-integer factor and adjusts actual rate accordingly", {
  # WHAT THIS TESTS:
  # After the non-integer factor warning, the function rounds to the nearest
  # integer factor and sets actual_target_rate = sampling_rate / rounded_factor.
  # For 512 Hz -> 200 Hz: factor = 2.56, rounded to 3,
  # actual_rate = 512 / 3 ≈ 170.67 Hz.
  # We verify the output sampling_rate matches the rounded computation.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1536)
  
  result <- suppressWarnings(
    downsample(eeg, target_rate = 200, verbose = FALSE)
  )
  
  rounded_factor  <- round(512 / 200)  # = 3
  expected_rate   <- 512 / rounded_factor  # ≈ 170.67
  
  expect_equal(result$sampling_rate, expected_rate)
})


# ============================================================================
# TEST SUITE 17: Preprocessing history
# ============================================================================

test_that("downsample() appends exactly one entry to preprocessing_history", {
  # WHAT THIS TESTS:
  # Every call to downsample() must log what it did in the preprocessing_history
  # list.  Starting from an empty history, the result must contain exactly one
  # entry.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(length(result$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
})

test_that("downsample() history entry contains old and new sampling rates", {
  # WHAT THIS TESTS:
  # The history string is constructed as:
  #   "Downsampled from <old> Hz to <new> Hz using '<method>' method ..."
  # We check that both Hz values appear in the logged string so the history
  # is informative enough for reproducibility audits.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  history_entry <- result$preprocessing_history[[length(result$preprocessing_history)]]
  
  expect_true(grepl("512", history_entry))
  expect_true(grepl("256", history_entry))
})

test_that("downsample() stacks history entries on a pre-existing history", {
  # WHAT THIS TESTS:
  # If the eeg object already has preprocessing history (e.g., from a prior
  # filtering step), downsample() must append to it, not overwrite it.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$preprocessing_history <- list("Prior step: band-pass filter applied")
  
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(length(result$preprocessing_history), 2)
  expect_equal(result$preprocessing_history[[1]], "Prior step: band-pass filter applied")
})


# ============================================================================
# TEST SUITE 18: Return value structure
# ============================================================================

test_that("downsample() returns an object of class 'eeg'", {
  # WHAT THIS TESTS:
  # The returned object must be an S3 instance of class 'eeg', compatible
  # with all other eeganalysis functions that dispatch on this class.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_s3_class(result, "eeg")
})

test_that("downsample() result has all required eeg fields", {
  # WHAT THIS TESTS:
  # The returned eeg object must contain all eight standard fields:
  # data, channels, sampling_rate, times, events, metadata, reference,
  # preprocessing_history.  Missing fields would break downstream functions.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expected_fields <- c("data", "channels", "sampling_rate", "times",
                       "events", "metadata", "reference", "preprocessing_history")
  
  expect_true(all(expected_fields %in% names(result)))
})

test_that("downsample() preserves metadata from the original eeg object", {
  # WHAT THIS TESTS:
  # Metadata (subject ID, session, etc.) is experiment-level information that
  # must be carried through unmodified.  Downsampling must not touch it.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$metadata <- list(subject = "S01", session = "A", date = "2026-02-18")
  
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(result$metadata, eeg$metadata)
})

test_that("downsample() preserves the reference field", {
  # WHAT THIS TESTS:
  # The reference scheme (e.g., "average", "CMS/DRL") describes how the EEG
  # was re-referenced.  It must be passed through unchanged.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$reference <- "average"
  
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(result$reference, "average")
})

test_that("downsample() preserves bads from the original eeg object", {
  # WHAT THIS TESTS:
  # The bad-channel list is a shared exclude list that later steps (re-
  # referencing, ICA, interpolation) read from eeg$bads.  Downsampling
  # rebuilds the eeg object, so it must hand the list back unchanged rather
  # than resetting it to "no bad channels".

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$bads <- "Ch2"

  result <- downsample(eeg, target_rate = 256, verbose = FALSE)

  expect_equal(result$bads, "Ch2")
})

test_that("downsample() preserves the montage from the original eeg object", {
  # WHAT THIS TESTS:
  # Electrode positions do not depend on the sampling rate, so an attached
  # montage must come back exactly as it went in.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$montage <- create_montage()

  result <- downsample(eeg, target_rate = 256, verbose = FALSE)

  expect_equal(result$montage, eeg$montage)
})

test_that("downsample() preserves annotations from the original eeg object", {
  # WHAT THIS TESTS:
  # Annotation onsets and durations are in seconds, not sample indices, so
  # (unlike events) they need no rescaling for the new sampling rate and must
  # come back exactly as they went in.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$annotations <- data.frame(
    onset       = c(0.5, 1.25),
    duration    = c(0.25, 0.5),
    description = c("BAD_muscle", "BAD_peak"),
    channel     = c(NA_character_, "Ch1"),
    stringsAsFactors = FALSE
  )

  result <- downsample(eeg, target_rate = 256, verbose = FALSE)

  expect_equal(result$annotations, eeg$annotations)
})

test_that("downsample() preserves the type of a derived channel", {
  # WHAT THIS TESTS:
  # channel_types is derived from the channel names when an eeg object is built,
  # but a channel made by set_bipolar_reference() is typed "eog" directly
  # (the name "VEOG" would not be classified that way).  Downsampling rebuilds
  # the eeg object, so it must put the original types back rather than turn
  # "VEOG" into an "eeg" channel.

  eeg <- make_eeg(n_channels = 4, sampling_rate = 512, n_timepoints = 1024)
  eeg <- set_bipolar_reference(eeg, anode = "Ch3", cathode = "Ch4", ch_name = "VEOG")

  result <- downsample(eeg, target_rate = 256, verbose = FALSE)

  expect_equal(result$channel_types, eeg$channel_types)
  expect_equal(result$channel_types[result$channels == "VEOG"], "eog")
})

test_that("downsample() still derives channel types when the input has none", {
  # WHAT THIS TESTS:
  # An object without a channel_types field (hand-built, or saved before the
  # field existed) must still come back with the types new_eeg() derives from
  # the channel names, not with the field deleted.

  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  eeg$channel_types <- NULL

  result <- downsample(eeg, target_rate = 256, verbose = FALSE)

  expect_equal(result$channel_types, rep("eeg", 3))
})

test_that("downsample() data matrix has correct storage type (numeric)", {
  # WHAT THIS TESTS:
  # The output data matrix must be a numeric (double) matrix, consistent with
  # the convention enforced by new_eeg().  Integer-type matrices could cause
  # subtle precision issues in downstream spectral analyses.
  
  eeg    <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_true(is.numeric(result$data))
  expect_true(is.matrix(result$data))
})


# ============================================================================
# TEST SUITE 19: Idempotence and reproducibility
# ============================================================================

test_that("downsample() produces identical results on repeated calls with same seed", {
  # WHAT THIS TESTS:
  # The function is deterministic given the same input data.  Calling it twice
  # on the same eeg object must produce bit-for-bit identical output.
  # (The signal package's filtfilt is deterministic for fixed inputs.)
  
  skip_if_not_installed("signal")
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024)
  
  result1 <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  result2 <- downsample(eeg, target_rate = 256, method = "decimate", verbose = FALSE)
  
  expect_equal(result1$data, result2$data)
  expect_equal(result1$times, result2$times)
})

test_that("downsample() does not modify the original eeg object (no side-effects)", {
  # WHAT THIS TESTS:
  # R passes objects by value, so the original eeg must be unchanged after
  # calling downsample().  The original sampling_rate, data dimensions, and
  # event table must be exactly as before.
  
  eeg     <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                      with_events = TRUE, event_onsets = c(100, 300))
  orig_sr <- eeg$sampling_rate
  orig_nc <- ncol(eeg$data)
  orig_ev <- nrow(eeg$events)
  
  downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(eeg$sampling_rate, orig_sr)
  expect_equal(ncol(eeg$data),    orig_nc)
  expect_equal(nrow(eeg$events),  orig_ev)
})


# ============================================================================
# TEST SUITE 20: Edge cases
# ============================================================================

test_that("downsample() handles single-channel EEG correctly", {
  # WHAT THIS TESTS:
  # A single-channel EEG (nrow=1) should work without issues.  Matrix
  # operations can behave differently for 1-row matrices (drop=FALSE is
  # needed), so we verify the output has exactly 1 channel.
  
  eeg <- make_eeg(n_channels = 1, sampling_rate = 512, n_timepoints = 1024)
  
  result <- downsample(eeg, target_rate = 256, verbose = FALSE)
  
  expect_equal(nrow(result$data), 1)
  expect_s3_class(result, "eeg")
})

test_that("downsample() handles large factor downsampling (factor = 8)", {
  # WHAT THIS TESTS:
  # A factor-8 reduction (512 -> 64 Hz) is aggressive but common for
  # slow-wave sleep research.  The output must have 1/8 the original samples
  # and report 64 Hz as the new rate.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 2048)
  
  result <- downsample(eeg, target_rate = 64, verbose = FALSE)
  
  expect_equal(result$sampling_rate, 64)
  expected_n <- length(seq(1, 2048, by = 8))
  expect_equal(ncol(result$data), expected_n)
})

test_that("downsample() works correctly with multiple events at different onsets", {
  # WHAT THIS TESTS:
  # A realistic scenario: many events at various positions throughout the
  # recording.  All events must be correctly rescaled and have valid onsets
  # within the new data range.
  
  eeg <- make_eeg(sampling_rate = 512, n_timepoints = 1024,
                  with_events = TRUE,
                  event_onsets = c(50, 100, 200, 300, 400, 512, 700, 900))
  
  result <- suppressWarnings(
    downsample(eeg, target_rate = 256, event_strategy = "round", verbose = FALSE)
  )
  
  # All retained events must have valid indices
  expect_true(all(result$events$onset >= 1))
  expect_true(all(result$events$onset <= ncol(result$data)))
})