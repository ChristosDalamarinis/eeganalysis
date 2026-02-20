# ============================================================================
#                         Test File for epoch2.R
# ============================================================================
#
# This test file provides comprehensive testing for the epoch2.R script, which
# contains the core EEG epoching pipeline for the eeganalysis package.
#
# Functions tested:
#   1. inspect_triggers()  - Diagnostic summary of event triggers in EEG data
#   2. epoch_eeg()         - Segments continuous EEG data into time-locked epochs
#   3. plot_epochs()       - Visualises epoched EEG data
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                         SHARED HELPER FUNCTIONS
# ============================================================================
#
# These helpers create minimal but realistic mock objects that are reused
# throughout the test suites, keeping individual tests focused and concise.
# ----------------------------------------------------------------------------

#' Build a minimal mock `eeg` object.
#'
#' @param n_channels  Number of EEG channels (rows in data matrix).
#' @param n_timepoints Number of time samples (columns in data matrix).
#' @param sr Sampling rate in Hz.
#' @param event_onsets Integer sample indices where events occur.
#' @param event_types  Corresponding trigger codes (stored as character,
#'   matching the real new_eeg() output convention).
#' @param data_values  Optional numeric matrix to use as the data payload.
#'   If NULL, random Gaussian noise is generated.
#' @return An S3 object of class "eeg".
make_mock_eeg <- function(n_channels   = 2,
                          n_timepoints = 1000,
                          sr           = 256,
                          event_onsets = c(200L, 500L, 800L),
                          event_types  = c("1", "2", "1"),
                          data_values  = NULL) {
  
  if (is.null(data_values)) {
    data_values <- matrix(rnorm(n_channels * n_timepoints),
                          nrow = n_channels, ncol = n_timepoints)
  }
  
  channels <- paste0("Ch", seq_len(n_channels))
  times    <- seq(0, (n_timepoints - 1) / sr, length.out = n_timepoints)
  
  # Guard against empty event inputs
  if (length(event_onsets) == 0) {
    events <- data.frame(
      onset        = integer(0),
      onset_time   = numeric(0),
      type         = character(0),
      description  = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    events <- data.frame(
      onset        = as.integer(event_onsets),
      onset_time   = event_onsets / sr,
      type         = as.character(event_types),
      description  = paste0("Trigger: ", event_types),
      stringsAsFactors = FALSE
    )
  }
  
  structure(
    list(
      data                  = data_values,
      channels              = channels,
      sampling_rate         = sr,
      times                 = times,
      events                = events,
      metadata              = list(),
      reference             = "original",
      preprocessing_history = list()
    ),
    class = "eeg"
  )
}

#' Build a minimal mock `eeg_epochs` object for plot_epochs() tests.
#'
#' @param n_channels   Number of channels.
#' @param n_samples    Number of time samples per epoch.
#' @param n_epochs     Number of epochs (trials).
#' @param tmin         Epoch start time (s).
#' @param tmax         Epoch end time (s).
#' @return An S3 object of class "eeg_epochs".
make_mock_epochs <- function(n_channels = 2,
                             n_samples  = 129,
                             n_epochs   = 5,
                             tmin       = -0.1,
                             tmax       =  0.4) {
  
  data_array <- array(rnorm(n_channels * n_samples * n_epochs),
                      dim = c(n_channels, n_samples, n_epochs))
  dimnames(data_array)[[1]] <- paste0("Ch", seq_len(n_channels))
  
  events <- data.frame(
    onset       = as.integer(seq(200, by = 300, length.out = n_epochs)),
    onset_time  = seq(200, by = 300, length.out = n_epochs) / 256,
    type        = rep(1L, n_epochs),
    description = rep("Trigger: 1", n_epochs),
    epoch_id    = seq_len(n_epochs),
    stringsAsFactors = FALSE
  )
  
  structure(
    list(
      data             = data_array,
      channels         = paste0("Ch", seq_len(n_channels)),
      times            = seq(tmin, tmax, length.out = n_samples),
      events           = events,
      sampling_rate    = 256,
      tmin             = tmin,
      tmax             = tmax,
      baseline         = c(-0.1, 0),
      baseline_method  = "mean",
      n_epochs         = n_epochs,
      rejected         = rep(FALSE, n_epochs),
      rejection_log    = data.frame(epoch_id   = integer(0),
                                    event_type = character(0),
                                    event_time = numeric(0),
                                    reason     = character(0),
                                    stringsAsFactors = FALSE),
      metadata         = list()
    ),
    class = "eeg_epochs"
  )
}


# ============================================================================
# TEST SUITE 1: inspect_triggers() — Trigger Inspection
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Input validation — wrong class
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: inspect_triggers() must only accept objects of class "eeg".
# Passing any other object (list, data.frame, etc.) must throw an error with
# a clear message before touching any downstream logic.
test_that("inspect_triggers rejects non-eeg objects", {
  expect_error(
    inspect_triggers(list(a = 1)),
    "Input must be an object of class 'eeg'"
  )
  expect_error(
    inspect_triggers(data.frame(onset = 1)),
    "Input must be an object of class 'eeg'"
  )
  expect_error(
    inspect_triggers("not_an_eeg"),
    "Input must be an object of class 'eeg'"
  )
})

# ----------------------------------------------------------------------------
# Test 1.2: No events in EEG object — graceful warning, not a crash
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When an "eeg" object has an empty events data frame,
# inspect_triggers() must return early with an invisible list containing
# n_events_raw = 0, rather than crashing or producing unintelligible output.
test_that("inspect_triggers handles eeg with no events gracefully", {
  eeg_no_events <- make_mock_eeg(event_onsets = integer(0),
                                 event_types  = character(0))
  
  result <- suppressWarnings(
    inspect_triggers(eeg_no_events, plot = FALSE)
  )
  
  expect_equal(result$n_events_raw, 0)
  expect_equal(result$n_events_cleaned, 0)
})

# ----------------------------------------------------------------------------
# Test 1.3: System trigger removal via 0xFFFF bitmask
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The core BioSemi-aware feature of inspect_triggers() is the
# bitwAnd(type, 0xFFFF) mask. Events whose trigger code is PURELY in bits 16-23
# (i.e., bitwAnd(type, 0xFFFF) == 0) must be counted as system triggers and
# removed from events_cleaned. Events with bits 0-15 set must be retained.
#
# 65536 = 2^16 → bitwAnd(65536, 0xFFFF) = 0 (system-only trigger → removed)
# 1     = 2^0  → bitwAnd(1,     0xFFFF) = 1 (experimental trigger → kept)
# 65537 = 2^16 + 2^0 → bitwAnd = 1           (both bits → kept, code becomes 1)
test_that("inspect_triggers removes system triggers via 0xFFFF bitmask", {
  eeg_mixed <- make_mock_eeg(
    n_timepoints = 2000,
    event_onsets = c(200L, 400L, 600L, 800L),
    event_types  = c("65536", "1", "65536", "2")   # two system, two experimental
  )
  
  result <- inspect_triggers(eeg_mixed, mode = "raw", plot = FALSE)
  
  expect_equal(result$n_events_raw, 4)
  expect_equal(result$n_system_removed, 2)
  expect_equal(result$n_events_cleaned, 2)
  # Only type 1 and type 2 remain
  expect_true(all(result$events_cleaned$type %in% c(1, 2)))
})

# ----------------------------------------------------------------------------
# Test 1.4: mode = "raw" returns all experimental triggers unfiltered
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: In "raw" mode, no auto-filtering is applied beyond the
# 0xFFFF bitmask. Consecutive identical events, brief events, and IEI
# violations are NOT removed. The n_events_cleaned count should equal the
# total number of events that pass the bitmask.
test_that("inspect_triggers mode='raw' skips all auto-filtering", {
  # 5 events: all experimental, some consecutive identical, tiny IEI
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(100L, 101L, 102L, 500L, 501L),   # consecutive; IEI < min_iei
    event_types  = c("1",  "1",  "1",  "2",  "2")
  )
  
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  # All 5 pass the mask (none are system triggers)
  expect_equal(result$n_events_cleaned, 5)
  expect_equal(result$n_system_removed, 0)
})

# ----------------------------------------------------------------------------
# Test 1.5: mode = "auto" merges consecutive identical events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: In "auto" mode, consecutive events of the same type
# (indicating a sustained trigger signal rather than distinct events) are
# collapsed to a single event. The function detects type changes and only
# keeps the first occurrence of each run.
test_that("inspect_triggers mode='auto' merges consecutive identical events", {
  # Three consecutive type-1 events, then one type-2 event = 2 unique events
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(200L, 202L, 204L, 800L),
    event_types  = c("1",  "1",  "1",  "2")
  )
  
  result <- inspect_triggers(eeg, mode = "auto", plot = FALSE)
  
  # After merging: only 1 + 1 = 2 events should remain
  expect_lte(result$n_events_cleaned, 2)
})

# ----------------------------------------------------------------------------
# Test 1.6: mode = "auto" removes brief events below trigger_threshold
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: After merging, "auto" mode removes events whose duration
# (defined as time to the next event) is less than trigger_threshold seconds.
# This filters hardware glitches that produce very short spurious triggers.
test_that("inspect_triggers mode='auto' removes brief events", {
  sr <- 1000   # 1 kHz for easy ms arithmetic
  # Event at 500 ms, then another at 501 ms (1 ms apart = 0.001 s < 0.002 s threshold)
  # Then a legitimate event at 600 ms (99 ms later)
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(500L, 501L, 600L),
    event_types  = c("1",  "2",  "1")
  )
  
  result_auto <- inspect_triggers(
    eeg, mode = "auto", trigger_threshold = 0.002, plot = FALSE
  )
  result_raw  <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  # Raw keeps all 3; auto should remove at least one brief event
  expect_lt(result_auto$n_events_cleaned, result_raw$n_events_cleaned)
})

# ----------------------------------------------------------------------------
# Test 1.7: mode = "auto" enforces minimum inter-event interval (min_iei)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: After brief-event removal, "auto" mode additionally drops
# events that occur within min_iei seconds of the preceding event. This guards
# against double-triggers or closely spaced noise bursts.
test_that("inspect_triggers mode='auto' enforces min_iei", {
  sr <- 1000
  # First two events are 5 ms apart (0.005 s) — below default min_iei of 0.01 s
  # Third event is 200 ms after the second — above min_iei
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(100L, 105L, 305L),
    event_types  = c("1",  "2",  "1")
  )
  
  result <- inspect_triggers(eeg, mode = "auto", min_iei = 0.01, plot = FALSE)
  
  # The event at 105 ms should be removed (IEI = 0.005 < 0.01)
  expect_lte(result$n_events_cleaned, 2)
})

# ----------------------------------------------------------------------------
# Test 1.8: exclude_types removes specified trigger codes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: In both "auto" and "manual" modes, the exclude_types
# argument explicitly removes events whose type matches. This allows the user
# to drop known nuisance triggers (e.g., response codes) from the analysis.
test_that("inspect_triggers exclude_types removes specified codes", {
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1",  "99", "1")    # trigger 99 should be excluded
  )
  
  result_with_exclude <- inspect_triggers(
    eeg, mode = "manual", exclude_types = 99, plot = FALSE
  )
  result_without_exclude <- inspect_triggers(
    eeg, mode = "raw", plot = FALSE
  )
  
  # With exclusion: 2 events; without: 3
  expect_equal(result_without_exclude$n_events_cleaned, 3)
  expect_equal(result_with_exclude$n_events_cleaned, 2)
  expect_false(99 %in% result_with_exclude$events_cleaned$type)
})

# ----------------------------------------------------------------------------
# Test 1.9: Return value has the expected structure
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: inspect_triggers() returns an invisible list. This test
# confirms all documented fields are present in the returned object so
# downstream code can rely on them.
test_that("inspect_triggers returns list with all documented fields", {
  eeg <- make_mock_eeg()
  
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  expected_fields <- c("n_events_raw", "n_events_cleaned", "n_system_removed",
                       "event_types", "event_counts",
                       "mean_iei", "median_iei", "event_rate", "events_cleaned")
  expect_true(all(expected_fields %in% names(result)))
})

# ----------------------------------------------------------------------------
# Test 1.10: Return value — events_cleaned is a data frame
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The events_cleaned component of the returned list must be a
# data frame (not a list or matrix) so it can be used directly in downstream
# analysis and CSV export.
test_that("inspect_triggers events_cleaned is a data frame", {
  eeg    <- make_mock_eeg()
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  expect_true(is.data.frame(result$events_cleaned))
})

# ----------------------------------------------------------------------------
# Test 1.11: CSV export writes a readable file
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When export_csv is set to a valid file path, inspect_triggers()
# must write the cleaned events data frame as a CSV. The file should exist and
# be readable as a data frame with the expected number of rows.
test_that("inspect_triggers exports cleaned events to CSV when requested", {
  eeg      <- make_mock_eeg()
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file), add = TRUE)
  
  inspect_triggers(eeg, mode = "raw", plot = FALSE, export_csv = tmp_file)
  
  expect_true(file.exists(tmp_file))
  exported <- read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(is.data.frame(exported))
  expect_gt(nrow(exported), 0)
})

# ----------------------------------------------------------------------------
# Test 1.12: Timing statistics are computed correctly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The returned mean_iei, median_iei, and event_rate values
# must be numeric (not NA) when more than one event is present, and must
# reflect the actual spacing between the cleaned events.
test_that("inspect_triggers computes timing statistics for multiple events", {
  sr  <- 256
  # Two events exactly 1 second apart
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(256L, 512L),    # 1 s apart at sr=256
    event_types  = c("1", "1")
  )
  
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  expect_false(is.na(result$mean_iei))
  expect_false(is.na(result$median_iei))
  expect_false(is.na(result$event_rate))
  expect_true(is.numeric(result$mean_iei))
})

# ----------------------------------------------------------------------------
# Test 1.13: Single event yields NA timing statistics
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With only one event there is no inter-event interval to
# compute, so mean_iei, median_iei, and event_rate must be NA. This prevents
# division-by-zero or misleading statistics.
test_that("inspect_triggers returns NA timing stats for single event", {
  eeg <- make_mock_eeg(
    n_timepoints = 1000,
    event_onsets = c(300L),
    event_types  = c("1")
  )
  
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  expect_true(is.na(result$mean_iei))
  expect_true(is.na(result$median_iei))
})

# ----------------------------------------------------------------------------
# Test 1.14: Mixed experimental and system triggers — counts are correct
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A realistic scenario where some events are pure system
# triggers (bits 16-23 only) and some are experimental (bits 0-15). The raw
# count must equal total input events, and n_system_removed must accurately
# reflect how many were masked out.
test_that("inspect_triggers correctly counts system vs experimental triggers", {
  eeg <- make_mock_eeg(
    n_timepoints = 3000,
    event_onsets = c(200L, 400L, 600L, 800L, 1000L),
    event_types  = c("65536", "65536", "65536", "1", "2")
  )
  
  result <- inspect_triggers(eeg, mode = "raw", plot = FALSE)
  
  expect_equal(result$n_events_raw,     5)
  expect_equal(result$n_system_removed, 3)
  expect_equal(result$n_events_cleaned, 2)
})


# ============================================================================
# TEST SUITE 2: epoch_eeg() — Core Epoching Function
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Input validation — non-eeg class is rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: epoch_eeg() must validate that its first argument is an
# object of class "eeg". Passing any other type must raise a descriptive error
# before any computation occurs.
test_that("epoch_eeg rejects non-eeg input", {
  expect_error(
    epoch_eeg(list(data = matrix(1:10, 2, 5))),
    "Input must be an object of class 'eeg'"
  )
  expect_error(
    epoch_eeg("not_an_eeg"),
    "Input must be an object of class 'eeg'"
  )
})

# ----------------------------------------------------------------------------
# Test 2.2: Input validation — eeg with no events is rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: epoch_eeg() requires events to define the locking points.
# An eeg object with an empty events table must produce an informative error
# telling the user to run inspect_triggers() first.
test_that("epoch_eeg rejects eeg object with no events", {
  eeg_no_events <- make_mock_eeg(event_onsets = integer(0),
                                 event_types  = character(0))
  expect_error(
    epoch_eeg(eeg_no_events, verbose = FALSE),
    "No events found"
  )
})

# ----------------------------------------------------------------------------
# Test 2.3: Input validation — tmin >= tmax is rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: An epoch window where the start time equals or exceeds the
# end time is physically meaningless. epoch_eeg() must detect this and error
# out with a message referencing both the offending tmin and tmax values.
test_that("epoch_eeg rejects tmin >= tmax", {
  eeg <- make_mock_eeg()
  expect_error(epoch_eeg(eeg, tmin = 0.5, tmax = 0.2, verbose = FALSE),
               "tmin")
  expect_error(epoch_eeg(eeg, tmin = 0.3, tmax = 0.3, verbose = FALSE),
               "tmin")
})

# ----------------------------------------------------------------------------
# Test 2.4: Input validation — invalid baseline is rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The baseline argument must be either NULL or a two-element
# numeric vector with start < end. Malformed inputs (wrong length, reversed
# order, or baseline extending outside the epoch window) must each be caught
# with a dedicated error message.
test_that("epoch_eeg rejects invalid baseline specifications", {
  eeg <- make_mock_eeg()
  
  # baseline[1] >= baseline[2]
  expect_error(
    epoch_eeg(eeg, tmin = -0.2, tmax = 0.5,
              baseline = c(0, -0.1), verbose = FALSE),
    "baseline"
  )
  
  # baseline extends before tmin
  expect_error(
    epoch_eeg(eeg, tmin = -0.1, tmax = 0.5,
              baseline = c(-0.3, 0), verbose = FALSE),
    "Baseline window must be within epoch window"
  )
  
  # baseline extends after tmax
  expect_error(
    epoch_eeg(eeg, tmin = -0.2, tmax = 0.5,
              baseline = c(-0.2, 0.8), verbose = FALSE),
    "Baseline window must be within epoch window"
  )
})

# ----------------------------------------------------------------------------
# Test 2.5: baseline_method = "none" disables baseline correction
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When baseline_method is "none", the function internally
# sets baseline to NULL, so no baseline correction is performed. The returned
# object should reflect this: baseline_method should be "none" and the baseline
# field should be NULL.
test_that("epoch_eeg disables baseline when baseline_method='none'", {
  eeg    <- make_mock_eeg()
  epochs <- epoch_eeg(eeg,
                      tmin             = -0.1, tmax = 0.4,
                      baseline         = c(-0.1, 0),
                      baseline_method  = "none",
                      verbose          = FALSE)
  
  expect_equal(epochs$baseline_method, "none")
  expect_null(epochs$baseline)
})

# ----------------------------------------------------------------------------
# Test 2.6: System triggers are filtered via 0xFFFF bitmask before epoching
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: epoch_eeg() applies the same BioSemi-aware bitmask as
# inspect_triggers(). Events whose type is purely in bits 16-23 (e.g., 65536)
# must not be used as epoch centres — only experimental events survive the mask.
# If ALL events are system triggers, the function must stop with an informative
# message.
test_that("epoch_eeg strips system triggers and errors when none remain", {
  eeg_system_only <- make_mock_eeg(
    event_onsets = c(200L, 500L),
    event_types  = c("65536", "65536")   # both purely system triggers
  )
  
  expect_error(
    epoch_eeg(eeg_system_only, events = "all", verbose = FALSE),
    "No experimental triggers found"
  )
})

# ----------------------------------------------------------------------------
# Test 2.7: events = "all" epochs around every experimental trigger
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When events = "all", epoch_eeg() must use every event that
# survives the 0xFFFF bitmask. The resulting n_epochs must equal the number of
# experimental events in the input (assuming no boundary rejections).
test_that("epoch_eeg with events='all' extracts epochs for every valid event", {
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1", "2", "3")
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(epochs$n_epochs, 3)
})

# ----------------------------------------------------------------------------
# Test 2.8: Selecting a specific event type
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When events is a numeric vector of trigger codes, only
# events matching those codes are used as epoch centres. The returned n_epochs
# must equal exactly the number of matching events (again, excluding boundary
# rejections).
test_that("epoch_eeg selects only specified event types", {
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1", "2", "1")
  )
  
  # Only epoch around trigger code 1 (appears at onsets 200 and 800)
  epochs <- epoch_eeg(eeg, events = 1,
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(epochs$n_epochs, 2)
})

# ----------------------------------------------------------------------------
# Test 2.9: No matching events — informative error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When the requested event codes don't match any experimental
# triggers in the data, epoch_eeg() must raise an error explaining which codes
# were requested but not found, rather than returning silently empty output.
test_that("epoch_eeg errors when requested event types are absent", {
  eeg <- make_mock_eeg(event_types = c("1", "2", "1"))
  
  expect_error(
    epoch_eeg(eeg, events = 99, verbose = FALSE),
    "No events matched"
  )
})

# ----------------------------------------------------------------------------
# Test 2.10: Output has class "eeg_epochs"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The returned object must carry the S3 class "eeg_epochs"
# so that downstream methods (plot_epochs(), average_epochs(), etc.) can
# perform class-based dispatch on it.
test_that("epoch_eeg returns an object of class 'eeg_epochs'", {
  eeg    <- make_mock_eeg(n_timepoints = 2000)
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  expect_s3_class(epochs, "eeg_epochs")
})

# ----------------------------------------------------------------------------
# Test 2.11: Output has all documented fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The "eeg_epochs" object must contain every field listed in
# the function's @return documentation. Missing fields would break downstream
# code that accesses these components by name.
test_that("epoch_eeg output contains all documented fields", {
  eeg    <- make_mock_eeg(n_timepoints = 2000)
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expected_fields <- c("data", "channels", "times", "events", "sampling_rate",
                       "tmin", "tmax", "baseline", "baseline_method",
                       "n_epochs", "rejected", "rejection_log", "metadata")
  expect_true(all(expected_fields %in% names(epochs)))
})

# ----------------------------------------------------------------------------
# Test 2.12: Data array has correct dimensions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The data component of the returned "eeg_epochs" object must
# be a 3-D array with dimensions [n_channels × n_samples_epoch × n_epochs].
# n_samples_epoch is derived from tmin/tmax and sampling_rate via:
#   smin = round(tmin * sr), smax = round(tmax * sr)
#   n_samples_epoch = smax - smin + 1
test_that("epoch_eeg data array has correct [channels × samples × epochs] dimensions", {
  sr           <- 256
  tmin         <- -0.1
  tmax         <-  0.4
  n_channels   <- 2
  n_epochs_exp <- 3
  
  smin             <- round(tmin * sr)    # -26
  smax             <- round(tmax * sr)    # 102
  n_samples_epoch  <- smax - smin + 1    # 129
  
  eeg <- make_mock_eeg(
    n_channels   = n_channels,
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1", "2", "1")
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(dim(epochs$data), c(n_channels, n_samples_epoch, n_epochs_exp))
})

# ----------------------------------------------------------------------------
# Test 2.13: Channel names are preserved in the data array
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: epoch_eeg() sets dimnames on the first dimension of the
# data array from eeg_obj$channels. This allows downstream code to index
# channels by name rather than position.
test_that("epoch_eeg preserves channel names as first-dimension dimnames", {
  eeg    <- make_mock_eeg(n_channels = 3, n_timepoints = 2000,
                          event_onsets = c(300L, 700L), event_types = c("1", "1"))
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(dimnames(epochs$data)[[1]], c("Ch1", "Ch2", "Ch3"))
})

# ----------------------------------------------------------------------------
# Test 2.14: Times vector has correct range and length
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The times vector in the returned object must span [tmin, tmax]
# and have length equal to n_samples_epoch, so it correctly maps sample indices
# to seconds for plotting and analysis.
test_that("epoch_eeg times vector spans [tmin, tmax] with correct length", {
  sr   <- 256
  tmin <- -0.1
  tmax <-  0.4
  
  smin            <- round(tmin * sr)
  smax            <- round(tmax * sr)
  n_samples_epoch <- smax - smin + 1
  
  eeg    <- make_mock_eeg(n_timepoints = 2000)
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(length(epochs$times), n_samples_epoch)
  # First and last time points must match tmin and tmax
  expect_equal(epochs$times[1],                   tmin)
  expect_equal(epochs$times[n_samples_epoch],      tmax)
})

# ----------------------------------------------------------------------------
# Test 2.15: Epochs extending beyond the data boundaries are rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When an event is too close to the start or end of the
# recording such that the requested epoch window would fall outside the data
# array, the epoch must be marked as rejected (valid_epochs[i] = FALSE) and
# logged in rejection_log with reason "Extends beyond data boundaries".
# The event must NOT appear in the returned epochs$events data frame.
test_that("epoch_eeg rejects epochs that extend beyond data boundaries", {
  sr   <- 256
  tmin <- -0.2   # needs 51 samples pre-stimulus
  tmax <-  0.4
  
  smin <- round(tmin * sr)   # -51
  
  # Place one event too early for the pre-stimulus window to fit
  # onset=10 → epoch_start = 10 + (-51) = -41 < 1 → must be rejected
  # Place another event safely in the middle
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(10L, 600L),
    event_types  = c("1", "1")
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL, verbose = FALSE)
  
  # Only the safe epoch should be extracted
  expect_equal(epochs$n_epochs, 1)
  
  # Rejection log must record the boundary violation
  expect_gt(nrow(epochs$rejection_log), 0)
  expect_true(any(grepl("boundaries", epochs$rejection_log$reason,
                        ignore.case = TRUE)))
})

# ----------------------------------------------------------------------------
# Test 2.16: Amplitude rejection removes epochs above the threshold
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When reject_threshold is set, any epoch whose peak
# absolute amplitude exceeds that value must be flagged as rejected. The
# rejection_log must record the reason and the offending peak value. Only
# clean epochs must appear in epochs$data.
test_that("epoch_eeg rejects epochs with amplitude exceeding threshold", {
  sr   <- 256
  tmin <- -0.1
  tmax <-  0.4
  
  smin <- round(tmin * sr)   # -26
  smax <- round(tmax * sr)   #  102
  
  n_samples <- smax - smin + 1   # 129
  
  # Build data so that the epoch centred at onset=300 contains an artifact
  data_clean    <- matrix(rnorm(2 * 2000, mean = 0, sd = 1), nrow = 2)
  # Inject a large spike into the window of the first event (onset=300)
  epoch_start   <- 300 + smin   # 300 - 26 = 274
  epoch_end     <- 300 + smax   # 300 + 102 = 402
  data_clean[1, epoch_start:epoch_end] <- 1000   # far above any threshold
  
  eeg <- make_mock_eeg(
    n_channels   = 2,
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(300L, 700L),
    event_types  = c("1", "1"),
    data_values  = data_clean
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL,
                      reject_threshold = 200,
                      verbose = FALSE)
  
  # Epoch 1 (the artifact) must be rejected; epoch 2 must survive
  expect_equal(epochs$n_epochs, 1)
  expect_true(any(grepl("Amplitude", epochs$rejection_log$reason)))
})

# ----------------------------------------------------------------------------
# Test 2.17: All epochs rejected raises an error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If every candidate epoch is rejected (either by boundary
# or amplitude criteria), the function must stop with a clear error message
# rather than returning an empty or malformed object.
test_that("epoch_eeg stops when all epochs are rejected", {
  sr   <- 256
  tmin <- -0.2
  tmax <-  0.4
  
  # Both events too close to the data boundaries
  eeg <- make_mock_eeg(
    n_timepoints = 500,
    sr           = sr,
    event_onsets = c(5L, 498L),
    event_types  = c("1", "1")
  )
  
  expect_error(
    epoch_eeg(eeg, events = "all",
              tmin = tmin, tmax = tmax,
              baseline = NULL, verbose = FALSE),
    "All epochs were rejected"
  )
})

# ----------------------------------------------------------------------------
# Test 2.18: Mean baseline correction makes the baseline window mean ~0
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The fundamental correctness of mean baseline correction.
# For each channel, epoch_eeg() subtracts the channel's mean amplitude in the
# baseline window from the entire epoch. After correction the mean of the
# baseline window for each channel should be numerically zero (within floating-
# point tolerance).
test_that("epoch_eeg mean baseline correction drives baseline mean to zero", {
  sr    <- 256
  tmin  <- -0.2
  tmax  <-  0.5
  n_ch  <- 2
  
  # Use a constant non-zero signal so the baseline effect is predictable
  # Channel 1: all 5.0; Channel 2: all -3.0
  n_tp <- 2000
  data_const        <- matrix(0, nrow = n_ch, ncol = n_tp)
  data_const[1, ]   <- 5.0
  data_const[2, ]   <- -3.0
  
  eeg <- make_mock_eeg(
    n_channels   = n_ch,
    n_timepoints = n_tp,
    sr           = sr,
    event_onsets = c(400L, 900L),
    event_types  = c("1", "1"),
    data_values  = data_const
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline        = c(-0.2, 0),
                      baseline_method = "mean",
                      verbose         = FALSE)
  
  # Identify baseline sample indices in the epochs$times vector
  bl_idx <- which(epochs$times >= -0.2 & epochs$times <= 0)
  
  for (ep in seq_len(epochs$n_epochs)) {
    for (ch in seq_len(n_ch)) {
      bl_mean <- mean(epochs$data[ch, bl_idx, ep])
      expect_lt(abs(bl_mean), 1e-9,
                label = sprintf("baseline mean for ch%d epoch%d", ch, ep))
    }
  }
})

# ----------------------------------------------------------------------------
# Test 2.19: Median baseline correction makes the baseline window median ~0
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Same logic as Test 2.18 but for "median" baseline_method.
# epoch_eeg() subtracts the channel's median in the baseline window; the
# median of the corrected baseline window must be numerically zero.
test_that("epoch_eeg median baseline correction drives baseline median to zero", {
  sr   <- 256
  tmin <- -0.2
  tmax <-  0.5
  n_ch <- 2
  n_tp <- 2000
  
  data_const        <- matrix(0, nrow = n_ch, ncol = n_tp)
  data_const[1, ]   <- 7.0
  data_const[2, ]   <- -2.5
  
  eeg <- make_mock_eeg(
    n_channels   = n_ch,
    n_timepoints = n_tp,
    sr           = sr,
    event_onsets = c(400L, 900L),
    event_types  = c("1", "1"),
    data_values  = data_const
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline        = c(-0.2, 0),
                      baseline_method = "median",
                      verbose         = FALSE)
  
  bl_idx <- which(epochs$times >= -0.2 & epochs$times <= 0)
  
  for (ep in seq_len(epochs$n_epochs)) {
    for (ch in seq_len(n_ch)) {
      bl_median <- median(epochs$data[ch, bl_idx, ep])
      expect_lt(abs(bl_median), 1e-9,
                label = sprintf("baseline median for ch%d epoch%d", ch, ep))
    }
  }
})

# ----------------------------------------------------------------------------
# Test 2.20: No baseline correction preserves original signal values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When baseline = NULL (no baseline correction), the extracted
# epoch data must exactly match the corresponding slice of the original EEG
# data matrix — no values should be modified.
test_that("epoch_eeg with no baseline leaves data values unchanged", {
  sr   <- 256
  tmin <- -0.1
  tmax <-  0.4
  
  smin <- round(tmin * sr)   # -26
  smax <- round(tmax * sr)   # 102
  
  n_tp <- 2000
  set.seed(42)
  raw_data <- matrix(rnorm(2 * n_tp), nrow = 2)
  
  eeg <- make_mock_eeg(
    n_channels   = 2,
    n_timepoints = n_tp,
    sr           = sr,
    event_onsets = c(400L),
    event_types  = c("1"),
    data_values  = raw_data
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL, verbose = FALSE)
  
  onset       <- 400L
  ep_start    <- onset + smin    # 374
  ep_end      <- onset + smax    # 502
  
  expected_ch1 <- raw_data[1, ep_start:ep_end]
  expect_equal(as.vector(epochs$data[1, , 1]), expected_ch1)
})

# ----------------------------------------------------------------------------
# Test 2.21: events metadata in output matches the selected events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The events field of the returned "eeg_epochs" object must
# contain exactly one row per accepted epoch, with epoch_id re-assigned
# consecutively from 1. The onset_time values must correspond to the original
# events that were used.
test_that("epoch_eeg events metadata matches accepted epochs", {
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1", "2", "1")
  )
  
  # Epoch around type 1 only (two events)
  epochs <- epoch_eeg(eeg, events = 1,
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(nrow(epochs$events), 2)
  expect_equal(epochs$events$epoch_id, c(1L, 2L))
})

# ----------------------------------------------------------------------------
# Test 2.22: rejected vector length matches original number of candidate epochs
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The "rejected" logical vector in the returned object records
# which of the ORIGINAL n_trials epochs were rejected, before subsetting to
# valid ones. Its length must equal the total number of events that matched the
# selection criteria, including rejected ones.
test_that("epoch_eeg rejected vector has correct length", {
  sr   <- 256
  tmin <- -0.2
  tmax <-  0.4
  
  # 3 events: onset=5 will be rejected (boundary), the other two accepted
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = sr,
    event_onsets = c(5L, 500L, 900L),
    event_types  = c("1", "1", "1")
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = tmin, tmax = tmax,
                      baseline = NULL, verbose = FALSE)
  
  # rejected must have one entry per original trial (3 total)
  expect_equal(length(epochs$rejected), 3)
  # Exactly one rejection (the boundary violation)
  expect_equal(sum(epochs$rejected), 1)
})

# ----------------------------------------------------------------------------
# Test 2.23: preload = FALSE returns NULL data
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When preload = FALSE, epoch data is NOT loaded into memory.
# The data field of the returned "eeg_epochs" object must be NULL, while all
# other metadata (times, events, n_epochs, etc.) must still be populated.
test_that("epoch_eeg with preload=FALSE returns NULL data", {
  eeg <- make_mock_eeg(n_timepoints = 2000)
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL,
                      preload  = FALSE,
                      verbose  = FALSE)
  
  expect_null(epochs$data)
  expect_equal(epochs$n_epochs, 3)   # metadata still populated
})

# ----------------------------------------------------------------------------
# Test 2.24: Multiple event types are handled independently
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When events contain multiple trigger codes and we request
# all of them, the returned epochs$events must correctly record the type for
# each accepted epoch so that downstream condition-based averaging works.
test_that("epoch_eeg preserves event type labels for multiple conditions", {
  eeg <- make_mock_eeg(
    n_timepoints = 2000,
    sr           = 256,
    event_onsets = c(200L, 500L, 800L),
    event_types  = c("1", "2", "1")
  )
  
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  # The type column in epochs$events must reflect original codes after bitmask
  expect_equal(sort(unique(epochs$events$type)), c(1L, 2L))
})

# ----------------------------------------------------------------------------
# Test 2.25: tmin and tmax are correctly stored in output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The returned object's tmin and tmax fields must exactly
# match the values passed by the user, preserving them for inspection and
# reproducibility.
test_that("epoch_eeg stores tmin and tmax in output object", {
  eeg    <- make_mock_eeg(n_timepoints = 2000)
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.15, tmax = 0.55,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(epochs$tmin, -0.15)
  expect_equal(epochs$tmax,  0.55)
})

# ----------------------------------------------------------------------------
# Test 2.26: sampling_rate is correctly carried through to output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The output object must record the same sampling rate as the
# input eeg object so downstream calculations (frequency-domain, plotting) use
# the correct rate.
test_that("epoch_eeg carries sampling_rate through to output", {
  eeg    <- make_mock_eeg(sr = 512, n_timepoints = 4000,
                          event_onsets = c(400L, 1200L), event_types = c("1", "1"))
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_equal(epochs$sampling_rate, 512)
})

# ----------------------------------------------------------------------------
# Test 2.27: rejection_log is a data frame with the expected columns
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The rejection_log must always be a data frame (even when
# empty) with columns epoch_id, event_type, event_time, and reason. Downstream
# quality-control code relies on this structure.
test_that("epoch_eeg rejection_log is a data frame with correct columns", {
  eeg    <- make_mock_eeg(n_timepoints = 2000)
  epochs <- epoch_eeg(eeg, events = "all",
                      tmin = -0.1, tmax = 0.4,
                      baseline = NULL, verbose = FALSE)
  
  expect_true(is.data.frame(epochs$rejection_log))
  expected_cols <- c("epoch_id", "event_type", "event_time", "reason")
  expect_true(all(expected_cols %in% names(epochs$rejection_log)))
})


# ============================================================================
# TEST SUITE 3: plot_epochs() — Epoch Visualisation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Input validation — non-eeg_epochs is rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: plot_epochs() requires an object of class "eeg_epochs".
# Any other class must trigger an informative error, preventing cryptic
# downstream failures.
test_that("plot_epochs rejects non-eeg_epochs input", {
  expect_error(
    plot_epochs(list(data = array(1:10, c(2, 5, 1)))),
    class = "error"
  )
  expect_error(
    plot_epochs("not_epochs"),
    class = "error"
  )
})

# ----------------------------------------------------------------------------
# Test 3.2: plot_type = "butterfly" produces a plot without error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The butterfly plot (default) overlays all channel traces for
# all (or up to max_trials) epochs. This test confirms the complete plotting
# pipeline executes without throwing an error for well-formed input.
test_that("plot_epochs butterfly plot runs without error", {
  epochs <- make_mock_epochs(n_channels = 3, n_epochs = 5)
  
  expect_no_error(
    plot_epochs(epochs, plot_type = "butterfly", add_legend = FALSE)
  )
})

# ----------------------------------------------------------------------------
# Test 3.3: plot_type = "image" produces a plot without error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The image plot averages across channels and displays a
# colour-coded time × trial matrix. This test confirms the averaging and
# image() call complete successfully.
test_that("plot_epochs image plot runs without error", {
  epochs <- make_mock_epochs(n_channels = 2, n_epochs = 4)
  
  expect_no_error(
    plot_epochs(epochs, plot_type = "image")
  )
})

# ----------------------------------------------------------------------------
# Test 3.4: plot_type = "evoked" produces a plot without error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The evoked plot averages across trials and draws the mean
# waveform with a 95% confidence interval shading. This test exercises the
# apply(…, mean) and polygon() calls.
test_that("plot_epochs evoked plot runs without error", {
  epochs <- make_mock_epochs(n_channels = 2, n_epochs = 10)
  
  expect_no_error(
    plot_epochs(epochs, plot_type = "evoked", add_legend = FALSE)
  )
})

# ----------------------------------------------------------------------------
# Test 3.5: plot_type = "trials" produces a plot without error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The trial-by-trial plot creates a grid of mini-panels,
# one per trial. This test confirms the grid layout (par(mfrow = ...)) and
# per-trial line drawing work correctly.
test_that("plot_epochs trials plot runs without error", {
  epochs <- make_mock_epochs(n_channels = 2, n_epochs = 6)
  
  expect_no_error(
    plot_epochs(epochs, plot_type = "trials", n_trials = 4)
  )
})

# ----------------------------------------------------------------------------
# Test 3.6: Subsetting channels by name
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When channels is a character vector, plot_epochs() must
# restrict the plotted data to those channels. This test verifies that passing
# a valid channel name does not raise an error, confirming the channel-lookup
# logic works correctly.
test_that("plot_epochs accepts a subset of channels by name", {
  epochs <- make_mock_epochs(n_channels = 4, n_epochs = 5)
  
  expect_no_error(
    plot_epochs(epochs, plot_type = "butterfly",
                channels = c("Ch1", "Ch3"), add_legend = FALSE)
  )
})

# ----------------------------------------------------------------------------
# Test 3.7: Subsetting events by type
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When events is specified, plot_epochs() must filter to only
# the epochs that belong to those event types. This test checks the filter path
# does not crash for well-formed input.
test_that("plot_epochs accepts event-type filter without error", {
  epochs <- make_mock_epochs(n_channels = 2, n_epochs = 5)
  # All mock epochs have type = 1 so this should run without error
  expect_no_error(
    plot_epochs(epochs, plot_type = "butterfly",
                events = 1, add_legend = FALSE)
  )
})