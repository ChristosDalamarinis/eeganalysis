# ============================================================================
#                       Test File for eeg_class.R
# ============================================================================
# 
# This test file provides comprehensive testing for the eeg_class.R script,
# which contains the core EEG S3 class definition and methods.
#
# Functions tested:
#   1. new_eeg()    - Constructor for creating EEG objects
#   2. print.eeg()  - Print method for displaying EEG objects
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# TEST SUITE 1: new_eeg() - EEG Object Constructor
# ============================================================================

#context("new_eeg() - EEG Object Constructor")

# ----------------------------------------------------------------------------
# Test 1.1: Basic object creation with minimal required parameters
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() can create a valid EEG object with
# only the three required parameters (data, channels, sampling_rate)
test_that("new_eeg creates valid object with minimal parameters", {
  # Create simple test data
  test_data <- matrix(rnorm(300), nrow = 3, ncol = 100)
  test_channels <- c("Cz", "Pz", "Oz")
  test_sr <- 256
  
  # Create EEG object
  eeg <- new_eeg(
    data = test_data,
    channels = test_channels,
    sampling_rate = test_sr
  )
  
  # Test class assignment
  expect_s3_class(eeg, "eeg")
  
  # Test data structure
  expect_true(is.list(eeg))
  expect_equal(length(eeg), 8)  # Should have 8 components
  
  # Test data dimensions and values
  expect_equal(dim(eeg$data), c(3, 100))
  expect_true(is.matrix(eeg$data))
  
  # Test channels
  expect_equal(eeg$channels, test_channels)
  expect_type(eeg$channels, "character")
  
  # Test sampling rate
  expect_equal(eeg$sampling_rate, test_sr)
  expect_type(eeg$sampling_rate, "double")
  
  # Test auto-generated times vector
  expect_equal(length(eeg$times), 100)
  expect_equal(eeg$times[1], 0)
  expect_equal(eeg$times[100], 99/256)
  
  # Test default values for optional parameters
  expect_true(is.data.frame(eeg$events))
  expect_equal(nrow(eeg$events), 0)
  expect_equal(ncol(eeg$events), 4)
  expect_equal(names(eeg$events), c("onset", "onset_time", "type", "description"))
  
  expect_true(is.list(eeg$metadata))
  expect_equal(length(eeg$metadata), 0)
  
  expect_equal(eeg$reference, "original")
  
  expect_true(is.list(eeg$preprocessing_history))
  expect_equal(length(eeg$preprocessing_history), 0)
})


# ----------------------------------------------------------------------------
# Test 1.2: Object creation with all parameters specified
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() correctly handles all optional
# parameters and stores them in the appropriate format
test_that("new_eeg creates valid object with all parameters", {
  # Create test data
  test_data <- matrix(rnorm(200), nrow = 2, ncol = 100)
  test_channels <- c("Fz", "Cz")
  test_sr <- 512
  test_times <- seq(0, 99/512, length.out = 100)
  
  # Create events dataframe
  test_events <- data.frame(
    onset = c(10, 50, 80),
    onset_time = c(10/512, 50/512, 80/512),
    type = c(1, 2, 1),
    description = c("stim_A", "stim_B", "stim_A")
  )
  
  # Create metadata list
  test_metadata <- list(
    subject_id = "S001",
    session = 1,
    date = "2026-02-16",
    experimenter = "JD"
  )
  
  # Create preprocessing history
  test_history <- list(
    "Imported from BDF file",
    "Applied 0.1-40 Hz bandpass filter"
  )
  
  # Create EEG object with all parameters
  eeg <- new_eeg(
    data = test_data,
    channels = test_channels,
    sampling_rate = test_sr,
    times = test_times,
    events = test_events,
    metadata = test_metadata,
    reference = "average",
    preprocessing_history = test_history
  )
  
  # Validate all components
  expect_s3_class(eeg, "eeg")
  expect_equal(eeg$sampling_rate, test_sr)
  expect_equal(eeg$channels, test_channels)
  expect_equal(dim(eeg$data), c(2, 100))
  
  # Test times
  expect_equal(length(eeg$times), 100)
  expect_equal(eeg$times, test_times)
  
  # Test events
  expect_equal(nrow(eeg$events), 3)
  expect_equal(eeg$events$type, c(1, 2, 1))
  
  # Test metadata
  expect_equal(length(eeg$metadata), 4)
  expect_equal(eeg$metadata$subject_id, "S001")
  expect_equal(eeg$metadata$session, 1)
  
  # Test reference
  expect_equal(eeg$reference, "average")
  
  # Test preprocessing history
  expect_equal(length(eeg$preprocessing_history), 2)
  expect_equal(eeg$preprocessing_history[[1]], "Imported from BDF file")
})


# ----------------------------------------------------------------------------
# Test 1.3: Data type conversion (data.frame to matrix)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() correctly converts data.frame input
# to matrix format, as the function should accept both but store as matrix
test_that("new_eeg converts data.frame to matrix", {
  # Create data.frame input
  test_df <- data.frame(
    ch1 = rnorm(50),
    ch2 = rnorm(50),
    ch3 = rnorm(50)
  )
  
  # Transpose to get channels as rows
  test_data <- t(test_df)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("ch1", "ch2", "ch3"),
    sampling_rate = 256
  )
  
  # Verify conversion to matrix
  expect_true(is.matrix(eeg$data))
  expect_equal(dim(eeg$data), c(3, 50))
})


# ----------------------------------------------------------------------------
# Test 1.4: Error handling - Channel count mismatch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() throws an appropriate error when
# the number of channels doesn't match the number of rows in the data matrix
test_that("new_eeg errors on channel count mismatch", {
  test_data <- matrix(rnorm(300), nrow = 3, ncol = 100)
  
  # Test with too few channels
  expect_error(
    new_eeg(
      data = test_data,
      channels = c("Cz", "Pz"),  # Only 2 channels, but data has 3 rows
      sampling_rate = 256
    ),
    regexp = "ERROR: Number of channels.*does not match.*number of columns in data"
  )
  
  # Test with too many channels
  expect_error(
    new_eeg(
      data = test_data,
      channels = c("Cz", "Pz", "Oz", "Fz"),  # 4 channels, but data has 3 rows
      sampling_rate = 256
    ),
    regexp = "ERROR: Number of channels.*does not match.*number of columns in data"
  )
})


# ----------------------------------------------------------------------------
# Test 1.5: Error handling - Time vector length mismatch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() throws an error when a custom
# times vector is provided but doesn't match the number of timepoints
test_that("new_eeg errors on time vector length mismatch", {
  test_data <- matrix(rnorm(300), nrow = 3, ncol = 100)
  test_channels <- c("Cz", "Pz", "Oz")
  
  # Create mismatched times vector
  wrong_times <- seq(0, 49/256, length.out = 50)  # Only 50 timepoints
  
  expect_error(
    new_eeg(
      data = test_data,
      channels = test_channels,
      sampling_rate = 256,
      times = wrong_times
    ),
    regexp = "ERROR: Length of times.*does not match.*number of rows in data"
  )
})


# ----------------------------------------------------------------------------
# Test 1.6: Auto-generated time vector accuracy
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when no times vector is provided, new_eeg()
# correctly generates one based on the sampling rate and number of timepoints
test_that("new_eeg generates correct time vector", {
  test_data <- matrix(rnorm(500), nrow = 5, ncol = 100)
  test_sr <- 1024
  
  eeg <- new_eeg(
    data = test_data,
    channels = paste0("Ch", 1:5),
    sampling_rate = test_sr
  )
  
  # Test time vector properties
  expect_equal(length(eeg$times), 100)
  expect_equal(eeg$times[1], 0)
  expect_equal(eeg$times[2], 1/1024)
  expect_equal(eeg$times[100], 99/1024)
  
  # Test that times are equally spaced
  time_diffs <- diff(eeg$times)
  expect_true(all(abs(time_diffs - 1/1024) < 1e-10))
})


# ----------------------------------------------------------------------------
# Test 1.7: Type coercion for all fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() properly coerces all input fields
# to their expected types (numeric, character, etc.)
test_that("new_eeg properly coerces field types", {
  test_data <- matrix(1:30, nrow = 3, ncol = 10)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("A", "B", "C"),
    sampling_rate = 100,
    reference = "linked_ears"
  )
  
  # Test type coercions
  expect_true(is.matrix(eeg$data))
  expect_type(eeg$data, "double")  # Should be coerced to numeric
  
  expect_type(eeg$channels, "character")
  expect_type(eeg$sampling_rate, "double")
  expect_type(eeg$times, "double")
  expect_type(eeg$reference, "character")
})


# ----------------------------------------------------------------------------
# Test 1.8: Edge case - Single channel
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() handles the edge case of a single
# channel correctly, maintaining proper matrix structure
test_that("new_eeg handles single channel correctly", {
  test_data <- matrix(rnorm(100), nrow = 1, ncol = 100)
  
  eeg <- new_eeg(
    data = test_data,
    channels = "Cz",
    sampling_rate = 256
  )
  
  expect_s3_class(eeg, "eeg")
  expect_equal(length(eeg$channels), 1)
  expect_equal(dim(eeg$data), c(1, 100))
  expect_equal(eeg$channels, "Cz")
})


# ----------------------------------------------------------------------------
# Test 1.9: Edge case - Single timepoint
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() handles the edge case of a single
# timepoint (though unusual for EEG data)
test_that("new_eeg handles single timepoint", {
  test_data <- matrix(c(1.5, 2.3, -0.8), nrow = 3, ncol = 1)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )
  
  expect_s3_class(eeg, "eeg")
  expect_equal(ncol(eeg$data), 1)
  expect_equal(length(eeg$times), 1)
  expect_equal(eeg$times[1], 0)
})


# ----------------------------------------------------------------------------
# Test 1.10: Empty events dataframe structure
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when no events are provided, new_eeg() creates
# an empty dataframe with the correct column structure
test_that("new_eeg creates correct empty events dataframe", {
  test_data <- matrix(rnorm(100), nrow = 2, ncol = 50)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 512
  )
  
  # Test events structure
  expect_true(is.data.frame(eeg$events))
  expect_equal(nrow(eeg$events), 0)
  expect_equal(ncol(eeg$events), 4)
  expect_equal(names(eeg$events), c("onset", "onset_time", "type", "description"))
  
  # Test column types
  expect_type(eeg$events$onset, "integer")
  expect_type(eeg$events$onset_time, "double")
  expect_type(eeg$events$type, "character")
  expect_type(eeg$events$description, "character")
})


# ----------------------------------------------------------------------------
# Test 1.11: Large data handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() can handle realistic EEG data sizes
# without errors or performance issues
test_that("new_eeg handles large realistic datasets", {
  # Simulate 64 channels, 10 minutes at 2048 Hz
  n_channels <- 64
  duration_sec <- 600  # 10 minutes
  sr <- 2048
  n_samples <- duration_sec * sr
  
  # Create large dataset
  large_data <- matrix(rnorm(n_channels * n_samples), 
                       nrow = n_channels, 
                       ncol = n_samples)
  
  # Create channel names
  large_channels <- paste0("Ch", 1:n_channels)
  
  # This should complete without error
  eeg <- new_eeg(
    data = large_data,
    channels = large_channels,
    sampling_rate = sr
  )
  
  expect_s3_class(eeg, "eeg")
  expect_equal(dim(eeg$data), c(64, n_samples))
  expect_equal(length(eeg$times), n_samples)
})


# ----------------------------------------------------------------------------
# Test 1.12: Numeric channel names handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that numeric channel names are properly converted
# to character strings
test_that("new_eeg converts numeric channel names to character", {
  test_data <- matrix(rnorm(150), nrow = 3, ncol = 50)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c(1, 2, 3),  # Numeric channel names
    sampling_rate = 256
  )
  
  expect_type(eeg$channels, "character")
  expect_equal(eeg$channels, c("1", "2", "3"))
})


# ----------------------------------------------------------------------------
# Test 1.13: Custom events with various types
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() correctly stores events with
# different data types and structures
test_that("new_eeg handles various event types", {
  test_data <- matrix(rnorm(200), nrow = 2, ncol = 100)
  
  # Events with mixed types
  test_events <- data.frame(
    onset = c(10L, 50L, 80L),
    onset_time = c(0.039, 0.195, 0.312),
    type = c("A", "B", "A"),
    description = c("Stimulus onset", "Response", "Stimulus onset")
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256,
    events = test_events
  )
  
  expect_equal(nrow(eeg$events), 3)
  expect_equal(eeg$events$type, c("A", "B", "A"))
  expect_equal(eeg$events$description[1], "Stimulus onset")
})


# ----------------------------------------------------------------------------
# Test 1.14: Metadata with complex structures
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_eeg() can store complex metadata including
# nested lists and various data types
test_that("new_eeg handles complex metadata structures", {
  test_data <- matrix(rnorm(100), nrow = 2, ncol = 50)
  
  complex_metadata <- list(
    subject_id = "S042",
    age = 25,
    gender = "F",
    session_info = list(
      date = "2026-02-16",
      time = "14:30",
      location = "Lab A"
    ),
    equipment = c("BioSemi ActiveTwo", "64-channel cap"),
    notes = "First session, participant alert"
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Fz", "Cz"),
    sampling_rate = 512,
    metadata = complex_metadata
  )
  
  expect_equal(length(eeg$metadata), 6)
  expect_equal(eeg$metadata$subject_id, "S042")
  expect_equal(eeg$metadata$age, 25)
  expect_true(is.list(eeg$metadata$session_info))
  expect_equal(eeg$metadata$session_info$date, "2026-02-16")
})


# ----------------------------------------------------------------------------
# Test 1.15: Preprocessing history tracking
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that preprocessing history is correctly stored and
# can track multiple processing steps
test_that("new_eeg maintains preprocessing history", {
  test_data <- matrix(rnorm(150), nrow = 3, ncol = 50)
  
  history <- list(
    "Step 1: Imported from data.bdf",
    "Step 2: High-pass filtered at 0.1 Hz",
    "Step 3: Notch filtered at 50 Hz",
    "Step 4: Re-referenced to average"
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Fz", "Cz", "Pz"),
    sampling_rate = 1024,
    preprocessing_history = history
  )
  
  expect_equal(length(eeg$preprocessing_history), 4)
  expect_equal(eeg$preprocessing_history[[3]], "Step 3: Notch filtered at 50 Hz")
})


# ============================================================================
#           TEST SUITE 2: print.eeg() - Print Method for EEG Objects
# ============================================================================
#
#
# ----------------------------------------------------------------------------
# Test 2.1: Basic print output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() produces output without errors
# and returns the object invisibly
test_that("print.eeg produces output and returns invisibly", {
  test_data <- matrix(rnorm(300), nrow = 3, ncol = 100)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )
  
  # Capture output
  output <- capture.output(result <- print(eeg))
  
  # Should return the object invisibly
  expect_identical(result, eeg)
  
  # Should produce some output
  expect_true(length(output) > 0)
  
  # Check for key sections in output
  expect_true(any(grepl("EEG Object Summary", output)))
  expect_true(any(grepl("RECORDING INFORMATION", output)))
  expect_true(any(grepl("Channels:", output)))
  expect_true(any(grepl("Sampling rate:", output)))
})


# ----------------------------------------------------------------------------
# Test 2.2: Channel list display formatting
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() correctly formats the channel list,
# showing first 5 channels and truncating with "... (+N more)" if there are more
test_that("print.eeg formats channel list correctly", {
  # Test with few channels (< 5)
  test_data_small <- matrix(rnorm(150), nrow = 3, ncol = 50)
  eeg_small <- new_eeg(
    data = test_data_small,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )
  
  output_small <- capture.output(print(eeg_small))
  channel_line <- output_small[grepl("List:", output_small)]
  
  # Should show all 3 channels
  expect_true(any(grepl("Cz, Pz, Oz", output_small)))
  expect_false(any(grepl("\\+.*more", output_small)))
  
  # Test with many channels (> 5)
  test_data_large <- matrix(rnorm(1000), nrow = 10, ncol = 100)
  eeg_large <- new_eeg(
    data = test_data_large,
    channels = paste0("Ch", 1:10),
    sampling_rate = 256
  )
  
  output_large <- capture.output(print(eeg_large))
  
  # Should show first 5 channels and truncate
  expect_true(any(grepl("Ch1, Ch2, Ch3, Ch4, Ch5", output_large)))
  expect_true(any(grepl("\\+5 more", output_large)))
})


# ----------------------------------------------------------------------------
# Test 2.3: Data statistics display
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() correctly displays amplitude
# statistics including range, mean, and standard deviation
test_that("print.eeg displays data statistics correctly", {
  # Create data with known statistics
  test_data <- matrix(seq(-10, 10, length.out = 300), nrow = 3, ncol = 100)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("A", "B", "C"),
    sampling_rate = 256
  )
  
  output <- capture.output(print(eeg))
  
  # Check for statistics section
  expect_true(any(grepl("DATA STATISTICS", output)))
  expect_true(any(grepl("Amplitude range:", output)))
  expect_true(any(grepl("Mean amplitude:", output)))
  expect_true(any(grepl("Std deviation:", output)))
  
  # Verify approximate values appear
  expect_true(any(grepl("-10", output)))
  expect_true(any(grepl("10", output)))
})


# ----------------------------------------------------------------------------
# Test 2.4: Events display
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() correctly displays event
# information including total count and event types
test_that("print.eeg displays events correctly", {
  test_data <- matrix(rnorm(200), nrow = 2, ncol = 100)
  
  # Test with events
  events_df <- data.frame(
    onset = c(10, 50, 80),
    onset_time = c(0.039, 0.195, 0.312),
    type = c(1, 2, 1),
    description = c("A", "B", "A")
  )
  
  eeg_with_events <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256,
    events = events_df
  )
  
  output_with <- capture.output(print(eeg_with_events))
  
  expect_true(any(grepl("EVENTS", output_with)))
  expect_true(any(grepl("Total events:.*3", output_with)))
  expect_true(any(grepl("Event types:.*1.*2", output_with)))
  
  # Test without events
  eeg_no_events <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256
  )
  
  output_no <- capture.output(print(eeg_no_events))
  expect_true(any(grepl("Total events:.*0", output_no)))
})


# ----------------------------------------------------------------------------
# Test 2.5: Metadata display
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() displays metadata when present
# and formats it appropriately
test_that("print.eeg displays metadata correctly", {
  test_data <- matrix(rnorm(150), nrow = 3, ncol = 50)
  
  metadata <- list(
    subject = "S001",
    age = 28,
    session = 1
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("A", "B", "C"),
    sampling_rate = 512,
    metadata = metadata
  )
  
  output <- capture.output(print(eeg))
  
  expect_true(any(grepl("METADATA", output)))
  expect_true(any(grepl("subject.*S001", output)))
  expect_true(any(grepl("age.*28", output)))
  expect_true(any(grepl("session.*1", output)))
})


# ----------------------------------------------------------------------------
# Test 2.6: Preprocessing history display
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() displays preprocessing history
# when present, with proper numbering
test_that("print.eeg displays preprocessing history correctly", {
  test_data <- matrix(rnorm(150), nrow = 3, ncol = 50)
  
  history <- list(
    "Imported from BDF",
    "Filtered 0.1-40 Hz",
    "Re-referenced to average"
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("A", "B", "C"),
    sampling_rate = 256,
    preprocessing_history = history
  )
  
  output <- capture.output(print(eeg))
  
  expect_true(any(grepl("PREPROCESSING HISTORY", output)))
  expect_true(any(grepl("1.*Imported from BDF", output)))
  expect_true(any(grepl("2.*Filtered 0.1-40 Hz", output)))
  expect_true(any(grepl("3.*Re-referenced to average", output)))
})


# ----------------------------------------------------------------------------
# Test 2.7: Reference scheme display
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() correctly displays the reference
# scheme used for the EEG data
test_that("print.eeg displays reference scheme", {
  test_data <- matrix(rnorm(100), nrow = 2, ncol = 50)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256,
    reference = "average"
  )
  
  output <- capture.output(print(eeg))
  
  expect_true(any(grepl("REFERENCE", output)))
  expect_true(any(grepl("Scheme:.*average", output)))
})


# ----------------------------------------------------------------------------
# Test 2.8: Duration calculation accuracy
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() correctly calculates and displays
# the recording duration
test_that("print.eeg calculates duration correctly", {
  # Create data with known duration
  sr <- 512
  duration <- 10  # 10 seconds
  n_samples <- sr * duration
  
  test_data <- matrix(rnorm(2 * n_samples), nrow = 2, ncol = n_samples)
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = sr
  )
  
  output <- capture.output(print(eeg))
  
  # Should show approximately 10 seconds (9.998 due to 0-indexing)
  expect_true(any(grepl("Duration:.*9\\.99|10\\.00", output)))
})


# ----------------------------------------------------------------------------
# Test 2.9: Long metadata value truncation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() truncates very long metadata
# values to prevent console overflow
test_that("print.eeg truncates long metadata values", {
  test_data <- matrix(rnorm(100), nrow = 2, ncol = 50)
  
  long_string <- paste(rep("A", 100), collapse = "")
  
  metadata <- list(
    normal = "short",
    very_long = long_string
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256,
    metadata = metadata
  )
  
  output <- capture.output(print(eeg))
  
  # Long values should be truncated with "..."
  long_line <- output[grepl("very_long", output)]
  expect_true(any(grepl("\\.\\.\\.", output)))
})


# ----------------------------------------------------------------------------
# Test 2.10: Empty object handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() handles an EEG object with minimal
# data (single channel, single timepoint) without errors
test_that("print.eeg handles minimal object", {
  test_data <- matrix(0, nrow = 1, ncol = 1)
  
  eeg <- new_eeg(
    data = test_data,
    channels = "Single",
    sampling_rate = 256
  )
  
  # Should not error
  expect_silent(output <- capture.output(print(eeg)))
  expect_true(length(output) > 0)
})


# ----------------------------------------------------------------------------
# Test 2.11: Multiple channel types in metadata
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that print.eeg() handles non-atomic metadata values
# appropriately (like nested lists or complex objects)
test_that("print.eeg handles complex metadata types", {
  test_data <- matrix(rnorm(100), nrow = 2, ncol = 50)
  
  metadata <- list(
    subject = "S001",
    nested_list = list(a = 1, b = 2),
    vector = c(1, 2, 3, 4, 5)
  )
  
  eeg <- new_eeg(
    data = test_data,
    channels = c("Ch1", "Ch2"),
    sampling_rate = 256,
    metadata = metadata
  )
  
  # Should not error when printing complex metadata
  expect_silent(output <- capture.output(print(eeg)))
  expect_true(any(grepl("nested_list.*list", output)))
})


# ============================================================================
#                     TEST SUITE 3: Integration Tests
# ============================================================================
#
#
# ----------------------------------------------------------------------------
# Test 3.1: Round-trip test
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that creating an object and printing it maintains
# data integrity
test_that("Creating and printing maintains data integrity", {
  original_data <- matrix(rnorm(200), nrow = 4, ncol = 50)
  original_channels <- c("F3", "F4", "C3", "C4")
  
  eeg <- new_eeg(
    data = original_data,
    channels = original_channels,
    sampling_rate = 512
  )
  
  # Print shouldn't modify the object
  output <- capture.output(result <- print(eeg))
  
  # Data should be unchanged
  expect_equal(result$data, original_data)
  expect_equal(result$channels, original_channels)
  expect_equal(result$sampling_rate, 512)
})


# ----------------------------------------------------------------------------
# Test 3.2: Realistic workflow simulation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Simulates a realistic workflow of creating an EEG object
# with various components and verifying it prints correctly
test_that("Realistic workflow - create, modify, print", {
  # Initial creation
  eeg <- new_eeg(
    data = matrix(rnorm(6400), nrow = 64, ncol = 100),
    channels = paste0("E", 1:64),
    sampling_rate = 2048,
    metadata = list(
      subject = "P001",
      date = "2026-02-16"
    )
  )
  
  # Verify initial state
  expect_s3_class(eeg, "eeg")
  expect_equal(length(eeg$channels), 64)
  
  # Print should work
  output <- capture.output(print(eeg))
  expect_true(any(grepl("64", output)))  # Channel count
  expect_true(any(grepl("2048", output)))  # Sampling rate
  expect_true(any(grepl("P001", output)))  # Metadata
})


# ============================================================================
#                             End of Test File
# ============================================================================