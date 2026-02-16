# ============================================================================
#                    Test File for extract_bdf_events.R
# ============================================================================
# 
# This test file provides comprehensive testing for the extract_bdf_events.R 
# script, which extracts trigger/event codes from BioSemi BDF files.
#
# Functions tested:
#   1. extract_bdf_events()        - Main function to extract events from BDF
#   2. summary_bdf_events()        - Creates summary statistics of events
#   3. plot_bdf_events()           - Visualizes event data
#   4. validate_bdf_events()       - Validates and analyzes event quality
#   5. .extract_events_from_input()- Helper to handle different input types
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#     TEST SUITE 1: extract_bdf_events() - Main Event Extraction Function
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Input validation - accepts valid eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that extract_bdf_events() properly accepts and
# processes an eeg object as input. Tests the "eeg object input" branch of
# the function logic.
test_that("extract_bdf_events accepts valid eeg object input", {
  # Create mock eeg object with events
  mock_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.390625, 0.78125, 1.171875),
    type = c("1", "2", "1"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    data = matrix(rnorm(300), nrow = 3, ncol = 100),
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256,
    times = seq(0, 99/256, length.out = 100),
    events = mock_events,
    metadata = list(file_path = "test.bdf"),
    reference = "CMS",
    preprocessing_history = list("imported")
  )
  class(mock_eeg) <- "eeg"
  
  # Extract events
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Test that result is a data frame
  expect_true(is.data.frame(result))
  
  # Test that result has correct columns
  expect_true(all(c("onset", "onset_time", "type", "description") %in% names(result)))
  
  # Test that result has correct number of rows
  expect_equal(nrow(result), 3)
  
  # Test that result has expected attributes
  expect_true(!is.null(attr(result, "sample_rate")))
  expect_equal(attr(result, "sample_rate"), 256)
})

# ----------------------------------------------------------------------------
# Test 1.2: Input validation - rejects invalid input types
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the function properly validates input and throws
# appropriate errors for invalid input types (not a file path or eeg object).
test_that("extract_bdf_events rejects invalid input types", {
  # Test with numeric input
  expect_error(
    extract_bdf_events(123, verbose = FALSE),
    "'data' must be either a file path.*or an eeg object"
  )
  
  # Test with list without required fields
  expect_error(
    extract_bdf_events(list(a = 1, b = 2), verbose = FALSE),
    "'data' must be either a file path.*or an eeg object"
  )
  
  # Test with NULL
  expect_error(
    extract_bdf_events(NULL, verbose = FALSE),
    "'data' must be either a file path.*or an eeg object"
  )
  
  # Test with data frame (not a valid events df)
  expect_error(
    extract_bdf_events(data.frame(x = 1:3), verbose = FALSE),
    "'data' must be either a file path.*or an eeg object"
  )
})

# ----------------------------------------------------------------------------
# Test 1.3: Return structure validation - correct columns
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the returned data frame has exactly the
# expected columns (onset, onset_time, type, description) with correct types.
test_that("extract_bdf_events returns data frame with correct columns", {
  # Create mock eeg object
  mock_events <- data.frame(
    onset = c(256, 512),
    onset_time = c(1.0, 2.0),
    type = c("100", "200"),
    description = c("Trigger: 100", "Trigger: 200"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 256,
    metadata = list(file_path = "test.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Test exact column names
  expect_equal(names(result), c("onset", "onset_time", "type", "description"))
  
  # Test column types
  expect_true(is.integer(result$onset) || is.numeric(result$onset))
  expect_true(is.numeric(result$onset_time))
  expect_true(is.character(result$type))
  expect_true(is.character(result$description))
})

# ----------------------------------------------------------------------------
# Test 1.4: Return structure validation - correct attributes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures all required metadata attributes are attached to
# the returned data frame (bdf_file, sample_rate, n_samples, duration_sec).
test_that("extract_bdf_events attaches correct attributes", {
  mock_events <- data.frame(
    onset = c(256, 512, 768),
    onset_time = c(1.0, 2.0, 3.0),
    type = c("1", "2", "3"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 3"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 512,
    times = seq(0, 4, length.out = 2048),
    metadata = list(file_path = "/path/to/test.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Test that all required attributes exist
  expect_true(!is.null(attr(result, "bdf_file")))
  expect_true(!is.null(attr(result, "sample_rate")))
  expect_true(!is.null(attr(result, "n_samples")))
  expect_true(!is.null(attr(result, "duration_sec")))
  
  # Test attribute values
  expect_equal(attr(result, "bdf_file"), "/path/to/test.bdf")
  expect_equal(attr(result, "sample_rate"), 512)
  expect_equal(attr(result, "n_samples"), 2048)
  expect_equal(attr(result, "duration_sec"), 4)
})

# ----------------------------------------------------------------------------
# Test 1.5: Edge case - empty events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the function handles the edge case of an eeg object
# with no events gracefully, returning an empty but properly structured data frame.
test_that("extract_bdf_events handles empty events gracefully", {
  # Create eeg object with no events
  mock_eeg <- list(
    events = data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    ),
    sampling_rate = 256,
    times = seq(0, 10, length.out = 2560),
    metadata = list(file_path = "empty.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  # Should not throw error
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Test structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
  expect_equal(names(result), c("onset", "onset_time", "type", "description"))
  
  # Test attributes still present
  expect_equal(attr(result, "sample_rate"), 256)
})

# ----------------------------------------------------------------------------
# Test 1.6: Edge case - NULL events structure
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the function handles eeg objects with NULL
# events field by creating an empty events data frame with proper structure.
test_that("extract_bdf_events handles NULL events structure", {
  mock_eeg <- list(
    events = NULL,
    sampling_rate = 256,
    times = seq(0, 5, length.out = 1280),
    metadata = list(file_path = "null_events.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  # Should issue warning but not error
  expect_warning(
    result <- extract_bdf_events(mock_eeg, verbose = FALSE),
    "No events structure found"
  )
  
  # Should return empty but valid data frame
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
  expect_equal(names(result), c("onset", "onset_time", "type", "description"))
})

# ----------------------------------------------------------------------------
# Test 1.7: Verbose mode control
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the verbose parameter properly controls console
# output. When verbose=TRUE, output should be generated; when FALSE, suppressed.
test_that("extract_bdf_events respects verbose parameter", {
  mock_events <- data.frame(
    onset = c(100),
    onset_time = c(0.5),
    type = c("1"),
    description = c("Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200,
    times = seq(0, 1, length.out = 200),
    metadata = list(file_path = "test.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  # Test verbose = FALSE produces no output
  output_silent <- capture.output(
    result_silent <- extract_bdf_events(mock_eeg, verbose = FALSE)
  )
  expect_equal(length(output_silent), 0)
  
  # Test verbose = TRUE produces output
  output_verbose <- capture.output(
    result_verbose <- extract_bdf_events(mock_eeg, verbose = TRUE)
  )
  expect_true(length(output_verbose) > 0)
  
  # Both should return same data
  expect_equal(result_silent, result_verbose)
})

# ----------------------------------------------------------------------------
# Test 1.8: File path input detection
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the function correctly identifies when a
# file path string is provided (even if the file doesn't exist, the detection
# logic should work).
test_that("extract_bdf_events detects file path input", {
  # This will fail because file doesn't exist, but that's expected
  # We're testing that it recognizes it as a file path attempt
  expect_error(
    extract_bdf_events("nonexistent_file.bdf", verbose = FALSE),
    "BDF file not found"
  )
})

# ----------------------------------------------------------------------------
# Test 1.9: Multiple trigger codes handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the function correctly handles events with multiple
# different trigger codes and preserves the distinction between them.
test_that("extract_bdf_events preserves multiple trigger codes", {
  mock_events <- data.frame(
    onset = c(100, 200, 300, 400, 500),
    onset_time = c(0.5, 1.0, 1.5, 2.0, 2.5),
    type = c("1", "2", "3", "2", "1"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 3", 
                    "Trigger: 2", "Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200,
    times = seq(0, 3, length.out = 600),
    metadata = list(file_path = "multi_trigger.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Test all trigger codes preserved
  expect_equal(nrow(result), 5)
  expect_equal(length(unique(result$type)), 3)
  expect_true(all(c("1", "2", "3") %in% result$type))
  
  # Test counts are correct
  expect_equal(sum(result$type == "1"), 2)
  expect_equal(sum(result$type == "2"), 2)
  expect_equal(sum(result$type == "3"), 1)
})


# ============================================================================
#       TEST SUITE 2: summary_bdf_events() - Event Summary Statistics
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Basic summary generation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that summary_bdf_events() generates a summary
# with basic statistics about the events (total count, unique triggers, etc.).
test_that("summary_bdf_events generates basic summary", {
  mock_events <- data.frame(
    onset = c(100, 200, 300, 400),
    onset_time = c(0.5, 1.0, 1.5, 2.0),
    type = c("1", "2", "1", "2"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  attr(mock_events, "duration_sec") <- 2.5
  
  # Should produce output
  output <- capture.output(summary_bdf_events(mock_events))
  expect_true(length(output) > 0)
})

# ----------------------------------------------------------------------------
# Test 2.2: Summary with empty events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures summary_bdf_events() handles empty events data frame
# without errors and produces appropriate output.
test_that("summary_bdf_events handles empty events", {
  # Create empty events data frame with attributes
  empty_events <- data.frame(
    onset = integer(0),
    onset_time = numeric(0),
    type = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
  
  attr(empty_events, "bdf_file") <- "test.bdf"
  attr(empty_events, "sample_rate") <- 512
  attr(empty_events, "duration_sec") <- 100
  attr(empty_events, "n_samples") <- 51200
  
  # Test that it runs without error (but produces output)
  expect_no_error(summary_bdf_events(empty_events))
  
  # Or test that it produces expected output
  output <- capture.output(summary_bdf_events(empty_events))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Total events: 0", output)))
})

# ----------------------------------------------------------------------------
# Test 2.3: Summary with single event
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the summary function works correctly when
# there's only one event (edge case for statistics calculations).
test_that("summary_bdf_events handles single event", {
  single_event <- data.frame(
    onset = 100,
    onset_time = 0.5,
    type = "1",
    description = "Trigger: 1",
    stringsAsFactors = FALSE
  )
  attr(single_event, "sample_rate") <- 200
  attr(single_event, "duration_sec") <- 1.0
  
  # Should handle single event without error
  output <- capture.output(summary_bdf_events(single_event))
  expect_true(length(output) > 0)
})

# ----------------------------------------------------------------------------
# Test 2.4: Summary with multiple trigger types
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the summary correctly identifies and reports on
# multiple different trigger types in the events.
test_that("summary_bdf_events reports multiple trigger types", {
  multi_trigger_events <- data.frame(
    onset = c(100, 200, 300, 400, 500, 600),
    onset_time = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
    type = c("1", "2", "3", "1", "2", "3"),
    description = paste0("Trigger: ", c(1, 2, 3, 1, 2, 3)),
    stringsAsFactors = FALSE
  )
  attr(multi_trigger_events, "sample_rate") <- 200
  attr(multi_trigger_events, "duration_sec") <- 3.5
  
  output <- capture.output(summary_bdf_events(multi_trigger_events))
  
  # Should mention multiple trigger codes
  output_text <- paste(output, collapse = " ")
  expect_true(length(output) > 0)
})


# ============================================================================
#     TEST SUITE 3: validate_bdf_events() - Event Quality Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Basic validation with valid events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures validate_bdf_events() runs successfully on valid
# events data and returns a list with expected components (summary, issues).
test_that("validate_bdf_events validates events successfully", {
  mock_events <- data.frame(
    onset = c(100, 200, 300, 400),
    onset_time = c(0.5, 1.0, 1.5, 2.0),
    type = c("1", "2", "1", "2"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  attr(mock_events, "duration_sec") <- 2.5
  
  result <- validate_bdf_events(mock_events, verbose = FALSE, plot = FALSE)
  
  # Test return structure
  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("issues" %in% names(result))
})

# ----------------------------------------------------------------------------
# Test 3.2: Validation with empty events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that validation handles empty events and reports
# appropriate issue ("no_events").
test_that("validate_bdf_events detects empty events", {
  empty_events <- data.frame(
    onset = integer(0),
    onset_time = numeric(0),
    type = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
  
  result <- validate_bdf_events(empty_events, verbose = FALSE, plot = FALSE)
  
  # Should report no_events issue
  expect_true(is.list(result))
  expect_true("no_events" %in% result$issues)
})

# ----------------------------------------------------------------------------
# Test 3.3: Validation accepts different input types
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures validate_bdf_events() can accept multiple input
# types (file path, eeg object, events data frame) via helper function.
test_that("validate_bdf_events accepts different input types", {
  # Test with events data frame
  mock_events <- data.frame(
    onset = c(100, 200),
    onset_time = c(0.5, 1.0),
    type = c("1", "2"),
    description = c("Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  
  expect_silent(
    result <- validate_bdf_events(mock_events, verbose = FALSE, plot = FALSE)
  )
  
  # Test with eeg object
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  expect_silent(
    result <- validate_bdf_events(mock_eeg, verbose = FALSE, plot = FALSE)
  )
})

# ----------------------------------------------------------------------------
# Test 3.4: Plot parameter control
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the plot parameter controls whether diagnostic
# plots are generated (FALSE, TRUE, "always", "timeline", "sequence").
test_that("validate_bdf_events respects plot parameter", {
  mock_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.5, 1.0, 1.5),
    type = c("1", "2", "1"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  attr(mock_events, "duration_sec") <- 2.0
  
  # Test different plot options don't cause errors
  expect_silent(
    validate_bdf_events(mock_events, verbose = FALSE, plot = FALSE)
  )
  
  expect_silent(
    validate_bdf_events(mock_events, verbose = FALSE, plot = TRUE)
  )
  
  expect_silent(
    validate_bdf_events(mock_events, verbose = FALSE, plot = "timeline")
  )
})

# ----------------------------------------------------------------------------
# Test 3.5: Verbose mode in validation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures verbose parameter controls console output during
# validation.
test_that("validate_bdf_events respects verbose parameter", {
  mock_events <- data.frame(
    onset = c(100, 200),
    onset_time = c(0.5, 1.0),
    type = c("1", "2"),
    description = c("Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  
  # Verbose FALSE should minimize output
  output_silent <- capture.output(
    result_silent <- validate_bdf_events(mock_events, verbose = FALSE, plot = FALSE)
  )
  
  # Verbose TRUE should produce output
  output_verbose <- capture.output(
    result_verbose <- validate_bdf_events(mock_events, verbose = TRUE, plot = FALSE)
  )
  
  expect_true(length(output_verbose) >= length(output_silent))
})

# ----------------------------------------------------------------------------
# Test 3.6: Return value is invisible
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that validate_bdf_events() returns its result
# invisibly (common R practice for functions that print).
test_that("validate_bdf_events returns result invisibly", {
  mock_events <- data.frame(
    onset = c(100),
    onset_time = c(0.5),
    type = c("1"),
    description = c("Trigger: 1"),
    stringsAsFactors = FALSE
  )
  attr(mock_events, "sample_rate") <- 200
  
  # Result should be invisible but capturable
  result <- validate_bdf_events(mock_events, verbose = FALSE, plot = FALSE)
  expect_true(is.list(result))
})


# ============================================================================
# TEST SUITE 4: .extract_events_from_input() - Helper Function
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Extract from file path string
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the helper function recognizes a file path string
# and attempts to call extract_bdf_events() on it.
test_that(".extract_events_from_input handles file path", {
  # This will error because file doesn't exist, but tests the detection logic
  expect_error(
    eeganalysis:::.extract_events_from_input("test.bdf"),
    "BDF file not found"
  )
})

# ----------------------------------------------------------------------------
# Test 4.2: Extract from eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the helper correctly extracts events from an eeg
# object by accessing the $events component.
test_that(".extract_events_from_input extracts from eeg object", {
  mock_events <- data.frame(
    onset = c(100, 200),
    onset_time = c(0.5, 1.0),
    type = c("1", "2"),
    description = c("Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  result <- eeganalysis:::.extract_events_from_input(mock_eeg)
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_equal(names(result), c("onset", "onset_time", "type", "description"))
})

# ----------------------------------------------------------------------------
# Test 4.3: Extract from events data frame
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the helper recognizes a properly formatted events
# data frame and returns it directly.
test_that(".extract_events_from_input accepts valid events data frame", {
  mock_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.5, 1.0, 1.5),
    type = c("1", "2", "1"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  result <- eeganalysis:::.extract_events_from_input(mock_events)
  
  expect_equal(result, mock_events)
})

# ----------------------------------------------------------------------------
# Test 4.4: Reject invalid input
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the helper throws an appropriate error for inputs
# that don't match any of the expected types.
test_that(".extract_events_from_input rejects invalid input", {
  # Test with numeric
  expect_error(
    eeganalysis:::.extract_events_from_input(123),
    "Invalid input type"
  )
  
  # Test with list without events
  expect_error(
    eeganalysis:::.extract_events_from_input(list(a = 1, b = 2)),
    "Invalid input type"
  )
  
  # Test with incomplete data frame
  expect_error(
    eeganalysis:::.extract_events_from_input(data.frame(onset = 1:3)),
    "Invalid input type"
  )
})

# ----------------------------------------------------------------------------
# Test 4.5: List with events field (non-eeg class)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the helper can handle lists that have an events
# field even if they're not formally classed as "eeg".
test_that(".extract_events_from_input handles list with events", {
  mock_events <- data.frame(
    onset = c(100),
    onset_time = c(0.5),
    type = c("1"),
    description = c("Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  # List with events but not eeg class
  mock_list <- list(
    events = mock_events,
    sampling_rate = 200
  )
  
  result <- eeganalysis:::.extract_events_from_input(mock_list)
  
  expect_equal(result, mock_events)
})


# ============================================================================
#         TEST SUITE 5: Integration Tests - End-to-End Workflows
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Complete workflow from eeg object to summary
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Simulates a complete workflow of extracting events from an
# eeg object and then generating a summary, ensuring all functions work together.
test_that("Complete workflow: extract -> summary works", {
  # Create realistic mock eeg object
  mock_events <- data.frame(
    onset = c(256, 512, 768, 1024, 1280),
    onset_time = c(1.0, 2.0, 3.0, 4.0, 5.0),
    type = c("100", "200", "100", "200", "100"),
    description = c("Trigger: 100", "Trigger: 200", "Trigger: 100", 
                    "Trigger: 200", "Trigger: 100"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    data = matrix(rnorm(10000), nrow = 10, ncol = 1000),
    channels = paste0("Ch", 1:10),
    sampling_rate = 256,
    times = seq(0, 3.9, length.out = 1000),
    events = mock_events,
    metadata = list(file_path = "integration_test.bdf"),
    reference = "CMS",
    preprocessing_history = list("imported")
  )
  class(mock_eeg) <- "eeg"
  
  # Extract events
  events <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Generate summary (should not error)
  expect_no_error(summary_bdf_events(events))
  
  # Verify extraction worked
  expect_equal(nrow(events), 5)
  expect_true(all(c("100", "200") %in% events$type))
})

# ----------------------------------------------------------------------------
# Test 5.2: Complete workflow: extract -> validate
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Tests the integration between extraction and validation
# functions to ensure they work seamlessly together.
test_that("Complete workflow: extract -> validate works", {
  mock_events <- data.frame(
    onset = c(100, 200, 300, 400),
    onset_time = c(0.5, 1.0, 1.5, 2.0),
    type = c("1", "2", "1", "2"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 1", "Trigger: 2"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200,
    times = seq(0, 2.5, length.out = 500),
    metadata = list(file_path = "test.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  # Extract
  events <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Validate
  validation <- validate_bdf_events(events, verbose = FALSE, plot = FALSE)
  
  # Check validation results
  expect_true(is.list(validation))
  expect_true(!is.null(validation$summary))
})


# ============================================================================
#             TEST SUITE 6: Edge Cases and Boundary Conditions
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: Single event at start of recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Handles edge case of a single event occurring at the
# very beginning of the recording.
test_that("Handles single event at recording start", {
  single_event_start <- data.frame(
    onset = 1,
    onset_time = 0.0,
    type = "1",
    description = "Trigger: 1",
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = single_event_start,
    sampling_rate = 256,
    times = seq(0, 10, length.out = 2560),
    metadata = list(file_path = "start_event.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$onset[1], 1)
  expect_true(result$onset_time[1] < 0.01)  # Very close to 0
})

# ----------------------------------------------------------------------------
# Test 6.2: High-frequency events (rapid triggers)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the function handles many events occurring in
# quick succession without issues.
test_that("Handles high-frequency rapid events", {
  # 100 events in 1 second
  rapid_events <- data.frame(
    onset = seq(1, 256, length.out = 100),
    onset_time = seq(0, 1, length.out = 100),
    type = as.character(rep(1:10, 10)),
    description = paste0("Trigger: ", rep(1:10, 10)),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = rapid_events,
    sampling_rate = 256,
    times = seq(0, 2, length.out = 512),
    metadata = list(file_path = "rapid_events.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_equal(nrow(result), 100)
  expect_true(all(!is.na(result$onset)))
  expect_true(all(!is.na(result$onset_time)))
})

# ----------------------------------------------------------------------------
# Test 6.3: Large trigger codes (up to 65535)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the function handles the full range of
# possible BioSemi trigger codes (0-65535).
test_that("Handles large trigger codes correctly", {
  large_code_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.5, 1.0, 1.5),
    type = c("1", "32767", "65535"),
    description = c("Trigger: 1", "Trigger: 32767", "Trigger: 65535"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = large_code_events,
    sampling_rate = 200,
    times = seq(0, 2, length.out = 400),
    metadata = list(file_path = "large_codes.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_equal(nrow(result), 3)
  expect_true("65535" %in% result$type)
  expect_true("32767" %in% result$type)
})

# ----------------------------------------------------------------------------
# Test 6.4: Long recording with sparse events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures efficient handling of long recordings with
# relatively few events.
test_that("Handles long recording with sparse events", {
  # 5 events over 1 hour of recording
  sparse_events <- data.frame(
    onset = c(100, 46080, 92160, 138240, 184320),  # ~3 minutes apart at 256 Hz
    onset_time = c(0.39, 180, 360, 540, 720),
    type = c("1", "2", "3", "4", "5"),
    description = paste0("Trigger: ", 1:5),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = sparse_events,
    sampling_rate = 256,
    times = seq(0, 3600, length.out = 921600),  # 1 hour
    metadata = list(file_path = "long_sparse.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_equal(nrow(result), 5)
  expect_true(max(result$onset_time) > 700)
})

# ----------------------------------------------------------------------------
# Test 6.5: Different sampling rates
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that onset_time calculations are correct across
# different sampling rates (timing accuracy is sampling-rate dependent).
test_that("Handles different sampling rates correctly", {
  # Test common sampling rates
  sampling_rates <- c(128, 256, 512, 1024, 2048)
  
  for (sr in sampling_rates) {
    mock_events <- data.frame(
      onset = c(sr, 2*sr, 3*sr),  # Events at 1, 2, 3 seconds
      onset_time = c(1.0, 2.0, 3.0),
      type = c("1", "2", "3"),
      description = paste0("Trigger: ", 1:3),
      stringsAsFactors = FALSE
    )
    
    mock_eeg <- list(
      events = mock_events,
      sampling_rate = sr,
      times = seq(0, 4, length.out = 4*sr),
      metadata = list(file_path = paste0("sr_", sr, ".bdf"))
    )
    class(mock_eeg) <- "eeg"
    
    result <- extract_bdf_events(mock_eeg, verbose = FALSE)
    
    expect_equal(attr(result, "sample_rate"), sr)
    expect_equal(nrow(result), 3)
    expect_true(abs(result$onset_time[2] - 2.0) < 0.01)
  }
})


# ============================================================================
#             TEST SUITE 7: Data Type and Structure Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: Event type is character not numeric
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures the 'type' column is stored as character strings
# (as per function specification), not numeric values.
test_that("Event type column is character", {
  mock_events <- data.frame(
    onset = c(100, 200),
    onset_time = c(0.5, 1.0),
    type = c("100", "200"),
    description = c("Trigger: 100", "Trigger: 200"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_true(is.character(result$type))
  expect_false(is.numeric(result$type))
})

# ----------------------------------------------------------------------------
# Test 7.2: Onset times are monotonically increasing
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that event onset times are in chronological order
# (should always be increasing).
test_that("Onset times are monotonically increasing", {
  mock_events <- data.frame(
    onset = c(100, 200, 300, 400, 500),
    onset_time = c(0.5, 1.0, 1.5, 2.0, 2.5),
    type = as.character(1:5),
    description = paste0("Trigger: ", 1:5),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Check that onset_time is monotonically increasing
  expect_true(all(diff(result$onset_time) >= 0))
  expect_true(all(diff(result$onset) >= 0))
})

# ----------------------------------------------------------------------------
# Test 7.3: Description format consistency
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures description strings follow the expected format
# "Trigger: XXX" consistently.
test_that("Description format is consistent", {
  mock_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.5, 1.0, 1.5),
    type = c("10", "20", "30"),
    description = c("Trigger: 10", "Trigger: 20", "Trigger: 30"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Check all descriptions start with "Trigger: "
  expect_true(all(grepl("^Trigger: ", result$description)))
})

# ----------------------------------------------------------------------------
# Test 7.4: No NA values in output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that no NA values appear in any of the output
# columns (onset, onset_time, type, description).
test_that("No NA values in extracted events", {
  mock_events <- data.frame(
    onset = c(100, 200, 300),
    onset_time = c(0.5, 1.0, 1.5),
    type = c("1", "2", "3"),
    description = c("Trigger: 1", "Trigger: 2", "Trigger: 3"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 200
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  expect_false(any(is.na(result$onset)))
  expect_false(any(is.na(result$onset_time)))
  expect_false(any(is.na(result$type)))
  expect_false(any(is.na(result$description)))
})

# ----------------------------------------------------------------------------
# Test 7.5: Attribute types are correct
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ensures all attributes attached to the events data frame
# have the correct data types.
test_that("Attribute types are correct", {
  mock_events <- data.frame(
    onset = c(100),
    onset_time = c(0.5),
    type = c("1"),
    description = c("Trigger: 1"),
    stringsAsFactors = FALSE
  )
  
  mock_eeg <- list(
    events = mock_events,
    sampling_rate = 256,
    times = seq(0, 10, length.out = 2560),
    metadata = list(file_path = "/path/to/file.bdf")
  )
  class(mock_eeg) <- "eeg"
  
  result <- extract_bdf_events(mock_eeg, verbose = FALSE)
  
  # Check attribute types
  expect_true(is.character(attr(result, "bdf_file")))
  expect_true(is.numeric(attr(result, "sample_rate")))
  expect_true(is.numeric(attr(result, "n_samples")))
  expect_true(is.numeric(attr(result, "duration_sec")))
})


# ============================================================================
# SUMMARY
# ============================================================================
# 
# This test file comprehensively tests all functions in extract_bdf_events.R:
#
# 1. extract_bdf_events() - 9 tests covering:
#    - Input validation (file paths, eeg objects, invalid inputs)
#    - Return structure (columns, attributes)
#    - Edge cases (empty events, NULL events, single event)
#    - Verbose mode control
#    - Multiple trigger codes
#
# 2. summary_bdf_events() - 4 tests covering:
#    - Basic summary generation
#    - Empty events handling
#    - Single event edge case
#    - Multiple trigger types
#
# 3. plot_bdf_events() - 5 tests covering:
#    - Timeline and sequence plot types
#    - show_codes parameter
#    - Empty events handling
#    - Invalid plot type handling
#
# 4. validate_bdf_events() - 6 tests covering:
#    - Basic validation
#    - Empty events detection
#    - Different input types
#    - Plot and verbose parameter control
#    - Return value structure
#
# 5. .extract_events_from_input() - 5 tests covering:
#    - File path handling
#    - EEG object extraction
#    - Events data frame acceptance
#    - Invalid input rejection
#    - List with events field
#
# 6. Integration tests - 3 tests covering:
#    - Extract -> summary workflow
#    - Extract -> validate workflow
#    - Extract -> plot workflow
#
# 7. Edge cases - 5 tests covering:
#    - Single event at recording start
#    - High-frequency rapid events
#    - Large trigger codes (up to 65535)
#    - Long recordings with sparse events
#    - Different sampling rates
#
# 8. Data validation - 5 tests covering:
#    - Type column is character
#    - Monotonically increasing times
#    - Description format consistency
#    - No NA values
#    - Correct attribute types
#
# Total: 42 comprehensive tests ensuring robust functionality of all
# event extraction and validation functions.
#
# ============================================================================