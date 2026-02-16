# ============================================================================
# Test File for read_bdf_native.R
# ============================================================================
#
# Tests for native BDF file reader functions
#
# Author: Generated test suite
# Date: Feb 2026
# ============================================================================

library(testthat)
library(eeganalysis)  # Adjust to your package name

# ============================================================================
# Tests for unpack_int24_native()
# ============================================================================

test_that("unpack_int24_native handles empty input", {
  result <- unpack_int24_native(raw(0))
  expect_equal(result, integer(0))
  expect_length(result, 0)
})

test_that("unpack_int24_native handles incomplete triplet", {
  # Only 2 bytes - should return empty
  result <- unpack_int24_native(as.raw(c(0x01, 0x02)))
  expect_equal(result, integer(0))
})

test_that("unpack_int24_native unpacks positive values correctly", {
  # Test value: 0x000001 (decimal 1)
  raw_vec <- as.raw(c(0x01, 0x00, 0x00))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, 1L)
  
  # Test value: 0x000100 (decimal 256)
  raw_vec <- as.raw(c(0x00, 0x01, 0x00))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, 256L)
  
  # Test value: 0x010000 (decimal 65536)
  raw_vec <- as.raw(c(0x00, 0x00, 0x01))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, 65536L)
})

test_that("unpack_int24_native unpacks negative values correctly (two's complement)", {
  # Test value: 0xFFFFFF (decimal -1 in 24-bit two's complement)
  raw_vec <- as.raw(c(0xFF, 0xFF, 0xFF))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, -1L)
  
  # Test value: 0x800000 (most negative: -8388608)
  raw_vec <- as.raw(c(0x00, 0x00, 0x80))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, -8388608L)
})

test_that("unpack_int24_native handles multiple values", {
  # Two values: 1 and 256
  raw_vec <- as.raw(c(0x01, 0x00, 0x00, 0x00, 0x01, 0x00))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, c(1L, 256L))
  expect_length(result, 2)
})

test_that("unpack_int24_native truncates incomplete final triplet", {
  # 7 bytes = 2 complete triplets + 1 incomplete byte
  raw_vec <- as.raw(c(0x01, 0x00, 0x00, 0x02, 0x00, 0x00, 0xFF))
  result <- unpack_int24_native(raw_vec)
  expect_equal(result, c(1L, 2L))
  expect_length(result, 2)
})

# ============================================================================
# Tests for digital_to_physical()
# ============================================================================

test_that("digital_to_physical converts correctly with standard parameters", {
  # Standard BioSemi range
  digital <- c(-8388608L, 0L, 8388607L)
  digital_min <- -8388608
  digital_max <- 8388607
  physical_min <- -262144  # microV
  physical_max <- 262144   # microV
  
  result <- digital_to_physical(digital, digital_min, digital_max, 
                                physical_min, physical_max)
  
  expect_equal(result[1], physical_min, tolerance = 1e-6)
  expect_equal(result[2], 0, tolerance = 1e-6)
  expect_equal(result[3], physical_max, tolerance = 1e-6)
})

test_that("digital_to_physical handles zero range", {
  # Edge case: constant signal
  digital <- c(100, 100, 100)
  digital_min <- 0
  digital_max <- 0  # Division by zero case
  physical_min <- 0
  physical_max <- 100
  
  # Should return Inf or NaN
  result <- digital_to_physical(digital, digital_min, digital_max, 
                                physical_min, physical_max)
  expect_true(all(is.infinite(result) | is.nan(result)))
})

test_that("digital_to_physical is linear", {
  digital <- seq(-1000, 1000, by = 100)
  digital_min <- -1000
  digital_max <- 1000
  physical_min <- -100
  physical_max <- 100
  
  result <- digital_to_physical(digital, digital_min, digital_max, 
                                physical_min, physical_max)
  
  # Check linearity: differences should be constant
  diffs <- diff(result)
  expect_true(all(abs(diffs - diffs[1]) < 1e-10))
})

test_that("digital_to_physical handles empty input", {
  result <- digital_to_physical(integer(0), -1000, 1000, -100, 100)
  expect_equal(result, numeric(0))
  expect_length(result, 0)
})

test_that("digital_to_physical preserves vector length", {
  digital <- 1:1000
  result <- digital_to_physical(digital, 0, 1000, 0, 100)
  expect_length(result, 1000)
})

# ============================================================================
#                     Tests for extract_events_native()
# ============================================================================

test_that("extract_events_native handles empty status signal", {
  result <- extract_events_native(numeric(0), 512)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_true(all(c("onset", "onset_time", "type", "description") %in% names(result)))
})

test_that("extract_events_native handles NULL status signal", {
  result <- extract_events_native(NULL, 512)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("extract_events_native detects single event", {
  # Create signal with one trigger: 0, 0, 0, 100, 100, 0, 0
  status_signal <- c(0, 0, 0, 100, 100, 0, 0)
  sampling_rate <- 512
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$onset[1], 4)  # Index where trigger appears
  expect_equal(result$type[1], 100)
  expect_equal(result$onset_time[1], 3 / sampling_rate)
})

test_that("extract_events_native detects multiple events", {
  # Create signal with three triggers
  status_signal <- c(0, 0, 10, 10, 0, 20, 20, 0, 30, 30, 0)
  sampling_rate <- 512
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  expect_equal(nrow(result), 3)
  expect_equal(result$type, c(10, 20, 30))
  expect_equal(result$onset, c(3, 6, 9))
})

test_that("extract_events_native masks to lower 16 bits", {
  # Create signal with high bits set (simulating BioSemi CMS/DRL bits)
  # 0x00010064 = 65636 = trigger 100 with high bit set
  status_signal <- c(0, 0, 65636, 65636, 0)
  sampling_rate <- 512
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  expect_equal(nrow(result), 1)
  # After masking to lower 16 bits: 65636 & 65535 = 100
  expect_equal(result$type[1], 100)
})

test_that("extract_events_native ignores zero triggers", {
  # Signal with transition to zero (should not be detected as event)
  status_signal <- c(0, 0, 100, 100, 0, 0)
  sampling_rate <- 512
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  # Should only detect onset of trigger 100, not offset to 0
  expect_equal(nrow(result), 1)
  expect_equal(result$type[1], 100)
})

test_that("extract_events_native calculates onset_time correctly", {
  status_signal <- c(0, 0, 0, 0, 100, 100, 0)
  sampling_rate <- 1000  # Use round number for easy checking
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  expect_equal(result$onset[1], 5)
  # onset_time = (onset - 1) / sampling_rate = 4 / 1000 = 0.004
  expect_equal(result$onset_time[1], 0.004)
})

test_that("extract_events_native creates correct description", {
  status_signal <- c(0, 0, 123, 123, 0)
  sampling_rate <- 512
  
  result <- extract_events_native(status_signal, sampling_rate)
  
  expect_equal(result$description[1], "Trigger: 123")
})

# ============================================================================
# Tests for read_bdf_header_native()
# ============================================================================

test_that("read_bdf_header_native fails on non-existent file", {
  expect_error(
    read_bdf_header_native("/path/to/nonexistent/file.bdf"),
    "cannot open the connection"
  )
})

# NOTE: Testing read_bdf_header_native with real files requires actual BDF data
# Consider adding tests with fixture files if available

# ============================================================================
# Tests for read_bdf_native()
# ============================================================================

test_that("read_bdf_native fails on non-existent file", {
  expect_error(
    read_bdf_native("/path/to/nonexistent/file.bdf", verbose = FALSE),
    "File not found"
  )
})

test_that("read_bdf_native warns on non-.bdf extension", {
  # Create temporary file with wrong extension
  temp_file <- tempfile(fileext = ".txt")
  writeLines("dummy", temp_file)
  on.exit(unlink(temp_file))
  
  expect_warning(
    read_bdf_native(temp_file, verbose = FALSE),
    "does not have .bdf extension"
  )
})

# ============================================================================
# Integration Tests (require real BDF files)
# ============================================================================

test_that("read_bdf_native works with real BDF file", {
  skip_if_not(file.exists("tests/testthat/fixtures/sample.bdf"), 
              "Sample BDF file not available")
  
  eeg_obj <- read_bdf_native("tests/testthat/fixtures/sample.bdf", verbose = FALSE)
  
  expect_s3_class(eeg_obj, "eeg")
  expect_true(is.matrix(eeg_obj$data))
  expect_true(nrow(eeg_obj$data) > 0)  # Has channels
  expect_true(ncol(eeg_obj$data) > 0)  # Has samples
  expect_true(length(eeg_obj$channels) > 0)
  expect_true(eeg_obj$sampling_rate > 0)
  expect_true(is.data.frame(eeg_obj$events))
})

test_that("read_bdf_native extracts correct number of channels", {
  skip_if_not(file.exists("tests/testthat/fixtures/sample.bdf"), 
              "Sample BDF file not available")
  
  eeg_obj <- read_bdf_native("tests/testthat/fixtures/sample.bdf", verbose = FALSE)
  
  # Number of channels in data should match length of channel labels
  expect_equal(nrow(eeg_obj$data), length(eeg_obj$channels))
})

test_that("read_bdf_native creates valid metadata", {
  skip_if_not(file.exists("tests/testthat/fixtures/sample.bdf"), 
              "Sample BDF file not available")
  
  eeg_obj <- read_bdf_native("tests/testthat/fixtures/sample.bdf", verbose = FALSE)
  
  expect_true(!is.null(eeg_obj$metadata))
  expect_true(is.list(eeg_obj$metadata))
  expect_true("file_path" %in% names(eeg_obj$metadata))
  expect_true("sampling_rate" %in% names(eeg_obj$metadata))
  expect_true("n_eeg_channels" %in% names(eeg_obj$metadata))
})

test_that("read_bdf_native handles different chunk sizes", {
  skip_if_not(file.exists("tests/testthat/fixtures/sample.bdf"), 
              "Sample BDF file not available")
  
  # Read with different chunk sizes
  eeg_obj_small <- read_bdf_native("tests/testthat/fixtures/sample.bdf", 
                                   verbose = FALSE, chunk_records = 10)
  eeg_obj_large <- read_bdf_native("tests/testthat/fixtures/sample.bdf", 
                                   verbose = FALSE, chunk_records = 500)
  
  # Results should be identical regardless of chunk size
  expect_equal(dim(eeg_obj_small$data), dim(eeg_obj_large$data))
  expect_equal(eeg_obj_small$channels, eeg_obj_large$channels)
  expect_equal(eeg_obj_small$sampling_rate, eeg_obj_large$sampling_rate)
})

# ============================================================================
# Performance and Memory Tests
# ============================================================================

test_that("unpack_int24_native is efficient with large inputs", {
  # Create large vector (1 million samples = 3 million bytes)
  n_samples <- 1e6
  raw_vec <- as.raw(sample(0:255, n_samples * 3, replace = TRUE))
  
  expect_silent({
    result <- unpack_int24_native(raw_vec)
  })
  
  expect_length(result, n_samples)
})

test_that("digital_to_physical handles large vectors efficiently", {
  digital <- sample(-8388608:8388607, 1e6, replace = TRUE)
  
  start_time <- Sys.time()
  result <- digital_to_physical(digital, -8388608, 8388607, -262144, 262144)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  expect_length(result, 1e6)
  expect_true(elapsed < 1.0)  # Should complete in less than 1 second
})

# ============================================================================
# Edge Cases and Error Handling
# ============================================================================

test_that("unpack_int24_native handles maximum and minimum values", {
  # Maximum positive: 0x7FFFFF (8388607)
  raw_max <- as.raw(c(0xFF, 0xFF, 0x7F))
  result_max <- unpack_int24_native(raw_max)
  expect_equal(result_max, 8388607L)
  
  # Minimum negative: 0x800000 (-8388608)
  raw_min <- as.raw(c(0x00, 0x00, 0x80))
  result_min <- unpack_int24_native(raw_min)
  expect_equal(result_min, -8388608L)
})

test_that("extract_events_native handles constant non-zero signal", {
  # All values the same (no events)
  status_signal <- rep(100, 100)
  result <- extract_events_native(status_signal, 512)
  
  # First sample might be detected as event if it's non-zero
  expect_true(nrow(result) <= 1)
})

test_that("extract_events_native handles very long signal", {
  # Simulate 1 hour of data at 512 Hz
  n_samples <- 512 * 3600
  status_signal <- numeric(n_samples)
  # Add a few events
  status_signal[c(1000, 50000, 100000)] <- c(10, 20, 30)
  
  result <- extract_events_native(status_signal, 512)
  
  expect_equal(nrow(result), 3)
  expect_true(all(result$onset %in% c(1000, 50000, 100000)))
})

# ============================================================================
# Documentation and Type Checking
# ============================================================================

test_that("functions return correct types", {
  # unpack_int24_native should return integer vector
  result1 <- unpack_int24_native(as.raw(c(0x01, 0x00, 0x00)))
  expect_type(result1, "integer")
  
  # digital_to_physical should return numeric vector
  result2 <- digital_to_physical(c(0, 100, 200), 0, 1000, 0, 100)
  expect_type(result2, "double")
  
  # extract_events_native should return data.frame
  result3 <- extract_events_native(c(0, 100, 0), 512)
  expect_s3_class(result3, "data.frame")
})

test_that("extract_events_native data frame has correct column types", {
  status_signal <- c(0, 100, 100, 0)
  result <- extract_events_native(status_signal, 512)
  
  if (nrow(result) > 0) {
    expect_type(result$onset, "integer")
    expect_type(result$onset_time, "double")
    expect_type(result$type, "integer")
    expect_type(result$description, "character")
  }
})