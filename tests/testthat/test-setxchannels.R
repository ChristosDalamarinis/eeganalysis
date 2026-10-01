# ============================================================================
#                   Test File for setexchannels.R
# ============================================================================
#
# This test file tests detect_external_channels() - a non-interactive,
# database-driven suggestion tool: does a channel name match a known
# BioSemi auxiliary port? It does not read or write channel_types, and
# does not rename anything.
#
# identify_external_channels() and apply_external_labels() used to be
# tested here too (an interactive-labeling-then-rename workflow built
# before channel_types existed on the eeg object). Both are retired: see
# R/setexchannels.R's own header, and set_channel_types() /
# read_bdf_native()'s eog=/misc= arguments for what replaced them.
#
# Function tested:
#   1. detect_external_channels() - Non-interactive detection of external channels
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)


# ============================================================================
#   TEST SUITE 1: detect_external_channels() - Non-Interactive Detection
# ============================================================================
#
# detect_external_channels() accepts a character vector, data frame, matrix,
# or list and returns a character vector of channel names whose position_type
# in the electrode database is one of:
# "External", "GSR", "Ergo/AUX", "Respiration", "Plethysmograph", "Temperature"
#
# It uses get_electrode_database() internally and performs case-insensitive
# lookups by converting each channel name to lowercase before matching.
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Character vector input - only EEG channels, no externals
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a character vector of standard EEG channel names is
# provided (e.g., "Cz", "Pz", "Fz"), detect_external_channels() should return
# an empty character vector because none of them are external channels in the
# electrode database. This verifies the base-case detection logic does not
# produce false positives.
test_that("detect_external_channels returns empty vector for pure EEG channels", {
  eeg_only_channels <- c("Fp1", "Fp2", "F3", "F4", "Fz", "C3", "C4", "Cz",
                         "P3", "P4", "Pz", "O1", "O2", "Oz", "T7", "T8")

  result <- detect_external_channels(eeg_only_channels)

  expect_true(is.character(result))
  expect_equal(length(result), 0)
})

# ----------------------------------------------------------------------------
# Test 1.2: Character vector input - mix of EEG and external channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a mix of standard EEG channels (e.g., "Cz") and known
# external BioSemi channels (e.g., "EXG1", "EXG2") is provided as a character
# vector, detect_external_channels() should return only the external channels.
# This confirms that the function correctly filters on position_type and does
# not include regular electrodes in its output.
test_that("detect_external_channels correctly isolates external channels from a mixed character vector", {
  mixed_channels <- c("Cz", "Pz", "Fz", "EXG1", "EXG2", "O1", "O2", "EXG3")

  result <- detect_external_channels(mixed_channels)

  expect_true(is.character(result))
  expect_equal(length(result), 3)
  expect_true("EXG1" %in% result)
  expect_true("EXG2" %in% result)
  expect_true("EXG3" %in% result)
  # EEG channels must NOT be included
  expect_false("Cz" %in% result)
  expect_false("Pz" %in% result)
  expect_false("Fz" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.3: Character vector input - all eight EXG electrodes detected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that all eight BioSemi external electrode slots
# (EXG1 through EXG8) are correctly identified as external channels. Their
# position_type in the database is "External", so each one must appear in
# the returned vector.
test_that("detect_external_channels detects all EXG1-EXG8 channels", {
  exg_channels <- c("EXG1", "EXG2", "EXG3", "EXG4",
                    "EXG5", "EXG6", "EXG7", "EXG8")

  result <- detect_external_channels(exg_channels)

  expect_equal(length(result), 8)
  for (ch in exg_channels) {
    expect_true(ch %in% result, info = paste(ch, "should be detected as external"))
  }
})

# ----------------------------------------------------------------------------
# Test 1.4: Character vector input - GSR, Plethysmograph, and Temperature
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies detection of BioSemi accessory channels that carry
# non-"External" position_types but are still treated as external: GSR1/GSR2
# (type "GSR"), Plet (type "Plethysmograph"), and Temp (type "Temperature").
# All four must be returned even though none are EXG-type electrodes.
test_that("detect_external_channels detects GSR, Plet, and Temp accessory channels", {
  accessory_channels <- c("GSR1", "GSR2", "Plet", "Temp", "Cz", "Pz")

  result <- detect_external_channels(accessory_channels)

  expect_true("GSR1" %in% result)
  expect_true("GSR2" %in% result)
  expect_true("Plet" %in% result)
  expect_true("Temp" %in% result)
  # EEG channels must not appear
  expect_false("Cz" %in% result)
  expect_false("Pz" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.5: Case-insensitive matching for character vector
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The function converts channel names to lowercase before
# looking them up in the electrode database. This test confirms that channel
# names provided in various cases (e.g., "exg1", "EXG1", "Exg1") are all
# detected as external channels, demonstrating that the tolower() logic works.
test_that("detect_external_channels performs case-insensitive matching", {
  # Provide channel names in different cases
  mixed_case <- c("exg1", "EXG2", "Exg3", "gsr1", "GSR2", "Cz")

  result <- detect_external_channels(mixed_case)

  # The function should match despite casing and return the original-cased names
  # All five non-EEG channels should be found
  expect_equal(length(result), 5)
  expect_false("Cz" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.6: Data frame input - uses column names by default
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a data frame is passed without specifying channel_col,
# detect_external_channels() falls back to reading colnames(data). This test
# verifies that a data frame whose column names include "EXG1" and "EXG2"
# yields those two channels in the result.
test_that("detect_external_channels reads column names from a data frame", {
  df <- data.frame(
    Cz   = rnorm(10),
    Pz   = rnorm(10),
    EXG1 = rnorm(10),
    EXG2 = rnorm(10),
    Fz   = rnorm(10)
  )

  result <- detect_external_channels(df)

  expect_true("EXG1" %in% result)
  expect_true("EXG2" %in% result)
  expect_false("Cz" %in% result)
  expect_false("Fz" %in% result)
  expect_equal(length(result), 2)
})

# ----------------------------------------------------------------------------
# Test 1.7: Data frame input - using the channel_col argument
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a data frame has a dedicated column containing channel
# names as strings (rather than using column names), channel_col can point to
# it. This test verifies that detect_external_channels() correctly reads that
# column, extracts unique values, and identifies any external channels among
# them. Only unique channel names should be considered (not every row).
test_that("detect_external_channels reads channel names from a specified channel_col", {
  df <- data.frame(
    channel = c("Cz", "Fz", "EXG1", "EXG2", "Cz", "EXG1"),
    amplitude = rnorm(6),
    stringsAsFactors = FALSE
  )

  result <- detect_external_channels(df, channel_col = "channel")

  expect_true("EXG1" %in% result)
  expect_true("EXG2" %in% result)
  expect_false("Cz" %in% result)
  # Duplicates in the column should not lead to duplicate entries in result
  expect_equal(length(result), 2)
})

# ----------------------------------------------------------------------------
# Test 1.8: Matrix input - uses column names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A numeric matrix with named columns is accepted by the
# function the same way a data frame is. This test confirms that colnames()
# are read correctly from a matrix and external channels are detected.
test_that("detect_external_channels reads column names from a matrix", {
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  colnames(mat) <- c("Cz", "Pz", "EXG1", "EXG4", "Oz")

  result <- detect_external_channels(mat)

  expect_true("EXG1" %in% result)
  expect_true("EXG4" %in% result)
  expect_false("Cz" %in% result)
  expect_false("Oz" %in% result)
  expect_equal(length(result), 2)
})

# ----------------------------------------------------------------------------
# Test 1.9: List input with $channels slot
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A list (such as an eeg object) can carry channel names in
# a $channels slot. detect_external_channels() checks for this slot first in
# list input. This test confirms those channels are read, trimmed, and
# processed correctly, and that external channels are returned.
test_that("detect_external_channels reads channel names from list$channels", {
  lst <- list(
    channels = c("Cz", "Pz", "EXG1", "EXG2", "Oz"),
    data = matrix(rnorm(50), nrow = 5)
  )

  result <- detect_external_channels(lst)

  expect_true("EXG1" %in% result)
  expect_true("EXG2" %in% result)
  expect_false("Cz" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.10: List input with $channel_names slot (fallback)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If a list does not have a $channels slot, the function falls
# back to $channel_names. This test verifies that fallback path works and
# correctly identifies external channels stored under the alternative slot name.
test_that("detect_external_channels falls back to list$channel_names", {
  lst <- list(
    channel_names = c("Fz", "GSR1", "Temp", "T7"),
    sampling_rate = 512
  )

  result <- detect_external_channels(lst)

  expect_true("GSR1" %in% result)
  expect_true("Temp" %in% result)
  expect_false("Fz" %in% result)
  expect_false("T7" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.11: eeg object (S3 class list with $channels) is handled correctly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: An eeg S3 object is a list with class "eeg" and a $channels
# slot. detect_external_channels() treats it the same as a plain list. This
# test creates a real eeg object using new_eeg() and verifies that external
# channels in $channels are correctly detected, confirming that the function
# works within a typical eeganalysis workflow.
test_that("detect_external_channels works correctly with an eeg object", {
  eeg_data <- new_eeg(
    data = matrix(rnorm(500), nrow = 5, ncol = 100),
    channels = c("Cz", "Pz", "EXG1", "EXG2", "Fz"),
    sampling_rate = 256
  )

  result <- detect_external_channels(eeg_data)

  expect_true("EXG1" %in% result)
  expect_true("EXG2" %in% result)
  expect_false("Cz" %in% result)
  expect_false("Pz" %in% result)
  expect_false("Fz" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.12: channel_col not found in data frame - stops with error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If the user specifies a channel_col that doesn't exist as a
# column in the data frame, detect_external_channels() should throw an error
# mentioning the missing column name. This prevents silent failures from
# mismatched column names.
test_that("detect_external_channels errors when channel_col is not in data frame", {
  df <- data.frame(channel = c("Cz", "EXG1"), value = c(1, 2))

  expect_error(
    detect_external_channels(df, channel_col = "nonexistent_col"),
    "nonexistent_col"
  )
})

# ----------------------------------------------------------------------------
# Test 1.13: List without $channels or $channel_names - stops with error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A list that contains neither $channels nor $channel_names
# provides no way for the function to extract channel names. It should throw
# a descriptive error rather than silently returning an empty result or crashing
# with an obscure message.
test_that("detect_external_channels errors on list missing channel slots", {
  bad_list <- list(data = matrix(rnorm(10)), sampling_rate = 256)

  expect_error(
    detect_external_channels(bad_list),
    regexp = "Cannot find channel names"
  )
})

# ----------------------------------------------------------------------------
# Test 1.14: Invalid input type - stops with error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When detect_external_channels() receives a type it cannot
# handle (e.g., a numeric scalar, a logical vector), it should stop with an
# informative error. This guards the function against misuse.
test_that("detect_external_channels errors on invalid input types", {
  expect_error(detect_external_channels(42))
  expect_error(detect_external_channels(TRUE))
})

# ----------------------------------------------------------------------------
# Test 1.15: Whitespace trimming in channel names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The function calls trimws() on extracted channel names.
# This test provides channels with leading/trailing spaces and verifies they
# are still correctly matched against the database. Without trimming, " EXG1"
# would not match "exg1" and would be silently missed.
test_that("detect_external_channels trims whitespace from channel names", {
  # Channels with extra whitespace - trimws() should handle these
  channels_with_spaces <- c("  EXG1  ", " EXG2", "Cz ", " Pz ")

  result <- detect_external_channels(channels_with_spaces)

  # EXG1 and EXG2 should be found after trimming
  expect_equal(length(result), 2)
  expect_false("Cz " %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.16: Unknown channel names not in the database are silently ignored
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Channel names that are not present in the electrode database
# at all (e.g., custom channel names like "MyCustomChannel") should simply be
# skipped - not error, and not returned as external. This is the correct
# behaviour since only database-known external channels should be flagged.
test_that("detect_external_channels silently ignores channels not in the database", {
  unknown_channels <- c("MyCustomChannel", "Sensor99", "UNKNOWN", "EXG1")

  result <- detect_external_channels(unknown_channels)

  # Only EXG1 is in the database with External type
  expect_equal(length(result), 1)
  expect_true("EXG1" %in% result)
})

# ----------------------------------------------------------------------------
# Test 1.17: Return type is always a character vector
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Regardless of whether any external channels are found,
# detect_external_channels() must always return a character vector (never NULL,
# a list, or a logical). This ensures downstream code can safely call
# length(), %in%, etc. without type guards.
test_that("detect_external_channels always returns a character vector", {
  # Case with no externals
  result_empty <- detect_external_channels(c("Cz", "Pz", "Fz"))
  expect_true(is.character(result_empty))

  # Case with externals
  result_found <- detect_external_channels(c("EXG1", "Cz"))
  expect_true(is.character(result_found))
})
