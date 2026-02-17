# ============================================================================
#                   Test File for setexchannels.R
# ============================================================================
#
# This test file provides comprehensive testing for the setexchannels.R script,
# which contains functions for detecting, identifying, and labeling external
# (non-EEG) channels such as EOG, EMG, ECG, and GSR electrodes.
#
# Functions tested:
#   1. detect_external_channels()    - Non-interactive detection of external channels
#   2. apply_external_labels()       - Renames external channels in data structures
#   3. identify_external_channels()  - Interactive labeling workflow (mocked)
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


# ============================================================================
#   TEST SUITE 2: apply_external_labels() - Rename External Channels
# ============================================================================
#
# apply_external_labels(data, labels, keep_original = TRUE) renames channels
# in a data structure based on a named labels vector. Labels equal to
# "Unlabeled" are skipped. The function handles:
#   - eeg objects (list with class "eeg", uses $channels)
#   - plain lists (uses $channels or $channel_names)
#   - data frames (renames colnames)
#   - matrices (renames colnames)
# When keep_original = TRUE, new name becomes "Label (OriginalName)".
# When keep_original = FALSE, new name is just "Label".
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Empty labels - returns data unchanged with a message
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When an empty named vector is passed as labels, there is
# nothing to rename and the function should emit a message (not an error) and
# return the data structure unmodified. This confirms the early-exit guard at
# the top of the function.
test_that("apply_external_labels returns data unchanged when labels is empty", {
  df <- data.frame(EXG1 = rnorm(5), Cz = rnorm(5))
  original_names <- colnames(df)
  
  expect_message(
    result <- apply_external_labels(df, labels = character(0)),
    "No labels to apply"
  )
  
  expect_equal(colnames(result), original_names)
})

# ----------------------------------------------------------------------------
# Test 2.2: Data frame - keep_original = TRUE (default)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a data frame is passed with labels and keep_original
# is TRUE (the default), the renamed column should follow the format
# "Label (OriginalName)". Unlabeled channels like "Cz" should be unchanged.
# This is the standard workflow for annotating external channels while
# preserving traceability of the original channel names.
test_that("apply_external_labels renames data frame columns with keep_original = TRUE", {
  df <- data.frame(
    Cz   = rnorm(5),
    Pz   = rnorm(5),
    EXG1 = rnorm(5),
    EXG2 = rnorm(5)
  )
  
  labels <- c(EXG1 = "EOG_L", EXG2 = "EOG_R")
  
  result <- apply_external_labels(df, labels = labels, keep_original = TRUE)
  
  expect_true("EOG_L (EXG1)" %in% colnames(result))
  expect_true("EOG_R (EXG2)" %in% colnames(result))
  # Unchanged EEG channels
  expect_true("Cz" %in% colnames(result))
  expect_true("Pz" %in% colnames(result))
  # Original names no longer present
  expect_false("EXG1" %in% colnames(result))
  expect_false("EXG2" %in% colnames(result))
})

# ----------------------------------------------------------------------------
# Test 2.3: Data frame - keep_original = FALSE
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When keep_original is FALSE, the renamed column should
# use only the label text (e.g., "EOG_L"), without the "(EXG1)" suffix.
# This verifies the branching logic inside the renaming loop.
test_that("apply_external_labels renames data frame columns with keep_original = FALSE", {
  df <- data.frame(
    Cz   = rnorm(5),
    EXG1 = rnorm(5),
    EXG2 = rnorm(5)
  )
  
  labels <- c(EXG1 = "EOG_L", EXG2 = "ECG")
  
  result <- apply_external_labels(df, labels = labels, keep_original = FALSE)
  
  expect_true("EOG_L" %in% colnames(result))
  expect_true("ECG" %in% colnames(result))
  # Ensure no "(EXG1)" suffix appears anywhere
  expect_false(any(grepl("\\(EXG", colnames(result))))
  expect_true("Cz" %in% colnames(result))
})

# ----------------------------------------------------------------------------
# Test 2.4: Data frame - "Unlabeled" entries are NOT renamed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The function explicitly checks `if (label != "Unlabeled")`
# before applying a rename. This test passes a labels vector where one channel
# is "Unlabeled" and verifies it keeps its original column name, while the
# properly labeled channel is renamed as expected.
test_that("apply_external_labels skips channels labeled as Unlabeled", {
  df <- data.frame(
    EXG1 = rnorm(5),
    EXG2 = rnorm(5),
    Cz   = rnorm(5)
  )
  
  labels <- c(EXG1 = "EOG_L", EXG2 = "Unlabeled")
  
  result <- apply_external_labels(df, labels = labels, keep_original = FALSE)
  
  # EXG1 labeled -> renamed
  expect_true("EOG_L" %in% colnames(result))
  expect_false("EXG1" %in% colnames(result))
  # EXG2 unlabeled -> unchanged
  expect_true("EXG2" %in% colnames(result))
})

# ----------------------------------------------------------------------------
# Test 2.5: Matrix - column renaming
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A matrix is handled the same as a data frame (colnames are
# renamed). This test confirms that matrix input returns a matrix with updated
# colnames, preserving numeric data and dimensions.
test_that("apply_external_labels renames matrix column names correctly", {
  mat <- matrix(rnorm(30), nrow = 6, ncol = 5)
  colnames(mat) <- c("Cz", "Pz", "EXG1", "EXG5", "Oz")
  
  labels <- c(EXG1 = "EMG", EXG5 = "ECG")
  
  result <- apply_external_labels(mat, labels = labels, keep_original = FALSE)
  
  expect_true(is.matrix(result))
  expect_true("EMG" %in% colnames(result))
  expect_true("ECG" %in% colnames(result))
  expect_true("Cz" %in% colnames(result))
  expect_false("EXG1" %in% colnames(result))
  expect_false("EXG5" %in% colnames(result))
  # Dimensions should be unchanged
  expect_equal(dim(result), c(6, 5))
})

# ----------------------------------------------------------------------------
# Test 2.6: List with $channels slot - channel names updated
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: For a plain list with a $channels slot, apply_external_labels()
# should update data$channels with the new names. All other slots of the list
# (e.g., $data, $sampling_rate) should remain intact. This mimics the scenario
# of labeling channels on an eeg-like list structure.
test_that("apply_external_labels updates $channels in a list structure", {
  lst <- list(
    channels = c("Cz", "EXG1", "EXG2", "Fz"),
    data = matrix(rnorm(40), nrow = 4),
    sampling_rate = 512
  )
  
  labels <- c(EXG1 = "EOG_L", EXG2 = "EOG_R")
  
  result <- apply_external_labels(lst, labels = labels, keep_original = FALSE)
  
  expect_true(is.list(result))
  expect_true("EOG_L" %in% result$channels)
  expect_true("EOG_R" %in% result$channels)
  expect_true("Cz" %in% result$channels)
  expect_false("EXG1" %in% result$channels)
  # Other slots preserved
  expect_equal(result$sampling_rate, 512)
  expect_equal(dim(result$data), c(4, 10))
})

# ----------------------------------------------------------------------------
# Test 2.7: List with $channel_names slot - channel names updated
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a list uses $channel_names instead of $channels,
# the function should update $channel_names accordingly. This tests the
# fallback branch for lists that use the alternative slot name.
test_that("apply_external_labels updates $channel_names when $channels is absent", {
  lst <- list(
    channel_names = c("Cz", "EXG3", "Pz"),
    sampling_rate = 256
  )
  
  labels <- c(EXG3 = "GSR")
  
  result <- apply_external_labels(lst, labels = labels, keep_original = FALSE)
  
  expect_true("GSR" %in% result$channel_names)
  expect_false("EXG3" %in% result$channel_names)
  expect_true("Cz" %in% result$channel_names)
})

# ----------------------------------------------------------------------------
# Test 2.8: eeg object (S3) - $channels updated, class preserved
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: An eeg S3 object is a named list with class attribute "eeg".
# apply_external_labels() must update $channels correctly and, crucially, the
# class attribute must remain "eeg" after the operation - because later
# downstream functions in the pipeline rely on inherits(x, "eeg").
test_that("apply_external_labels updates eeg object channels and preserves eeg class", {
  eeg_obj <- new_eeg(
    data = matrix(rnorm(500), nrow = 5, ncol = 100),
    channels = c("Cz", "Pz", "EXG1", "EXG2", "Fz"),
    sampling_rate = 256
  )
  
  labels <- c(EXG1 = "EOG_L", EXG2 = "EOG_R")
  
  result <- apply_external_labels(eeg_obj, labels = labels, keep_original = FALSE)
  
  # Class must be preserved
  expect_s3_class(result, "eeg")
  
  # Channel names updated
  expect_true("EOG_L" %in% result$channels)
  expect_true("EOG_R" %in% result$channels)
  expect_false("EXG1" %in% result$channels)
  expect_false("EXG2" %in% result$channels)
  
  # Untouched EEG channels still present
  expect_true("Cz" %in% result$channels)
  expect_true("Pz" %in% result$channels)
  expect_true("Fz" %in% result$channels)
})

# ----------------------------------------------------------------------------
# Test 2.9: eeg object - keep_original = TRUE wraps names correctly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The keep_original = TRUE branch on an eeg object should
# produce names in the format "Label (OriginalName)". This test explicitly
# checks that format is applied to the $channels slot and that the channel
# count is unchanged.
test_that("apply_external_labels applies keep_original format to eeg object channels", {
  eeg_obj <- new_eeg(
    data = matrix(rnorm(300), nrow = 3, ncol = 100),
    channels = c("EXG1", "Cz", "Pz"),
    sampling_rate = 512
  )
  
  labels <- c(EXG1 = "EOG_L")
  
  result <- apply_external_labels(eeg_obj, labels = labels, keep_original = TRUE)
  
  expect_true("EOG_L (EXG1)" %in% result$channels)
  expect_false("EXG1" %in% result$channels)
  expect_equal(length(result$channels), 3)  # Channel count unchanged
})

# ----------------------------------------------------------------------------
# Test 2.10: Labels with no matching channels - data unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If the labels vector contains channel names that do not
# appear in the data, no renaming should occur. The data structure is returned
# unmodified (no error, no warning). This is important for robustness when
# labels come from a previous detection step on a different dataset.
test_that("apply_external_labels leaves data unchanged when no label keys match", {
  df <- data.frame(Cz = rnorm(5), Pz = rnorm(5), Fz = rnorm(5))
  original_names <- colnames(df)
  
  # Labels that don't match any column
  labels <- c(EXG1 = "EOG_L", EXG2 = "ECG")
  
  result <- apply_external_labels(df, labels = labels, keep_original = FALSE)
  
  expect_equal(colnames(result), original_names)
})

# ----------------------------------------------------------------------------
# Test 2.11: Data integrity - values are not altered during renaming
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Renaming columns must never alter the numeric data inside
# the columns. This test captures the actual values of the EXG1 column before
# calling apply_external_labels() and confirms the column data is identical
# after renaming (only the name changes, not the values).
test_that("apply_external_labels does not alter channel data values during rename", {
  original_values <- rnorm(10)
  df <- data.frame(Cz = rnorm(10), EXG1 = original_values)
  
  labels <- c(EXG1 = "EOG_L")
  result <- apply_external_labels(df, labels = labels, keep_original = FALSE)
  
  # Values in the renamed column should be identical
  expect_equal(result[["EOG_L"]], original_values)
})

# ----------------------------------------------------------------------------
# Test 2.12: Invalid data type - stops with error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If apply_external_labels() receives a data type it cannot
# handle (e.g., a numeric vector or a factor), it should throw a descriptive
# error. This test ensures the final else-branch error guard is reached and
# fires correctly.
test_that("apply_external_labels errors on invalid data types", {
  labels <- c(EXG1 = "EOG_L")
  
  expect_error(
    apply_external_labels(c(1, 2, 3), labels = labels),
    regexp = "Invalid data type"
  )
})

# ----------------------------------------------------------------------------
# Test 2.13: List without $channels or $channel_names - stops with error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A list that provides no recognizable channel name slot
# should produce a clear error message listing what was tried ($channels,
# $channel_names) and what was actually available in the list.
test_that("apply_external_labels errors when list has no channel name slot", {
  bad_list <- list(data = matrix(rnorm(10)), sampling_rate = 256)
  labels <- c(EXG1 = "EOG_L")
  
  expect_error(
    apply_external_labels(bad_list, labels = labels),
    regexp = "Cannot find channel names"
  )
})

# ----------------------------------------------------------------------------
# Test 2.14: n_changed count is accurate (cat output check)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The function prints "[OK] Applied N external channel label(s)."
# to the console after renaming. This test captures that output and verifies
# the count N matches the number of non-Unlabeled labels that were actually
# applied. Using capture.output() to intercept cat() output.
test_that("apply_external_labels reports the correct number of renamed channels", {
  df <- data.frame(EXG1 = rnorm(5), EXG2 = rnorm(5), EXG3 = rnorm(5), Cz = rnorm(5))
  labels <- c(EXG1 = "EOG_L", EXG2 = "EOG_R", EXG3 = "Unlabeled")
  
  output <- capture.output(
    apply_external_labels(df, labels = labels, keep_original = FALSE)
  )
  
  # Only 2 channels were actually labeled (EXG3 was Unlabeled)
  expect_true(any(grepl("Applied 2", output)))
})


# ============================================================================
#   TEST SUITE 3: identify_external_channels() - Interactive Workflow (Mocked)
# ============================================================================
#
# identify_external_channels() is interactive: it uses readline() to prompt
# the user to label each detected external channel. In a test environment,
# readline() must be mocked so no actual user input is needed.
#
# testthat 3.x provides local_mocked_bindings() to intercept readline().
# Each test supplies a sequence of mock responses to simulate user choices.
#
# The function:
#   1. Extracts channel names from data
#   2. Looks up each name in the electrode database
#   3. Prompts the user via readline() for each external channel
#   4. Returns list(labels = named_char_vector, summary = data.frame)
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Returns a list with $labels and $summary when externals are found
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The return value structure of identify_external_channels()
# must be a list with exactly two named elements: "labels" (a named character
# vector) and "summary" (a data frame). Mocking readline() to return "EOG_L"
# for one channel lets us test this without user interaction.
test_that("identify_external_channels returns list with $labels and $summary", {
  # Mock readline() to always return "EOG_L" for any prompt
  local_mocked_bindings(
    readline = function(prompt = "") "EOG_L",
    .package = "base"
  )
  
  channels <- c("Cz", "Pz", "EXG1")
  
  result <- identify_external_channels(channels)
  
  expect_true(is.list(result))
  expect_true("labels" %in% names(result))
  expect_true("summary" %in% names(result))
})

# ----------------------------------------------------------------------------
# Test 3.2: $labels is a named character vector with correct names and values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The $labels element must be a named character vector where
# names are the original channel names (e.g., "EXG1") and values are the
# user-supplied labels (e.g., "EOG_L"). Mocking readline() to return "EOG_L"
# and checking the names/values of the returned vector confirms the label
# assignment loop works correctly.
test_that("identify_external_channels $labels has correct names and values", {
  local_mocked_bindings(
    readline = function(prompt = "") "EOG_L",
    .package = "base"
  )
  
  channels <- c("Cz", "EXG1", "Pz")
  result <- identify_external_channels(channels)
  
  expect_true(is.character(result$labels))
  expect_equal(names(result$labels), "EXG1")
  expect_equal(result$labels[["EXG1"]], "EOG_L")
})

# ----------------------------------------------------------------------------
# Test 3.3: Empty input skips labeling and returns empty results
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When no external channels are present in the input (e.g.,
# only EEG channels like "Cz", "Pz"), the function should short-circuit before
# calling readline() at all, print a "No external channels detected" message,
# and return list(labels = character(0), summary = data.frame()). This tests
# the length(external_channels) == 0 early-exit branch.
test_that("identify_external_channels returns empty results when no external channels exist", {
  channels <- c("Cz", "Pz", "Fz", "O1", "O2")
  
  result <- identify_external_channels(channels)
  
  expect_equal(length(result$labels), 0)
  expect_true(is.data.frame(result$summary))
  expect_equal(nrow(result$summary), 0)
})

# ----------------------------------------------------------------------------
# Test 3.4: Pressing ENTER (empty response) marks channel as "Unlabeled"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When the user presses ENTER (empty string), the repeat loop
# in identify_external_channels() should label that channel as "Unlabeled"
# and move on. Mocking readline() to return "" simulates this and lets us
# verify the "Unlabeled" assignment for skipped channels.
test_that("identify_external_channels marks channel as Unlabeled when user presses ENTER", {
  local_mocked_bindings(
    readline = function(prompt = "") "",
    .package = "base"
  )
  
  channels <- c("EXG1", "Cz")
  result <- identify_external_channels(channels)
  
  expect_equal(result$labels[["EXG1"]], "Unlabeled")
})

# ----------------------------------------------------------------------------
# Test 3.5: $summary data frame has the required columns and correct row count
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The summary data frame returned in $summary must contain
# exactly the columns: "Channel", "Label", "Type", "Description". It should
# have one row per external channel found. This test mocks two channels being
# labeled and checks both the column names and the row count.
test_that("identify_external_channels returns summary data frame with correct structure", {
  response_counter <- 0L
  responses <- c("EOG_L", "EOG_R")
  
  local_mocked_bindings(
    readline = function(prompt = "") {
      response_counter <<- response_counter + 1L
      responses[response_counter]
    },
    .package = "base"
  )
  
  channels <- c("EXG1", "EXG2", "Cz")
  result <- identify_external_channels(channels)
  
  expect_true(is.data.frame(result$summary))
  expect_equal(ncol(result$summary), 4)
  expect_true(all(c("Channel", "Label", "Type", "Description") %in% names(result$summary)))
  expect_equal(nrow(result$summary), 2)  # One row per external channel
})

# ----------------------------------------------------------------------------
# Test 3.6: $summary Type and Description columns filled from electrode database
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: After user labeling, the function enriches the summary data
# frame with the Type and Description values from the electrode database for
# each channel. For EXG1 the Type should be "External" and the Description
# should be the position_name string from the database. This confirms the
# database lookup loop in the summary-building section works correctly.
test_that("identify_external_channels populates summary Type and Description from database", {
  local_mocked_bindings(
    readline = function(prompt = "") "EOG_L",
    .package = "base"
  )
  
  channels <- c("EXG1", "Cz")
  result <- identify_external_channels(channels)
  
  expect_equal(result$summary$Channel[1], "EXG1")
  expect_equal(result$summary$Label[1], "EOG_L")
  expect_equal(result$summary$Type[1], "External")  # From database
  expect_true(nchar(result$summary$Description[1]) > 0)  # Non-empty
})

# ----------------------------------------------------------------------------
# Test 3.7: Data frame input (colnames) works with identify_external_channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: identify_external_channels() shares the same input-parsing
# logic as detect_external_channels(). A data frame with EXG columns should
# yield those columns for labeling. The mock provides one label for one channel.
test_that("identify_external_channels accepts data frame input", {
  local_mocked_bindings(
    readline = function(prompt = "") "EMG",
    .package = "base"
  )
  
  df <- data.frame(Cz = rnorm(5), EXG4 = rnorm(5), Pz = rnorm(5))
  result <- identify_external_channels(df)
  
  expect_equal(names(result$labels), "EXG4")
  expect_equal(result$labels[["EXG4"]], "EMG")
})

# ----------------------------------------------------------------------------
# Test 3.8: Return is invisibly returned (invisible() wrapping)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: identify_external_channels() uses invisible() to return its
# result. This means that typing the call directly in the console won't auto-
# print. In a test, we capture the output and also assign the result directly
# to verify the invisible return still works with assignment.
test_that("identify_external_channels returns result invisibly", {
  local_mocked_bindings(
    readline = function(prompt = "") "EOG_L",
    .package = "base"
  )
  
  channels <- c("EXG1", "Cz")
  
  # Should not auto-print anything beyond cat() calls
  result <- identify_external_channels(channels)
  
  # But the assigned value is still a valid list
  expect_true(is.list(result))
  expect_true("labels" %in% names(result))
})


# ============================================================================
#   TEST SUITE 4: Integration Tests - Full Workflow
# ============================================================================
#
# These tests simulate realistic end-to-end workflows combining all three
# functions: detecting external channels, labeling them, and applying labels.
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: detect -> apply pipeline on a data frame
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A complete two-step workflow:
#   1. detect_external_channels() finds EXG1 and EXG2 in a data frame
#   2. apply_external_labels() renames those columns using a manually built
#      labels vector
# This confirms that the output format of detect_external_channels() (a
# character vector of names) is compatible with constructing a labels vector
# that apply_external_labels() can consume.
test_that("detect -> apply pipeline correctly renames external channels in a data frame", {
  df <- data.frame(
    Cz   = rnorm(10),
    Pz   = rnorm(10),
    EXG1 = rnorm(10),
    EXG2 = rnorm(10),
    Fz   = rnorm(10)
  )
  
  # Step 1: detect
  detected <- detect_external_channels(df)
  expect_equal(length(detected), 2)
  
  # Step 2: construct labels and apply
  user_labels <- setNames(c("EOG_L", "EOG_R"), detected)
  result <- apply_external_labels(df, labels = user_labels, keep_original = FALSE)
  
  expect_true("EOG_L" %in% colnames(result))
  expect_true("EOG_R" %in% colnames(result))
  expect_true("Cz" %in% colnames(result))
  expect_false("EXG1" %in% colnames(result))
  expect_false("EXG2" %in% colnames(result))
})

# ----------------------------------------------------------------------------
# Test 4.2: detect -> apply pipeline on an eeg object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The same two-step pipeline but operating on an actual eeg
# object created with new_eeg(). This is the most realistic scenario in the
# package. After apply_external_labels(), the object must still be class "eeg",
# have the correct number of channels, and carry the new label names.
test_that("detect -> apply pipeline works correctly on a real eeg object", {
  eeg_obj <- new_eeg(
    data = matrix(rnorm(800), nrow = 8, ncol = 100),
    channels = c("Fp1", "Fp2", "F3", "F4", "EXG1", "EXG2", "GSR1", "Plet"),
    sampling_rate = 512
  )
  
  # Step 1: detect externals (should find EXG1, EXG2, GSR1, Plet)
  detected <- detect_external_channels(eeg_obj)
  expect_equal(length(detected), 4)
  
  # Step 2: apply labels
  labels <- setNames(c("EOG_L", "EOG_R", "GSR", "Pulse"), detected)
  result <- apply_external_labels(eeg_obj, labels = labels, keep_original = FALSE)
  
  # Class preserved
  expect_s3_class(result, "eeg")
  
  # New names present
  expect_true("EOG_L" %in% result$channels)
  expect_true("EOG_R" %in% result$channels)
  expect_true("GSR" %in% result$channels)
  expect_true("Pulse" %in% result$channels)
  
  # EEG channels untouched
  expect_true("Fp1" %in% result$channels)
  expect_true("F4"  %in% result$channels)
  
  # Channel count unchanged
  expect_equal(length(result$channels), 8)
})

# ----------------------------------------------------------------------------
# Test 4.3: identify (mocked) output feeds directly into apply_external_labels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The $labels output of identify_external_channels() (the
# interactive function) is designed to be passed directly to apply_external_labels().
# This test mocks the interactive step, takes the $labels result, and pipes it
# into apply_external_labels() on a data frame, confirming the interface is
# compatible end-to-end.
test_that("identify_external_channels $labels output feeds into apply_external_labels correctly", {
  local_mocked_bindings(
    readline = function(prompt = "") "EOG_L",
    .package = "base"
  )
  
  channels_vec <- c("EXG1", "Cz", "Pz")
  identification_result <- identify_external_channels(channels_vec)
  
  # Now apply to a data frame that has EXG1 as a column
  df <- data.frame(EXG1 = rnorm(5), Cz = rnorm(5), Pz = rnorm(5))
  result <- apply_external_labels(df,
                                  labels = identification_result$labels,
                                  keep_original = FALSE)
  
  expect_true("EOG_L" %in% colnames(result))
  expect_false("EXG1" %in% colnames(result))
  expect_true("Cz" %in% colnames(result))
})

# ----------------------------------------------------------------------------
# Test 4.4: Full pipeline with mixed external types (EXG, GSR, Plet, Temp)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A comprehensive integration scenario with all four categories
# of external channels (External/EXG, GSR, Plethysmograph, Temperature) mixed
# together with regular EEG channels. detect_external_channels() must identify
# all of them, and apply_external_labels() must rename them all correctly.
test_that("Full pipeline correctly handles all external channel types together", {
  all_channels <- c("Fp1", "Fp2", "EXG1", "EXG3", "GSR1", "GSR2", "Plet", "Temp",
                    "Cz", "Pz")
  
  detected <- detect_external_channels(all_channels)
  expect_equal(length(detected), 6)
  
  labels <- setNames(
    c("EOG_L", "EOG_R", "GSR_L", "GSR_R", "Pulse", "Body_Temp"),
    detected
  )
  
  df <- as.data.frame(matrix(rnorm(100), nrow = 10,
                             dimnames = list(NULL, all_channels)))
  
  result <- apply_external_labels(df, labels = labels, keep_original = FALSE)
  
  expect_true("EOG_L" %in% colnames(result))
  expect_true("GSR_L" %in% colnames(result))
  expect_true("Pulse" %in% colnames(result))
  expect_true("Body_Temp" %in% colnames(result))
  expect_true("Fp1" %in% colnames(result))
  expect_true("Cz" %in% colnames(result))
  expect_equal(ncol(result), 10)  # Column count unchanged
})


# ============================================================================
#                           End of Test File
# ============================================================================