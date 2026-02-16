# ============================================================================
#                    Test File for channel_info2.R
# ============================================================================
# 
# This test file provides comprehensive testing for the channel_info2.R script,
# which contains functions for BioSemi channel information and inspection.
#
# Functions tested:
#   1. get_electrode_database()           - Returns comprehensive electrode database
#   2. get_electrode_position()           - Looks up position info for an electrode
#   3. scan_biosemi_channels()            - Extracts channel info from BDF header
#   4. inspect_biosemi_file()             - Prints formatted channel summary
#   5. detect_electrode_naming_system()   - Detects electrode naming convention
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#       TEST SUITE 1: get_electrode_database() - Electrode Database
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Database structure and completeness
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that get_electrode_database() returns a properly
# structured database with all expected BioSemi 64-channel electrodes plus
# external channels
test_that("get_electrode_database returns complete database structure", {
  db <- get_electrode_database()
  
  # Test that database is a list
  expect_true(is.list(db))
  
  # Test that database is not empty
  expect_true(length(db) > 0)
  
  # Test that database contains expected standard electrodes
  standard_electrodes <- c("Fp1", "Fp2", "Cz", "Pz", "Oz", "Fz", "C3", "C4")
  for (electrode in standard_electrodes) {
    expect_true(electrode %in% names(db), 
                info = paste("Standard electrode", electrode, "should be in database"))
  }
  
  # Test that database contains BioSemi naming
  biosemi_electrodes <- c("A1", "A16", "A32", "B1", "B16", "B32")
  for (electrode in biosemi_electrodes) {
    expect_true(electrode %in% names(db),
                info = paste("BioSemi electrode", electrode, "should be in database"))
  }
  
  # Test that database contains external channels
  external_channels <- c("EXG1", "EXG2", "Status", "GSR1", "GSR2")
  for (channel in external_channels) {
    expect_true(channel %in% names(db),
                info = paste("External channel", channel, "should be in database"))
  }
})

# ----------------------------------------------------------------------------
# Test 1.2: Electrode entry structure
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that each electrode entry has all required fields
# with correct types and valid values
test_that("electrode entries have correct structure and fields", {
  db <- get_electrode_database()
  
  # Test a standard electrode (Cz)
  cz_entry <- db$Cz
  
  # Test required fields exist
  expect_true("position_name" %in% names(cz_entry))
  expect_true("position_type" %in% names(cz_entry))
  expect_true("region" %in% names(cz_entry))
  expect_true("standard_system" %in% names(cz_entry))
  expect_true("biosemi_name" %in% names(cz_entry))
  expect_true("standard_name" %in% names(cz_entry))
  expect_true("cartesian_coords" %in% names(cz_entry))
  expect_true("spherical_coords" %in% names(cz_entry))
  
  # Test field types
  expect_type(cz_entry$position_name, "character")
  expect_type(cz_entry$position_type, "character")
  expect_type(cz_entry$region, "character")
  expect_type(cz_entry$standard_system, "character")
  
  # Test coordinate structures
  expect_true(is.list(cz_entry$cartesian_coords))
  expect_true(is.list(cz_entry$spherical_coords))
  
  # Test Cartesian coordinates
  expect_true("x" %in% names(cz_entry$cartesian_coords))
  expect_true("y" %in% names(cz_entry$cartesian_coords))
  expect_true("z" %in% names(cz_entry$cartesian_coords))
  
  # Test Spherical coordinates
  expect_true("theta" %in% names(cz_entry$spherical_coords))
  expect_true("phi" %in% names(cz_entry$spherical_coords))
  expect_true("radius" %in% names(cz_entry$spherical_coords))
  
  # Test coordinate values are numeric
  expect_type(cz_entry$cartesian_coords$x, "double")
  expect_type(cz_entry$spherical_coords$theta, "double")
})

# ----------------------------------------------------------------------------
# Test 1.3: Dual naming system (standard and BioSemi)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that electrodes can be accessed by both standard
# 10-20/10-10 names AND BioSemi A/B names, and that they point to the same data
test_that("dual naming system works correctly", {
  db <- get_electrode_database()
  
  # Test that Cz and B16 are the same electrode
  cz_entry <- db$Cz
  b16_entry <- db$B16
  
  expect_equal(cz_entry, b16_entry)
  
  # Test standard name matches
  expect_equal(cz_entry$standard_name, "Cz")
  expect_equal(b16_entry$standard_name, "Cz")
  
  # Test BioSemi name matches
  expect_equal(cz_entry$biosemi_name, "B16")
  expect_equal(b16_entry$biosemi_name, "B16")
  
  # Test another pair: Fp1 and A1
  fp1_entry <- db$Fp1
  a1_entry <- db$A1
  
  expect_equal(fp1_entry, a1_entry)
  expect_equal(fp1_entry$standard_name, "Fp1")
  expect_equal(a1_entry$biosemi_name, "A1")
})

# ----------------------------------------------------------------------------
# Test 1.4: Case-insensitive lookup
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that electrode lookup is case-insensitive,
# allowing users to use "cz", "CZ", or "Cz" interchangeably
test_that("case-insensitive electrode lookup works", {
  db <- get_electrode_database()
  
  # Test lowercase
  expect_true("cz" %in% names(db))
  expect_true("fp1" %in% names(db))
  expect_true("a1" %in% names(db))
  
  # Test that lowercase and uppercase return the same data
  cz_lower <- db$cz
  cz_upper <- db$Cz
  
  expect_equal(cz_lower, cz_upper)
})

# ----------------------------------------------------------------------------
# Test 1.5: External channel entries
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that external channels (EOG, ECG, etc.) are
# properly included with appropriate metadata
test_that("external channels are properly defined", {
  db <- get_electrode_database()
  
  # Test EXG1 (external electrode)
  exg1_entry <- db$EXG1
  
  expect_equal(exg1_entry$region, "External")
  expect_equal(exg1_entry$standard_system, "BioSemi")
  expect_type(exg1_entry$position_name, "character")
  
  # Test that coordinates are NA for external channels
  expect_true(is.na(exg1_entry$cartesian_coords$x))
  expect_true(is.na(exg1_entry$spherical_coords$theta))
  
  # Test Status channel
  status_entry <- db$Status
  expect_equal(status_entry$region, "Metadata")
  expect_equal(status_entry$position_type, "Status")
})


# ============================================================================
#       TEST SUITE 2: get_electrode_position() - Position Lookup
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Valid electrode lookup
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that get_electrode_position() correctly returns
# position information for known electrodes
test_that("get_electrode_position returns correct info for valid electrodes", {
  # Test standard electrode
  cz_info <- get_electrode_position("Cz")
  
  expect_true(is.list(cz_info))
  expect_equal(cz_info$position_type, "Central")
  expect_equal(cz_info$region, "Central")
  expect_true(cz_info$standard_system %in% c("10-20", "10-10"))
  
  # Test frontal electrode
  fp1_info <- get_electrode_position("Fp1")
  
  expect_equal(fp1_info$region, "Frontal")
  expect_true(grepl("Frontal", fp1_info$position_type, ignore.case = TRUE))
})

# ----------------------------------------------------------------------------
# Test 2.2: BioSemi naming lookup
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that BioSemi A/B naming can be used to look up
# electrode positions
test_that("get_electrode_position works with BioSemi naming", {
  # Test A-series
  a1_info <- get_electrode_position("A1")
  
  expect_true(is.list(a1_info))
  expect_equal(a1_info$standard_name, "Fp1")
  expect_equal(a1_info$biosemi_name, "A1")
  
  # Test B-series
  b16_info <- get_electrode_position("B16")
  
  expect_equal(b16_info$standard_name, "Cz")
  expect_equal(b16_info$biosemi_name, "B16")
})

# ----------------------------------------------------------------------------
# Test 2.3: Unknown electrode handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that unknown electrodes return appropriate
# "Unknown" status rather than throwing errors
test_that("get_electrode_position handles unknown electrodes gracefully", {
  unknown_info <- get_electrode_position("XYZ999")
  
  expect_true(is.list(unknown_info))
  expect_equal(unknown_info$position_type, "Unknown")
  expect_equal(unknown_info$region, "Unknown")
  expect_equal(unknown_info$standard_system, "Unknown")
  expect_true(grepl("Unknown electrode", unknown_info$position_name))
})

# ----------------------------------------------------------------------------
# Test 2.4: External channel lookup
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that external channels return correct metadata
test_that("get_electrode_position works for external channels", {
  exg1_info <- get_electrode_position("EXG1")
  
  expect_equal(exg1_info$region, "External")
  expect_equal(exg1_info$standard_system, "BioSemi")
  
  status_info <- get_electrode_position("Status")
  
  expect_equal(status_info$region, "Metadata")
  expect_equal(status_info$position_type, "Status")
})


# ============================================================================
#       TEST SUITE 3: scan_biosemi_channels() - BDF Channel Scanning
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: File validation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that scan_biosemi_channels() properly validates
# file existence and provides helpful error messages
test_that("scan_biosemi_channels validates file existence", {
  # Test with non-existent file
  expect_error(
    scan_biosemi_channels("/path/to/nonexistent/file.bdf"),
    "File not found"
  )
})

# ----------------------------------------------------------------------------
# Test 3.2: File extension warning
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that non-.bdf files trigger a warning but still
# attempt to process
test_that("scan_biosemi_channels warns for non-.bdf extensions", {
  # Create a temporary file with wrong extension
  temp_file <- tempfile(fileext = ".txt")
  file.create(temp_file)
  
  # Should warn but not error (will error on actual read, which is expected)
  expect_warning(
    expect_error(scan_biosemi_channels(temp_file)),
    "does not have .bdf extension"
  )
  
  # Cleanup
  unlink(temp_file)
})

# ----------------------------------------------------------------------------
# Test 3.3: Return structure (when valid BDF available)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that scan_biosemi_channels() returns a properly
# structured data frame with expected columns
# NOTE: This test requires a valid BDF file to be available
test_that("scan_biosemi_channels returns correct data frame structure", {
  skip_if_not(exists("test_bdf_file") && file.exists(test_bdf_file),
              message = "Valid test BDF file not available")
  
  channel_info <- scan_biosemi_channels(test_bdf_file)
  
  # Test it returns a data frame
  expect_true(is.data.frame(channel_info))
  
  # Test required columns exist
  expect_true("channel_number" %in% names(channel_info))
  expect_true("channel_name" %in% names(channel_info))
  expect_true("position_type" %in% names(channel_info))
  expect_true("region" %in% names(channel_info))
  expect_true("standard_system" %in% names(channel_info))
  
  # Test channel numbers are sequential
  expect_equal(channel_info$channel_number, 1:nrow(channel_info))
  
  # Test column types
  expect_type(channel_info$channel_number, "integer")
  expect_type(channel_info$channel_name, "character")
  expect_type(channel_info$position_type, "character")
})


# ============================================================================
#     TEST SUITE 4: inspect_biosemi_file() - Formatted File Inspection
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: File validation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that inspect_biosemi_file() properly validates
# file existence before attempting inspection
test_that("inspect_biosemi_file validates file existence", {
  expect_error(
    inspect_biosemi_file("/path/to/nonexistent/file.bdf"),
    "File not found"
  )
})

# ----------------------------------------------------------------------------
# Test 4.2: Console output and invisible return
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that inspect_biosemi_file() prints to console
# and invisibly returns the channel info data frame
# NOTE: This test requires a valid BDF file
test_that("inspect_biosemi_file prints and returns correctly", {
  skip_if_not(exists("test_bdf_file") && file.exists(test_bdf_file),
              message = "Valid test BDF file not available")
  
  # Capture output
  output <- capture.output(result <- inspect_biosemi_file(test_bdf_file))
  
  # Test that output was generated
  expect_true(length(output) > 0)
  
  # Test that output contains expected headers
  output_text <- paste(output, collapse = "\n")
  expect_true(grepl("BioSemi BDF File Channel Inspection", output_text))
  expect_true(grepl("FILE INFORMATION", output_text))
  expect_true(grepl("CHANNEL SUMMARY", output_text))
  
  # Test that result is invisibly returned as data frame
  expect_true(is.data.frame(result))
  expect_true("channel_name" %in% names(result))
})


# ============================================================================
# TEST SUITE 5: detect_electrode_naming_system() - Naming Convention Detection
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Input type detection - file path
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that detect_electrode_naming_system() can detect
# and process file path inputs
test_that("detect_electrode_naming_system handles file path input", {
  skip_if_not(exists("test_bdf_file") && file.exists(test_bdf_file),
              message = "Valid test BDF file not available")
  
  result <- detect_electrode_naming_system(test_bdf_file)
  
  # Test return structure
  expect_true(is.list(result))
  expect_true("naming_system" %in% names(result))
  expect_true("confidence" %in% names(result))
  expect_true("channels_detected" %in% names(result))
  expect_true("summary" %in% names(result))
})

# ----------------------------------------------------------------------------
# Test 5.2: Input type detection - character vector
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that detect_electrode_naming_system() correctly
# processes character vectors of channel names
test_that("detect_electrode_naming_system handles character vector input", {
  # Test with standard 10-20 naming
  standard_channels <- c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", 
                         "O1", "O2", "Fz", "Cz", "Pz")
  
  result <- detect_electrode_naming_system(standard_channels)
  
  expect_equal(result$naming_system, "10-20/10-10")
  expect_true(result$confidence >= 0.9)
  expect_true(is.data.frame(result$channels_detected))
  expect_equal(nrow(result$channels_detected), length(standard_channels))
})

# ----------------------------------------------------------------------------
# Test 5.3: Standard 10-20/10-10 detection
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies accurate detection of standard international
# electrode naming convention
test_that("detect_electrode_naming_system correctly identifies standard naming", {
  standard_channels <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                         "T7", "C3", "Cz", "C4", "T8",
                         "P7", "P3", "Pz", "P4", "P8",
                         "O1", "Oz", "O2")
  
  result <- detect_electrode_naming_system(standard_channels)
  
  expect_equal(result$naming_system, "10-20/10-10")
  expect_equal(result$confidence, 1.0)
  
  # Check channels_detected data frame
  expect_true(all(result$channels_detected$system_category == "Standard_10-20/10-10"))
})

# ----------------------------------------------------------------------------
# Test 5.4: BioSemi A/B naming detection
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies accurate detection of BioSemi internal naming
# convention (A1-A32, B1-B32)
test_that("detect_electrode_naming_system correctly identifies BioSemi naming", {
  biosemi_channels <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
                        "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")
  
  result <- detect_electrode_naming_system(biosemi_channels)
  
  expect_equal(result$naming_system, "BioSemi_AB")
  expect_equal(result$confidence, 1.0)
  
  # Check that standard equivalents are provided
  expect_true(all(!is.na(result$channels_detected$standard_equivalent)))
  
  # Verify A1 maps to Fp1
  a1_row <- result$channels_detected[result$channels_detected$channel_name == "A1", ]
  expect_equal(a1_row$standard_equivalent, "Fp1")
})

# ----------------------------------------------------------------------------
# Test 5.5: Mixed naming detection
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies detection of mixed naming conventions (both
# standard and BioSemi naming present)
test_that("detect_electrode_naming_system correctly identifies mixed naming", {
  mixed_channels <- c("Fp1", "Fp2", "A3", "A4", "B1", "B2", "Cz", "Pz")
  
  result <- detect_electrode_naming_system(mixed_channels)
  
  expect_equal(result$naming_system, "Mixed")
  expect_true(result$confidence > 0 && result$confidence <= 1)
})

# ----------------------------------------------------------------------------
# Test 5.6: Unknown naming detection
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that unrecognized electrode names result in
# "Unknown" system classification
test_that("detect_electrode_naming_system handles unknown naming", {
  unknown_channels <- c("Ch1", "Ch2", "Ch3", "Ch4", "Ch5")
  
  result <- detect_electrode_naming_system(unknown_channels)
  
  expect_equal(result$naming_system, "Unknown")
  expect_equal(result$confidence, 0)
})

# ----------------------------------------------------------------------------
# Test 5.7: External channel filtering
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that external channels (EXG, Status, etc.) are
# properly excluded from naming system detection
test_that("detect_electrode_naming_system excludes external channels from detection", {
  channels_with_external <- c("Fp1", "Fp2", "Cz", "Pz", "EXG1", "EXG2", "Status")
  
  result <- detect_electrode_naming_system(channels_with_external)
  
  # Should detect as standard naming, ignoring external channels
  expect_equal(result$naming_system, "10-20/10-10")
  
  # Confidence should be based on non-external channels only
  # 4 standard channels out of 4 non-external = 1.0
  expect_equal(result$confidence, 1.0)
})

# ----------------------------------------------------------------------------
# Test 5.8: Input from data frame (column names)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that detect_electrode_naming_system() can extract
# channel names from data frame column names
test_that("detect_electrode_naming_system handles data frame input", {
  # Create dummy data frame with channel names as columns
  test_df <- data.frame(
    Fp1 = rnorm(100),
    Fp2 = rnorm(100),
    Cz = rnorm(100),
    Pz = rnorm(100)
  )
  
  result <- detect_electrode_naming_system(test_df)
  
  expect_equal(result$naming_system, "10-20/10-10")
  expect_equal(nrow(result$channels_detected), 4)
})

# ----------------------------------------------------------------------------
# Test 5.9: Input from matrix (column names)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that detect_electrode_naming_system() can extract
# channel names from matrix column names
test_that("detect_electrode_naming_system handles matrix input", {
  # Create matrix with channel names
  test_matrix <- matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(test_matrix) <- c("A1", "A2", "B1", "B2")
  
  result <- detect_electrode_naming_system(test_matrix)
  
  expect_equal(result$naming_system, "BioSemi_AB")
})

# ----------------------------------------------------------------------------
# Test 5.10: Input from list (various channel locations)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that detect_electrode_naming_system() can extract
# channel names from list structures in common locations
test_that("detect_electrode_naming_system handles list input with $channels", {
  # Test list with $channels component
  test_list <- list(
    channels = c("Fp1", "Fp2", "Cz", "Pz"),
    data = matrix(rnorm(400), nrow = 100, ncol = 4)
  )
  
  result <- detect_electrode_naming_system(test_list)
  
  expect_equal(result$naming_system, "10-20/10-10")
  expect_equal(nrow(result$channels_detected), 4)
})

# ----------------------------------------------------------------------------
# Test 5.11: Summary output format
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the summary component is properly formatted
# as a character vector with expected content
test_that("detect_electrode_naming_system generates proper summary", {
  channels <- c("Fp1", "Fp2", "Cz", "Pz")
  
  result <- detect_electrode_naming_system(channels)
  
  expect_type(result$summary, "character")
  expect_true(length(result$summary) > 0)
  
  # Check summary contains key information
  summary_text <- paste(result$summary, collapse = "\n")
  expect_true(grepl("PRIMARY NAMING SYSTEM", summary_text))
})

# ----------------------------------------------------------------------------
# Test 5.12: Confidence calculation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that confidence scores are calculated correctly
# based on the percentage of recognized electrodes
test_that("detect_electrode_naming_system calculates confidence correctly", {
  # All recognized
  all_recognized <- c("Fp1", "Fp2", "Cz", "Pz")
  result_all <- detect_electrode_naming_system(all_recognized)
  expect_equal(result_all$confidence, 1.0)
  
  # Half recognized, half unknown
  half_recognized <- c("Fp1", "Fp2", "Unknown1", "Unknown2")
  result_half <- detect_electrode_naming_system(half_recognized)
  expect_equal(result_half$confidence, 0.5)
  
  # None recognized
  none_recognized <- c("Unknown1", "Unknown2", "Unknown3")
  result_none <- detect_electrode_naming_system(none_recognized)
  expect_equal(result_none$confidence, 0)
})

# ----------------------------------------------------------------------------
# Test 5.13: Invalid input handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that invalid input types generate appropriate
# error messages
test_that("detect_electrode_naming_system handles invalid input types", {
  # Test with numeric input
  expect_error(
    detect_electrode_naming_system(123),
    "Invalid input type"
  )
  
  # Test with list without channel information
  empty_list <- list(foo = "bar", baz = 123)
  expect_error(
    detect_electrode_naming_system(empty_list),
    "could not find channel names"
  )
  
  # Test with data frame without column names
  df_no_names <- data.frame(matrix(rnorm(100), nrow = 10))
  names(df_no_names) <- NULL
  expect_error(
    detect_electrode_naming_system(df_no_names),
    "has no column names"
  )
})

# ----------------------------------------------------------------------------
# Test 5.14: Whitespace handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that channel names with leading/trailing whitespace
# are properly trimmed and processed
test_that("detect_electrode_naming_system handles whitespace correctly", {
  channels_with_spaces <- c(" Fp1 ", "  Fp2", "Cz  ", " Pz ")
  
  result <- detect_electrode_naming_system(channels_with_spaces)
  
  expect_equal(result$naming_system, "10-20/10-10")
  
  # Check that trimmed names are in channels_detected
  expect_true(all(c("Fp1", "Fp2", "Cz", "Pz") %in% 
                    result$channels_detected$channel_name))
})


# ============================================================================
#                         TEST SUITE 6: Integration Tests
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: Database and position lookup integration
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that get_electrode_position() correctly retrieves
# data from get_electrode_database()
test_that("get_electrode_position integrates correctly with database", {
  db <- get_electrode_database()
  
  # Test multiple electrodes
  test_electrodes <- c("Cz", "Fp1", "O2", "A1", "B32", "EXG1")
  
  for (electrode in test_electrodes) {
    position_info <- get_electrode_position(electrode)
    db_entry <- db[[electrode]]
    
    # Verify they contain the same information
    expect_equal(position_info$position_type, db_entry$position_type)
    expect_equal(position_info$region, db_entry$region)
    expect_equal(position_info$standard_system, db_entry$standard_system)
  }
})

# ----------------------------------------------------------------------------
# Test 6.2: Naming detection and position lookup integration
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that channels detected by detect_electrode_naming_system()
# can be successfully looked up with get_electrode_position()
test_that("naming detection integrates with position lookup", {
  channels <- c("Fp1", "A2", "Cz", "EXG1")
  
  result <- detect_electrode_naming_system(channels)
  
  # For each detected channel, verify position lookup works
  for (i in 1:nrow(result$channels_detected)) {
    channel_name <- result$channels_detected$channel_name[i]
    position_info <- get_electrode_position(channel_name)
    
    expect_true(!is.null(position_info))
    expect_type(position_info$position_type, "character")
  }
})


# ============================================================================
#               TEST SUITE 7: Edge Cases and Error Handling
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: Empty input handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies graceful handling of empty inputs
test_that("functions handle empty inputs appropriately", {
  # Empty character vector
  expect_error(
    detect_electrode_naming_system(character(0)),
    "No channel names found"
  )
  
  # Empty data frame
  empty_df <- data.frame()
  expect_error(
    detect_electrode_naming_system(empty_df),
    "has no column names"
  )
})

# ----------------------------------------------------------------------------
# Test 7.2: Special characters in channel names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies handling of channel names with special characters
test_that("functions handle special characters in channel names", {
  special_channels <- c("EOG_L", "EOG_R", "ECG-1", "Ref-Cz")
  
  # Should not crash, though may classify as unknown
  result <- detect_electrode_naming_system(special_channels)
  
  expect_true(is.list(result))
  expect_type(result$naming_system, "character")
})

# ----------------------------------------------------------------------------
# Test 7.3: Very large channel sets
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that functions can handle large numbers of channels
# without performance issues or errors
test_that("functions handle large channel sets", {
  # Create 128 channel names (standard high-density EEG setup)
  large_channel_set <- paste0("Ch", 1:128)
  
  result <- detect_electrode_naming_system(large_channel_set)
  
  expect_equal(nrow(result$channels_detected), 128)
  expect_equal(result$naming_system, "Unknown")  # Custom numbering
})

# ----------------------------------------------------------------------------
# Test 7.4: Duplicate channel names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies handling of duplicate channel names in input
test_that("functions handle duplicate channel names", {
  duplicate_channels <- c("Fp1", "Fp1", "Cz", "Cz", "Pz")
  
  # Should process without error (duplicates are kept)
  result <- detect_electrode_naming_system(duplicate_channels)
  
  expect_equal(nrow(result$channels_detected), 5)
})

# ----------------------------------------------------------------------------
# Test 7.5: NULL and NA handling
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies appropriate error handling for NULL and NA inputs
test_that("functions handle NULL and NA inputs appropriately", {
  # NULL input
  expect_error(
    get_electrode_position(NULL)
  )
  
  # NA input
  na_result <- get_electrode_position(NA_character_)
  expect_equal(na_result$position_type, "Unknown")
})


# ============================================================================
#                         SUMMARY OF TEST COVERAGE
# ============================================================================
# 
# This test file provides comprehensive coverage for channel_info2.R:
#
# 1. get_electrode_database():
#    - Database structure and completeness
#    - Electrode entry field validation
#    - Dual naming system (standard + BioSemi)
#    - Case-insensitive lookup
#    - External channel entries
#
# 2. get_electrode_position():
#    - Valid electrode lookup
#    - BioSemi naming support
#    - Unknown electrode handling
#    - External channel lookup
#
# 3. scan_biosemi_channels():
#    - File validation
#    - File extension warnings
#    - Return structure validation
#    (Note: Full testing requires valid BDF files)
#
# 4. inspect_biosemi_file():
#    - File validation
#    - Console output and invisible return
#    (Note: Full testing requires valid BDF files)
#
# 5. detect_electrode_naming_system():
#    - Multiple input type detection (file, vector, df, matrix, list)
#    - Standard 10-20/10-10 naming detection
#    - BioSemi A/B naming detection
#    - Mixed naming detection
#    - Unknown naming handling
#    - External channel filtering
#    - Confidence calculation
#    - Summary generation
#    - Invalid input handling
#    - Whitespace trimming
#
# 6. Integration Tests:
#    - Database and position lookup integration
#    - Naming detection and position lookup integration
#
# 7. Edge Cases:
#    - Empty inputs
#    - Special characters
#    - Large channel sets
#    - Duplicate names
#    - NULL and NA handling
#
# Total: 50+ individual test cases covering all major functions and use cases
#
# ============================================================================