# ============================================================================
#                       Test File for montage.R
# ============================================================================
#
# This test file provides comprehensive testing for the montage.R script,
# which defines the 'montage' S3 class (channel scalp positions) and the
# functions to build one from the built-in BioSemi electrode template and
# attach it to an 'eeg' object.
#
# Functions tested:
#   1. new_montage()     - Low-level montage constructor
#   2. create_montage()  - Builds a montage from the BioSemi 64 template
#   3. set_montage()     - Attaches a montage to an eeg object
#   4. print.montage()   - Print method for montage objects
#
# Author: Christos Dalamarinis
# Date: July 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#          TEST SUITE 1: new_montage() - Low-Level Constructor
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Basic object creation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_montage() builds a valid 'montage'
# object from arbitrary channel names and a positions data frame.
test_that("new_montage creates a valid montage object", {
  positions <- data.frame(
    channel = c("A", "B"),
    x = c(1, 2), y = c(3, 4), z = c(5, 6),
    theta = c(10, 20), phi = c(30, 40), radius = c(87.54, 87.54),
    stringsAsFactors = FALSE
  )

  m <- new_montage(channels = c("A", "B"), positions = positions,
                    coord_frame = "custom")

  expect_s3_class(m, "montage")
  expect_equal(m$channels, c("A", "B"))
  expect_equal(m$coord_frame, "custom")
  expect_equal(nrow(m$positions), 2)
})

# ----------------------------------------------------------------------------
# Test 1.2: Default coord_frame
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the default coord_frame is "biosemi64".
test_that("new_montage defaults coord_frame to biosemi64", {
  m <- new_montage(channels = "Cz",
                    positions = data.frame(channel = "Cz", x = 0, y = 0, z = 88,
                                            theta = 0, phi = 0, radius = 87.54))
  expect_equal(m$coord_frame, "biosemi64")
})

# ============================================================================
#     TEST SUITE 2: create_montage() - Build From BioSemi Template
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Default call returns all 64 standard electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that create_montage() with no arguments returns
# a montage with all 64 standard 10-20/10-10 BioSemi electrodes, and that
# external/status channels are excluded from the default set.
test_that("create_montage() with no arguments returns all 64 standard electrodes", {
  m <- create_montage()

  expect_s3_class(m, "montage")
  expect_equal(length(m$channels), 64)
  expect_equal(nrow(m$positions), 64)
  expect_true(all(c("Cz", "Fz", "Pz", "Oz", "Fp1", "Fp2") %in% m$channels))
  expect_false(any(grepl("^EXG|^GSR|^Status$", m$channels, ignore.case = TRUE)))
})

# ----------------------------------------------------------------------------
# Test 2.2: Subset selection, mixed standard/BioSemi naming, case-insensitivity
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that create_montage() accepts a subset of channels
# mixing standard 10-20 names, BioSemi A/B names, and non-standard casing -
# all should resolve to the same underlying electrode positions.
test_that("create_montage() resolves standard names, BioSemi names, and mixed case", {
  m <- create_montage(c("Cz", "cz", "CZ", "B16"))

  expect_equal(nrow(m$positions), 4)

  # Cz and B16 are the same physical electrode (see get_electrode_database())
  cz_rows <- m$positions[m$positions$channel %in% c("Cz", "cz", "CZ", "B16"), ]
  expect_equal(length(unique(cz_rows$x)), 1)
  expect_equal(length(unique(cz_rows$y)), 1)
  expect_equal(length(unique(cz_rows$z)), 1)
  expect_equal(unique(cz_rows$x), 0)
  expect_equal(unique(cz_rows$z), 88)
})

# ----------------------------------------------------------------------------
# Test 2.3: Unmatched channels are dropped with a warning
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that channels not present in the BioSemi template
# (typos, or external channels with no coordinates) are dropped with a
# warning rather than causing an error, as long as at least one channel
# resolves.
test_that("create_montage() warns and drops unmatched/external channels", {
  expect_warning(
    m <- create_montage(c("Cz", "NOT_A_CHANNEL", "EXG1")),
    "NOT_A_CHANNEL"
  )
  expect_equal(m$channels, "Cz")
})

# ----------------------------------------------------------------------------
# Test 2.4: All channels unmatched raises an error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() errors informatively when none
# of the requested channels can be resolved.
test_that("create_montage() errors when no channels resolve", {
  expect_error(
    suppressWarnings(create_montage(c("NOT_REAL_1", "NOT_REAL_2"))),
    "None of the requested channels"
  )
})

# ----------------------------------------------------------------------------
# Test 2.5: Unknown template errors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() rejects unsupported templates.
test_that("create_montage() errors on an unknown template", {
  expect_error(create_montage(template = "not_a_real_template"), "Unknown montage template")
})

# ============================================================================
#          TEST SUITE 3: set_montage() - Attach Montage to eeg Object
# ============================================================================

# Helper: build a small synthetic eeg object with real electrode names plus
# one external channel, so channel_types has both "eeg" and "external".
make_test_eeg <- function() {
  chans <- c("Cz", "Fz", "Pz", "Oz", "T7", "T8", "EXG1")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  new_eeg(data = mat, channels = chans, sampling_rate = 100)
}

# ----------------------------------------------------------------------------
# Test 3.1: Successful attachment
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() attaches a montage restricted to
# the eeg object's actual "eeg"-type channels (excluding the external
# channel), appends to preprocessing_history, and returns an 'eeg' object.
test_that("set_montage() attaches montage scoped to eeg-type channels only", {
  eeg <- make_test_eeg()
  eeg2 <- suppressWarnings(set_montage(eeg, create_montage()))

  expect_s3_class(eeg2, "eeg")
  expect_s3_class(eeg2$montage, "montage")

  # EXG1 is classified "external" and must not receive a scalp position
  expect_false("EXG1" %in% eeg2$montage$channels)
  expect_true(all(c("Cz", "Fz", "Pz", "Oz", "T7", "T8") %in% eeg2$montage$channels))

  # preprocessing_history grew by one entry
  expect_equal(length(eeg2$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
  expect_true(grepl("Montage attached", eeg2$preprocessing_history[[length(eeg2$preprocessing_history)]]))
})

# ----------------------------------------------------------------------------
# Test 3.2: Input validation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() requires an 'eeg' object and a
# 'montage' object, in the correct argument positions.
test_that("set_montage() validates input classes", {
  eeg <- make_test_eeg()
  m <- create_montage()

  expect_error(set_montage(list(), m), "class 'eeg'")
  expect_error(set_montage(eeg, list()), "class 'montage'")
})

# ----------------------------------------------------------------------------
# Test 3.3: Warnings for mismatched channel sets
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() warns (not errors) when the eeg
# object has EEG channels missing from the montage, or the montage has
# channels not present in the data.
test_that("set_montage() warns on partial channel overlap", {
  chans <- c("Cz", "Fz", "NotInMontage")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100)

  small_montage <- create_montage(c("Cz", "Fz", "Pz"))

  expect_warning(
    expect_warning(
      eeg2 <- set_montage(eeg, small_montage),
      "NotInMontage"
    ),
    "Pz"
  )
  expect_setequal(eeg2$montage$channels, c("Cz", "Fz"))
})

# ----------------------------------------------------------------------------
# Test 3.4: No overlap at all errors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() errors when there is no overlap
# at all between the montage and the eeg object's channels.
test_that("set_montage() errors when there is no channel overlap", {
  chans <- c("NotAnElectrode1", "NotAnElectrode2")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100)

  m <- create_montage(c("Cz", "Fz"))

  expect_error(suppressWarnings(set_montage(eeg, m)), "No overlap")
})

# ============================================================================
#              TEST SUITE 4: print.montage() - Print Method
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Print output contains key summary fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies print.montage() runs without error and reports
# the coordinate template and channel count.
test_that("print.montage displays a formatted summary", {
  m <- create_montage(c("Cz", "Fz", "Pz"))

  output <- capture.output(print(m))
  expect_true(any(grepl("Montage Object Summary", output)))
  expect_true(any(grepl("biosemi64", output)))
  expect_true(any(grepl("3", output)))
})

# ----------------------------------------------------------------------------
# Test 4.2: Print returns its argument invisibly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies print.montage() follows the standard R print
# method convention of invisibly returning x.
test_that("print.montage returns x invisibly", {
  m <- create_montage(c("Cz", "Fz"))
  ret <- withVisible(print(m))
  expect_false(ret$visible)
  expect_identical(ret$value, m)
})

# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - new_montage(): basic construction, default coord_frame
# - create_montage(): default 64-electrode set, mixed/case-insensitive naming,
#   unmatched-channel warnings, all-unmatched error, unknown template error
# - set_montage(): scoping to eeg-type channels, input validation, partial
#   overlap warnings, no-overlap error, preprocessing_history update
# - print.montage(): formatted output, invisible return
# ============================================================================
