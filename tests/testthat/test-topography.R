# ============================================================================
#                       Test File for topography.R
# ============================================================================
#
# This test file provides comprehensive testing for the topography.R script,
# which renders a 2D interpolated scalp heatmap ("topoplot") of one value per
# channel using a montage's electrode positions.
#
# Functions tested:
#   1. plot_topography() - Draws the interpolated scalp topography
#
# Author: Christos Dalamarinis
# Date: July 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# Helper: build a small synthetic eeg object with a montage already attached,
# and a corresponding named vector of per-channel values (mimicking the
# output of eeg_band_power()).
make_topo_fixture <- function() {
  chans <- c("Cz", "Fz", "Pz", "Oz", "T7", "T8")
  mat <- matrix(rnorm(length(chans) * 200), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100)
  eeg <- set_montage(eeg, create_montage(chans))
  values <- setNames(runif(length(chans)), chans)
  list(eeg = eeg, values = values)
}

# Runs `expr` against a throwaway PNG device so no plot window is displayed
# during the test run, and always cleans the device up afterwards.
with_null_device <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  force(expr)
}

# ============================================================================
#           TEST SUITE 1: Input Validation and Error Handling
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Requires an 'eeg' object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() rejects a non-'eeg' first
# argument.
test_that("plot_topography errors when eeg_obj is not an 'eeg' object", {
  fx <- make_topo_fixture()
  expect_error(plot_topography(list(), fx$values), "class 'eeg'")
})

# ----------------------------------------------------------------------------
# Test 1.2: Requires a montage
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() errors clearly when neither
# eeg_obj$montage nor an explicit montage argument is available.
test_that("plot_topography errors when no montage is available", {
  chans <- c("Cz", "Fz", "Pz")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  eeg_no_montage <- new_eeg(data = mat, channels = chans, sampling_rate = 100)
  values <- setNames(runif(length(chans)), chans)

  expect_error(plot_topography(eeg_no_montage, values), "No montage available")
})

# ----------------------------------------------------------------------------
# Test 1.3: Requires named values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() errors when 'values' has no
# names to match against montage channels.
test_that("plot_topography errors when values is not named", {
  fx <- make_topo_fixture()
  expect_error(plot_topography(fx$eeg, unname(fx$values)), "named numeric vector")
})

# ----------------------------------------------------------------------------
# Test 1.4: Warns on channel mismatches between values and montage
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() warns (not errors) when
# 'values' includes channels absent from the montage, or the montage has
# channels with no supplied value - as long as enough channels remain.
test_that("plot_topography warns on partial value/montage overlap", {
  fx <- make_topo_fixture()
  values <- fx$values
  names(values)[1] <- "NotAMontageChannel"

  with_null_device({
    expect_warning(
      expect_warning(
        res <- plot_topography(fx$eeg, values),
        "no matching montage channel"
      ),
      "no value supplied"
    )
  })
  expect_equal(nrow(res$channel_positions), length(fx$values) - 1)
})

# ----------------------------------------------------------------------------
# Test 1.5: Errors when fewer than 3 channels can be plotted
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() refuses to interpolate a
# surface from fewer than 3 positioned channels.
test_that("plot_topography errors with fewer than 3 usable channels", {
  fx <- make_topo_fixture()
  values <- fx$values[1:2]

  expect_error(
    suppressWarnings(plot_topography(fx$eeg, values)),
    "At least 3 channels"
  )
})

# ============================================================================
#                  TEST SUITE 2: Successful Rendering
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Returns the expected invisible structure
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a successful call returns (invisibly) a list with
# grid_x, grid_y, grid_z, and channel_positions, using the montage attached
# to the eeg object.
test_that("plot_topography returns grid and channel position data", {
  fx <- make_topo_fixture()

  res <- with_null_device(plot_topography(fx$eeg, fx$values))

  expect_type(res, "list")
  expect_named(res, c("grid_x", "grid_y", "grid_z", "channel_positions"))
  expect_type(res$grid_x, "double")
  expect_type(res$grid_y, "double")
  expect_true(is.matrix(res$grid_z))
  expect_equal(dim(res$grid_z), c(length(res$grid_x), length(res$grid_y)))
  expect_s3_class(res$channel_positions, "data.frame")
  expect_equal(nrow(res$channel_positions), length(fx$values))
})

# ----------------------------------------------------------------------------
# Test 2.2: interpolate_res controls grid resolution
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the interpolate_res parameter controls the size
# of the returned interpolation grid.
test_that("plot_topography respects interpolate_res", {
  fx <- make_topo_fixture()

  res <- with_null_device(
    plot_topography(fx$eeg, fx$values, interpolate_res = 25)
  )

  expect_equal(length(res$grid_x), 25)
  expect_equal(length(res$grid_y), 25)
})

# ----------------------------------------------------------------------------
# Test 2.3: An explicit montage argument overrides eeg_obj$montage
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() accepts a montage passed
# directly, even when the eeg object has no montage attached at all.
test_that("plot_topography accepts an explicit montage argument", {
  chans <- c("Cz", "Fz", "Pz", "T7")
  mat <- matrix(rnorm(length(chans) * 100), nrow = length(chans))
  eeg_no_montage <- new_eeg(data = mat, channels = chans, sampling_rate = 100)
  values <- setNames(runif(length(chans)), chans)

  res <- with_null_device(
    plot_topography(eeg_no_montage, values, montage = create_montage(chans))
  )

  expect_equal(nrow(res$channel_positions), length(chans))
})

# ----------------------------------------------------------------------------
# Test 2.4: Collinear channels are rejected with a clear error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() gives a clear error (rather
# than a cryptic akima::interp failure) when the selected channels have no
# spread in one spatial dimension, e.g. a purely midline selection.
test_that("plot_topography errors clearly on collinear channels", {
  chans <- c("Cz", "Fz", "Pz", "Oz")
  mat <- matrix(rnorm(length(chans) * 100), nrow = length(chans))
  eeg_midline <- new_eeg(data = mat, channels = chans, sampling_rate = 100)
  eeg_midline <- set_montage(eeg_midline, create_montage(chans))
  values <- setNames(runif(length(chans)), chans)

  expect_error(
    with_null_device(plot_topography(eeg_midline, values)),
    "collinear"
  )
})

# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - Input validation: non-'eeg' object, missing montage, unnamed values
# - Warnings: partial overlap between supplied values and montage channels
# - Error: fewer than 3 usable channels for interpolation
# - Successful rendering: return structure, interpolate_res, explicit
#   montage argument overriding/substituting for eeg_obj$montage
# ============================================================================
