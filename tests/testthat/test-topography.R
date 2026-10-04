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

# Helper: a smooth, non-random map on the full 64-channel layout - a broad
# bump that is 1 at Cz and falls towards 0 at the edge of the head - with an
# eeg object that carries the montage. Unlike make_topo_fixture() it has no
# random values, so tests can rely on its exact shape, and it has enough
# electrodes for a real interpolation: akima::interp() returns a flat
# all-zero surface when given fewer than 10 points, so the 6-channel fixture
# above cannot show anything about colours or extrapolation.
make_blob_fixture <- function() {
  montage   <- create_montage()
  positions <- montage$positions
  xyz       <- as.matrix(positions[, c("x", "y", "z")])
  xyz       <- xyz / sqrt(rowSums(xyz^2))
  cosang    <- as.vector(xyz %*% xyz[positions$channel == "Cz", ])
  angle     <- acos(pmin(pmax(cosang, -1), 1))
  values    <- setNames(exp(-(angle / 0.6)^2), positions$channel)
  eeg <- new_eeg(data = matrix(0, nrow = nrow(positions), ncol = 4),
                 channels = positions$channel, sampling_rate = 100,
                 montage = montage)
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

# Draws `expr` on a throwaway PNG device and returns the bytes of the image,
# so that two drawings can be compared.
png_bytes <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit(unlink(tmp), add = TRUE)
  tryCatch(force(expr), finally = grDevices::dev.off())
  readBin(tmp, "raw", file.size(tmp))
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
# grid_x, grid_y, grid_z, zlim, and channel_positions, using the montage
# attached to the eeg object.
test_that("plot_topography returns grid and channel position data", {
  fx <- make_topo_fixture()

  res <- with_null_device(plot_topography(fx$eeg, fx$values))

  expect_type(res, "list")
  expect_named(
    res,
    c("grid_x", "grid_y", "grid_z", "zlim", "channel_positions")
  )
  expect_type(res$grid_x, "double")
  expect_type(res$grid_y, "double")
  expect_type(res$zlim, "double")
  expect_length(res$zlim, 2)
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
#            TEST SUITE 3: eeg_obj$bads Exclusion
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Channels marked bad are excluded with a warning
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() reads eeg_obj$bads (the shared
# exclude list, see new_eeg()) and drops bad channels from the topography
# with a warning, instead of feeding their raw value into the interpolation.
test_that("plot_topography excludes eeg_obj$bads with a warning", {
  chans <- c("Cz", "Fz", "Pz", "Oz", "T7", "T8")
  mat <- matrix(rnorm(length(chans) * 200), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100, bads = "T7")
  eeg <- set_montage(eeg, create_montage(chans))
  values <- setNames(runif(length(chans)), chans)

  with_null_device({
    expect_warning(res <- plot_topography(eeg, values), "marked bad.*T7")
  })

  expect_false("T7" %in% res$channel_positions$channel)
  expect_equal(nrow(res$channel_positions), length(chans) - 1)
})

# ----------------------------------------------------------------------------
# Test 3.2: No bads means no exclusion and no warning
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies plot_topography() does not warn or drop any
# channel when eeg_obj$bads is empty (the default).
test_that("plot_topography does not warn when there are no bad channels", {
  fx <- make_topo_fixture()

  expect_no_warning(res <- with_null_device(plot_topography(fx$eeg, fx$values)))
  expect_equal(nrow(res$channel_positions), length(fx$values))
})

# ----------------------------------------------------------------------------
# Test 3.3: Excluding bads can drop below the 3-channel minimum
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that if excluding bad channels leaves fewer than
# 3 usable channels, plot_topography() still errors clearly (after warning
# about the exclusion).
test_that("plot_topography errors when bads exclusion drops below 3 channels", {
  chans <- c("Cz", "Fz", "Pz")
  mat <- matrix(rnorm(length(chans) * 100), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100, bads = "Pz")
  eeg <- set_montage(eeg, create_montage(chans))
  values <- setNames(runif(length(chans)), chans)

  expect_error(
    suppressWarnings(with_null_device(plot_topography(eeg, values))),
    "At least 3 channels"
  )
})

# ============================================================================
#          TEST SUITE 4: zero_centered and extrapolate Arguments
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: The defaults keep the behaviour from before the new arguments
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that leaving zero_centered and extrapolate out
# gives the same result as passing their defaults (FALSE and TRUE), and that
# the colour limits (zlim) are then simply the range of the plotted values, as
# before the arguments existed.
test_that("zero_centered and extrapolate default to the old behaviour", {
  fx <- make_blob_fixture()

  by_default  <- with_null_device(plot_topography(fx$eeg, fx$values))
  spelled_out <- with_null_device(
    plot_topography(fx$eeg, fx$values,
                    zero_centered = FALSE, extrapolate = TRUE)
  )

  expect_equal(by_default, spelled_out)
  expect_equal(by_default$zlim, range(by_default$grid_z, na.rm = TRUE))
})

# ----------------------------------------------------------------------------
# Test 4.2: zero_centered makes the colour limits symmetric around zero
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies zero_centered = TRUE leaves the plotted values
# alone and only changes the colour limits: they become symmetric around zero
# and just wide enough to hold the strongest plotted value. The bump's values
# are all positive (0 to 1), so the default limits are not symmetric.
test_that("zero_centered makes the colour limits symmetric around zero", {
  fx <- make_blob_fixture()

  default_map  <- with_null_device(plot_topography(fx$eeg, fx$values))
  centered_map <- with_null_device(
    plot_topography(fx$eeg, fx$values, zero_centered = TRUE)
  )

  strongest <- max(abs(range(default_map$grid_z, na.rm = TRUE)))

  expect_equal(centered_map$grid_z, default_map$grid_z)
  expect_equal(centered_map$zlim, c(-strongest, strongest))
  expect_false(isTRUE(all.equal(centered_map$zlim, default_map$zlim)))
})

# ----------------------------------------------------------------------------
# Test 4.3: extrapolate = FALSE empties the area outside the electrodes only
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies extrapolate = FALSE only blanks out cells (NA)
# outside the area the electrodes span: more cells are empty than by default,
# every cell that was already empty stays empty, and every cell that is still
# painted holds exactly the value it had before. No new warnings.
test_that("extrapolate = FALSE empties the area outside the electrodes only", {
  fx <- make_blob_fixture()

  default_map <- with_null_device(plot_topography(fx$eeg, fx$values))
  expect_no_warning(
    inside_map <- with_null_device(
      plot_topography(fx$eeg, fx$values, extrapolate = FALSE)
    )
  )

  empty_before <- is.na(default_map$grid_z)
  empty_now    <- is.na(inside_map$grid_z)

  expect_gt(sum(empty_now), sum(empty_before))
  expect_true(all(empty_now[empty_before]))
  expect_equal(inside_map$grid_z[!empty_now], default_map$grid_z[!empty_now])
})

# ----------------------------------------------------------------------------
# Test 4.4: Painting outside the electrodes can invent values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies why extrapolate = FALSE exists. On a smooth bump
# centred on Cz (1 at Cz, falling to about 0 at the edge) the default map
# extrapolates beyond the electrodes to well over the largest value any
# electrode has (the test asks for more than 1.5 times it), while with
# extrapolate = FALSE the largest plotted value stays at the electrodes' own
# maximum.
test_that("extrapolate = FALSE removes values invented beyond the electrodes", {
  fx <- make_blob_fixture()
  electrode_max <- max(fx$values)

  default_map <- with_null_device(plot_topography(fx$eeg, fx$values))
  inside_map  <- with_null_device(
    plot_topography(fx$eeg, fx$values, extrapolate = FALSE)
  )

  expect_gt(max(default_map$grid_z, na.rm = TRUE), 1.5 * electrode_max)
  expect_lte(max(inside_map$grid_z, na.rm = TRUE), 1.05 * electrode_max)
})

# ----------------------------------------------------------------------------
# Test 4.5: The two arguments work together
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that with zero_centered = TRUE and
# extrapolate = FALSE together the colour limits are symmetric around zero AND
# come from the map painted only inside the electrodes, so the values invented
# outside them no longer stretch the colour scale. The map is the bump from
# Test 4.4 shifted down by 0.5 so that it has both signs.
test_that("zero_centered and extrapolate = FALSE work together", {
  fx <- make_blob_fixture()
  signed <- fx$values - 0.5

  centered_only <- with_null_device(
    plot_topography(fx$eeg, signed, zero_centered = TRUE)
  )
  both <- with_null_device(
    plot_topography(fx$eeg, signed,
                    zero_centered = TRUE, extrapolate = FALSE)
  )

  strongest_inside <- max(abs(range(both$grid_z, na.rm = TRUE)))

  expect_equal(both$zlim, c(-strongest_inside, strongest_inside))
  expect_lt(both$zlim[2], centered_only$zlim[2])
})

# ----------------------------------------------------------------------------
# Test 4.6: zero_centered changes the drawn picture, not just zlim
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the colour limits reach the heatmap itself: the
# same map drawn with and without zero_centered gives different images, while
# drawing the same thing twice gives identical ones (so a difference means the
# switch, not noise). Without this, zlim could be computed correctly and still
# never be used for drawing.
test_that("zero_centered changes the drawn picture", {
  fx <- make_blob_fixture()

  default_png  <- png_bytes(plot_topography(fx$eeg, fx$values))
  repeat_png   <- png_bytes(plot_topography(fx$eeg, fx$values))
  centered_png <- png_bytes(
    plot_topography(fx$eeg, fx$values, zero_centered = TRUE)
  )

  expect_identical(default_png, repeat_png)
  expect_false(identical(default_png, centered_png))
})

# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - Input validation: non-'eeg' object, missing montage, unnamed values
# - Warnings: partial overlap between supplied values and montage channels
# - Error: fewer than 3 usable channels for interpolation
# - Successful rendering: return structure, interpolate_res, explicit
#   montage argument overriding/substituting for eeg_obj$montage
# - eeg_obj$bads exclusion: warns and drops bad channels, no-op when empty,
#   can trigger the 3-channel minimum error
# - Colour scale and extrapolation: the defaults keep the old behaviour
#   (picture and zlim unchanged), zero_centered makes the colour limits
#   symmetric around zero, extrapolate = FALSE leaves the area outside the
#   electrodes empty without changing the rest, painting outside the
#   electrodes can invent values that extrapolate = FALSE removes, the two
#   arguments work together, and zero_centered changes the drawn picture
#   (not just the returned zlim)
# ============================================================================
