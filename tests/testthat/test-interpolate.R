# ============================================================================
#                       Test File for interpolate.R
# ============================================================================
#
# Functions tested:
#   1. interpolate_bads()             - Exported orchestrator
#   2. make_interpolation_matrix()    - Internal: spherical-spline weights
#   3. calc_g()                       - Internal: Legendre-series kernel
#   4. calc_h()                       - Internal: companion kernel of calc_g()
#   5. .legendre_series_eval()        - Internal: Legendre series evaluation
#
# Test suites:
#   1. .legendre_series_eval() against closed-form Legendre polynomials
#   2. calc_g() shape/symmetry, calc_h() relation to calc_g()
#   3. make_interpolation_matrix() shape + row-sums-to-1 property
#   4. interpolate_bads() input validation
#   5. interpolate_bads() no-op / partial montage coverage
#   6. interpolate_bads() core reconstruction (single bad channel)
#   7. interpolate_bads() multiple simultaneous bad channels
#   8. interpolate_bads() reset_bads = FALSE
#   9. interpolate_bads() exclude
#   10. interpolate_bads() no good channels left
#
# Author: Christos Dalamarinis
# Date: Sep 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

.standard_channels <- c("Fp1", "Fp2", "F3", "F4", "Fz",
                         "C3", "C4", "Cz", "P3", "P4", "Pz", "O1", "O2")

# Every channel carries the SAME shared smooth waveform (a spatially
# constant field) plus small independent per-channel noise, montage
# attached. A spatially constant field is the simplest case any correct
# spatial interpolator must reconstruct well - gives a real, checkable
# ground truth instead of just a "did it run" test. Returns both the eeg
# object and the true shared waveform to compare reconstructions against.
.make_constant_field_eeg <- function(n_tp = 300, seed = 1, noise_sd = 0.5) {
  set.seed(seed)
  shared <- sin(seq(0, 4 * pi, length.out = n_tp)) * 20
  data <- t(vapply(seq_along(.standard_channels), function(i) {
    shared + rnorm(n_tp, sd = noise_sd)
  }, numeric(n_tp)))
  eeg <- new_eeg(data = data, channels = .standard_channels, sampling_rate = 256)
  eeg <- set_montage(eeg, create_montage(.standard_channels))
  list(eeg = eeg, shared = shared)
}

# ============================================================================
# TEST SUITE 1: .legendre_series_eval()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: matches closed-form Legendre polynomials P0-P3
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: the 3-term recurrence used internally is correct, checked
# against the textbook closed-form polynomials it's supposed to reproduce.
test_that(".legendre_series_eval matches closed-form Legendre polynomials", {
  x <- c(-1, -0.5, 0, 0.5, 1)

  expect_equal(eeganalysis:::.legendre_series_eval(x, c(1)), rep(1, length(x)))
  expect_equal(eeganalysis:::.legendre_series_eval(x, c(0, 1)), x)
  expect_equal(eeganalysis:::.legendre_series_eval(x, c(0, 0, 1)),
               (3 * x^2 - 1) / 2)
  expect_equal(eeganalysis:::.legendre_series_eval(x, c(0, 0, 0, 1)),
               (5 * x^3 - 3 * x) / 2, tolerance = 1e-10)

  # Linear combination: 2*P0 - 1*P2
  expect_equal(eeganalysis:::.legendre_series_eval(x, c(2, 0, -1)),
               2 - (3 * x^2 - 1) / 2)
})

# ----------------------------------------------------------------------------
# Test 1.2: preserves matrix shape
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: calc_g() calls this on full cosang matrices, not just
# vectors - shape must survive elementwise evaluation.
test_that(".legendre_series_eval preserves matrix shape", {
  m <- matrix(c(-1, -0.5, 0, 0.5, 1, 0.25), nrow = 2)
  out <- eeganalysis:::.legendre_series_eval(m, c(0, 1))  # P1(x) = x
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(m))
  expect_equal(out, m)
})

# ============================================================================
# TEST SUITE 2: calc_g() and calc_h()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: symmetric cosang produces a symmetric G
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: cosang_from = pos_from %*% t(pos_from) is always
# symmetric in make_interpolation_matrix() - calc_g() must preserve that.
test_that("calc_g returns a symmetric matrix for symmetric cosang", {
  set.seed(2)
  pos <- matrix(rnorm(5 * 3), nrow = 5)
  pos <- pos / sqrt(rowSums(pos^2))
  cosang <- pos %*% t(pos)

  g <- eeganalysis:::calc_g(cosang)
  expect_equal(dim(g), c(5, 5))
  expect_equal(g, t(g))
})

# ----------------------------------------------------------------------------
# Test 2.2: non-square input shape is preserved
# ----------------------------------------------------------------------------
test_that("calc_g output shape matches a non-square cosang input", {
  cosang <- matrix(runif(3 * 5, -1, 1), nrow = 3, ncol = 5)
  g <- eeganalysis:::calc_g(cosang)
  expect_equal(dim(g), c(3, 5))
})

# ----------------------------------------------------------------------------
# Test 2.3: calc_h is calc_g with the stiffness lowered by one
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: the two kernels differ only by one power in the weight of
# each term, so calc_h(x, m) must give the same numbers as calc_g(x, m - 1).
# Guards the "(stiffness - 1)" brackets in calc_h(): without them R still
# runs, it just silently returns wrong numbers.
test_that("calc_h equals calc_g with the stiffness lowered by one", {
  x <- c(-0.9, -0.5, 0, 0.3, 0.7, 0.95)

  # Default stiffness is 4, so the default calc_h() matches calc_g() at 3.
  expect_equal(eeganalysis:::calc_h(x),
               eeganalysis:::calc_g(x, stiffness = 3))

  # Holds for another stiffness too, so the exponent is not hard-coded.
  expect_equal(eeganalysis:::calc_h(x, stiffness = 3),
               eeganalysis:::calc_g(x, stiffness = 2))
})

# ----------------------------------------------------------------------------
# Test 2.4: calc_h is minus the surface Laplacian of calc_g
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: that calc_h() is the right function mathematically, not
# just a copy of its own formula. For a function f of x = cos(angle) on a
# sphere, the surface Laplacian is d/dx[(1 - x^2) * df/dx], and the "h"
# kernel is minus that applied to the "g" kernel. The derivatives are taken
# numerically (central differences), so this check does not depend on the
# formula inside calc_h(). Measured agreement is about 1e-6 (relative); the
# tolerance leaves a wide margin.
test_that("calc_h is minus the surface Laplacian of calc_g", {
  calc_g <- eeganalysis:::calc_g
  x      <- c(-0.9, -0.5, 0, 0.3, 0.7, 0.95)
  eps    <- 1e-3

  dg    <- function(u) (calc_g(u + eps) - calc_g(u - eps)) / (2 * eps)
  flux  <- function(u) (1 - u^2) * dg(u)
  lap_g <- (flux(x + eps) - flux(x - eps)) / (2 * eps)

  expect_equal(-lap_g, eeganalysis:::calc_h(x), tolerance = 1e-4)
})

# ============================================================================
# TEST SUITE 3: make_interpolation_matrix()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1-3.2: output shape is n_to x n_from, including a single target
# ----------------------------------------------------------------------------
test_that("make_interpolation_matrix returns an n_to x n_from matrix", {
  set.seed(3)
  pos_from <- matrix(rnorm(6 * 3), nrow = 6)
  pos_to   <- matrix(rnorm(2 * 3), nrow = 2)

  m <- eeganalysis:::make_interpolation_matrix(pos_from, pos_to)
  expect_equal(dim(m), c(2, 6))
})

test_that("make_interpolation_matrix handles a single target position", {
  set.seed(4)
  pos_from <- matrix(rnorm(6 * 3), nrow = 6)
  pos_to   <- matrix(rnorm(1 * 3), nrow = 1)

  m <- eeganalysis:::make_interpolation_matrix(pos_from, pos_to)
  expect_equal(dim(m), c(1, 6))
})

# ----------------------------------------------------------------------------
# Test 3.3: rows sum to 1 - reproduces a spatially constant field exactly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: a provable property of the Perrin et al. (1989)
# construction - the sum(spline coefficients) = 0 constraint guarantees
# that applying the weight matrix to a constant vector returns that same
# constant, i.e. every row sums to 1. Uses real 10-20 positions (via
# create_montage()), not random points, so G_from is well-conditioned.
test_that("make_interpolation_matrix rows sum to 1 (reproduces a constant field)", {
  montage <- create_montage(.standard_channels)
  pos <- as.matrix(montage$positions[, c("x", "y", "z")])

  pos_from <- pos                               # all 13 real channels
  pos_to   <- pos[c(1, 6, 11), , drop = FALSE]  # a few arbitrary targets

  m <- eeganalysis:::make_interpolation_matrix(pos_from, pos_to)
  expect_equal(rowSums(m), rep(1, 3), tolerance = 1e-3)
})

# ============================================================================
# TEST SUITE 4: interpolate_bads() input validation
# ============================================================================

test_that("interpolate_bads validates its inputs", {
  eeg <- .make_constant_field_eeg()$eeg

  expect_error(interpolate_bads(list()), "class 'eeg'")

  eeg_no_montage <- new_eeg(data = eeg$data, channels = eeg$channels,
                             sampling_rate = eeg$sampling_rate, bads = "Fp1")
  expect_error(interpolate_bads(eeg_no_montage), "No montage attached")

  expect_error(interpolate_bads(eeg, origin = c(0, 0)), "length 3")
})

# ============================================================================
# TEST SUITE 5: no-op / partial montage coverage
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: no bad channels at all - message, object unchanged
# ----------------------------------------------------------------------------
test_that("interpolate_bads is a no-op (with message) when there are no bad channels", {
  eeg <- .make_constant_field_eeg()$eeg
  expect_message(result <- interpolate_bads(eeg), "no bad EEG channel")
  expect_identical(result$data, eeg$data)
  expect_equal(result$bads, eeg$bads)
})

# ----------------------------------------------------------------------------
# Test 5.2: a bad channel with no montage position is warned about and skipped
# ----------------------------------------------------------------------------
test_that("interpolate_bads warns and skips a bad channel with no montage position", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg

  # Re-attach a montage that has no position for Fp1, then mark Fp1 bad by
  # hand - it's a real EEG channel, just one the montage doesn't cover.
  partial_montage <- create_montage(setdiff(.standard_channels, "Fp1"))
  eeg <- suppressWarnings(set_montage(eeg, partial_montage))
  eeg$bads <- "Fp1"

  expect_warning(
    expect_message(result <- interpolate_bads(eeg), "no bad EEG channel"),
    "no position in the attached montage.*Fp1"
  )
  expect_equal(result$bads, "Fp1")
})

# ============================================================================
# TEST SUITE 6: core reconstruction (single bad channel)
# ============================================================================

test_that("interpolate_bads reconstructs a single bad channel from a spatially constant field", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg
  shared <- fixture$shared

  o1_idx <- match("O1", eeg$channels)
  eeg$data[o1_idx, ] <- rnorm(ncol(eeg$data), sd = 500)  # corrupt O1
  corrupted <- eeg$data[o1_idx, ]
  eeg$bads <- "O1"

  result <- interpolate_bads(eeg)

  expect_false(isTRUE(all.equal(result$data[o1_idx, ], corrupted)))
  expect_gt(cor(result$data[o1_idx, ], shared), 0.9)

  other_idx <- setdiff(seq_along(result$channels), o1_idx)
  expect_equal(result$data[other_idx, ], eeg$data[other_idx, ])

  expect_false("O1" %in% result$bads)
  expect_equal(length(result$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
})

# ============================================================================
# TEST SUITE 7: multiple simultaneous bad channels
# ============================================================================

test_that("interpolate_bads reconstructs multiple simultaneous bad channels", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg
  shared <- fixture$shared

  bad_names <- c("O1", "O2", "Fp1")
  bad_idx <- match(bad_names, eeg$channels)
  for (i in bad_idx) {
    eeg$data[i, ] <- rnorm(ncol(eeg$data), sd = 500)
  }
  eeg$bads <- bad_names

  result <- interpolate_bads(eeg)

  for (i in bad_idx) {
    expect_gt(cor(result$data[i, ], shared), 0.9)
  }
  expect_equal(result$bads, character(0))
})

# ============================================================================
# TEST SUITE 8: reset_bads = FALSE
# ============================================================================

test_that("interpolate_bads keeps channels in bads when reset_bads = FALSE", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg
  shared <- fixture$shared

  o1_idx <- match("O1", eeg$channels)
  eeg$data[o1_idx, ] <- rnorm(ncol(eeg$data), sd = 500)
  eeg$bads <- "O1"

  result <- interpolate_bads(eeg, reset_bads = FALSE)

  expect_gt(cor(result$data[o1_idx, ], shared), 0.9)
  expect_equal(result$bads, "O1")
})

# ============================================================================
# TEST SUITE 9: exclude
# ============================================================================

test_that("interpolate_bads leaves excluded channels untouched, even if bad", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg

  o1_idx  <- match("O1", eeg$channels)
  fp1_idx <- match("Fp1", eeg$channels)
  eeg$data[o1_idx, ]  <- rnorm(ncol(eeg$data), sd = 500)
  eeg$data[fp1_idx, ] <- rnorm(ncol(eeg$data), sd = 500)
  corrupted_fp1 <- eeg$data[fp1_idx, ]
  eeg$bads <- c("O1", "Fp1")

  result <- interpolate_bads(eeg, exclude = "Fp1")

  expect_equal(result$data[fp1_idx, ], corrupted_fp1)  # untouched
  expect_true("Fp1" %in% result$bads)                  # still marked bad
  expect_false("O1" %in% result$bads)                  # O1 still fixed normally
})

# ============================================================================
# TEST SUITE 10: no good channels left
# ============================================================================

test_that("interpolate_bads errors when no good channels with a position remain", {
  fixture <- .make_constant_field_eeg()
  eeg <- fixture$eeg
  eeg$bads <- .standard_channels  # every montage-eligible channel is bad

  expect_error(interpolate_bads(eeg), "No good EEG channels")
})
