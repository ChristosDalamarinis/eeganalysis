# ============================================================================
#                          Test File for csd.R
# ============================================================================
#
# Functions tested:
#   1. .csd_matrix()         - Internal: CSD weight table
#   2. .csd_apply_matrix()   - Internal: applies the table to eeg / eeg_epochs
#   3. compute_csd()         - Exported orchestrator
#
# Test suites:
#   1. .csd_matrix() against known reference values
#   2. .csd_matrix() structural properties
#   3. .csd_apply_matrix() eeg and eeg_epochs
#   4. compute_csd() input validation
#   5. compute_csd() data checks (bads, positions, NA, sphericity)
#   6. compute_csd() core behaviour on eeg objects
#   7. compute_csd() non-EEG channels
#   8. compute_csd() head radius and parameters
#   9. compute_csd() on eeg_epochs
#   10. compute_csd() physical behaviour (sharpening, reference-free)
#
# Author: Christos Dalamarinis
# Date: Oct 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

# Full 64-channel BioSemi template: the table values in suite 1 are for this
# exact montage.
.csd_mont <- create_montage()

.csd_pos <- function(mont = .csd_mont) {
  pos <- as.matrix(mont$positions[, c("x", "y", "z")])
  rownames(pos) <- mont$channels
  pos
}

# Random 64-channel continuous object with the montage attached.
.make_csd_eeg <- function(n_tp = 500, seed = 1) {
  set.seed(seed)
  new_eeg(data = matrix(rnorm(64 * n_tp), 64, n_tp),
          channels = .csd_mont$channels,
          sampling_rate = 256,
          montage = .csd_mont)
}

# Three epochs cut from an object with three events. Epoched objects do not
# carry a montage, so the montage is passed to compute_csd() explicitly.
.make_csd_epochs <- function(seed = 1) {
  eeg <- .make_csd_eeg(n_tp = 1000, seed = seed)
  ev <- data.frame(onset = c(300L, 500L, 700L),
                   onset_time = (c(300, 500, 700) - 1) / 256,
                   type = 1L,
                   description = "Trigger: 1")
  eeg_ev <- new_eeg(data = eeg$data, channels = .csd_mont$channels,
                    sampling_rate = 256, events = ev, montage = .csd_mont)
  epoch_eeg(eeg_ev, events = "all", tmin = -0.1, tmax = 0.3,
            baseline = NULL, baseline_method = "none", verbose = FALSE)
}

# ============================================================================
# TEST SUITE 1: .csd_matrix() against known reference values
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: weights match values computed independently (radius = 1)
# ----------------------------------------------------------------------------
test_that(".csd_matrix reproduces independently computed weights", {
  pos <- .csd_pos()
  X <- eeganalysis:::.csd_matrix(pos)
  dimnames(X) <- list(.csd_mont$channels, .csd_mont$channels)

  expect_equal(X["Cz", "Cz"],   7.61238, tolerance = 1e-4)
  expect_equal(X["Cz", "C1"],   1.20699, tolerance = 1e-4)
  expect_equal(X["Cz", "FC1"], -0.81932, tolerance = 1e-4)
  expect_equal(X["Cz", "Iz"],  -0.10804, tolerance = 1e-4)
  expect_equal(X["Fp1", "Fp1"], 5.02264, tolerance = 1e-4)
  expect_equal(X["Oz", "POz"],  1.27884, tolerance = 1e-4)
})

# ----------------------------------------------------------------------------
# Test 1.2: every electrode's own weight is positive and in the known range
# ----------------------------------------------------------------------------
test_that(".csd_matrix own-weights are positive and within 4.41 to 8.60", {
  X <- eeganalysis:::.csd_matrix(.csd_pos())
  own <- diag(X)

  expect_true(all(own > 0))
  expect_gte(min(own), 4.40)
  expect_lte(max(own), 8.61)
})

# ----------------------------------------------------------------------------
# Test 1.3: the weights scale with 1 / radius^2
# ----------------------------------------------------------------------------
test_that(".csd_matrix divides the table by radius squared", {
  pos <- .csd_pos()
  X1 <- eeganalysis:::.csd_matrix(pos, radius = 1)
  X2 <- eeganalysis:::.csd_matrix(pos, radius = 0.5)

  expect_equal(X2, 4 * X1)
})

# ============================================================================
# TEST SUITE 2: .csd_matrix() structural properties
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: square, with one row and column per channel
# ----------------------------------------------------------------------------
test_that(".csd_matrix returns an n x n matrix", {
  pos <- .csd_pos()
  X <- eeganalysis:::.csd_matrix(pos)
  expect_equal(dim(X), c(64L, 64L))

  X10 <- eeganalysis:::.csd_matrix(pos[1:10, ])
  expect_equal(dim(X10), c(10L, 10L))
})

# ----------------------------------------------------------------------------
# Test 2.2: every row sums to zero (why a constant has no effect)
# ----------------------------------------------------------------------------
test_that(".csd_matrix rows sum to zero", {
  X <- eeganalysis:::.csd_matrix(.csd_pos())
  expect_lt(max(abs(rowSums(X))), 1e-8)
})

# ----------------------------------------------------------------------------
# Test 2.3: scaling the positions does not change the table
# ----------------------------------------------------------------------------
test_that(".csd_matrix only uses the direction of each position", {
  pos <- .csd_pos()
  expect_equal(eeganalysis:::.csd_matrix(pos),
               eeganalysis:::.csd_matrix(pos * 3))
})

# ----------------------------------------------------------------------------
# Test 2.4: two electrodes at the same position cannot be inverted
# ----------------------------------------------------------------------------
test_that(".csd_matrix errors when the kernel matrix is singular", {
  pos <- .csd_pos()[1:6, ]
  pos[2, ] <- pos[1, ]
  expect_error(eeganalysis:::.csd_matrix(pos, lambda2 = 0),
               "could not be inverted")
})

# ============================================================================
# TEST SUITE 3: .csd_apply_matrix()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: eeg object - one matrix product, other rows untouched
# ----------------------------------------------------------------------------
test_that(".csd_apply_matrix multiplies the picked rows of an eeg object", {
  eeg <- .make_csd_eeg(n_tp = 50)
  picks <- 1:10
  X <- matrix(rnorm(100), 10, 10)

  out <- eeganalysis:::.csd_apply_matrix(eeg, picks, X)

  expect_equal(out$data[picks, ], X %*% eeg$data[picks, ])
  expect_identical(out$data[-picks, ], eeg$data[-picks, ])
})

# ----------------------------------------------------------------------------
# Test 3.2: eeg_epochs - same as applying the table epoch by epoch
# ----------------------------------------------------------------------------
test_that(".csd_apply_matrix on epochs equals epoch-by-epoch application", {
  ep <- .make_csd_epochs()
  picks <- 1:10
  X <- matrix(rnorm(100), 10, 10)

  out <- eeganalysis:::.csd_apply_matrix(ep, picks, X)

  for (i in seq_len(dim(ep$data)[3])) {
    expect_equal(out$data[picks, , i], X %*% ep$data[picks, , i],
                 ignore_attr = TRUE)
  }
  expect_identical(out$data[-picks, , ], ep$data[-picks, , ])
  expect_equal(dim(out$data), dim(ep$data))
})

# ----------------------------------------------------------------------------
# Test 3.3: the picks do not have to be contiguous
# ----------------------------------------------------------------------------
test_that(".csd_apply_matrix handles non-contiguous picks", {
  eeg <- .make_csd_eeg(n_tp = 50)
  picks <- c(2, 5, 9, 30)
  X <- matrix(rnorm(16), 4, 4)

  out <- eeganalysis:::.csd_apply_matrix(eeg, picks, X)

  expect_equal(out$data[picks, ], X %*% eeg$data[picks, ])
  expect_identical(out$data[-picks, ], eeg$data[-picks, ])
})

# ============================================================================
# TEST SUITE 4: compute_csd() input validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: object class and data shape
# ----------------------------------------------------------------------------
test_that("compute_csd rejects bad objects", {
  expect_error(compute_csd(list()), "class 'eeg' or 'eeg_epochs'")
  expect_error(compute_csd(matrix(1, 2, 2)), "class 'eeg' or 'eeg_epochs'")

  eeg <- .make_csd_eeg()
  eeg$data <- as.vector(eeg$data)
  expect_error(compute_csd(eeg), "numeric matrix")

  ep <- .make_csd_epochs()
  ep_null <- ep
  ep_null$data <- NULL
  expect_error(compute_csd(ep_null, montage = .csd_mont),
               "epoch data not loaded")

  ep_2d <- ep
  ep_2d$data <- ep$data[, , 1]
  expect_error(compute_csd(ep_2d, montage = .csd_mont), "3D array")
})

# ----------------------------------------------------------------------------
# Test 4.2: numeric parameters
# ----------------------------------------------------------------------------
test_that("compute_csd validates its numeric parameters", {
  eeg <- .make_csd_eeg()

  expect_error(compute_csd(eeg, lambda2 = -1), "'lambda2'")
  expect_error(compute_csd(eeg, lambda2 = 1), "'lambda2'")
  expect_error(compute_csd(eeg, lambda2 = c(0.1, 0.2)), "'lambda2'")
  expect_error(compute_csd(eeg, lambda2 = "a"), "'lambda2'")
  expect_error(compute_csd(eeg, lambda2 = NA_real_), "'lambda2'")

  expect_error(compute_csd(eeg, stiffness = -1), "'stiffness'")
  expect_error(compute_csd(eeg, stiffness = c(1, 2)), "'stiffness'")

  expect_error(compute_csd(eeg, n_legendre_terms = 0), "'n_legendre_terms'")
  expect_error(compute_csd(eeg, n_legendre_terms = 2.5), "'n_legendre_terms'")
  expect_error(compute_csd(eeg, n_legendre_terms = c(5, 6)),
               "'n_legendre_terms'")

  expect_error(compute_csd(eeg, origin = c(0, 0)), "'origin'")
  expect_error(compute_csd(eeg, origin = c(0, 0, NA)), "'origin'")
  expect_error(compute_csd(eeg, origin = "a"), "'origin'")

  expect_error(compute_csd(eeg, head_radius = 0), "'head_radius'")
  expect_error(compute_csd(eeg, head_radius = -5), "'head_radius'")
  expect_error(compute_csd(eeg, head_radius = c(80, 90)), "'head_radius'")
  expect_error(compute_csd(eeg, head_radius = "a"), "'head_radius'")

  expect_error(compute_csd(eeg, verbose = "yes"), "'verbose'")
  expect_error(compute_csd(eeg, verbose = NA), "'verbose'")
})

# ----------------------------------------------------------------------------
# Test 4.3: a montage is required
# ----------------------------------------------------------------------------
test_that("compute_csd needs a montage", {
  eeg <- .make_csd_eeg()
  eeg$montage <- NULL
  expect_error(compute_csd(eeg), "No montage available")
  expect_error(compute_csd(eeg, montage = "not a montage"),
               "No montage available")

  # An explicit montage works even when the object carries none
  expect_no_error(compute_csd(eeg, montage = .csd_mont, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 4.4: epochs carry no montage, so one must be passed
# ----------------------------------------------------------------------------
test_that("compute_csd on epochs without a montage errors", {
  ep <- .make_csd_epochs()
  expect_error(compute_csd(ep), "No montage available")
})

# ----------------------------------------------------------------------------
# Test 4.5: CSD cannot be applied twice
# ----------------------------------------------------------------------------
test_that("compute_csd refuses to run twice", {
  eeg <- .make_csd_eeg()
  out <- compute_csd(eeg, verbose = FALSE)
  expect_error(compute_csd(out), "already been applied")
})

# ----------------------------------------------------------------------------
# Test 4.6: no EEG channels at all
# ----------------------------------------------------------------------------
test_that("compute_csd errors when there are no EEG channels", {
  eeg <- .make_csd_eeg()
  eeg$channel_types <- rep("eog", 64)
  expect_error(compute_csd(eeg), "No EEG channels found")
})

# ============================================================================
# TEST SUITE 5: compute_csd() data checks
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: bad EEG channels are refused and named
# ----------------------------------------------------------------------------
test_that("compute_csd refuses bad EEG channels", {
  eeg <- .make_csd_eeg()
  eeg$bads <- c("Fz", "Cz")
  expect_error(compute_csd(eeg), "bad EEG channels: Fz, Cz")
  expect_error(compute_csd(eeg), "interpolate_bads")
})

# ----------------------------------------------------------------------------
# Test 5.2: a bad channel that is not an EEG channel does not block CSD
# ----------------------------------------------------------------------------
test_that("compute_csd ignores bads that are not EEG channels", {
  eeg <- set_channel_types(.make_csd_eeg(), c(Fp1 = "eog"))
  eeg$bads <- "Fp1"
  expect_no_error(compute_csd(eeg, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 5.3: an EEG channel missing from the montage
# ----------------------------------------------------------------------------
test_that("compute_csd errors for EEG channels with no montage position", {
  eeg <- .make_csd_eeg()
  keep <- setdiff(.csd_mont$channels, c("Fz", "Cz"))
  small <- create_montage(keep)
  eeg$montage <- small
  expect_error(compute_csd(eeg), "no position in the montage: Fz, Cz")
})

# ----------------------------------------------------------------------------
# Test 5.4: NA and infinite data
# ----------------------------------------------------------------------------
test_that("compute_csd errors on NA or non-finite EEG data", {
  eeg <- .make_csd_eeg()
  eeg$data[3, 10] <- NA
  expect_error(compute_csd(eeg), "NA or non-finite")

  eeg <- .make_csd_eeg()
  eeg$data[3, 10] <- Inf
  expect_error(compute_csd(eeg), "NA or non-finite")

  ep <- .make_csd_epochs()
  ep$data[3, 10, 2] <- NaN
  expect_error(compute_csd(ep, montage = .csd_mont), "NA or non-finite")
})

# ----------------------------------------------------------------------------
# Test 5.5: NA in a non-EEG channel is allowed (it is never touched)
# ----------------------------------------------------------------------------
test_that("compute_csd ignores NA in non-EEG channels", {
  eeg <- set_channel_types(.make_csd_eeg(), c(Fp1 = "eog"))
  eeg$data[1, 5] <- NA
  expect_no_error(compute_csd(eeg, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 5.6: two channels at the same position are named in the error
# ----------------------------------------------------------------------------
test_that("compute_csd errors when two channels share a position", {
  eeg <- .make_csd_eeg()
  mont <- .csd_mont
  mont$positions[2, c("x", "y", "z")] <- mont$positions[1, c("x", "y", "z")]
  eeg$montage <- mont
  expect_error(compute_csd(eeg), "share the same position.*Fp1, AF7")
})

# ----------------------------------------------------------------------------
# Test 5.7: a channel at the centre of the sphere
# ----------------------------------------------------------------------------
test_that("compute_csd errors for a channel sitting at the origin", {
  eeg <- .make_csd_eeg()
  mont <- .csd_mont
  mont$positions[5, c("x", "y", "z")] <- c(0, 0, 0)
  eeg$montage <- mont
  expect_error(compute_csd(eeg), "Zero or non-finite channel position")
})

# ----------------------------------------------------------------------------
# Test 5.8: electrodes far from a sphere trigger the warning
# ----------------------------------------------------------------------------
test_that("compute_csd warns when positions are not spherical", {
  eeg <- .make_csd_eeg()

  # No warning for the template (largest deviation is about 0.6 %)
  expect_no_warning(compute_csd(eeg, verbose = FALSE))

  mont <- .csd_mont
  mont$positions[10, c("x", "y", "z")] <-
    1.5 * mont$positions[10, c("x", "y", "z")]
  eeg$montage <- mont
  expect_warning(compute_csd(eeg, verbose = FALSE),
                 "not close to spherical")
})

# ============================================================================
# TEST SUITE 6: compute_csd() core behaviour on eeg objects
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: output shape, class, and the CSD marks
# ----------------------------------------------------------------------------
test_that("compute_csd returns a marked eeg object of the same shape", {
  eeg <- .make_csd_eeg()
  out <- compute_csd(eeg, verbose = FALSE)

  expect_s3_class(out, "eeg")
  expect_equal(dim(out$data), dim(eeg$data))
  expect_equal(out$reference, "CSD")
  expect_equal(out$metadata$reference_scheme, "CSD")
  expect_false(isTRUE(all.equal(out$data, eeg$data)))
})

# ----------------------------------------------------------------------------
# Test 6.2: data equals the weight table times the voltages
# ----------------------------------------------------------------------------
test_that("compute_csd output is the weight table applied to the data", {
  eeg <- .make_csd_eeg()
  out <- compute_csd(eeg, verbose = FALSE)

  pos <- .csd_pos()
  radius_m <- mean(sqrt(rowSums(pos^2))) / 1000
  X <- eeganalysis:::.csd_matrix(pos, radius = radius_m)

  expect_equal(out$data, X %*% eeg$data, ignore_attr = TRUE)
})

# ----------------------------------------------------------------------------
# Test 6.3: history entry
# ----------------------------------------------------------------------------
test_that("compute_csd appends one history entry", {
  eeg <- .make_csd_eeg()
  out <- compute_csd(eeg, verbose = FALSE)

  n_before <- length(eeg$preprocessing_history)
  expect_length(out$preprocessing_history, n_before + 1)

  entry <- out$preprocessing_history[[n_before + 1]]
  expect_match(entry, "^compute_csd\\(\\)")
  expect_match(entry, "64 EEG channel")
  expect_match(entry, "estimated from the montage")
  expect_match(entry, "uV/m\\^2")
})

# ----------------------------------------------------------------------------
# Test 6.4: all other fields are carried through unchanged
# ----------------------------------------------------------------------------
test_that("compute_csd leaves every other field unchanged", {
  eeg <- .make_csd_eeg()
  eeg$bads <- character(0)
  out <- compute_csd(eeg, verbose = FALSE)

  for (f in c("channels", "channel_types", "bads", "annotations",
              "sampling_rate", "times", "events", "montage")) {
    expect_identical(out[[f]], eeg[[f]], info = f)
  }
})

# ----------------------------------------------------------------------------
# Test 6.5: the input object is not modified
# ----------------------------------------------------------------------------
test_that("compute_csd does not modify its input", {
  eeg <- .make_csd_eeg()
  eeg_copy <- eeg
  invisible(compute_csd(eeg, verbose = FALSE))
  expect_identical(eeg, eeg_copy)
})

# ----------------------------------------------------------------------------
# Test 6.6: verbose controls the console message
# ----------------------------------------------------------------------------
test_that("compute_csd prints a summary only when verbose", {
  eeg <- .make_csd_eeg()
  expect_message(compute_csd(eeg), "CSD applied to 64 EEG channel")
  expect_no_message(compute_csd(eeg, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 6.7: a channel-wise rescale of the montage does not change the result
# ----------------------------------------------------------------------------
test_that("compute_csd is unchanged by the units of the montage positions", {
  eeg <- .make_csd_eeg()
  mont <- .csd_mont
  mont$positions[, c("x", "y", "z")] <- 2 * mont$positions[, c("x", "y", "z")]

  # Same head radius given explicitly: positions in a different scale give
  # the same directions, hence the same table
  a <- compute_csd(eeg, head_radius = 87.5, verbose = FALSE)
  b <- compute_csd(eeg, montage = mont, head_radius = 87.5, verbose = FALSE)
  expect_equal(a$data, b$data)
})

# ============================================================================
# TEST SUITE 7: compute_csd() non-EEG channels
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: non-EEG channels are untouched, the rest uses only EEG positions
# ----------------------------------------------------------------------------
test_that("compute_csd leaves non-EEG channels exactly as they are", {
  eeg <- set_channel_types(.make_csd_eeg(), c(Fp1 = "eog"))
  out <- compute_csd(eeg, verbose = FALSE)

  expect_identical(out$data[1, ], eeg$data[1, ])
  expect_false(isTRUE(all.equal(out$data[-1, ], eeg$data[-1, ])))

  # The other 63 channels use only their own 63 positions
  pos <- .csd_pos()[-1, ]
  radius_m <- mean(sqrt(rowSums(pos^2))) / 1000
  X <- eeganalysis:::.csd_matrix(pos, radius = radius_m)
  expect_equal(out$data[-1, ], X %*% eeg$data[-1, ], ignore_attr = TRUE)

  expect_identical(out$channel_types, eeg$channel_types)
})

# ============================================================================
# TEST SUITE 8: compute_csd() head radius and parameters
# ============================================================================

# ----------------------------------------------------------------------------
# Test 8.1: halving the head radius multiplies the output by 4
# ----------------------------------------------------------------------------
test_that("compute_csd output scales with 1 / head_radius^2", {
  eeg <- .make_csd_eeg()
  a <- compute_csd(eeg, head_radius = 87.49, verbose = FALSE)
  b <- compute_csd(eeg, head_radius = 43.745, verbose = FALSE)

  expect_equal(b$data, 4 * a$data)
})

# ----------------------------------------------------------------------------
# Test 8.2: the automatic radius is the mean electrode distance
# ----------------------------------------------------------------------------
test_that("compute_csd estimates the head radius from the montage", {
  eeg <- .make_csd_eeg()
  auto <- compute_csd(eeg, verbose = FALSE)

  pos <- .csd_pos()
  r <- mean(sqrt(rowSums(pos^2)))
  expect_equal(r, 87.49, tolerance = 1e-3)

  explicit <- compute_csd(eeg, head_radius = r, verbose = FALSE)
  expect_equal(auto$data, explicit$data)

  entry <- auto$preprocessing_history[[length(auto$preprocessing_history)]]
  expect_match(entry, "87.49 mm \\(estimated from the montage\\)")
})

# ----------------------------------------------------------------------------
# Test 8.3: a user-supplied radius is reported as such
# ----------------------------------------------------------------------------
test_that("compute_csd records a user-supplied head radius", {
  eeg <- .make_csd_eeg()
  out <- compute_csd(eeg, head_radius = 90, verbose = FALSE)
  entry <- out$preprocessing_history[[length(out$preprocessing_history)]]
  expect_match(entry, "90.00 mm \\(user supplied\\)")
})

# ----------------------------------------------------------------------------
# Test 8.4: regularisation and stiffness change the result and are recorded
# ----------------------------------------------------------------------------
test_that("compute_csd responds to lambda2, stiffness and n_legendre_terms", {
  eeg <- .make_csd_eeg()
  base <- compute_csd(eeg, verbose = FALSE)

  smoother <- compute_csd(eeg, lambda2 = 1e-2, verbose = FALSE)
  stiffer  <- compute_csd(eeg, stiffness = 3, verbose = FALSE)
  fewer    <- compute_csd(eeg, n_legendre_terms = 30, verbose = FALSE)

  expect_false(isTRUE(all.equal(base$data, smoother$data)))
  expect_false(isTRUE(all.equal(base$data, stiffer$data)))
  expect_false(isTRUE(all.equal(base$data, fewer$data)))

  entry <- smoother$preprocessing_history[[
    length(smoother$preprocessing_history)]]
  expect_match(entry, "lambda2 = 0.01")
  expect_match(entry, "stiffness = 4")
  expect_match(entry, "n_legendre_terms = 50")
})

# ----------------------------------------------------------------------------
# Test 8.5: more smoothing means a smaller response to channel-wise noise
# ----------------------------------------------------------------------------
test_that("a larger lambda2 gives a smoother (smaller) result on noise", {
  eeg <- .make_csd_eeg()
  low  <- compute_csd(eeg, lambda2 = 1e-5, verbose = FALSE)
  high <- compute_csd(eeg, lambda2 = 1e-1, verbose = FALSE)

  expect_lt(sd(high$data), sd(low$data))
})

# ============================================================================
# TEST SUITE 9: compute_csd() on eeg_epochs
# ============================================================================

# ----------------------------------------------------------------------------
# Test 9.1: shape, class and marks
# ----------------------------------------------------------------------------
test_that("compute_csd works on epochs given a montage", {
  ep <- .make_csd_epochs()
  out <- compute_csd(ep, montage = .csd_mont, verbose = FALSE)

  expect_s3_class(out, "eeg_epochs")
  expect_equal(dim(out$data), c(64, 104, 3))
  expect_equal(out$reference, "CSD")
  expect_length(out$preprocessing_history,
                length(ep$preprocessing_history) + 1)
})

# ----------------------------------------------------------------------------
# Test 9.2: equals applying the continuous result epoch by epoch
# ----------------------------------------------------------------------------
test_that("compute_csd on epochs equals processing each epoch separately", {
  ep <- .make_csd_epochs()
  out <- compute_csd(ep, montage = .csd_mont, verbose = FALSE)

  pos <- .csd_pos()
  radius_m <- mean(sqrt(rowSums(pos^2))) / 1000
  X <- eeganalysis:::.csd_matrix(pos, radius = radius_m)

  for (i in seq_len(dim(ep$data)[3])) {
    expect_equal(out$data[, , i], X %*% ep$data[, , i], ignore_attr = TRUE)
  }
})

# ----------------------------------------------------------------------------
# Test 9.3: linear, so CSD then average equals average then CSD
# ----------------------------------------------------------------------------
test_that("CSD commutes with averaging over epochs", {
  ep <- .make_csd_epochs()
  out <- compute_csd(ep, montage = .csd_mont, verbose = FALSE)

  mean_of_csd <- apply(out$data, c(1, 2), mean)

  avg_data <- apply(ep$data, c(1, 2), mean)
  avg_obj <- .make_csd_eeg(n_tp = ncol(avg_data))
  avg_obj$data <- avg_data
  csd_of_mean <- compute_csd(avg_obj, verbose = FALSE)$data

  expect_equal(mean_of_csd, csd_of_mean, ignore_attr = TRUE)
})

# ----------------------------------------------------------------------------
# Test 9.4: every epoch field is carried through, input untouched
# ----------------------------------------------------------------------------
test_that("compute_csd keeps all eeg_epochs fields and the input intact", {
  ep <- .make_csd_epochs()
  ep_copy <- ep
  out <- compute_csd(ep, montage = .csd_mont, verbose = FALSE)

  expect_identical(ep, ep_copy)

  for (f in c("channels", "times", "events", "sampling_rate", "tmin", "tmax",
              "baseline", "baseline_method", "n_epochs", "rejected",
              "rejection_log", "channel_types", "bads")) {
    expect_identical(out[[f]], ep[[f]], info = f)
  }
  expect_identical(names(out), names(ep))
})

# ----------------------------------------------------------------------------
# Test 9.5: non-EEG channels in epochs are untouched
# ----------------------------------------------------------------------------
test_that("compute_csd leaves non-EEG channels of epochs untouched", {
  ep <- .make_csd_epochs()
  ep$channel_types[1] <- "eog"
  out <- compute_csd(ep, montage = .csd_mont, verbose = FALSE)

  expect_identical(out$data[1, , ], ep$data[1, , ])
  expect_false(isTRUE(all.equal(out$data[-1, , ], ep$data[-1, , ])))
})

# ============================================================================
# TEST SUITE 10: compute_csd() physical behaviour
# ============================================================================

# ----------------------------------------------------------------------------
# Test 10.1: reference-free - a constant added to every channel disappears
# ----------------------------------------------------------------------------
test_that("compute_csd is unaffected by a constant added to all channels", {
  eeg <- .make_csd_eeg()
  shifted <- eeg
  shifted$data <- shifted$data + 50

  a <- compute_csd(eeg, verbose = FALSE)
  b <- compute_csd(shifted, verbose = FALSE)

  # Tolerance is relative to the size of the output (the table is divided by
  # the head radius in metres squared, so values are large)
  expect_lt(max(abs(a$data - b$data)), 1e-6)
})

# ----------------------------------------------------------------------------
# Test 10.2: a different reference (average reference) gives the same CSD
# ----------------------------------------------------------------------------
test_that("compute_csd gives the same result after re-referencing", {
  eeg <- .make_csd_eeg()
  rereferenced <- eeg
  rereferenced$data <- sweep(eeg$data, 2, colMeans(eeg$data), "-")

  a <- compute_csd(eeg, verbose = FALSE)
  b <- compute_csd(rereferenced, verbose = FALSE)

  expect_lt(max(abs(a$data - b$data)), 1e-6)
})

# ----------------------------------------------------------------------------
# Test 10.3: sharpening - a broad blob gives a focal peak with the right sign
# ----------------------------------------------------------------------------
test_that("compute_csd turns a broad blob into a positive peak at its centre", {
  pos <- .csd_pos()
  unit <- pos / sqrt(rowSums(pos^2))
  cz <- unit["Cz", ]
  angle <- acos(pmin(pmax(as.vector(unit %*% cz), -1), 1))

  blob <- exp(-(angle / 0.6)^2)               # broad bump centred on Cz
  eeg <- .make_csd_eeg(n_tp = 4)
  eeg$data <- matrix(blob, 64, 4)

  out <- compute_csd(eeg, verbose = FALSE)
  csd <- out$data[, 1]
  names(csd) <- .csd_mont$channels

  expect_equal(names(which.max(csd)), "Cz")
  expect_gt(csd["Cz"], 0)

  # More focal: fewer electrodes stay near the peak than in the voltage map
  near_volt <- sum(blob >= 0.5 * max(blob))
  near_csd  <- sum(csd >= 0.5 * max(csd))
  expect_lt(near_csd, near_volt)
})

# ----------------------------------------------------------------------------
# Test 10.4: a spatially constant field gives exactly zero
# ----------------------------------------------------------------------------
test_that("a spatially flat pattern gives zero CSD", {
  eeg <- .make_csd_eeg(n_tp = 10)
  eeg$data <- matrix(rep(seq_len(10), each = 64), 64, 10)

  out <- compute_csd(eeg, verbose = FALSE)
  expect_lt(max(abs(out$data)), 1e-6)
})

# ----------------------------------------------------------------------------
# Test 10.5: linearity in the data
# ----------------------------------------------------------------------------
test_that("compute_csd is linear in the data", {
  a <- .make_csd_eeg(seed = 1)
  b <- .make_csd_eeg(seed = 2)
  combo <- a
  combo$data <- 2 * a$data + 3 * b$data

  ca <- compute_csd(a, verbose = FALSE)$data
  cb <- compute_csd(b, verbose = FALSE)$data
  cc <- compute_csd(combo, verbose = FALSE)$data

  expect_equal(cc, 2 * ca + 3 * cb)
})
