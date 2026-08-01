# ============================================================================
#                       Test File for ica1.R
# ============================================================================
#
# This test file provides comprehensive testing for ica1.R, which implements
# ICA-based artifact removal for EEG data (construction, fitting, source
# inspection, and reconstruction), following the mne.preprocessing.ICA design
# and a native port of scikit-learn's FastICA rotation.
#
# There is no vendored MNE/sklearn reference CSV for ICA (unlike filter1.R /
# fourier.R), so numerical correctness here is established two ways instead:
#   1. Exact hand-derived values for the small, purely arithmetic private
#      helpers (G-functions, decorrelation, pre-whitening, n_components
#      resolution).
#   2. Property-based / analytic-invariant checks for the linear-algebra
#      pipeline (orthonormality of the unmixing matrix, PCA whitening to
#      unit variance, the documented round-trip identity when
#      n_components_ equals the number of fitted channels) plus a classic
#      blind-source-separation sanity check (mixing two known independent
#      non-Gaussian signals and confirming fit_ica() recovers them).
#
# Functions tested:
#   new_ica(), print.eeg_ica(), fit_ica(), get_sources(),
#   get_component_topography(), plot_ica_topography(), ica_component_summary(),
#   set_exclude(), apply_ica(), plot_ica_sources(), plot_ica_overlay(),
#   and the internal (non-exported) helpers used to build them.
#
# Author: Christos Dalamarinis
# Date: July 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

# Runs `expr` against a throwaway PNG device so no plot window is displayed
# during the test run, and always cleans the device up afterwards (same
# pattern as test-topography.R).
with_null_device <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  force(expr)
}

# General-purpose ICA fixture: 3 "eeg" channels (Fp1, Fp2, Cz) that are a
# fixed linear mixture of 2 independent non-Gaussian sources, an "external"
# channel (EXG1) that tracks source 1 like an EOG channel would, a "status"
# channel of random bits (must always be excluded from default fit_ica()
# picks), and a small amount of independent noise added to every fitted
# channel so the fitted covariance is always strictly full rank (no exact
# zero eigenvalues -> no Inf/NaN in .pca_whiten(), and the documented
# "reconstructs to floating-point precision when n_components_ equals the
# number of fitted channels" round trip actually holds).
make_ica_fixture <- function(n_samples = 3000, sampling_rate = 256, seed = 123) {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate
  s1 <- sin(2 * pi * 5 * t)
  s2 <- sign(sin(2 * pi * 11 * t))
  S <- rbind(s1, s2)

  A <- matrix(c(1.0,  0.5,
                0.3,  1.0,
                0.8, -0.6), nrow = 3, byrow = TRUE)

  eeg_signal <- A %*% S
  eeg_noise  <- matrix(rnorm(3 * n_samples, sd = 0.02), nrow = 3)
  eeg_data   <- eeg_signal + eeg_noise

  exg1   <- s1 + rnorm(n_samples, sd = 0.05)
  status <- sample(0:1, n_samples, replace = TRUE)

  data     <- rbind(eeg_data, exg1, status)
  channels <- c("Fp1", "Fp2", "Cz", "EXG1", "Status")

  eeg <- new_eeg(data = data, channels = channels, sampling_rate = sampling_rate)
  list(eeg = eeg, sources = S, mixing = A)
}

# Clean blind-source-separation fixture: exactly 2 known independent
# non-Gaussian sources (a sine wave and a square wave, the classic FastICA
# textbook pair) linearly mixed into 3 "eeg" channels with only a trace of
# noise, so fit_ica(n_components = 2) should recover both sources with high
# fidelity (up to the usual ICA sign/permutation ambiguity).
make_source_separation_fixture <- function(n_samples = 4000, sampling_rate = 256, seed = 7) {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate
  s1 <- sin(2 * pi * 5 * t)
  s2 <- sign(sin(2 * pi * 11 * t))
  S <- rbind(s1, s2)

  A <- matrix(c(1.0,  0.5,
                0.3,  1.0,
                0.8, -0.6), nrow = 3, byrow = TRUE)

  noise <- matrix(rnorm(3 * n_samples, sd = 0.01), nrow = 3)
  data  <- A %*% S + noise

  eeg <- new_eeg(data = data, channels = c("Fp1", "Fp2", "Cz"),
                 sampling_rate = sampling_rate)
  list(eeg = eeg, sources = S)
}

# ============================================================================
#                    TEST SUITE 1: new_ica() - Constructor
# ============================================================================

test_that("new_ica returns an unfitted eeg_ica object with documented defaults", {
  ica <- new_ica()

  expect_s3_class(ica, "eeg_ica")
  expect_equal(ica$current_fit, "unfitted")
  expect_null(ica$n_components)
  expect_null(ica$noise_cov)
  expect_null(ica$random_state)
  expect_equal(ica$method, "fastica")
  expect_equal(ica$max_iter, 1000L)
  expect_true(is.integer(ica$max_iter))
  expect_null(ica$ch_names)

  # fitted (trailing-underscore) fields stay NULL until fit_ica()
  for (fld in c("n_components_", "pre_whitener_", "pca_components_", "pca_mean_",
                "pca_explained_variance_", "mixing_matrix_", "unmixing_matrix_",
                "n_samples_", "n_iter_")) {
    expect_null(ica[[fld]], label = fld)
  }

  expect_equal(ica$exclude, integer(0))
  expect_equal(ica$labels_, list())
})

test_that("new_ica resolves default fit_params", {
  ica <- new_ica()
  expect_equal(ica$fit_params$algorithm, "parallel")
  expect_equal(ica$fit_params$fun, "logcosh")
  expect_equal(ica$fit_params$fun_args, list(alpha = 1))
  expect_equal(ica$fit_params$tol, 1e-04)
})

test_that("new_ica fit_params partial override keeps other defaults (modifyList)", {
  ica <- new_ica(fit_params = list(tol = 1e-6))
  expect_equal(ica$fit_params$tol, 1e-6)
  expect_equal(ica$fit_params$algorithm, "parallel")
  expect_equal(ica$fit_params$fun, "logcosh")
})

test_that("new_ica accepts an integer n_components > 1", {
  ica <- new_ica(n_components = 15)
  expect_equal(ica$n_components, 15)
})

test_that("new_ica rejects n_components == 1", {
  expect_error(new_ica(n_components = 1), "greater than 1")
})

test_that("new_ica accepts a float n_components strictly in (0, 1)", {
  ica <- new_ica(n_components = 0.95)
  expect_equal(ica$n_components, 0.95)
})

test_that("new_ica rejects a float n_components outside (0, 1)", {
  expect_error(new_ica(n_components = 1.5),  "between 0.0 and 1.0")
  expect_error(new_ica(n_components = -0.5), "between 0.0 and 1.0")
})

test_that("new_ica rejects non-scalar / non-numeric n_components", {
  expect_error(new_ica(n_components = c(1, 2)), "single number")
  expect_error(new_ica(n_components = "5"),     "single number")
})

test_that("new_ica rejects unimplemented methods", {
  expect_error(new_ica(method = "infomax"), "not yet implemented")
  expect_error(new_ica(method = c("fastica", "infomax")), "single character string")
})

test_that("new_ica resolves max_iter = 'auto' to 1000L", {
  ica <- new_ica(max_iter = "auto")
  expect_identical(ica$max_iter, 1000L)
})

test_that("new_ica accepts and coerces a custom positive max_iter", {
  ica <- new_ica(max_iter = 250)
  expect_identical(ica$max_iter, 250L)
})

test_that("new_ica rejects non-positive or non-integer max_iter", {
  expect_error(new_ica(max_iter = 0),     "positive integer")
  expect_error(new_ica(max_iter = -10),   "positive integer")
  expect_error(new_ica(max_iter = 100.5), "positive integer")
})

test_that("new_ica rejects invalid fit_params$algorithm / fun / tol", {
  expect_error(new_ica(fit_params = list(algorithm = "bogus")), "algorithm")
  expect_error(new_ica(fit_params = list(fun = "bogus")),       "fun")
  expect_error(new_ica(fit_params = list(tol = -1)),            "tol")
  expect_error(new_ica(fit_params = list(fun_args = "not a list")), "fun_args")
})

test_that("new_ica rejects a non-list fit_params", {
  expect_error(new_ica(fit_params = "not a list"), "named list")
})

test_that("new_ica stores noise_cov and random_state as supplied", {
  cov_mat <- diag(3)
  ica <- new_ica(noise_cov = cov_mat, random_state = 42)
  expect_equal(ica$noise_cov, cov_mat)
  expect_equal(ica$random_state, 42)
})

# ============================================================================
#              TEST SUITE 2: print.eeg_ica() - Print Method
# ============================================================================

test_that("print.eeg_ica shows 'not fitted' for a fresh object and returns invisibly", {
  ica <- new_ica(n_components = 10)
  output <- capture.output(result <- print(ica))

  expect_identical(result, ica)
  expect_true(any(grepl("ICA Object Summary", output)))
  expect_true(any(grepl("Not fitted", output)))
  expect_true(any(grepl("n_components:.*10", output)))
  expect_true(any(grepl("method:.*fastica", output)))
})

test_that("print.eeg_ica shows fit statistics once fitted", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  output <- capture.output(print(ica))

  expect_true(any(grepl("Fitted on:.*raw", output)))
  expect_true(any(grepl("n_components_:.*2", output)))
  expect_true(any(grepl("n_samples_:", output)))
  expect_true(any(grepl("n_iter_:", output)))
  expect_true(any(grepl("Excluded ICs:.*0", output)))
})

test_that("print.eeg_ica shows a Labels line when labels_ is populated", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  ica$labels_ <- list(eog = 1L, ecg = c(2L, 3L))

  output <- capture.output(print(ica))
  expect_true(any(grepl("Labels:.*eog = 1", output)))
})

# ============================================================================
#     TEST SUITE 3: .population_sd() - population standard deviation
# ============================================================================

test_that(".population_sd matches numpy's ddof = 0 convention", {
  x <- c(1, 2, 3, 4)
  expect_equal(.population_sd(x), sqrt(mean((x - mean(x))^2)))
  # Independent derivation: population sd = sample sd(n-1) * sqrt((n-1)/n)
  n <- length(x)
  expect_equal(.population_sd(x), sd(x) * sqrt((n - 1) / n))
})

test_that(".population_sd flattens matrices before computing (axis = None)", {
  m <- matrix(1:6, nrow = 2)
  expect_equal(.population_sd(m), .population_sd(as.vector(m)))
})

test_that(".population_sd of a constant vector is 0", {
  expect_equal(.population_sd(rep(5, 10)), 0)
})

# ============================================================================
#   TEST SUITE 4: .compute_pre_whitener() / .pre_whiten() - private helpers
# ============================================================================

test_that(".compute_pre_whitener pools one std across each channel *type*", {
  set.seed(1)
  data <- rbind(
    rnorm(500, sd = 10),  # eeg
    rnorm(500, sd = 20),  # eeg
    rnorm(500, sd = 2),   # external
    rnorm(500, sd = 4)    # external
  )
  channels      <- c("A", "B", "C", "D")
  channel_types <- c("eeg", "eeg", "external", "external")

  pw <- .compute_pre_whitener(data, channels, channel_types)

  expect_named(pw, channels)
  expect_equal(pw[["A"]], pw[["B"]])
  expect_equal(pw[["C"]], pw[["D"]])
  expect_false(isTRUE(all.equal(pw[["A"]], pw[["C"]])))

  # pooled value equals the population sd of the *flattened* group
  expect_equal(pw[["A"]], .population_sd(c(data[1, ], data[2, ])))
  expect_equal(pw[["C"]], .population_sd(c(data[3, ], data[4, ])))
})

test_that(".compute_pre_whitener errors on noise_cov (not yet implemented)", {
  data <- matrix(rnorm(20), nrow = 2)
  expect_error(
    .compute_pre_whitener(data, c("A", "B"), c("eeg", "eeg"), noise_cov = diag(2)),
    "not yet implemented"
  )
})

test_that(".compute_pre_whitener errors on mismatched channels/channel_types length", {
  data <- matrix(rnorm(20), nrow = 2)
  expect_error(.compute_pre_whitener(data, c("A"), c("eeg", "eeg")), "channels")
  expect_error(.compute_pre_whitener(data, c("A", "B"), c("eeg")),   "channel_types")
})

test_that(".pre_whiten divides each row by its pre-whitener entry", {
  data <- rbind(c(2, 4, 6), c(10, 20, 30))
  pw   <- c(2, 10)
  expect_equal(.pre_whiten(data, pw), rbind(c(1, 2, 3), c(1, 2, 3)))
})

# ============================================================================
#            TEST SUITE 5: .pca_whiten() - PCA whitening helper
# ============================================================================

test_that(".pca_whiten mean-centers, whitens to unit variance, and orders eigenvalues", {
  set.seed(42)
  # 4 channels with real cross-correlation (not independent), plenty of samples
  raw <- matrix(rnorm(4 * 500), nrow = 4)
  mixing <- matrix(rnorm(16), nrow = 4)
  data <- mixing %*% raw

  pca <- .pca_whiten(data)

  expect_equal(dim(pca$data), dim(data))
  expect_equal(pca$pca_mean_, rowMeans(data))

  # each whitened PC row has unit sample variance (ddof = 1, same as the
  # eigenvalues were computed with)
  row_vars <- apply(pca$data, 1, var)
  expect_equal(row_vars, rep(1, nrow(data)), tolerance = 1e-8)

  # decreasing, non-negative eigenvalues
  expect_true(all(diff(pca$pca_explained_variance_) <= 1e-10))
  expect_true(all(pca$pca_explained_variance_ >= 0))

  # components are orthonormal rows
  I_hat <- pca$pca_components_ %*% t(pca$pca_components_)
  expect_equal(I_hat, diag(nrow(data)), tolerance = 1e-8)
})

test_that(".pca_whiten applies a deterministic sign convention (largest-|entry| positive)", {
  set.seed(99)
  data <- matrix(rnorm(4 * 300), nrow = 4)
  pca  <- .pca_whiten(data)

  for (i in seq_len(nrow(pca$pca_components_))) {
    row  <- pca$pca_components_[i, ]
    peak <- which.max(abs(row))
    expect_gt(row[peak], 0)
  }
})

# ============================================================================
#          TEST SUITE 6: neg-entropy G-functions (.g_logcosh/.g_exp/.g_cube)
# ============================================================================

test_that(".axis_mean reduces vectors to a scalar and matrices row-wise", {
  expect_equal(.axis_mean(c(1, 2, 3)), 2)
  m <- matrix(1:6, nrow = 2)
  expect_equal(.axis_mean(m), rowMeans(m))
})

test_that(".g_logcosh on a vector returns the elementwise, un-reduced derivative", {
  x <- c(0, 1, -1)
  out <- .g_logcosh(x)
  expect_equal(out$gx, tanh(x))
  expect_equal(out$g_x, 1 - tanh(x)^2)     # un-reduced: same length as x
  expect_length(out$g_x, length(x))
})

test_that(".g_logcosh on a matrix reduces g_x to one value per row", {
  x   <- matrix(c(0, 1, -1, 2), nrow = 2)
  out <- .g_logcosh(x)
  expect_equal(out$gx, tanh(x))
  expect_equal(out$g_x, rowMeans(1 - tanh(x)^2))
  expect_length(out$g_x, nrow(x))
})

test_that(".g_logcosh honors a custom alpha", {
  out <- .g_logcosh(1, fun_args = list(alpha = 2))
  expect_equal(out$gx, tanh(2))
  expect_equal(out$g_x, 2 * (1 - tanh(2)^2))
})

test_that(".g_exp always reduces g_x regardless of vector/matrix input", {
  x <- c(0, 1, -1)
  out <- .g_exp(x)
  expect_equal(out$gx, x * exp(-(x^2) / 2))
  expect_equal(out$g_x, mean((1 - x^2) * exp(-(x^2) / 2)))
  expect_length(out$g_x, 1)

  m <- matrix(c(0, 1, -1, 2), nrow = 2)
  out_m <- .g_exp(m)
  expect_equal(out_m$g_x, rowMeans((1 - m^2) * exp(-(m^2) / 2)))
  expect_length(out_m$g_x, 2)
})

test_that(".g_cube matches x^3 / mean(3 x^2)", {
  x   <- c(1, 2, 3)
  out <- .g_cube(x)
  expect_equal(out$gx, x^3)
  expect_equal(out$g_x, mean(3 * x^2))
})

# ============================================================================
#     TEST SUITE 7: decorrelation helpers (.sym_decorrelation / .gs_decorrelation)
# ============================================================================

test_that(".sym_decorrelation returns an orthonormal matrix (W W^T == I)", {
  set.seed(5)
  W <- matrix(rnorm(9), nrow = 3)
  W_dec <- .sym_decorrelation(W)
  expect_equal(W_dec %*% t(W_dec), diag(3), tolerance = 1e-10)
})

test_that(".sym_decorrelation leaves an already-orthonormal matrix unchanged", {
  W <- diag(3)
  expect_equal(.sym_decorrelation(W), W, tolerance = 1e-10)
})

test_that(".gs_decorrelation with n_prior = 0 is a no-op", {
  w <- c(0.5, 0.5, 0.5)
  W <- matrix(0, 3, 3)
  expect_equal(.gs_decorrelation(w, W, 0), w)
})

test_that(".gs_decorrelation projects out overlap with prior orthonormal rows", {
  W <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0))
  w <- c(0.5, 0.5, 0.5)

  result_1 <- .gs_decorrelation(w, W, 1)
  expect_equal(result_1, c(0, 0.5, 0.5))

  result_2 <- .gs_decorrelation(w, W, 2)
  expect_equal(result_2, c(0, 0, 0.5))
  # orthogonal to both prior rows
  expect_equal(as.numeric(result_2 %*% W[1, ]), 0)
  expect_equal(as.numeric(result_2 %*% W[2, ]), 0)
})

# ============================================================================
#    TEST SUITE 8: FastICA rotation algorithms (.ica_par() / .ica_def())
# ============================================================================

# Small, deterministic PCA-whitened fixture shared by the rotation tests:
# 3 components x 800 samples of independent non-Gaussian data, already
# whitened via .pca_whiten() for consistency with how fit_ica() calls these.
make_whitened_fixture <- function(seed = 11) {
  set.seed(seed)
  raw <- rbind(runif(800, -1, 1), rexp(800) - 1, sign(rnorm(800)) * runif(800))
  .pca_whiten(raw)$data
}

test_that(".ica_par returns an orthonormal unmixing matrix", {
  X <- make_whitened_fixture()
  set.seed(1)
  w_init <- matrix(rnorm(9), 3, 3)
  out <- .ica_par(X, w_init, .g_logcosh, list(alpha = 1), max_iter = 500, tol = 1e-6)

  expect_equal(dim(out$W), c(3, 3))
  expect_equal(out$W %*% t(out$W), diag(3), tolerance = 1e-6)
  expect_true(out$n_iter >= 1 && out$n_iter <= 500)
})

test_that(".ica_par warns when it fails to converge within max_iter", {
  X <- make_whitened_fixture()
  set.seed(1)
  w_init <- matrix(rnorm(9), 3, 3)
  expect_warning(
    .ica_par(X, w_init, .g_logcosh, list(alpha = 1), max_iter = 1, tol = 0),
    "did not converge"
  )
})

test_that(".ica_def returns rows that are each unit-norm and mutually orthogonal", {
  X <- make_whitened_fixture()
  set.seed(1)
  w_init <- matrix(rnorm(9), 3, 3)
  out <- .ica_def(X, w_init, .g_logcosh, list(alpha = 1), max_iter = 500, tol = 1e-6)

  expect_equal(dim(out$W), c(3, 3))
  expect_equal(out$W %*% t(out$W), diag(3), tolerance = 1e-6)
  expect_true(out$n_iter >= 1 && out$n_iter <= 500)
})

test_that(".ica_par and .g_exp / .g_cube run without error and stay orthonormal", {
  X <- make_whitened_fixture()
  set.seed(2)
  w_init <- matrix(rnorm(9), 3, 3)

  out_exp  <- .ica_par(X, w_init, .g_exp,  list(), max_iter = 500, tol = 1e-6)
  out_cube <- .ica_par(X, w_init, .g_cube, list(), max_iter = 500, tol = 1e-6)

  expect_equal(out_exp$W  %*% t(out_exp$W),  diag(3), tolerance = 1e-6)
  expect_equal(out_cube$W %*% t(out_cube$W), diag(3), tolerance = 1e-6)
})

# ============================================================================
#        TEST SUITE 9: .resolve_n_components() - n_components resolution
# ============================================================================

test_that(".resolve_n_components with an integer just passes it through", {
  expect_equal(.resolve_n_components(3, c(10, 5, 3, 2, 1)), 3L)
})

test_that(".resolve_n_components errors when the integer is out of range", {
  expect_error(.resolve_n_components(1, c(10, 5, 3)),  "greater than 1")
  expect_error(.resolve_n_components(10, c(10, 5, 3)), "less than or equal")
})

test_that(".resolve_n_components with a float selects smallest n crossing the threshold", {
  ev <- c(0.6, 0.3, 0.05, 0.05)   # sums to 1.0; cumsum ~= .6, .9, .95, 1.0
  # thresholds deliberately kept clear of the cumulative-sum values
  # themselves, since e.g. 0.6 + 0.3 == 0.8999999999999999 in floating point
  # -- comparing against a threshold of exactly 0.9 would be a coin flip.
  expect_equal(.resolve_n_components(0.85, ev), 2L)
  expect_equal(.resolve_n_components(0.5,  ev), 1L)
  expect_equal(.resolve_n_components(0.96, ev), 4L)
})

test_that(".resolve_n_components with NULL uses the 0.999999 default threshold", {
  ev <- c(0.9999995, 0.0000005)   # first component alone already clears 0.999999
  expect_equal(.resolve_n_components(NULL, ev), 1L)
})

# ============================================================================
#                  TEST SUITE 10: fit_ica() - public orchestration
# ============================================================================

test_that("fit_ica rejects non-'eeg_ica' / non-'eeg' inputs", {
  fx <- make_ica_fixture()
  expect_error(fit_ica(list(), fx$eeg), "class 'eeg_ica'")
  expect_error(fit_ica(new_ica(), list()), "class 'eeg'")
})

test_that("fit_ica default picks exclude the status channel", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg))
  expect_false("Status" %in% ica$ch_names)
  expect_setequal(ica$ch_names, c("Fp1", "Fp2", "Cz", "EXG1"))
})

test_that("fit_ica accepts character picks and restricts ch_names accordingly", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(
    fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
            picks = c("Fp1", "Fp2", "Cz"))
  )
  expect_equal(ica$ch_names, c("Fp1", "Fp2", "Cz"))
})

test_that("fit_ica accepts numeric index picks", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(
    fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg, picks = 1:3)
  )
  expect_equal(ica$ch_names, fx$eeg$channels[1:3])
})

test_that("fit_ica errors on unknown channel name in picks", {
  fx <- make_ica_fixture()
  expect_error(
    fit_ica(new_ica(n_components = 2), fx$eeg, picks = c("Fp1", "BOGUS")),
    "not found in eeg\\$channels"
  )
})

test_that("fit_ica errors when fewer than 2 channels are selected", {
  fx <- make_ica_fixture()
  expect_error(fit_ica(new_ica(), fx$eeg, picks = "Fp1"), "at least 2 channels")
})

test_that("fit_ica populates all fitted (trailing-underscore) fields and flips current_fit", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))

  expect_equal(ica$current_fit, "raw")
  expect_equal(ica$n_components_, 3L)
  expect_length(ica$pre_whitener_, 3L)
  expect_equal(dim(ica$pca_components_), c(3L, 3L))
  expect_length(ica$pca_mean_, 3L)
  expect_length(ica$pca_explained_variance_, 3L)
  expect_equal(dim(ica$unmixing_matrix_), c(3L, 3L))
  expect_equal(dim(ica$mixing_matrix_), c(3L, 3L))
  expect_equal(ica$n_samples_, ncol(fx$eeg$data))
  expect_true(ica$n_iter_ >= 1)

  # mixing_matrix_ is the actual inverse of unmixing_matrix_
  expect_equal(ica$mixing_matrix_ %*% ica$unmixing_matrix_, diag(3), tolerance = 1e-8)
})

test_that("fit_ica with a fixed random_state is exactly reproducible", {
  fx <- make_ica_fixture()
  ica1 <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 42), fx$eeg))
  ica2 <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 42), fx$eeg))
  expect_identical(ica1$unmixing_matrix_, ica2$unmixing_matrix_)
})

test_that("fit_ica's random_state does not perturb the caller's global RNG stream", {
  fx <- make_ica_fixture()

  set.seed(123)
  expected_next <- { set.seed(123); runif(3) }

  set.seed(123)
  ica <- new_ica(n_components = 2, random_state = 999)
  invisible(suppressWarnings(fit_ica(ica, fx$eeg)))
  actual_next <- runif(3)

  expect_equal(actual_next, expected_next)
})

test_that("fit_ica (parallel, logcosh) recovers two known independent sources", {
  fx  <- make_source_separation_fixture()
  ica <- fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg)
  src_hat <- get_sources(ica, fx$eeg)

  cor_mat <- cor(t(fx$sources), t(src_hat))
  best_per_true_source <- apply(abs(cor_mat), 1, max)
  expect_true(all(best_per_true_source > 0.9))

  best_match_idx <- apply(abs(cor_mat), 1, which.max)
  expect_equal(length(unique(best_match_idx)), 2L)
})

test_that("fit_ica algorithm = 'deflation' also recovers the two known sources", {
  fx  <- make_source_separation_fixture()
  ica <- fit_ica(
    new_ica(n_components = 2, random_state = 1,
            fit_params = list(algorithm = "deflation")),
    fx$eeg
  )
  src_hat <- get_sources(ica, fx$eeg)
  cor_mat <- cor(t(fx$sources), t(src_hat))
  expect_true(all(apply(abs(cor_mat), 1, max) > 0.9))
})

# ============================================================================
#              TEST SUITE 11: get_sources() - component time-courses
# ============================================================================

test_that("get_sources errors on an unfitted ica or non-'eeg' input", {
  fx <- make_ica_fixture()
  expect_error(get_sources(new_ica(), fx$eeg), "has not been fit")
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  expect_error(get_sources(ica, list()), "class 'eeg'")
})

test_that("get_sources errors when eeg is missing a fitted channel", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  eeg_missing <- new_eeg(data = fx$eeg$data[1:2, ], channels = c("Fp1", "Fp2"),
                          sampling_rate = fx$eeg$sampling_rate)
  expect_error(get_sources(ica, eeg_missing), "missing channel")
})

test_that("get_sources returns n_components_ x n_samples with IC row names", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  src <- get_sources(ica, fx$eeg)

  expect_equal(dim(src), c(3L, ncol(fx$eeg$data)))
  expect_equal(rownames(src), c("IC1", "IC2", "IC3"))
})

test_that("get_sources replays the fitted transform exactly by hand", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))

  data     <- fx$eeg$data[match(ica$ch_names, fx$eeg$channels), , drop = FALSE]
  data_pw  <- .pre_whiten(data, ica$pre_whitener_)
  centered <- sweep(data_pw, 1, ica$pca_mean_, "-")
  whitened <- sweep(ica$pca_components_ %*% centered, 1,
                     sqrt(ica$pca_explained_variance_[1:3]), "/")
  expected <- ica$unmixing_matrix_ %*% whitened

  expect_equal(unname(get_sources(ica, fx$eeg)), unname(expected), tolerance = 1e-10)
})

# ============================================================================
#     TEST SUITE 12: get_component_topography() - scalp weights per IC
# ============================================================================

test_that("get_component_topography errors on unfitted ica or bad component index", {
  fx <- make_ica_fixture()
  expect_error(get_component_topography(new_ica(), 1), "has not been fit")

  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  expect_error(get_component_topography(ica, 0), "must be an integer in")
  expect_error(get_component_topography(ica, 3), "must be an integer in")
})

test_that("get_component_topography returns one named weight per fitted channel", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  topo <- get_component_topography(ica, 1)
  expect_named(topo, ica$ch_names)
  expect_length(topo, 3L)
})

test_that("summing topography_k outer source_k over all components reconstructs the raw data", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  sources <- get_sources(ica, fx$eeg)

  recon <- matrix(0, nrow = 3, ncol = ncol(sources))
  for (k in 1:3) {
    recon <- recon + outer(get_component_topography(ica, k), sources[k, ])
  }
  # get_component_topography() deliberately drops the constant pca_mean_
  # offset (a DC baseline, not a per-component contribution - see its
  # derivation comment above .population_sd() in R/ica1.R) - add it back in
  # original channel units, exactly as apply_ica()'s own remix step does.
  recon <- recon + ica$pca_mean_ * ica$pre_whitener_

  original <- fx$eeg$data[match(ica$ch_names, fx$eeg$channels), , drop = FALSE]
  expect_equal(unname(recon), unname(original), tolerance = 1e-6)
})

# ============================================================================
#           TEST SUITE 13: ica_component_summary() - inspection table
# ============================================================================

test_that("ica_component_summary errors on unfitted ica or non-'eeg' input", {
  fx <- make_ica_fixture()
  expect_error(ica_component_summary(new_ica(), fx$eeg), "has not been fit")
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  expect_error(ica_component_summary(ica, list()), "class 'eeg'")
})

test_that("ica_component_summary variance_explained sums to ~1 at full rank", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  summ <- ica_component_summary(ica, fx$eeg)

  expect_equal(summ$component, 1:3)
  expect_equal(sum(summ$variance_explained), 1, tolerance = 1e-6)
})

test_that("ica_component_summary adds one corr_<channel> column per external channel", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 4, random_state = 1), fx$eeg))
  summ <- ica_component_summary(ica, fx$eeg)

  expect_true("corr_EXG1" %in% names(summ))
  expect_true(all(abs(summ$corr_EXG1) <= 1 + 1e-8))
})

test_that("ica_component_summary omits corr_ columns when there are no external channels", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  summ <- ica_component_summary(ica, fx$eeg)
  expect_equal(names(summ), c("component", "variance_explained"))
})

# ============================================================================
#           TEST SUITE 14: set_exclude() - record exclusion choices
# ============================================================================

test_that("set_exclude errors on an unfitted ica", {
  expect_error(set_exclude(new_ica(), 1), "has not been fit")
})

test_that("set_exclude validates component indices are in range", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  expect_error(set_exclude(ica, 0), "must be integer indices")
  expect_error(set_exclude(ica, 4), "must be integer indices")
})

test_that("set_exclude stores a sorted, de-duplicated index set", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, c(3, 1, 1))
  expect_equal(ica$exclude, c(1L, 3L))
})

test_that("set_exclude replaces rather than appends", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, c(1, 2))
  ica <- set_exclude(ica, 3)
  expect_equal(ica$exclude, 3L)
})

test_that("set_exclude(ica, NULL) clears the exclude list", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, c(1, 2))
  ica <- set_exclude(ica, NULL)
  expect_equal(ica$exclude, integer(0))
})

# ============================================================================
#         TEST SUITE 15: apply_ica() - reconstruct cleaned EEG data
# ============================================================================

test_that("apply_ica errors on an unfitted ica or non-'eeg' input", {
  fx <- make_ica_fixture()
  expect_error(apply_ica(new_ica(), fx$eeg), "has not been fit")
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  expect_error(apply_ica(ica, list()), "class 'eeg'")
})

test_that("apply_ica warns when ica$exclude is empty", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  expect_warning(apply_ica(ica, fx$eeg), "exclude.*is empty")
})

test_that("apply_ica with nothing excluded round-trips to floating-point precision at full rank", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  cleaned <- suppressWarnings(apply_ica(ica, fx$eeg))

  idx <- match(ica$ch_names, fx$eeg$channels)
  expect_equal(cleaned$data[idx, ], fx$eeg$data[idx, ], tolerance = 1e-6)
})

test_that("apply_ica leaves channels outside ch_names untouched", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, 1)
  cleaned <- apply_ica(ica, fx$eeg)

  exg_idx <- which(fx$eeg$channels == "EXG1")
  status_idx <- which(fx$eeg$channels == "Status")
  expect_equal(cleaned$data[exg_idx, ], fx$eeg$data[exg_idx, ])
  expect_equal(cleaned$data[status_idx, ], fx$eeg$data[status_idx, ])
})

test_that("apply_ica returns class eeg and appends a descriptive history entry", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, c(1, 2))
  cleaned <- apply_ica(ica, fx$eeg)

  expect_s3_class(cleaned, "eeg")
  last_entry <- cleaned$preprocessing_history[[length(cleaned$preprocessing_history)]]
  expect_true(grepl("ICA applied: 2 component\\(s\\) removed \\(1, 2\\)", last_entry))
})

test_that("apply_ica removes an artifact-correlated component and reduces its footprint elsewhere", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 4, random_state = 1), fx$eeg))
  summ <- ica_component_summary(ica, fx$eeg)
  bad_comp <- summ$component[which.max(abs(summ$corr_EXG1))]

  fp1_idx  <- which(fx$eeg$channels == "Fp1")
  exg1_idx <- which(fx$eeg$channels == "EXG1")

  before_cor <- cor(fx$eeg$data[fp1_idx, ], fx$eeg$data[exg1_idx, ])

  ica_excl <- set_exclude(ica, bad_comp)
  cleaned  <- apply_ica(ica_excl, fx$eeg)
  after_cor <- cor(cleaned$data[fp1_idx, ], fx$eeg$data[exg1_idx, ])

  expect_gt(abs(before_cor), 0.3)
  expect_lt(abs(after_cor), abs(before_cor) - 0.2)
})

# ============================================================================
#     TEST SUITE 16: plot_ica_sources() - stacked component time-courses
# ============================================================================

test_that("plot_ica_sources validates its inputs", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))

  expect_error(plot_ica_sources(list(), fx$eeg), "class 'eeg_ica'")
  expect_error(plot_ica_sources(new_ica(), fx$eeg), "has not been fit")
  expect_error(plot_ica_sources(ica, list()), "class 'eeg'")
  with_null_device({
    expect_error(plot_ica_sources(ica, fx$eeg, components = 5), "components")
    expect_error(plot_ica_sources(ica, fx$eeg, tmin = 100, tmax = 50), "less than")
  })
})

test_that("plot_ica_sources draws without error and returns the expected invisible summary", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))

  with_null_device({
    result <- plot_ica_sources(ica, fx$eeg, components = c(1, 2))
    expect_equal(result$components, c(1, 2))
    expect_equal(result$n_samples, ncol(fx$eeg$data))
    expect_true(result$tmin_ms < result$tmax_ms)
  })
})

test_that("plot_ica_sources defaults to plotting every component", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  with_null_device({
    result <- plot_ica_sources(ica, fx$eeg)
    expect_equal(result$components, 1:3)
  })
})

# ============================================================================
#      TEST SUITE 17: plot_ica_topography() - scalp map for one component
# ============================================================================

test_that("plot_ica_topography draws a component's scalp map without error", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 3, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  eeg_m <- set_montage(fx$eeg, create_montage(c("Fp1", "Fp2", "Cz")))

  with_null_device({
    expect_no_error(plot_ica_topography(ica, eeg_m, component = 1))
  })
})

# ============================================================================
#       TEST SUITE 18: plot_ica_overlay() - before/after channel preview
# ============================================================================

test_that("plot_ica_overlay validates its inputs", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, 1)

  expect_error(plot_ica_overlay(list(), fx$eeg), "class 'eeg_ica'")
  expect_error(plot_ica_overlay(new_ica(), fx$eeg), "has not been fit")
  with_null_device({
    expect_error(plot_ica_overlay(ica, fx$eeg, channel = "Status"), "not one of the channels")
    expect_error(plot_ica_overlay(ica, fx$eeg, channel = "BOGUS"), "not found in 'eeg'")
  })
})

test_that("plot_ica_overlay draws without error and reports the excluded components", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, 1)

  with_null_device({
    result <- plot_ica_overlay(ica, fx$eeg, channel = "Fp1")
    expect_equal(result$channel, "Fp1")
    expect_equal(result$excluded, 1L)
    expect_equal(result$n_samples, ncol(fx$eeg$data))
  })
})

test_that("plot_ica_overlay defaults to the first fitted channel", {
  fx  <- make_ica_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg,
                                   picks = c("Fp1", "Fp2", "Cz")))
  ica <- set_exclude(ica, 1)

  with_null_device({
    result <- plot_ica_overlay(ica, fx$eeg)
    expect_equal(result$channel, ica$ch_names[1])
  })
})

# ============================================================================
#                             End of Test File
# ============================================================================
