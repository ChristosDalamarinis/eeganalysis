# ============================================================================
#                     Test File for regression.R  (Phase 1)
# ============================================================================
#
# Tests EOG regression on continuous data:
#   new_eog_regression(), fit_eog_regression(), apply_eog_regression(),
#   print.eeg_eog_regression()
#
# HOW THE MATHS IS CHECKED
# ------------------------
# The fixture plants KNOWN weights: every EEG channel = independent noise +
# a fixed DC offset + beta_true x (VEOG, HEOG). With the noise switched off,
# fit_eog_regression() must recover beta_true exactly; with noise on, to a
# small tolerance. On top of that:
#
#   - parity: the fit must equal a line-by-line transcription of MNE's
#     per-channel loop (mne/preprocessing/_regress.py, EOGRegression.fit)
#   - orthogonality: after fitting AND applying on the same data, the cleaned
#     channels have zero covariance with the mean-removed EOG. That is a
#     property of ordinary least squares, so any indexing or sign slip in
#     fit/apply would break it
#   - DC offsets, EOG channels and every non-target channel must come back
#     untouched (only the eye-related fluctuation is removed)
#
# The Phase 1 scope is continuous data only, so epoched input must fail with a
# clear message rather than silently doing the wrong thing.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

EEG_NAMES <- c("Fp1", "Fp2", "Fz", "Cz", "Pz", "Oz")
EOG_NAMES <- c("VEOG (EXG1)", "HEOG (EXG2)")

# Row layout of the fixture's data matrix:
#   1-6 EEG channels | 7 VEOG | 8 HEOG | 9 Status

# Continuous recording with planted eye artifacts. VEOG = blink-like bumps at
# irregular times, HEOG = slow oscillations (independent of the blinks). The
# EOG channels carry the sources exactly, so beta_true is the exact target of
# the fit. noise_sd = 0 gives a noise-free recording (exact recovery).
make_eog_fixture <- function(n_samples = 4000, sampling_rate = 250, seed = 1,
                             noise_sd = 1, reference = "Common Average") {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate

  veog <- rep(0, n_samples)
  centers <- sort(sample(seq(100, n_samples - 100), 25))
  for (b in centers) {
    win <- max(1, b - 60):min(n_samples, b + 60)
    veog[win] <- veog[win] + 120 * exp(-((win - b)^2) / (2 * 15^2))
  }
  heog <- 40 * sin(2 * pi * 0.30 * t + 1) + 25 * sin(2 * pi * 0.11 * t)

  beta_true <- rbind(
    Fp1 = c(0.60,  0.20),
    Fp2 = c(0.60, -0.20),
    Fz  = c(0.35,  0.05),
    Cz  = c(0.15,  0.00),
    Pz  = c(0.08,  0.00),
    Oz  = c(0.03,  0.00)
  )
  colnames(beta_true) <- EOG_NAMES
  dc <- c(30, -20, 10, 5, -15, 25)

  brain <- matrix(rnorm(length(EEG_NAMES) * n_samples, sd = noise_sd),
                  nrow = length(EEG_NAMES))
  eeg_part <- brain + dc + beta_true %*% rbind(veog, heog)

  status <- rep(0, n_samples)
  status[seq(200, n_samples, by = 400)] <- 65536

  data <- unname(rbind(eeg_part, veog, heog, status))
  eeg <- new_eeg(data = data,
                 channels = c(EEG_NAMES, EOG_NAMES, "Status"),
                 sampling_rate = sampling_rate,
                 reference = reference)
  list(eeg = eeg, beta_true = beta_true, dc = dc)
}

demean_rows <- function(M) M - rowMeans(M)

# MNE's EOGRegression.fit() loop, transcribed one-to-one (_regress.py:185-199):
# ref_data mean-removed, cov_ref = ref_data @ ref_data.T, then per channel
# coef[pi] = solve(cov_ref, ref_data @ cov_data.T).
naive_mne_fit <- function(data, tgt_idx, art_idx) {
  ref_data <- demean_rows(data[art_idx, , drop = FALSE])
  cov_ref  <- ref_data %*% t(ref_data)
  coef <- matrix(0, nrow = length(tgt_idx), ncol = length(art_idx))
  for (pi in seq_along(tgt_idx)) {
    this_data <- data[tgt_idx[pi], ]
    cov_data  <- this_data - mean(this_data)
    coef[pi, ] <- t(solve(cov_ref, ref_data %*% cov_data))
  }
  coef
}

permute_channels <- function(eeg, perm) {
  out <- eeg
  out$data <- eeg$data[perm, , drop = FALSE]
  out$channels <- eeg$channels[perm]
  out$channel_types <- eeg$channel_types[perm]
  out
}

drop_channel <- function(eeg, name) {
  keep <- eeg$channels != name
  out <- eeg
  out$data <- eeg$data[keep, , drop = FALSE]
  out$channels <- eeg$channels[keep]
  out$channel_types <- eeg$channel_types[keep]
  out
}

add_external_channel <- function(eeg, name, values) {
  out <- eeg
  out$data <- rbind(eeg$data, values)
  out$channels <- c(eeg$channels, name)
  out$channel_types <- c(eeg$channel_types, "external")
  out
}

# ============================================================================
# SUITE 1 - FIXTURE SANITY (guards every assumption the rest relies on)
# ============================================================================

test_that("fixture channels are typed as intended", {
  fx <- make_eog_fixture()
  expect_equal(unname(fx$eeg$channel_types),
               c(rep("eeg", 6), "external", "external", "status"))
  expect_equal(fx$eeg$reference, "Common Average")
})

# ============================================================================
# SUITE 2 - new_eog_regression()
# ============================================================================

test_that("new_eog_regression() builds a model and sets the dimnames", {
  m <- new_eog_regression(matrix(c(1, 2, 3, 4), nrow = 2),
                          ch_names = c("A", "B"),
                          ch_names_artifact = c("E1", "E2"),
                          reference = "Common Average", n_samples = 100)
  expect_s3_class(m, "eeg_eog_regression")
  expect_equal(dimnames(m$coef_), list(c("A", "B"), c("E1", "E2")))
  expect_equal(m$reference, "Common Average")
  expect_equal(m$fit_on, "continuous")
  expect_equal(m$n_samples_, 100)
})

test_that("new_eog_regression() defaults: unknown reference and n_samples", {
  m <- new_eog_regression(matrix(0.5, 1, 1), "A", "E1")
  expect_true(is.na(m$reference))
  expect_true(is.na(m$n_samples_))
})

test_that("new_eog_regression() accepts fit_on = 'epochs' and rejects others", {
  ok <- matrix(0.5, 1, 1)
  expect_equal(new_eog_regression(ok, "A", "E1", fit_on = "epochs")$fit_on,
               "epochs")
  expect_error(new_eog_regression(ok, "A", "E1", fit_on = "raw"))
})

test_that("new_eog_regression() rejects bad input", {
  ok <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  expect_error(new_eog_regression(1:4, c("A", "B"), c("E1", "E2")),
               "numeric matrix")
  bad <- ok; bad[1, 1] <- NA
  expect_error(new_eog_regression(bad, c("A", "B"), c("E1", "E2")),
               "NA, NaN or Inf")
  expect_error(new_eog_regression(ok, c("A", "B", "C"), c("E1", "E2")),
               "coef' has")
  expect_error(new_eog_regression(ok, c("A", "A"), c("E1", "E2")),
               "unique, non-missing")
  expect_error(new_eog_regression(ok, c("A", "B"), c("E1", NA)),
               "unique, non-missing")
  expect_error(new_eog_regression(ok, c("A", "B"), c("B", "E2")),
               "appear in both")
  expect_error(new_eog_regression(ok, c("A", "B"), c("E1", "E2"),
                                  reference = c("x", "y")),
               "'reference' must")
  expect_error(new_eog_regression(ok, c("A", "B"), c("E1", "E2"),
                                  n_samples = "many"),
               "'n_samples' must")
})

# ============================================================================
# SUITE 3 - fit_eog_regression(): the numbers
# ============================================================================

test_that("noise-free data: planted weights are recovered exactly", {
  fx <- make_eog_fixture(noise_sd = 0)
  m <- fit_eog_regression(fx$eeg)
  expect_equal(m$coef_, fx$beta_true, tolerance = 1e-8)
})

test_that("noisy data: planted weights are recovered to a small tolerance", {
  fx <- make_eog_fixture(noise_sd = 1)
  m <- fit_eog_regression(fx$eeg)
  expect_lt(max(abs(m$coef_ - fx$beta_true)), 0.01)
})

test_that("fit equals a transcription of MNE's per-channel loop", {
  fx <- make_eog_fixture(noise_sd = 1, seed = 5)
  m <- fit_eog_regression(fx$eeg)
  expected <- naive_mne_fit(fx$eeg$data, tgt_idx = 1:6, art_idx = 7:8)
  expect_equal(unname(m$coef_), expected, tolerance = 1e-9)
})

test_that("a channel's weights do not depend on which other channels are fit", {
  fx <- make_eog_fixture()
  full <- fit_eog_regression(fx$eeg)
  part <- fit_eog_regression(fx$eeg, picks = c("Cz", "Fp1"))
  expect_equal(part$coef_, full$coef_[c("Cz", "Fp1"), , drop = FALSE])
})

test_that("a constant (DC) offset on a target does not change its weights", {
  fx <- make_eog_fixture(seed = 3)
  shifted <- fx$eeg
  shifted$data[1:6, ] <- shifted$data[1:6, ] + 5000
  expect_equal(fit_eog_regression(shifted)$coef_,
               fit_eog_regression(fx$eeg)$coef_, tolerance = 1e-7)
})

test_that("the model records what it was fit on", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  expect_s3_class(m, "eeg_eog_regression")
  expect_equal(m$reference, "Common Average")
  expect_equal(m$fit_on, "continuous")
  expect_equal(m$n_samples_, 4000)
  expect_equal(dim(m$coef_), c(6, 2))
})

# ============================================================================
# SUITE 4 - fit_eog_regression(): channel selection
# ============================================================================

test_that("default targets are the good EEG channels; EOG and Status are not", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  expect_equal(m$ch_names, EEG_NAMES)
  expect_equal(m$ch_names_artifact, EOG_NAMES)
})

test_that("bads are left out of the default targets", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- "Pz"
  m <- fit_eog_regression(eeg)
  expect_equal(m$ch_names, setdiff(EEG_NAMES, "Pz"))
})

test_that("an explicit picks is taken as-is, bad channels included", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- "Pz"
  m <- fit_eog_regression(eeg, picks = c("Pz", "Cz"))
  expect_equal(m$ch_names, c("Pz", "Cz"))
})

test_that("picks and picks_artifact accept channel indices", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg, picks = 1:3, picks_artifact = 7:8)
  expect_equal(m$ch_names, EEG_NAMES[1:3])
  expect_equal(m$ch_names_artifact, EOG_NAMES)
})

test_that("an EOG channel marked bad is skipped when auto-detecting", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- "HEOG (EXG2)"
  m <- fit_eog_regression(eeg)
  expect_equal(m$ch_names_artifact, "VEOG (EXG1)")
  expect_equal(ncol(m$coef_), 1)
})

test_that("if every auto-detected EOG channel is bad, it errors; explicit works", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- EOG_NAMES
  expect_error(fit_eog_regression(eeg), "marked bad")
  m <- fit_eog_regression(eeg, picks_artifact = EOG_NAMES)
  expect_equal(m$ch_names_artifact, EOG_NAMES)
})

test_that("no EOG channel in the recording gives an actionable error", {
  eeg <- make_eog_fixture()$eeg
  eeg <- drop_channel(drop_channel(eeg, EOG_NAMES[1]), EOG_NAMES[2])
  expect_error(fit_eog_regression(eeg), "no EOG channel found")
  expect_error(fit_eog_regression(eeg), "picks_artifact")
})

test_that("explicit picks_artifact can name any channel", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg, picks_artifact = "VEOG (EXG1)")
  expect_equal(m$ch_names_artifact, "VEOG (EXG1)")
  expect_equal(ncol(m$coef_), 1)
})

test_that("unknown, overlapping, empty or duplicated picks are refused", {
  eeg <- make_eog_fixture()$eeg
  expect_error(fit_eog_regression(eeg, picks = "Nope"), "not found in eeg\\$channels")
  expect_error(fit_eog_regression(eeg, picks_artifact = "Nope"), "not found")
  expect_error(fit_eog_regression(eeg, picks = c("Fp1", "VEOG (EXG1)")),
               "cannot be regressed on itself")
  expect_error(fit_eog_regression(eeg, picks = 99), "valid channel indices")
  expect_error(fit_eog_regression(eeg, picks = 1.5), "valid channel indices")
  expect_error(fit_eog_regression(eeg, picks = character(0)), "is empty")
  expect_error(fit_eog_regression(eeg, picks = c("Fp1", "Fp1")), "more than once")
  expect_error(fit_eog_regression(eeg, picks = TRUE), "channel names or")
})

test_that("with no good EEG channel left, it errors", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- EEG_NAMES
  expect_error(fit_eog_regression(eeg), "no target channels")
})

# ============================================================================
# SUITE 5 - fit_eog_regression(): guards and errors
# ============================================================================

test_that("data with no reference applied is refused", {
  for (ref in c("original", "Biosemi CMS/DRL")) {
    eeg <- make_eog_fixture(reference = ref)$eeg
    expect_error(fit_eog_regression(eeg), "no EEG reference")
    expect_error(fit_eog_regression(eeg), "eeg_rereference")
  }
  eeg <- make_eog_fixture()$eeg
  eeg$reference <- NULL
  expect_error(fit_eog_regression(eeg), "no EEG reference")
})

test_that("other reference schemes are accepted", {
  for (ref in c("Common Average", "M1+M2")) {
    eeg <- make_eog_fixture(reference = ref)$eeg
    expect_s3_class(fit_eog_regression(eeg), "eeg_eog_regression")
  }
})

test_that("wrong input types give clear errors", {
  expect_error(fit_eog_regression(1:10), "class 'eeg'")
  epochs <- structure(list(), class = "eeg_epochs")
  expect_error(fit_eog_regression(epochs), "not supported yet")
})

test_that("non-finite values in a target or an EOG channel are refused", {
  eeg <- make_eog_fixture()$eeg
  eeg_t <- eeg
  eeg_t$data[4, 10] <- NaN                         # Cz
  expect_error(fit_eog_regression(eeg_t), "target channel.*Cz")
  eeg_e <- eeg
  eeg_e$data[7, 10] <- Inf                         # VEOG
  expect_error(fit_eog_regression(eeg_e), "EOG channel.*VEOG")
})

test_that("exactly collinear EOG channels give an informative error", {
  fx <- make_eog_fixture()
  eeg <- add_external_channel(fx$eeg, "VEOG2 (EXG3)", fx$eeg$data[7, ])
  expect_error(fit_eog_regression(eeg), "collinear")
})

test_that("nearly collinear EOG channels give a warning, not an error", {
  fx <- make_eog_fixture()
  veog <- fx$eeg$data[7, ]
  set.seed(11)
  near_copy <- veog + rnorm(length(veog), sd = stats::sd(veog) * 1e-6)
  eeg <- add_external_channel(fx$eeg, "VEOG2 (EXG3)", near_copy)
  expect_warning(m <- fit_eog_regression(eeg), "nearly collinear")
  expect_s3_class(m, "eeg_eog_regression")
})

# ============================================================================
# SUITE 6 - apply_eog_regression(): what it changes and what it must not
# ============================================================================

test_that("noise-free: cleaned channels are flat at their original level", {
  fx <- make_eog_fixture(noise_sd = 0)
  out <- apply_eog_regression(fit_eog_regression(fx$eeg), fx$eeg)
  n <- ncol(fx$eeg$data)
  expected <- matrix(rowMeans(fx$eeg$data[1:6, ]), nrow = 6, ncol = n)
  expect_equal(out$data[1:6, ], expected, tolerance = 1e-6)
})

test_that("after an in-sample fit the cleaned channels are orthogonal to the EOG", {
  fx <- make_eog_fixture(noise_sd = 1)
  out <- apply_eog_regression(fit_eog_regression(fx$eeg), fx$eeg)
  R <- demean_rows(fx$eeg$data[7:8, ])
  before <- tcrossprod(R, demean_rows(fx$eeg$data[1:6, ]))
  after  <- tcrossprod(R, demean_rows(out$data[1:6, ]))
  expect_gt(max(abs(before)), 1e4)                 # the artifact was really there
  expect_lt(max(abs(after)), 1e-8 * max(abs(before)))
})

test_that("each channel's overall level (DC offset) is preserved", {
  fx <- make_eog_fixture(noise_sd = 1)
  out <- apply_eog_regression(fit_eog_regression(fx$eeg), fx$eeg)
  expect_equal(rowMeans(out$data[1:6, ]), rowMeans(fx$eeg$data[1:6, ]))
})

test_that("EOG, Status and channels outside the model come back untouched", {
  eeg <- make_eog_fixture()$eeg
  eeg$bads <- "Pz"                                 # left out of the fit
  m <- fit_eog_regression(eeg)
  out <- apply_eog_regression(m, eeg)
  expect_identical(out$data[7:9, ], eeg$data[7:9, ])          # EOG + Status
  expect_identical(out$data[5, ], eeg$data[5, ])              # Pz (bad)
  expect_false(isTRUE(all.equal(out$data[1, ], eeg$data[1, ])))  # Fp1 changed
})

test_that("apply leaves the input alone and every other field intact", {
  fx <- make_eog_fixture()
  before <- fx$eeg
  out <- apply_eog_regression(fit_eog_regression(fx$eeg), fx$eeg)
  expect_identical(fx$eeg, before)
  keep <- setdiff(names(out), c("data", "preprocessing_history"))
  expect_identical(out[keep], fx$eeg[keep])
  expect_s3_class(out, "eeg")
  expect_equal(dim(out$data), dim(fx$eeg$data))
})

test_that("apply logs the step in preprocessing_history", {
  fx <- make_eog_fixture()
  out <- apply_eog_regression(fit_eog_regression(fx$eeg), fx$eeg)
  expect_length(out$preprocessing_history,
                length(fx$eeg$preprocessing_history) + 1)
  entry <- out$preprocessing_history[[length(out$preprocessing_history)]]
  expect_match(entry, "EOG regression applied: 6 channel\\(s\\)")
  expect_match(entry, "VEOG \\(EXG1\\), HEOG \\(EXG2\\)")
  expect_match(entry, "Common Average")
})

test_that("a model fit on one recording cleans another (same channel names)", {
  a <- make_eog_fixture(seed = 1)
  b <- make_eog_fixture(seed = 2)
  m <- fit_eog_regression(a$eeg)
  out_b <- apply_eog_regression(m, b$eeg)
  r_before <- cor(b$eeg$data[1, ], b$eeg$data[7, ])      # Fp1 vs VEOG
  r_after  <- cor(out_b$data[1, ], out_b$data[7, ])
  expect_gt(abs(r_before), 0.9)
  expect_lt(abs(r_after), 0.1)
})

test_that("long recordings are handled block by block without seams", {
  # 100,000 samples = three full 32,768-sample blocks plus a partial last one
  fx <- make_eog_fixture(n_samples = 100000, seed = 4)
  m <- fit_eog_regression(fx$eeg)
  out <- apply_eog_regression(m, fx$eeg)
  R <- demean_rows(fx$eeg$data[7:8, ])
  expected <- fx$eeg$data[1:6, ] - m$coef_ %*% R
  expect_equal(unname(out$data[1:6, ]), unname(expected))
})

test_that("a NaN in the EOG stays local to its own sample", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  eeg <- fx$eeg
  eeg$data[7, 100] <- NaN
  out <- apply_eog_regression(m, eeg)
  expect_true(all(is.nan(out$data[1:6, 100])))
  expect_true(all(is.finite(out$data[1:6, -100])))
})

# ============================================================================
# SUITE 7 - apply_eog_regression(): matching, guards and errors
# ============================================================================

test_that("channels are matched by name, so channel order does not matter", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  perm <- c(9, 7, 3, 1, 8, 5, 2, 6, 4)
  out_orig <- apply_eog_regression(m, fx$eeg)
  out_perm <- apply_eog_regression(m, permute_channels(fx$eeg, perm))
  expect_equal(out_perm$channels, fx$eeg$channels[perm])
  expect_equal(out_perm$data, out_orig$data[perm, ])
})

test_that("a missing target or EOG channel is named in the error", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  expect_error(apply_eog_regression(m, drop_channel(fx$eeg, "Cz")),
               "missing target channel.*Cz")
  expect_error(apply_eog_regression(m, drop_channel(fx$eeg, "VEOG (EXG1)")),
               "missing EOG channel.*VEOG")
})

test_that("a different reference than at fit time warns", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  other <- fx$eeg
  other$reference <- "M1+M2"
  expect_warning(apply_eog_regression(m, other), "depend on the reference")
})

test_that("an unknown model reference does not warn", {
  fx <- make_eog_fixture()
  m0 <- fit_eog_regression(fx$eeg)
  m <- new_eog_regression(m0$coef_, m0$ch_names, m0$ch_names_artifact)
  expect_warning(apply_eog_regression(m, fx$eeg), NA)
})

test_that("applying to data with no reference applied is refused", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  cms <- fx$eeg
  cms$reference <- "Biosemi CMS/DRL"
  expect_error(apply_eog_regression(m, cms), "no EEG reference")
})

test_that("wrong input types give clear errors", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  expect_error(apply_eog_regression(list(), fx$eeg), "eeg_eog_regression")
  expect_error(apply_eog_regression(m, 1:10), "class 'eeg'")
  epochs <- structure(list(), class = "eeg_epochs")
  expect_error(apply_eog_regression(m, epochs), "not supported yet")
})

test_that("weights brought in via new_eog_regression() work like fitted ones", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  external <- new_eog_regression(unname(m$coef_), m$ch_names,
                                 m$ch_names_artifact,
                                 reference = "Common Average")
  expect_equal(apply_eog_regression(external, fx$eeg)$data,
               apply_eog_regression(m, fx$eeg)$data)
})

# ============================================================================
# SUITE 8 - print() and saving
# ============================================================================

test_that("print shows the essentials and returns the model invisibly", {
  fx <- make_eog_fixture()
  m <- fit_eog_regression(fx$eeg)
  expect_output(print(m), "EOG Regression Model")
  expect_output(print(m), "VEOG \\(EXG1\\), HEOG \\(EXG2\\)")
  expect_output(print(m), "Common Average")
  expect_output(print(m), "continuous data \\(4000 samples\\)")
  res <- NULL
  capture.output(res <- withVisible(print(m)))
  expect_false(res$visible)
  expect_identical(res$value, m)
})

test_that("print abbreviates a long target list", {
  m <- new_eog_regression(matrix(0.1, nrow = 8, ncol = 1),
                          ch_names = paste0("Ch", 1:8),
                          ch_names_artifact = "EOG")
  expect_output(print(m), "\\+3 more")
  expect_output(print(m), "unknown")
})

test_that("a model survives saveRDS() / readRDS() unchanged", {
  m <- fit_eog_regression(make_eog_fixture()$eeg)
  f <- tempfile(fileext = ".rds")
  saveRDS(m, f)
  reloaded <- readRDS(f)
  unlink(f)
  expect_identical(reloaded, m)
})
