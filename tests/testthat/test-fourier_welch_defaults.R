# ============================================================================
#           Test File for the default segments of eeg_psd_welch()
# ============================================================================
#
# eeg_psd_welch() used to cut the signal into 2 s segments with 50 % overlap
# by default. The defaults are now:
#   - segments of 2048 samples, or the whole signal when it is shorter
#     (2048 is a sample count: 4 s at 512 Hz, 2 s at 1024 Hz)
#   - no overlap: segments side by side, no sample used twice
#   - Hann window (unchanged)
#
# The main test rebuilds the estimate by hand: cut the signal into 2048-sample
# pieces side by side, remove each piece's mean, multiply by a Hann window,
# take |FFT|^2, scale to uV^2/Hz, double every bin that has a mirror twin, and
# average the pieces.
#
# test-fourier.R is skipped as a whole for now, so these tests live in their
# own file to make sure they run.
#
# Test suites:
#   1. Default segment length and overlap
#   2. Same numbers as a hand-made Welch estimate
#   3. Argument checks
#
# Author: Christos Dalamarinis
# Date: Nov 2026
# ============================================================================

library(testthat)
library(eeganalysis)

SR_D <- 512   # sampling rate (Hz)

# A small 'eeg' object with one channel
welch_eeg <- function(x, sr = SR_D) {
  structure(list(data = matrix(x, nrow = 1L), channels = "Cz",
                 sampling_rate = sr),
            class = "eeg")
}

# Welch by hand: Hann pieces of seg_len samples side by side, mean removed,
# |FFT|^2 / (sr * sum(w^2)), mirror bins doubled (seg_len is even here, so
# not the first bin and not the last one), averaged over the pieces.
hand_welch <- function(x, sr, seg_len) {
  n_seg  <- length(x) %/% seg_len
  w      <- 0.5 * (1 - cos(2 * pi * (0:(seg_len - 1)) / seg_len))
  n_half <- seg_len %/% 2 + 1
  pieces <- sapply(seq_len(n_seg), function(j) {
    seg <- x[(j - 1) * seg_len + seq_len(seg_len)]
    seg <- (seg - mean(seg)) * w
    p   <- Mod(fft(seg)[seq_len(n_half)])^2 / (sr * sum(w^2))
    p[2:(n_half - 1)] <- 2 * p[2:(n_half - 1)]
    p
  })
  rowMeans(pieces)
}


# ============================================================================
# Suite 1: default segment length and overlap
# ============================================================================

test_that("default segments are 2048 samples, side by side, Hann", {
  set.seed(1L)
  s <- eeg_psd_welch(welch_eeg(rnorm(10 * SR_D)), verbose = FALSE)
  expect_equal(s$n_fft, 2048)
  expect_equal(s$window_length, 2048 / SR_D)          # 4 s at 512 Hz
  expect_equal(s$overlap, 0)
  expect_equal(s$n_segments, 2)                       # 5120 %/% 2048
  expect_equal(s$window_type, "hann")
  expect_equal(s$frequencies, (0:1024) * SR_D / 2048, tolerance = 1e-12)
})

test_that("2048 is a sample count, not a time", {
  set.seed(2L)
  s <- eeg_psd_welch(welch_eeg(rnorm(20 * 256), sr = 256), verbose = FALSE)
  expect_equal(s$n_fft, 2048)
  expect_equal(s$window_length, 8)                    # 8 s at 256 Hz
})

test_that("a signal shorter than 2048 samples is one segment, no warning", {
  set.seed(3L)
  expect_no_warning(
    s <- eeg_psd_welch(welch_eeg(rnorm(1000)), verbose = FALSE)
  )
  expect_equal(s$n_fft, 1000)
  expect_equal(s$n_segments, 1)
})

test_that("each 513-sample epoch is one segment by default", {
  set.seed(4L)
  ep <- structure(list(data = array(rnorm(513 * 3), c(1L, 513L, 3L)),
                       channels = "Cz", sampling_rate = SR_D),
                  class = "eeg_epochs")
  expect_no_warning(s <- eeg_psd_welch(ep, per_epoch = TRUE, verbose = FALSE))
  expect_equal(s$n_fft, 513)
  expect_equal(s$frequencies, (0:256) * SR_D / 513, tolerance = 1e-12)
})


# ============================================================================
# Suite 2: same numbers as a hand-made Welch estimate
# ============================================================================

test_that("the default estimate equals Welch done by hand", {
  set.seed(5L)
  x <- 20 * sin(2 * pi * 10 * (0:(10 * SR_D - 1)) / SR_D) + rnorm(10 * SR_D, sd = 5)
  s <- eeg_psd_welch(welch_eeg(x), verbose = FALSE)
  expect_equal(unname(s$power[1L, ]), hand_welch(x, SR_D, 2048),
               tolerance = 1e-12)
})

test_that("the old settings are still available by name", {
  set.seed(6L)
  x <- rnorm(10 * SR_D)
  s <- eeg_psd_welch(welch_eeg(x), window_length = 2, overlap = 0.5,
                     verbose = FALSE)
  expect_equal(s$n_fft, 1024)
  expect_equal(s$n_segments, 9)                       # starts every 512 samples
})


# ============================================================================
# Suite 3: argument checks
# ============================================================================

test_that("window_length must be NULL or a positive number", {
  eeg <- welch_eeg(rnorm(4096))
  expect_no_error(eeg_psd_welch(eeg, window_length = NULL, verbose = FALSE))
  expect_error(eeg_psd_welch(eeg, window_length = -1, verbose = FALSE),
               "window_length")
  expect_error(eeg_psd_welch(eeg, window_length = "2", verbose = FALSE),
               "window_length")
})

test_that("a given window_length longer than the signal still warns", {
  expect_warning(eeg_psd_welch(welch_eeg(rnorm(1000)), window_length = 10,
                               verbose = FALSE),
                 "exceeds signal length")
})
