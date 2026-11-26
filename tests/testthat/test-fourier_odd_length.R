# ============================================================================
#         Test File for fourier.R spectra of an ODD number of samples
# ============================================================================
#
# The FFT of n samples at sr Hz gives bins spaced sr / n apart, and the
# one-sided spectrum keeps bins 0 .. floor(n / 2). Every kept bin except 0 Hz
# has a mirror twin in the dropped half, so its power is doubled to put the
# twin's share back. The only other bin without a twin is the Nyquist bin
# (exactly sr / 2), and it exists only when n is EVEN.
#
# fourier.R used to assume n is always even:
#   1. Frequency labels: it spread the labels evenly from 0 to sr / 2, so for
#      odd n every label was too high by n / (n - 1) (513 samples at 512 Hz:
#      last label 256 Hz instead of 255.5 Hz; 5 samples at 10 Hz: a 2 Hz
#      wave was labelled 2.5 Hz).
#   2. Doubling: it never doubled the last kept bin, so for odd n the value
#      at the top frequency came out at half its true size.
#
# Odd n is common: an epoch from -0.2 to 0.8 s at 512 Hz has 513 samples,
# because epoch_eeg() keeps both ends.
#
# test-fourier.R is skipped as a whole for now, so these tests live in their
# own file to make sure they run.
#
# Test suites:
#   1. Frequency labels (eeg_fft, eeg_psd_welch, eeg_multitaper)
#   2. The top bin is doubled when n is odd, not when n is even
#
# Author: Christos Dalamarinis
# Date: Nov 2026
# ============================================================================

library(testthat)
library(eeganalysis)

SR_O <- 512   # sampling rate (Hz) for the epoch-sized tests

# A small 'eeg' object with one channel
odd_eeg <- function(x, sr = SR_O) {
  structure(list(data = matrix(x, nrow = 1L), channels = "Cz",
                 sampling_rate = sr), class = "eeg")
}

# A small 'eeg_epochs' object: one channel, every column of m is one epoch
odd_epochs <- function(m, sr = SR_O) {
  structure(list(data = array(m, dim = c(1L, nrow(m), ncol(m))),
                 channels = "Cz", sampling_rate = sr), class = "eeg_epochs")
}

# The true frequency of each kept bin: k * sr / n, k = 0 .. floor(n / 2)
true_freqs <- function(n, sr) (0:floor(n / 2)) * sr / n


# ============================================================================
# TEST SUITE 1: Frequency labels
# ============================================================================

test_that("eeg_fft(): a 2 Hz wave in 5 samples at 10 Hz is labelled 2 Hz", {
  x    <- cos(2 * pi * 2 * (0:4) / 10)
  spec <- eeg_fft(odd_eeg(x, sr = 10), window = "none", verbose = FALSE)
  expect_equal(spec$frequencies, c(0, 2, 4), tolerance = 1e-12)
  expect_equal(spec$frequencies[which.max(spec$power[1L, ])], 2,
               tolerance = 1e-12)
})

test_that("all three spectra label a 513-sample epoch sr / 513 apart", {
  set.seed(3L)
  ep    <- odd_epochs(matrix(rnorm(513 * 4), nrow = 513))
  truth <- true_freqs(513, SR_O)        # 0, 0.998, ..., 255.501 Hz

  expect_equal(eeg_fft(ep, verbose = FALSE)$frequencies, truth,
               tolerance = 1e-12)
  expect_equal(eeg_psd_welch(ep, window_length = 513 / SR_O,
                             verbose = FALSE)$frequencies, truth,
               tolerance = 1e-12)
  expect_equal(eeg_multitaper(ep, verbose = FALSE)$frequencies, truth,
               tolerance = 1e-12)
})

test_that("odd segments and odd n_fft of continuous data are labelled right", {
  set.seed(4L)
  eeg <- odd_eeg(rnorm(10 * SR_O))
  # eeg_psd_welch(): 2 s segments (1024 samples) with an odd 1025-point FFT
  expect_equal(eeg_psd_welch(eeg, n_fft = 1025, verbose = FALSE)$frequencies,
               true_freqs(1025, SR_O), tolerance = 1e-12)
  # eeg_multitaper(): 513-sample segments
  expect_equal(eeg_multitaper(eeg, segment_length = 513 / SR_O,
                              verbose = FALSE)$frequencies,
               true_freqs(513, SR_O), tolerance = 1e-12)
})

test_that("even lengths keep ending exactly at sr / 2", {
  set.seed(5L)
  ep <- odd_epochs(matrix(rnorm(512 * 4), nrow = 512))
  expect_equal(eeg_fft(ep, verbose = FALSE)$frequencies,
               true_freqs(512, SR_O), tolerance = 1e-12)
  expect_equal(tail(eeg_multitaper(ep, verbose = FALSE)$frequencies, 1L), 256)
})


# ============================================================================
# TEST SUITE 2: The top bin is doubled when n is odd, not when n is even
# ============================================================================

test_that("eeg_fft(): a wave on the top bin keeps its full height", {
  # odd: 5 samples at 10 Hz, a 4 Hz wave of height 1 sits on the last bin,
  # and that bin has a mirror twin -> must be doubled to give 1
  odd  <- eeg_fft(odd_eeg(cos(2 * pi * 4 * (0:4) / 10), sr = 10),
                  window = "none", verbose = FALSE)
  # even: 4 samples at 10 Hz, a 5 Hz wave of height 1 is the Nyquist bin,
  # which has no twin -> must NOT be doubled (stays 1, not 2)
  even <- eeg_fft(odd_eeg(cos(2 * pi * 5 * (0:3) / 10), sr = 10),
                  window = "none", verbose = FALSE)
  expect_equal(tail(odd$amplitude[1L, ], 1L),  1, tolerance = 1e-12)
  expect_equal(tail(even$amplitude[1L, ], 1L), 1, tolerance = 1e-12)
})

test_that("eeg_psd_welch(): the PSD of one odd segment adds up to mean(x^2)", {
  # One segment, no window, no mean removal: the PSD added up over all
  # frequencies (times the true step sr / n) must give back mean(x^2) exactly.
  set.seed(6L)
  for (n in c(513L, 512L)) {
    x    <- rnorm(n, sd = 5)
    spec <- eeg_psd_welch(odd_eeg(x), window_length = n / SR_O,
                          window = "none", remove_dc = FALSE, verbose = FALSE)
    expect_equal(sum(spec$power[1L, ]) * SR_O / n, mean(x^2),
                 tolerance = 1e-10, info = paste("n =", n))
  }
})

test_that("eeg_multitaper(): the top bin of odd epochs is not halved", {
  # White noise has a flat PSD, so the top bin must sit at the same level as
  # the rest (it used to come out at about half).
  set.seed(7L)
  ep   <- odd_epochs(matrix(rnorm(513 * 100, sd = 5), nrow = 513))
  spec <- eeg_multitaper(ep, per_epoch = TRUE, verbose = FALSE)
  p    <- spec$power[1L, ]
  expect_equal(tail(p, 1L) / mean(p[20:236]), 1, tolerance = 0.2)
})
