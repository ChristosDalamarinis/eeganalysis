# ============================================================================
#                 Test File for the units of fourier.R spectra
# ============================================================================
#
# A power spectral density (PSD) is in uV^2/Hz. Two facts let us check that
# without any other program:
#
#   1. A sine wave of amplitude A has a mean power of A^2 / 2. Adding up the
#      PSD over all frequencies (each value times the width of one frequency
#      step) must give back that number: 20 uV -> 200 uV^2.
#   2. White noise with standard deviation s, sampled at sr Hz, has a flat
#      PSD of 2 * s^2 / sr: s = 5 uV at 512 Hz -> 0.0977 uV^2/Hz.
#
# eeg_multitaper() used to divide by the segment length one time too many,
# so its PSD was n times too small (1024x for 2 s at 512 Hz). These tests
# catch that, for every kind of input and for eeg_band_power() too.
#
# test-fourier.R is skipped as a whole for now, so these tests live in their
# own file to make sure they run.
#
# Test suites:
#   1. Sine wave: the PSD adds up to A^2 / 2 (eeg, eeg_epochs, eeg_evoked)
#   2. White noise: the PSD level is 2 * s^2 / sr, same as Welch
#   3. eeg_band_power(method = "multitaper")
#
# Author: Christos Dalamarinis
# Date: Nov 2026
# ============================================================================

library(testthat)
library(eeganalysis)

SR_U <- 512   # sampling rate (Hz)
AMP  <- 20    # sine amplitude (uV): mean power AMP^2 / 2 = 200 uV^2

# A small 'eeg' object with one channel
unit_eeg <- function(x) {
  structure(list(data = matrix(x, nrow = 1L), channels = "Cz",
                 sampling_rate = SR_U), class = "eeg")
}

# A small 'eeg_epochs' object: one channel, every column of m is one epoch
unit_epochs <- function(m) {
  structure(list(data = array(m, dim = c(1L, nrow(m), ncol(m))),
                 channels = "Cz", sampling_rate = SR_U), class = "eeg_epochs")
}

# A small 'eeg_evoked' object: one channel, one condition
unit_evoked <- function(x) {
  structure(list(evoked = list(cond1 = matrix(x, nrow = 1L)),
                 conditions = "cond1", channels = "Cz",
                 sampling_rate = SR_U), class = "eeg_evoked")
}

# PSD added up over all frequencies: sum(PSD) * frequency step
total_power <- function(psd, freqs) sum(psd) * (freqs[2L] - freqs[1L])

sine_wave <- function(n_samp, freq = 10) {
  AMP * sin(2 * pi * freq * (seq_len(n_samp) - 1) / SR_U)
}


# ============================================================================
# TEST SUITE 1: Sine wave - the multitaper PSD adds up to A^2 / 2
# ============================================================================

test_that("eeg_multitaper() on 'eeg': a 20 uV sine adds up to 200 uV^2", {
  spec <- eeg_multitaper(unit_eeg(sine_wave(10 * SR_U)), verbose = FALSE)
  expect_equal(total_power(spec$power[1L, ], spec$frequencies), 200,
               tolerance = 0.01)
})

test_that("eeg_multitaper() on 'eeg_epochs': a 20 uV sine adds up to 200 uV^2", {
  m <- cbind(sine_wave(1024), sine_wave(1024))
  for (pe in c(FALSE, TRUE)) {
    spec <- eeg_multitaper(unit_epochs(m), per_epoch = pe, verbose = FALSE)
    expect_equal(total_power(spec$power[1L, ], spec$frequencies), 200,
                 tolerance = 0.01, info = paste("per_epoch =", pe))
  }
})

test_that("eeg_multitaper() on 'eeg_evoked': a 20 uV sine adds up to 200 uV^2", {
  spec <- eeg_multitaper(unit_evoked(sine_wave(1024)), verbose = FALSE)
  expect_equal(total_power(spec$power$cond1[1L, ], spec$frequencies), 200,
               tolerance = 0.01)
})


# ============================================================================
# TEST SUITE 2: White noise - the PSD level is 2 * s^2 / sr
# ============================================================================

test_that("eeg_multitaper() and eeg_psd_welch() give the white-noise level", {
  set.seed(1L)
  eeg   <- unit_eeg(rnorm(60 * SR_U, sd = 5))   # 60 s, s = 5 uV
  truth <- 2 * 5^2 / SR_U                       # 0.0977 uV^2/Hz

  level <- function(spec) {
    keep <- spec$frequencies >= 5 & spec$frequencies <= 200
    mean(spec$power[1L, keep])
  }
  mt <- level(eeg_multitaper(eeg, verbose = FALSE))
  w  <- level(eeg_psd_welch(eeg, verbose = FALSE))

  expect_equal(mt, truth, tolerance = 0.05)
  expect_equal(w,  truth, tolerance = 0.05)
  expect_equal(mt, w,     tolerance = 0.05)
})


# ============================================================================
# TEST SUITE 3: eeg_band_power(method = "multitaper")
# ============================================================================

test_that("eeg_band_power(method = 'multitaper') holds the sine's 200 uV^2", {
  eeg  <- unit_eeg(sine_wave(10 * SR_U))
  band <- list(around_10 = c(5, 15))
  mt <- eeg_band_power(eeg, bands = band, method = "multitaper",
                       verbose = FALSE)
  w  <- eeg_band_power(eeg, bands = band, method = "welch", verbose = FALSE)

  expect_equal(mt$around_10, 200, tolerance = 0.02)
  expect_equal(mt$around_10, w$around_10, tolerance = 0.02)
})
