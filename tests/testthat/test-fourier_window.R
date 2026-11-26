# ============================================================================
#              Test File for the spectral windows of fourier.R
# ============================================================================
#
# Before the FFT, eeg_fft() and eeg_psd_welch() multiply each piece of signal
# by a window (a hill that fades the edges to zero). The FFT treats the piece
# as if it repeats forever, end joined to start, so the window must be
# PERIODIC: exactly one cycle of a cosine over the n samples, so that two
# copies side by side join into one smooth wave. The periodic window of n
# samples is the symmetric window of n + 1 samples with its last zero dropped.
#
# .make_window() used to build SYMMETRIC windows (dividing by n - 1, the form
# used for filter design). These tests check three facts that only a
# periodic window has:
#
#   1. Its FFT has nothing except its cosine terms: Hann and Hamming use
#      bins 0 and +-1 only, Blackman bins 0, +-1 and +-2. Every other bin is 0.
#   2. Two Hann windows shifted by half a window add up to exactly 1.
#   3. A sine that sits exactly on a frequency bin spreads into its two
#      neighbouring bins only, each with exactly 1/4 of the peak power
#      (Hann), and into no other bin.
#
# test-fourier.R is skipped as a whole for now, so these tests live in their
# own file to make sure they run.
#
# Test suites:
#   1. .make_window(): the periodic shapes
#   2. eeg_psd_welch() and eeg_fft(): a sine on a bin
#
# Author: Christos Dalamarinis
# Date: Nov 2026
# ============================================================================

library(testthat)
library(eeganalysis)

SR_W <- 512   # sampling rate (Hz)

# Largest FFT bin outside the cosine terms, divided by n (0 for a periodic window)
off_cosine <- function(w, n_terms) {
  n    <- length(w)
  keep <- c(1L, 1L + seq_len(n_terms), n + 1L - seq_len(n_terms))
  max(Mod(fft(w))[-keep]) / n
}

# A small 'eeg' object with one channel
win_eeg <- function(x) {
  structure(list(data = matrix(x, nrow = 1L), channels = "Cz",
                 sampling_rate = SR_W), class = "eeg")
}

# 20 uV sine at 10 Hz; with 2 s pieces (0.5 Hz steps) 10 Hz sits on a bin
sine_10hz <- function(n_samp) {
  20 * sin(2 * pi * 10 * (seq_len(n_samp) - 1) / SR_W)
}

# Power in the two bins next to 10 Hz, and the largest power farther away,
# both relative to the power at 10 Hz
neighbours <- function(power, freqs) {
  k    <- which(abs(freqs - 10) < 1e-9)
  step <- freqs[2L] - freqs[1L]
  far  <- abs(freqs - 10) > 1.5 * step
  list(below = power[k - 1L] / power[k],
       above = power[k + 1L] / power[k],
       far   = max(power[far]) / power[k])
}


# ============================================================================
# TEST SUITE 1: .make_window() - the periodic shapes
# ============================================================================

test_that(".make_window() gives the periodic Hann of 4 samples", {
  # symmetric Hann of 5 samples is 0, 0.5, 1, 0.5, 0: drop the last zero
  expect_equal(eeganalysis:::.make_window(4L, "hann"), c(0, 0.5, 1, 0.5),
               tolerance = 1e-12)
})

test_that(".make_window(type = 'none') is all ones", {
  expect_identical(eeganalysis:::.make_window(7L, "none"), rep(1.0, 7L))
})

test_that("the FFT of each window has only its cosine terms (even and odd n)", {
  for (n in c(16L, 15L)) {
    expect_lt(off_cosine(eeganalysis:::.make_window(n, "hann"), 1L), 1e-12)
    expect_lt(off_cosine(eeganalysis:::.make_window(n, "hamming"), 1L), 1e-12)
    expect_lt(off_cosine(eeganalysis:::.make_window(n, "blackman"), 2L), 1e-12)
  }
})

test_that("two Hann windows half a window apart add up to exactly 1", {
  w <- eeganalysis:::.make_window(16L, "hann")
  expect_equal(w[1:8] + w[9:16], rep(1, 8L), tolerance = 1e-12)
})


# ============================================================================
# TEST SUITE 2: eeg_psd_welch() and eeg_fft() - a sine on a bin
# ============================================================================

test_that("eeg_psd_welch(): a sine on a bin spreads into its 2 neighbours only", {
  spec <- eeg_psd_welch(win_eeg(sine_10hz(10 * SR_W)), window_length = 2,
                        window = "hann", verbose = FALSE)
  nb <- neighbours(spec$power[1L, ], spec$frequencies)
  expect_equal(nb$below, 0.25, tolerance = 1e-9)
  expect_equal(nb$above, 0.25, tolerance = 1e-9)
  expect_lt(nb$far, 1e-12)
})

test_that("eeg_fft(): a sine on a bin spreads into its 2 neighbours only", {
  spec <- eeg_fft(win_eeg(sine_10hz(2 * SR_W)), window = "hann",
                  verbose = FALSE)
  nb <- neighbours(spec$power[1L, ], spec$frequencies)
  expect_equal(nb$below, 0.25, tolerance = 1e-9)
  expect_equal(nb$above, 0.25, tolerance = 1e-9)
  expect_lt(nb$far, 1e-12)
})
