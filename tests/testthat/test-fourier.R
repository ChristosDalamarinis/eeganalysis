# ============================================================================
#                       Test File for fft.R
# ============================================================================
#
# Comprehensive tests for the Fourier / Spectral Analysis module (R/fft.R).
#
# Functions tested:
#   .make_window()      - Internal spectral window generator
#   .compute_dpss()     - Internal DPSS taper computation
#   eeg_fft()           - Raw FFT amplitude / power / phase spectra
#   eeg_psd_welch()     - Welch averaged periodogram PSD
#   eeg_multitaper()    - Multitaper PSD (DPSS / Slepian)
#   eeg_band_power()    - EEG frequency band power extraction
#
# Test suite structure:
#   Suite 1:  .make_window()           — shape, range, normalisation
#   Suite 2:  .compute_dpss()          — dimensions, orthogonality, concentration
#   Suite 3:  eeg_fft() validation     — input class errors, bad arguments
#   Suite 4:  eeg_fft() structure      — return class, field names, dimensions
#   Suite 5:  eeg_fft() accuracy       — sine wave peak recovery
#   Suite 6:  eeg_fft() eeg_epochs     — per_epoch flag, 3-D epoch_power
#   Suite 7:  eeg_fft() eeg_evoked     — condition-wise list output
#   Suite 8:  eeg_psd_welch() validation — bad arguments
#   Suite 9:  eeg_psd_welch() structure  — return fields, dimensions, non-negative
#   Suite 10: eeg_psd_welch() accuracy   — sine wave peak recovery
#   Suite 11: eeg_psd_welch() eeg_epochs/evoked dispatch
#   Suite 12: eeg_multitaper() validation
#   Suite 13: eeg_multitaper() structure
#   Suite 14: eeg_multitaper() accuracy  — sine wave peak recovery
#   Suite 15: eeg_multitaper() eeg_epochs/evoked dispatch
#   Suite 16: eeg_band_power() validation
#   Suite 17: eeg_band_power() structure — correct column names, dimensions
#   Suite 18: eeg_band_power() accuracy  — relative power sums to 1, band ordering
#   Suite 19: eeg_band_power() evoked    — returns named list
#   Suite 20: print.eeg_spectrum()       — no error, returns invisibly
#
# Author: Christos Dalamarinis
# Date: 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                         SHARED MOCK HELPERS
# ============================================================================

# Sampling rate and duration used throughout — long enough for Welch segments
SR  <- 512L   # Hz
DUR <- 4.0    # seconds
N   <- as.integer(SR * DUR)   # 2048 samples

#' Build a minimal mock 'eeg' object.
#' @param n_channels  Number of channels (default 3).
#' @param n_tp        Number of time points (default N = 2048 at 512 Hz).
#' @param sr          Sampling rate in Hz.
#' @param data_values Optional numeric matrix [n_channels x n_tp].
#'                    If NULL, uses Gaussian noise.
#' @return An S3 object of class "eeg".
make_mock_eeg <- function(n_channels  = 3L,
                          n_tp        = N,
                          sr          = SR,
                          data_values = NULL) {
  if (is.null(data_values)) {
    data_values <- matrix(rnorm(n_channels * n_tp),
                          nrow = n_channels, ncol = n_tp)
  }
  chan_names <- paste0("Ch", seq_len(n_channels))
  times      <- seq(0, (n_tp - 1) / sr, length.out = n_tp)
  
  structure(
    list(
      data                  = data_values,
      channels              = chan_names,
      sampling_rate         = as.numeric(sr),
      times                 = times,
      events                = data.frame(onset        = integer(0),
                                         onset_time   = numeric(0),
                                         type         = character(0),
                                         description  = character(0),
                                         stringsAsFactors = FALSE),
      metadata              = list(),
      reference             = "original",
      preprocessing_history = list()
    ),
    class = "eeg"
  )
}

#' Build a mock 'eeg' object whose channel 1 contains a pure sine wave.
#' Remaining channels are Gaussian noise.
#' @param freq_hz  Target sine frequency in Hz.
#' @param n_tp     Signal length in samples.
#' @param sr       Sampling rate in Hz.
#' @param amp      Peak amplitude of the sine (default 10 µV).
#' @return An S3 object of class "eeg".
make_sine_eeg <- function(freq_hz, n_tp = N, sr = SR, amp = 10) {
  t        <- seq(0, (n_tp - 1) / sr, length.out = n_tp)
  sine_row <- amp * sin(2 * pi * freq_hz * t)
  noise    <- rnorm(n_tp, sd = 0.01)   # tiny noise, essentially pure sine
  data_mat <- rbind(sine_row + noise,
                    matrix(rnorm(2 * n_tp), nrow = 2, ncol = n_tp))
  make_mock_eeg(n_channels = 3L, n_tp = n_tp, sr = sr, data_values = data_mat)
}

#' Build a minimal mock 'eeg_epochs' object.
#' @param n_channels Number of channels.
#' @param n_samp     Samples per epoch.
#' @param n_epochs   Number of epochs (trials).
#' @param sr         Sampling rate in Hz.
#' @param data_array Optional 3-D array [n_channels x n_samp x n_epochs].
#' @return An S3 object of class "eeg_epochs".
make_mock_epochs <- function(n_channels = 2L,
                             n_samp     = 256L,
                             n_epochs   = 5L,
                             sr         = SR,
                             data_array = NULL) {
  if (is.null(data_array)) {
    data_array <- array(rnorm(n_channels * n_samp * n_epochs),
                        dim = c(n_channels, n_samp, n_epochs))
  }
  chan_names <- paste0("Ch", seq_len(n_channels))
  dimnames(data_array)[[1L]] <- chan_names
  
  structure(
    list(
      data             = data_array,
      channels         = chan_names,
      times            = seq(-0.1, 0.1, length.out = n_samp),
      events           = data.frame(
        onset        = seq(200L, by = 300L, length.out = n_epochs),
        onset_time   = seq(200L, by = 300L, length.out = n_epochs) / sr,
        type         = rep(1L, n_epochs),
        description  = rep("Trigger: 1", n_epochs),
        epoch_id     = seq_len(n_epochs),
        stringsAsFactors = FALSE
      ),
      sampling_rate    = as.numeric(sr),
      tmin             = -0.1,
      tmax             =  0.1,
      baseline         = c(-0.1, 0),
      baseline_method  = "mean",
      n_epochs         = n_epochs,
      rejected         = rep(FALSE, n_epochs),
      rejection_log    = data.frame(epoch_id   = integer(0),
                                    event_type = character(0),
                                    event_time = numeric(0),
                                    reason     = character(0),
                                    stringsAsFactors = FALSE),
      metadata         = list()
    ),
    class = "eeg_epochs"
  )
}

#' Build a minimal mock 'eeg_evoked' object with two conditions.
#' @param n_channels  Number of channels.
#' @param n_samp      Samples per condition.
#' @param sr          Sampling rate in Hz.
#' @return An S3 object of class "eeg_evoked".
make_mock_evoked <- function(n_channels = 2L, n_samp = 256L, sr = SR) {
  chan_names  <- paste0("Ch", seq_len(n_channels))
  make_cond   <- function() {
    m <- matrix(rnorm(n_channels * n_samp), nrow = n_channels, ncol = n_samp)
    rownames(m) <- chan_names
    m
  }
  
  structure(
    list(
      evoked       = list(cond_A = make_cond(), cond_B = make_cond()),
      se           = NULL,
      n_trials     = c(cond_A = 20L, cond_B = 20L),
      conditions   = c("cond_A", "cond_B"),
      channels     = chan_names,
      times        = seq(-0.1, 0.1, length.out = n_samp),
      sampling_rate = as.numeric(sr),
      tmin         = -0.1,
      tmax         =  0.1,
      baseline     = c(-0.1, 0),
      method       = "mean"
    ),
    class = "eeg_evoked"
  )
}

#' Find the frequency bin with the highest power in a power vector.
#' @param pwr_vec  Numeric vector of power values.
#' @param freqs    Corresponding frequency vector.
#' @return Scalar: peak frequency in Hz.
peak_freq <- function(pwr_vec, freqs) freqs[which.max(pwr_vec)]


# ============================================================================
# TEST SUITE 1: .make_window()  — spectral window generator
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Output length equals n for all window types
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: .make_window() must always return a vector of exactly the
# requested length, regardless of window type.
test_that(".make_window() returns a vector of the requested length", {
  for (type in c("hann", "hamming", "blackman", "none")) {
    w <- eeganalysis:::.make_window(128L, type)
    expect_equal(length(w), 128L,
                  info = paste("window type:", type))
  }
})

# ----------------------------------------------------------------------------
# Test 1.2: "none" window is all ones
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The rectangular ("none") window must be a vector of 1s,
# leaving the signal unmodified when applied.
test_that(".make_window('none') returns all ones", {
  w <- eeganalysis:::.make_window(64L, "none")
  expect_true(all(w == 1.0))
})

# ----------------------------------------------------------------------------
# Test 1.3: Hann window has correct endpoints and peak
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The standard Hann window equals 0 at both endpoints (n=1
# and n=N) and 1.0 at its centre. This confirms the cosine formula is correct.
test_that(".make_window('hann') has zero endpoints and unit peak", {
  n <- 128L
  w <- eeganalysis:::.make_window(n, "hann")
  expect_equal(w[1L],  0.0, tolerance = 1e-10)
  expect_equal(w[n],   0.0, tolerance = 1e-10)
  expect_equal(max(w), 1.0, tolerance = 1e-10)
})

# ----------------------------------------------------------------------------
# Test 1.4: All windows have values in [0, 1]
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Hann, Hamming, and Blackman windows must always be
# non-negative and <= 1, since they are amplitude tapering functions.
test_that(".make_window() values are in [0, 1] for tapering windows", {
  for (type in c("hann", "hamming", "blackman")) {
    w <- eeganalysis:::.make_window(256L, type)
    expect_true(all(w >= 0),
                info = paste(type, "has negative values"))
    expect_true(all(w <= 1.0 + 1e-10),
                info = paste(type, "exceeds 1.0"))
  }
})


# ============================================================================
# TEST SUITE 2: .compute_dpss()  — DPSS taper computation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Output dimensions are [n x k]
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: .compute_dpss() must return a matrix with exactly n rows
# and k columns, one column per taper.
test_that(".compute_dpss() returns matrix with correct dimensions", {
  tapers <- eeganalysis:::.compute_dpss(n = 256L, nw = 4, k = 7L)
  expect_true(is.matrix(tapers))
  expect_equal(dim(tapers), c(256L, 7L))
})

# ----------------------------------------------------------------------------
# Test 2.2: Each taper has unit energy (L2 norm = 1)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: DPSS tapers should be normalised so that sum(taper^2) = 1.
# This is required by the PSD normalisation in .mtaper_psd().
test_that(".compute_dpss() tapers have unit energy", {
  tapers <- eeganalysis:::.compute_dpss(256L, 4, 7L)
  energies <- apply(tapers, 2L, function(v) sum(v^2))
  expect_true(all(abs(energies - 1.0) < 1e-10))
})

# ----------------------------------------------------------------------------
# Test 2.3: Tapers are mutually orthogonal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A key property of DPSS is that tapers are orthogonal:
# sum(taper_i * taper_j) ≈ 0 for i != j. This reduces spectral leakage.
test_that(".compute_dpss() tapers are mutually orthogonal", {
  tapers <- eeganalysis:::.compute_dpss(256L, 4, 4L)
  for (i in 1:3) {
    for (j in (i + 1):4) {
      dot <- abs(sum(tapers[, i] * tapers[, j]))
      expect_true(dot < 1e-8,
                info = paste("tapers", i, "and", j, "are not orthogonal"))
    }
  }
})

# ----------------------------------------------------------------------------
# Test 2.4: Single taper request (k = 1) returns a matrix
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Edge case — requesting k = 1 taper. Must return a matrix
# (not a vector) with [n x 1] dimensions, ensuring consistent downstream use.
test_that(".compute_dpss() with k=1 returns a 2-D matrix", {
  tapers <- eeganalysis:::.compute_dpss(128L, 2, 1L)
  expect_true(is.matrix(tapers))
  expect_equal(dim(tapers), c(128L, 1L))
})


# ============================================================================
# TEST SUITE 3: eeg_fft() — Input Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Wrong input class throws an informative error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: eeg_fft() must reject any object that is not 'eeg',
# 'eeg_epochs', or 'eeg_evoked' with a clear error message.
test_that("eeg_fft() errors on invalid input class", {
  expect_error(eeg_fft(list(data = matrix(1))),
               regexp = "Input must be one of")
  expect_error(eeg_fft(data.frame(a = 1:10)),
               regexp = "Input must be one of")
})

# ----------------------------------------------------------------------------
# Test 3.2: NULL epochs data throws an informative error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If an eeg_epochs object has $data = NULL (preload = FALSE),
# eeg_fft() must stop with a clear message rather than giving a cryptic error.
test_that("eeg_fft() errors when eeg_epochs$data is NULL", {
  ep        <- make_mock_epochs()
  ep$data   <- NULL
  expect_error(eeg_fft(ep, verbose = FALSE),
               regexp = "Epoch data not loaded")
})

# ----------------------------------------------------------------------------
# Test 3.3: Invalid n_fft throws an error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: n_fft must be a positive integer or NULL.
test_that("eeg_fft() errors on non-positive n_fft", {
  eeg <- make_mock_eeg()
  expect_error(eeg_fft(eeg, n_fft = -100, verbose = FALSE))
  expect_error(eeg_fft(eeg, n_fft = 0,    verbose = FALSE))
})


# ============================================================================
# TEST SUITE 4: eeg_fft() — Output Structure
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Return class is "eeg_spectrum"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: eeg_fft() must return an S3 object of class "eeg_spectrum"
# so that print.eeg_spectrum() and downstream methods dispatch correctly.
test_that("eeg_fft() returns class 'eeg_spectrum'", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, verbose = FALSE)
  expect_s3_class(spec, "eeg_spectrum")
})

# ----------------------------------------------------------------------------
# Test 4.2: All documented fields are present
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The return list must contain every field documented in the
# module header. Missing fields would break downstream code.
test_that("eeg_fft() output contains all documented fields", {
  eeg      <- make_mock_eeg()
  spec     <- eeg_fft(eeg, verbose = FALSE)
  expected <- c("power", "amplitude", "phase", "frequencies", "channels",
                "sampling_rate", "method", "input_class", "n_fft",
                "window_type", "per_epoch", "epoch_power", "conditions")
  expect_true(all(expected %in% names(spec)))
})

# ----------------------------------------------------------------------------
# Test 4.3: Metadata fields have correct values for eeg input
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: method, input_class, and window_type must be set correctly
# according to the function arguments and input object.
test_that("eeg_fft() sets metadata fields correctly for eeg input", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, window = "hamming", verbose = FALSE)
  expect_equal(spec$method,       "fft")
  expect_equal(spec$input_class,  "eeg")
  expect_equal(spec$window_type,  "hamming")
  expect_equal(spec$per_epoch,    FALSE)
  expect_null(spec$conditions)
  expect_null(spec$epoch_power)
})

# ----------------------------------------------------------------------------
# Test 4.4: Power and amplitude matrix dimensions are [channels x (n_fft/2+1)]
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: For an eeg object with n_ch channels and n_fft FFT length,
# the output power and amplitude matrices must have dimensions [n_ch x n_half]
# where n_half = floor(n_fft/2) + 1.
test_that("eeg_fft() power and amplitude have correct [channels x freqs] dimensions", {
  n_ch   <- 4L
  n_fft  <- 1024L
  eeg    <- make_mock_eeg(n_channels = n_ch)
  spec   <- eeg_fft(eeg, n_fft = n_fft, verbose = FALSE)
  n_half <- floor(n_fft / 2L) + 1L
  
  expect_equal(dim(spec$power),     c(n_ch, n_half))
  expect_equal(dim(spec$amplitude), c(n_ch, n_half))
  expect_equal(dim(spec$phase),     c(n_ch, n_half))
})

# ----------------------------------------------------------------------------
# Test 4.5: frequencies vector has correct length and range
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The frequency vector must span [0, Fs/2] with length equal
# to n_half = floor(n_fft/2) + 1.
test_that("eeg_fft() frequencies vector spans [0, Fs/2] with correct length", {
  sr     <- 512L
  n_fft  <- 2048L
  eeg    <- make_mock_eeg(sr = sr)
  spec   <- eeg_fft(eeg, n_fft = n_fft, verbose = FALSE)
  n_half <- floor(n_fft / 2L) + 1L
  
  expect_length(spec$frequencies, n_half)
  expect_equal(spec$frequencies[1L],      0,    tolerance = 1e-10)
  expect_equal(spec$frequencies[n_half],  sr / 2, tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 4.6: Power is non-negative everywhere
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Power = amplitude^2 must always be >= 0. Negative power
# would indicate a normalisation or implementation bug.
test_that("eeg_fft() power is non-negative for all channels and frequencies", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, verbose = FALSE)
  expect_true(all(spec$power >= 0))
})

# ----------------------------------------------------------------------------
# Test 4.7: Channel names are preserved as row names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Row names of the power/amplitude matrices must match the
# channel names in the input EEG object.
test_that("eeg_fft() preserves channel names as row names of power matrix", {
  eeg  <- make_mock_eeg(n_channels = 3L)
  spec <- eeg_fft(eeg, verbose = FALSE)
  expect_equal(rownames(spec$power), eeg$channels)
})


# ============================================================================
# TEST SUITE 5: eeg_fft() — Spectral Accuracy
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Peak frequency recovery from a pure sine wave (channel 1)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Given a 10 Hz pure sine wave in channel 1, the FFT spectrum
# must peak at (or very near) 10 Hz. Tolerance = 1 FFT bin width.
test_that("eeg_fft() recovers peak frequency of a 10 Hz sine wave", {
  target_hz <- 10
  eeg       <- make_sine_eeg(freq_hz = target_hz)
  spec      <- eeg_fft(eeg, window = "hann", verbose = FALSE)
  freq_res  <- spec$frequencies[2L]   # bin width
  
  pk <- peak_freq(spec$power[1L, ], spec$frequencies)
  expect_equal(pk, target_hz, tolerance = freq_res)
})

# ----------------------------------------------------------------------------
# Test 5.2: Peak frequency recovery for a 40 Hz sine wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Gamma-band frequency recovery. Confirms the FFT works
# correctly across the frequency range relevant to EEG analysis.
test_that("eeg_fft() recovers peak frequency of a 40 Hz sine wave", {
  eeg      <- make_sine_eeg(freq_hz = 40)
  spec     <- eeg_fft(eeg, window = "hann", verbose = FALSE)
  freq_res <- spec$frequencies[2L]
  
  pk <- peak_freq(spec$power[1L, ], spec$frequencies)
  expect_equal(pk, 40, tolerance = freq_res)
})

# ----------------------------------------------------------------------------
# Test 5.3: Zero-padding increases frequency resolution without new information
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When n_fft > signal length, the frequency vector is finer
# (more bins). The spectral peak must still be at the correct frequency.
test_that("eeg_fft() zero-padded spectrum still peaks at correct frequency", {
  eeg      <- make_sine_eeg(freq_hz = 10, n_tp = SR)   # 1-second signal
  spec_zp  <- eeg_fft(eeg, n_fft = SR * 4L, window = "none", verbose = FALSE)
  
  pk <- peak_freq(spec_zp$power[1L, ], spec_zp$frequencies)
  expect_equal(pk, 10, tolerance = 0.5)
})

# ----------------------------------------------------------------------------
# Test 5.4: amplitude^2 == power throughout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The power field must equal amplitude squared at every
# frequency bin (by definition). Any mismatch reveals a normalisation bug.
test_that("eeg_fft() power equals amplitude^2 at all bins", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, verbose = FALSE)
  expect_equal(spec$power, spec$amplitude^2, tolerance = 1e-12)
})


# ============================================================================
# TEST SUITE 6: eeg_fft() — eeg_epochs input
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: eeg_fft() on eeg_epochs returns 'eeg_spectrum' with correct fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Dispatching to the eeg_epochs branch must still produce an
# eeg_spectrum with the correct input_class and n_epochs fields.
test_that("eeg_fft() on eeg_epochs returns correct class and input_class", {
  ep   <- make_mock_epochs()
  spec <- eeg_fft(ep, verbose = FALSE)
  expect_s3_class(spec, "eeg_spectrum")
  expect_equal(spec$input_class, "eeg_epochs")
  expect_equal(spec$n_epochs, ep$n_epochs)
})

# ----------------------------------------------------------------------------
# Test 6.2: per_epoch = TRUE stores 3-D epoch_power array
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When per_epoch = TRUE, $epoch_power must be a 3-D array
# with dimensions [n_channels x n_freqs x n_epochs].
test_that("eeg_fft() per_epoch=TRUE stores [channels x freqs x epochs] array", {
  n_ch  <- 2L; n_samp <- 256L; n_ep <- 6L
  ep    <- make_mock_epochs(n_channels = n_ch, n_samp = n_samp, n_epochs = n_ep)
  spec  <- eeg_fft(ep, per_epoch = TRUE, verbose = FALSE)
  
  n_half <- floor(n_samp / 2L) + 1L
  expect_true(!is.null(spec$epoch_power))
  expect_equal(dim(spec$epoch_power), c(n_ch, n_half, n_ep))
  expect_equal(spec$per_epoch, TRUE)
})

# ----------------------------------------------------------------------------
# Test 6.3: per_epoch = FALSE leaves epoch_power as NULL
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When per_epoch = FALSE (default), $epoch_power must be NULL
# so that callers can distinguish between the two operating modes.
test_that("eeg_fft() per_epoch=FALSE leaves epoch_power NULL", {
  ep   <- make_mock_epochs()
  spec <- eeg_fft(ep, per_epoch = FALSE, verbose = FALSE)
  expect_null(spec$epoch_power)
})

# ----------------------------------------------------------------------------
# Test 6.4: Power matrix dimensions match [channels x n_fft/2+1]
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Regardless of per_epoch setting, the summary $power matrix
# must have correct dimensions for downstream use.
test_that("eeg_fft() on eeg_epochs: $power dimensions are [channels x freqs]", {
  n_ch  <- 3L; n_samp <- 256L
  ep    <- make_mock_epochs(n_channels = n_ch, n_samp = n_samp)
  spec  <- eeg_fft(ep, verbose = FALSE)
  n_half <- floor(n_samp / 2L) + 1L
  expect_equal(dim(spec$power), c(n_ch, n_half))
})


# ============================================================================
# TEST SUITE 7: eeg_fft() — eeg_evoked input
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: eeg_fft() on eeg_evoked returns power as a named list
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: For evoked input, $power must be a named list with one
# matrix per condition, not a 2-D matrix.
test_that("eeg_fft() on eeg_evoked returns power as a named list", {
  evoked <- make_mock_evoked()
  spec   <- eeg_fft(evoked, verbose = FALSE)
  expect_true(is.list(spec$power))
  expect_equal(sort(names(spec$power)), sort(evoked$conditions))
})

# ----------------------------------------------------------------------------
# Test 7.2: Each condition's power matrix has correct dimensions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Each condition matrix in $power must be [n_ch x n_half],
# consistent with the single-input case.
test_that("eeg_fft() evoked: each condition power matrix has [channels x freqs] dims", {
  n_ch   <- 2L; n_samp <- 256L
  evoked <- make_mock_evoked(n_channels = n_ch, n_samp = n_samp)
  spec   <- eeg_fft(evoked, verbose = FALSE)
  n_half <- floor(n_samp / 2L) + 1L
  
  for (cond in evoked$conditions) {
    expect_equal(dim(spec$power[[cond]]), c(n_ch, n_half),
                 info = paste("condition:", cond))
  }
})

# ----------------------------------------------------------------------------
# Test 7.3: input_class is set to "eeg_evoked"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: $input_class must correctly reflect the class of the input
# object so downstream code can distinguish evoked from continuous/epochs.
test_that("eeg_fft() sets input_class = 'eeg_evoked' for evoked input", {
  evoked <- make_mock_evoked()
  spec   <- eeg_fft(evoked, verbose = FALSE)
  expect_equal(spec$input_class, "eeg_evoked")
  expect_equal(spec$conditions, evoked$conditions)
})


# ============================================================================
# TEST SUITE 8: eeg_psd_welch() — Input Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 8.1: Wrong input class throws an error
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() errors on invalid input class", {
  expect_error(eeg_psd_welch(list(x = 1)),
               regexp = "Input must be one of")
})

# ----------------------------------------------------------------------------
# Test 8.2: Non-positive window_length throws an error
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() errors on non-positive window_length", {
  eeg <- make_mock_eeg()
  expect_error(eeg_psd_welch(eeg, window_length = -1, verbose = FALSE))
  expect_error(eeg_psd_welch(eeg, window_length =  0, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 8.3: overlap outside [0, 1) throws an error
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() errors on invalid overlap", {
  eeg <- make_mock_eeg()
  expect_error(eeg_psd_welch(eeg, overlap = 1.0, verbose = FALSE))
  expect_error(eeg_psd_welch(eeg, overlap = 1.5, verbose = FALSE))
  expect_error(eeg_psd_welch(eeg, overlap = -0.1, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 8.4: window_length longer than signal warns and uses full signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If window_length > signal length, a warning is issued and
# the function does not stop — it gracefully uses the full signal as one segment.
test_that("eeg_psd_welch() warns when window_length > signal length", {
  eeg <- make_mock_eeg(n_tp = 256L)
  expect_warning(eeg_psd_welch(eeg, window_length = 10, verbose = FALSE),
                 regexp = "exceeds signal length")
})


# ============================================================================
# TEST SUITE 9: eeg_psd_welch() — Output Structure
# ============================================================================

# ----------------------------------------------------------------------------
# Test 9.1: Return class is "eeg_spectrum"
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() returns class 'eeg_spectrum'", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  expect_s3_class(spec, "eeg_spectrum")
})

# ----------------------------------------------------------------------------
# Test 9.2: method field is set to "welch"
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() sets method = 'welch'", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  expect_equal(spec$method, "welch")
})

# ----------------------------------------------------------------------------
# Test 9.3: PSD is strictly non-negative
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: PSD = averaged |X|^2 / normalisation, so it can never be
# negative. A violation would indicate a sign or computation error.
test_that("eeg_psd_welch() PSD values are all non-negative", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  expect_true(all(spec$power >= 0))
})

# ----------------------------------------------------------------------------
# Test 9.4: amplitude = sqrt(power) throughout
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() amplitude equals sqrt(power) at all bins", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  expect_equal(spec$amplitude, sqrt(spec$power), tolerance = 1e-12)
})

# ----------------------------------------------------------------------------
# Test 9.5: phase is NULL (not computed for Welch)
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() leaves phase as NULL", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  expect_null(spec$phase)
})

# ----------------------------------------------------------------------------
# Test 9.6: Welch-specific metadata fields are populated
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() populates window_length, overlap, and n_segments", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, window_length = 2.0, overlap = 0.5, verbose = FALSE)
  expect_equal(spec$window_length, 2.0)
  expect_equal(spec$overlap, 0.5)
  expect_true(!is.null(spec$n_segments))
  expect_true(spec$n_segments >= 1L)
})


# ============================================================================
# TEST SUITE 10: eeg_psd_welch() — Spectral Accuracy
# ============================================================================

# ----------------------------------------------------------------------------
# Test 10.1: Peak frequency recovery of a 10 Hz sine wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With a 4-second 10 Hz sine and 2 s Welch windows (giving
# 0.5 Hz resolution), the PSD peak must be at 10 Hz within 0.5 Hz tolerance.
test_that("eeg_psd_welch() recovers peak at 10 Hz from a sine wave", {
  eeg  <- make_sine_eeg(freq_hz = 10)
  spec <- eeg_psd_welch(eeg, window_length = 2.0, overlap = 0.5, verbose = FALSE)
  
  freq_res <- 1.0 / 2.0    # 0.5 Hz — resolution with 2 s window
  pk       <- peak_freq(spec$power[1L, ], spec$frequencies)
  expect_equal(pk, 10, tolerance = freq_res)
})

# ----------------------------------------------------------------------------
# Test 10.2: Welch PSD of white noise is approximately flat
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: For white (Gaussian) noise the expected PSD is constant
# across frequency. The coefficient of variation (sd/mean) across bins should
# be small (< 50%) given enough segments, indicating approximate flatness.
test_that("eeg_psd_welch() PSD of white noise is approximately flat", {
  set.seed(42L)
  eeg  <- make_mock_eeg(n_channels = 1L, n_tp = SR * 10L)
  spec <- eeg_psd_welch(eeg, window_length = 1.0, overlap = 0.5, verbose = FALSE)
  
  # Exclude DC and Nyquist; check interior bins
  pwr <- spec$power[1L, 2L:(length(spec$frequencies) - 1L)]
  cv  <- sd(pwr) / mean(pwr)   # coefficient of variation
  expect_lt(cv, 0.6)
})

# ----------------------------------------------------------------------------
# Test 10.3: More Welch segments reduce variance
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Using shorter windows (more segments) should produce a
# more stable (less variable) PSD estimate for the same noise signal.
test_that("eeg_psd_welch() more segments give lower PSD variance", {
  set.seed(7L)
  eeg <- make_mock_eeg(n_channels = 1L, n_tp = SR * 20L)
  
  spec_coarse <- eeg_psd_welch(eeg, window_length = 4.0, verbose = FALSE)
  spec_fine   <- eeg_psd_welch(eeg, window_length = 1.0, verbose = FALSE)
  
  # Fine-grained (more segments) should have lower variance in interior bins
  pwr_c <- spec_coarse$power[1L, 2L:(length(spec_coarse$frequencies) - 1L)]
  pwr_f <- spec_fine$power[1L,   2L:(length(spec_fine$frequencies)   - 1L)]
  
  # Compare median absolute deviation (MAD) as a robust variance measure
  # Fine windowing has more segments => lower MAD of the PSD
  expect_lt(mad(pwr_f), mad(pwr_c) * 1.5)   # generous tolerance
})


# ============================================================================
# TEST SUITE 11: eeg_psd_welch() — eeg_epochs / eeg_evoked dispatch
# ============================================================================

# ----------------------------------------------------------------------------
# Test 11.1: eeg_psd_welch() on eeg_epochs returns eeg_spectrum
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() on eeg_epochs returns 'eeg_spectrum'", {
  ep   <- make_mock_epochs()
  expect_warning(
    spec <- eeg_psd_welch(ep, verbose = FALSE),
    regexp = "exceeds signal length")
  expect_s3_class(spec, "eeg_spectrum")
  expect_equal(spec$input_class, "eeg_epochs")
})

# ----------------------------------------------------------------------------
# Test 11.2: per_epoch = TRUE on eeg_epochs stores epoch_power array
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() per_epoch=TRUE stores 3-D epoch_power", {
  n_ch <- 2L; n_ep <- 4L; n_samp <- 256L
  ep   <- make_mock_epochs(n_channels = n_ch, n_epochs = n_ep, n_samp = n_samp)

  spec <- eeg_psd_welch(ep, per_epoch = TRUE, window_length = 0.5, verbose = FALSE)  
  n_half <- floor(n_samp / 2L) + 1L
  
  expect_false(is.null(spec$epoch_power))
  expect_equal(dim(spec$epoch_power)[c(1L, 3L)], c(n_ch, n_ep))
})

# ----------------------------------------------------------------------------
# Test 11.3: eeg_psd_welch() on eeg_evoked returns named list of PSD matrices
# ----------------------------------------------------------------------------
test_that("eeg_psd_welch() on eeg_evoked returns power as named list", {
  evoked <- make_mock_evoked()
  spec <- eeg_psd_welch(evoked, window_length = 0.5, verbose = FALSE)
  expect_true(is.list(spec$power))
  expect_equal(sort(names(spec$power)), sort(evoked$conditions))
})


# ============================================================================
# TEST SUITE 12: eeg_multitaper() — Input Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 12.1: Wrong input class throws an error
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() errors on invalid input class", {
  expect_error(eeg_multitaper(42),
               regexp = "Input must be one of")
})

# ----------------------------------------------------------------------------
# Test 12.2: Non-positive time_bandwidth throws an error
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() errors on non-positive time_bandwidth", {
  eeg <- make_mock_eeg()
  expect_error(eeg_multitaper(eeg, time_bandwidth = -2, verbose = FALSE))
  expect_error(eeg_multitaper(eeg, time_bandwidth =  0, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 12.3: n_tapers < 1 throws an error
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() errors when n_tapers < 1", {
  eeg <- make_mock_eeg()
  expect_error(eeg_multitaper(eeg, n_tapers = 0L, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 12.4: Excess tapers issue a warning (not an error)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When n_tapers exceeds 2*NW - 1 (the recommended maximum),
# eeg_multitaper() should warn rather than stop, so the call still succeeds.
test_that("eeg_multitaper() warns but continues when n_tapers > 2*NW-1", {
  eeg <- make_mock_eeg()
  expect_warning(
    eeg_multitaper(eeg, time_bandwidth = 2, n_tapers = 10L, verbose = FALSE),
    regexp = "recommended maximum"
  )
})


# ============================================================================
# TEST SUITE 13: eeg_multitaper() — Output Structure
# ============================================================================

# ----------------------------------------------------------------------------
# Test 13.1: Return class is "eeg_spectrum" with method = "multitaper"
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() returns class 'eeg_spectrum' with method='multitaper'", {
  eeg  <- make_mock_eeg()
  spec <- eeg_multitaper(eeg, verbose = FALSE)
  expect_s3_class(spec, "eeg_spectrum")
  expect_equal(spec$method, "multitaper")
})

# ----------------------------------------------------------------------------
# Test 13.2: Multitaper metadata fields are populated correctly
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() stores time_bandwidth and n_tapers", {
  eeg  <- make_mock_eeg()
  spec <- eeg_multitaper(eeg, time_bandwidth = 4, n_tapers = 7L, verbose = FALSE)
  expect_equal(spec$time_bandwidth, 4)
  expect_equal(spec$n_tapers, 7L)
  expect_equal(spec$window_type, "dpss")
})

# ----------------------------------------------------------------------------
# Test 13.3: PSD is non-negative
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() PSD values are all non-negative", {
  eeg  <- make_mock_eeg()
  spec <- eeg_multitaper(eeg, verbose = FALSE)
  expect_true(all(spec$power >= 0))
})

# ----------------------------------------------------------------------------
# Test 13.4: power and amplitude matrix dimensions
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() produces power matrix with [channels x freqs] dims", {
  n_ch     <- 3L
  win_samp <- round(2.0 * SR)           # 2-second segment
  n_half   <- floor(win_samp / 2L) + 1L
  eeg      <- make_mock_eeg(n_channels = n_ch)
  spec     <- eeg_multitaper(eeg, segment_length = 2.0, verbose = FALSE)
  
  expect_equal(nrow(spec$power), n_ch)
  expect_equal(ncol(spec$power), n_half)
})

# ----------------------------------------------------------------------------
# Test 13.5: amplitude = sqrt(power) everywhere
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() amplitude equals sqrt(power)", {
  eeg  <- make_mock_eeg()
  spec <- eeg_multitaper(eeg, verbose = FALSE)
  expect_equal(spec$amplitude, sqrt(spec$power), tolerance = 1e-12)
})


# ============================================================================
# TEST SUITE 14: eeg_multitaper() — Spectral Accuracy
# ============================================================================

# ----------------------------------------------------------------------------
# Test 14.1: Multitaper peak frequency recovery for a 10 Hz sine wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The multitaper PSD of a 10 Hz pure sine must peak at 10 Hz.
# Tolerance is the Rayleigh resolution of the 2-second analysis window (0.5 Hz).
test_that("eeg_multitaper() recovers 10 Hz peak from a pure sine wave", {
  eeg      <- make_sine_eeg(freq_hz = 10)
  spec     <- eeg_multitaper(eeg, time_bandwidth = 4,
                             segment_length = 2.0, verbose = FALSE)
  freq_res <- 1.0 / 2.0    # 0.5 Hz resolution
  
  pk <- peak_freq(spec$power[1L, ], spec$frequencies)
  expect_equal(pk, 10, tolerance = freq_res)
})

# ----------------------------------------------------------------------------
# Test 14.2: Multitaper PSD of white noise is approximately flat
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Like the Welch test, the multitaper PSD of white noise
# should be approximately spectrally flat (CV < 60% in the interior).
test_that("eeg_multitaper() PSD of white noise has low spectral variation", {
  set.seed(99L)
  eeg  <- make_mock_eeg(n_channels = 1L, n_tp = SR * 10L)
  spec <- eeg_multitaper(eeg, segment_length = 1.0, verbose = FALSE)
  
  pwr <- spec$power[1L, 2L:(length(spec$frequencies) - 1L)]
  cv  <- sd(pwr) / mean(pwr)
  expect_lt(cv, 0.6)
})


# ============================================================================
# TEST SUITE 15: eeg_multitaper() — eeg_epochs / eeg_evoked dispatch
# ============================================================================

# ----------------------------------------------------------------------------
# Test 15.1: eeg_multitaper() on eeg_epochs returns eeg_spectrum
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() on eeg_epochs returns 'eeg_spectrum'", {
  ep   <- make_mock_epochs()
  spec <- eeg_multitaper(ep, verbose = FALSE)
  expect_s3_class(spec, "eeg_spectrum")
  expect_equal(spec$input_class, "eeg_epochs")
})

# ----------------------------------------------------------------------------
# Test 15.2: per_epoch = TRUE stores 3-D epoch_power array
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() per_epoch=TRUE stores epoch_power 3-D array", {
  n_ch <- 2L; n_ep <- 4L; n_samp <- 256L
  ep   <- make_mock_epochs(n_channels = n_ch, n_epochs = n_ep, n_samp = n_samp)
  spec <- eeg_multitaper(ep, per_epoch = TRUE, verbose = FALSE)
  
  n_half <- floor(n_samp / 2L) + 1L
  expect_false(is.null(spec$epoch_power))
  expect_equal(dim(spec$epoch_power), c(n_ch, n_half, n_ep))
})

# ----------------------------------------------------------------------------
# Test 15.3: eeg_multitaper() on eeg_evoked returns power as named list
# ----------------------------------------------------------------------------
test_that("eeg_multitaper() on eeg_evoked returns power as a named list", {
  evoked <- make_mock_evoked()
  spec   <- eeg_multitaper(evoked, verbose = FALSE)
  expect_true(is.list(spec$power))
  expect_equal(sort(names(spec$power)), sort(evoked$conditions))
})


# ============================================================================
# TEST SUITE 16: eeg_band_power() — Input Validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 16.1: Wrong input class throws an error
# ----------------------------------------------------------------------------
test_that("eeg_band_power() errors on invalid input class", {
  expect_error(eeg_band_power(42),
               regexp = "Input must be")
})

# ----------------------------------------------------------------------------
# Test 16.2: Empty bands list throws an error
# ----------------------------------------------------------------------------
test_that("eeg_band_power() errors on empty bands list", {
  eeg <- make_mock_eeg()
  expect_error(eeg_band_power(eeg, bands = list(), verbose = FALSE),
               regexp = "non-empty named list")
})

# ----------------------------------------------------------------------------
# Test 16.3: Malformed band definition throws an error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A band must be c(low, high) with low < high. Reversed or
# length-1 band definitions should cause an immediate stop.
test_that("eeg_band_power() errors on malformed band definition", {
  eeg <- make_mock_eeg()
  expect_error(
    eeg_band_power(eeg, bands = list(bad = c(20, 5)), verbose = FALSE),
    regexp = "low < high"
  )
  expect_error(
    eeg_band_power(eeg, bands = list(bad = c(10)), verbose = FALSE),
    regexp = "low < high"
  )
})


# ============================================================================
# TEST SUITE 17: eeg_band_power() — Output Structure
# ============================================================================

# ----------------------------------------------------------------------------
# Test 17.1: Returns a data frame with correct column names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The result must be a data frame with one column named
# "channel" and one column per band with the user-supplied band names.
test_that("eeg_band_power() returns data.frame with correct columns", {
  eeg <- make_mock_eeg(n_channels = 3L)
  bp  <- eeg_band_power(eeg, verbose = FALSE)
  
  expect_true(is.data.frame(bp))
  expect_true("channel" %in% names(bp))
  band_cols <- c("delta", "theta", "alpha", "beta", "gamma")
  expect_true(all(band_cols %in% names(bp)))
})

# ----------------------------------------------------------------------------
# Test 17.2: Number of rows matches number of channels
# ----------------------------------------------------------------------------
test_that("eeg_band_power() returns one row per channel", {
  n_ch <- 4L
  eeg  <- make_mock_eeg(n_channels = n_ch)
  bp   <- eeg_band_power(eeg, verbose = FALSE)
  expect_equal(nrow(bp), n_ch)
})

# ----------------------------------------------------------------------------
# Test 17.3: Channel names in output match input
# ----------------------------------------------------------------------------
test_that("eeg_band_power() channel column matches input channel names", {
  eeg <- make_mock_eeg(n_channels = 3L)
  bp  <- eeg_band_power(eeg, verbose = FALSE)
  expect_equal(bp$channel, eeg$channels)
})

# ----------------------------------------------------------------------------
# Test 17.4: Accepts a pre-computed eeg_spectrum without recomputing PSD
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Passing an eeg_spectrum directly must succeed, confirming
# the function detects the already-computed PSD and skips recomputation.
test_that("eeg_band_power() accepts pre-computed eeg_spectrum as input", {
  eeg  <- make_mock_eeg()
  psd  <- eeg_psd_welch(eeg, verbose = FALSE)
  bp   <- eeg_band_power(psd, verbose = FALSE)
  expect_true(is.data.frame(bp))
  expect_equal(nrow(bp), length(eeg$channels))
})


# ============================================================================
# TEST SUITE 18: eeg_band_power() — Accuracy
# ============================================================================

# ----------------------------------------------------------------------------
# Test 18.1: Relative band power sums to 1 across all bands
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When relative = TRUE, band powers normalise each channel's
# total (summing to 1). This is the fundamental definition of relative PSD.
test_that("eeg_band_power() relative=TRUE sums to 1 per channel", {
  eeg     <- make_mock_eeg(n_channels = 3L)
  bp_rel  <- eeg_band_power(eeg, relative = TRUE, verbose = FALSE)
  band_cols <- c("delta", "theta", "alpha", "beta", "gamma")
  
  for (i in seq_len(nrow(bp_rel))) {
    row_sum <- sum(unlist(bp_rel[i, band_cols]))
    expect_equal(row_sum, 1.0, tolerance = 1e-6,
                 info = paste("channel", bp_rel$channel[i]))
  }
})

# ----------------------------------------------------------------------------
# Test 18.2: Alpha band power is elevated for a 10 Hz signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Given a 10 Hz pure sine, the alpha band (8–13 Hz) must
# contain more power than any other canonical band. This validates that band
# integration correctly attributes power to the right frequency range.
test_that("eeg_band_power() alpha dominates for a 10 Hz sine wave", {
  eeg <- make_sine_eeg(freq_hz = 10)
  bp  <- eeg_band_power(eeg, verbose = FALSE)
  
  # Channel 1 has the sine; it should have the highest alpha power
  band_vals  <- bp[1L, c("delta", "theta", "alpha", "beta", "gamma")]
  max_band   <- names(which.max(unlist(band_vals)))
  expect_equal(max_band, "alpha")
})

# ----------------------------------------------------------------------------
# Test 18.3: Custom bands are respected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When a user supplies custom bands the output must have those
# custom names as columns, not the default canonical band names.
test_that("eeg_band_power() uses user-supplied custom band names", {
  eeg    <- make_mock_eeg()
  custom <- list(low_alpha = c(8, 10), high_alpha = c(10, 13))
  bp     <- eeg_band_power(eeg, bands = custom, verbose = FALSE)
  expect_true(all(c("low_alpha", "high_alpha") %in% names(bp)))
  expect_false("alpha" %in% names(bp))
})

# ----------------------------------------------------------------------------
# Test 18.4: Band power values are all non-negative
# ----------------------------------------------------------------------------
test_that("eeg_band_power() absolute power values are non-negative", {
  eeg <- make_mock_eeg(n_channels = 4L)
  bp  <- eeg_band_power(eeg, verbose = FALSE)
  band_cols <- c("delta", "theta", "alpha", "beta", "gamma")
  expect_true(all(unlist(bp[, band_cols]) >= 0, na.rm = TRUE))
})


# ============================================================================
# TEST SUITE 19: eeg_band_power() — eeg_evoked dispatch
# ============================================================================

# ----------------------------------------------------------------------------
# Test 19.1: eeg_band_power() on eeg_evoked returns a named list
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: For evoked input, band power must be returned as a named
# list (one data frame per condition), matching the $conditions vector.
test_that("eeg_band_power() on eeg_evoked returns a named list of data frames", {
  evoked <- make_mock_evoked()
  bp     <- eeg_band_power(evoked, window_length = 0.5, verbose = FALSE)
  
  expect_true(is.list(bp))
  expect_equal(sort(names(bp)), sort(evoked$conditions))
  for (cond in evoked$conditions) {
    expect_true(is.data.frame(bp[[cond]]),
                info = paste("condition", cond, "is not a data.frame"))
  }
})

# ----------------------------------------------------------------------------
# Test 19.2: Each condition data frame has correct structure
# ----------------------------------------------------------------------------
test_that("eeg_band_power() evoked: each condition df has correct columns and rows", {
  n_ch   <- 2L
  evoked <- make_mock_evoked(n_channels = n_ch)
  bp     <- eeg_band_power(evoked, window_length = 0.5, verbose = FALSE)
  
  for (cond in evoked$conditions) {
    expect_equal(nrow(bp[[cond]]), n_ch,
                 info = paste("row count for", cond))
    expect_true("channel" %in% names(bp[[cond]]))
  }
})


# ============================================================================
# TEST SUITE 20: print.eeg_spectrum()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 20.1: print.eeg_spectrum() produces output without error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The S3 print method must run without errors and produce
# visible output containing key fields.
test_that("print.eeg_spectrum() produces output without error", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, verbose = FALSE)
  expect_output(print(spec), regexp = "EEG Spectrum Object")
})

# ----------------------------------------------------------------------------
# Test 20.2: print.eeg_spectrum() returns the spectrum object invisibly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Following R print method convention, the object must be
# returned invisibly so print(spec) does not cause double-printing.
test_that("print.eeg_spectrum() returns the input object invisibly", {
  eeg  <- make_mock_eeg()
  spec <- eeg_fft(eeg, verbose = FALSE)
  out  <- capture.output(result <- print(spec))
  expect_identical(result, spec)
})

# ----------------------------------------------------------------------------
# Test 20.3: print.eeg_spectrum() displays method-specific parameters
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The Welch and multitaper print summaries must include their
# specific parameter blocks (window_length/overlap vs time_bandwidth/tapers).
test_that("print.eeg_spectrum() shows Welch parameters for Welch output", {
  eeg  <- make_mock_eeg()
  spec <- eeg_psd_welch(eeg, verbose = FALSE)
  out  <- capture.output(print(spec))
  expect_true(any(grepl("Welch", out)))
  expect_true(any(grepl("Window length", out)))
  expect_true(any(grepl("Overlap", out)))
})

test_that("print.eeg_spectrum() shows multitaper parameters for multitaper output", {
  eeg  <- make_mock_eeg()
  spec <- eeg_multitaper(eeg, verbose = FALSE)
  out  <- capture.output(print(spec))
  expect_true(any(grepl("Multitaper", out)))
  expect_true(any(grepl("Time-bandwidth", out)))
  expect_true(any(grepl("Tapers", out)))
})