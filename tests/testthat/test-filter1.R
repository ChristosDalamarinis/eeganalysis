# test_filter1.R
# Unit tests for filter1.R — FIR filtering pipeline
#
# Run with:
#   testthat::test_file("test_filter1.R")
#
# Requires: testthat, mne_reference/ reference CSVs

library(testthat)

# ─── Minimal eeg constructor (eeg_class.R not present here) ──────────────────

new_eeg <- function(data, channels, sampling_rate, times = NULL,
                    events = NULL, metadata = NULL, reference = "unknown",
                    preprocessing_history = list()) {
  structure(
    list(
      data                  = data,
      channels              = channels,
      sampling_rate         = sampling_rate,
      times                 = times,
      events                = events,
      metadata              = metadata,
      reference             = reference,
      preprocessing_history = preprocessing_history
    ),
    class = "eeg"
  )
}

library(eeganalysis)

#REF      <- "C:/Users/dalam/OneDrive/R Projects/R Package/translating/mne_reference"
load_ref <- function(name) scan(file.path(REF, name), quiet = TRUE)

raw_signal <- load_ref("raw_signal.csv")
eeg1ch <- new_eeg(
  data          = matrix(raw_signal, nrow = 1L),
  channels      = "CH1",
  sampling_rate = 256L
)

# ─── .next_fast_len() ────────────────────────────────────────────────────────

test_that(".next_fast_len returns target unchanged when <= 6", {
  expect_equal(.next_fast_len(1L), 1L)
  expect_equal(.next_fast_len(6L), 6L)
})

test_that(".next_fast_len returns exact powers of 2", {
  for (exp in c(3L, 4L, 8L, 10L)) {
    n <- 2L^exp
    expect_equal(.next_fast_len(n), n,
                 label = paste0(".next_fast_len(", n, ")"))
  }
})

test_that(".next_fast_len finds next 5-smooth number (lookup table range)", {
  expect_equal(.next_fast_len(7L),   8L)
  expect_equal(.next_fast_len(11L),  12L)
  expect_equal(.next_fast_len(13L),  15L)
  expect_equal(.next_fast_len(101L), 108L)
  expect_equal(.next_fast_len(1000L), 1000L)
  expect_equal(.next_fast_len(1001L), 1024L)
})

test_that(".next_fast_len result for large target has only prime factors 2,3,5", {
  result <- .next_fast_len(10001L)
  expect_gte(result, 10001L)
  n <- result
  for (p in c(2L, 3L, 5L)) while (n %% p == 0L) n <- n %/% p
  expect_equal(n, 1L)
})

# ─── .auto_trans_bandwidth() ─────────────────────────────────────────────────

test_that(".auto_trans_bandwidth highpass-only: correct formula and NULL h", {
  r <- .auto_trans_bandwidth(l_freq = 1.0, sfreq = 256)
  expect_equal(r$l_trans_bw, min(max(0.25 * 1.0, 2.0), 1.0))
  expect_null(r$h_trans_bw)
})

test_that(".auto_trans_bandwidth lowpass-only: correct formula and NULL l", {
  r <- .auto_trans_bandwidth(h_freq = 40.0, sfreq = 256)
  expect_null(r$l_trans_bw)
  expect_equal(r$h_trans_bw, min(max(0.25 * 40.0, 2.0), 128.0 - 40.0))
})

test_that(".auto_trans_bandwidth bandpass: both computed independently", {
  r <- .auto_trans_bandwidth(l_freq = 1.0, h_freq = 40.0, sfreq = 256)
  expect_equal(r$l_trans_bw, min(max(0.25 * 1.0,  2.0), 1.0))
  expect_equal(r$h_trans_bw, min(max(0.25 * 40.0, 2.0), 128.0 - 40.0))
})

test_that(".auto_trans_bandwidth 2 Hz floor kicks in for small l_freq", {
  # 0.25 * 1.0 = 0.25, floor raises it to 2.0, then clamped by l_freq=1 -> 1.0
  r <- .auto_trans_bandwidth(l_freq = 1.0, sfreq = 256)
  expect_equal(r$l_trans_bw, 1.0)

  # 0.25 * 10 = 2.5 > 2.0 floor, clamped by l_freq=10 -> 2.5
  r2 <- .auto_trans_bandwidth(l_freq = 10.0, sfreq = 256)
  expect_equal(r2$l_trans_bw, 2.5)
})

# ─── .auto_filter_length() ───────────────────────────────────────────────────

test_that(".auto_filter_length returns an odd integer", {
  N <- .auto_filter_length(l_trans_bw = 1.0, sfreq = 256)
  expect_true(is.integer(N))
  expect_equal(N %% 2L, 1L)
})

test_that(".auto_filter_length: narrower bandwidth -> longer filter", {
  N_narrow <- .auto_filter_length(l_trans_bw = 1.0, sfreq = 256)
  N_wide   <- .auto_filter_length(l_trans_bw = 10.0, sfreq = 256)
  expect_gt(N_narrow, N_wide)
})

test_that(".auto_filter_length uses minimum of two bandwidths", {
  # min(1, 10) = 1 -> same as l_only = 1
  N_both <- .auto_filter_length(l_trans_bw = 1.0, h_trans_bw = 10.0, sfreq = 256)
  N_l    <- .auto_filter_length(l_trans_bw = 1.0, sfreq = 256)
  expect_equal(N_both, N_l)
})

test_that(".auto_filter_length length factors differ by window", {
  N_ham <- .auto_filter_length(l_trans_bw = 2.0, sfreq = 256, window = "hamming")
  N_han <- .auto_filter_length(l_trans_bw = 2.0, sfreq = 256, window = "hann")
  N_bla <- .auto_filter_length(l_trans_bw = 2.0, sfreq = 256, window = "blackman")
  expect_gt(N_bla, N_ham)  # blackman factor 5.0 > hamming 3.3
  expect_gt(N_ham, N_han)  # hamming 3.3 > hann 3.1
})

test_that(".auto_filter_length errors without any bandwidth argument", {
  expect_error(.auto_filter_length(sfreq = 256),
               "At least one of l_trans_bw or h_trans_bw must be non-NULL")
})

test_that(".auto_filter_length errors on unknown window name", {
  expect_error(.auto_filter_length(l_trans_bw = 1.0, sfreq = 256, window = "kaiser"),
               "Unknown window")
})

# ─── .sinc_norm() ────────────────────────────────────────────────────────────

test_that(".sinc_norm returns 1 at x = 0", {
  expect_equal(.sinc_norm(0), 1.0)
})

test_that(".sinc_norm returns 0 at nonzero integers", {
  expect_equal(.sinc_norm(1),  0.0, tolerance = 1e-15)
  expect_equal(.sinc_norm(-2), 0.0, tolerance = 1e-15)
  expect_equal(.sinc_norm(3),  0.0, tolerance = 1e-15)
})

test_that(".sinc_norm matches sin(pi*x)/(pi*x) for non-integer x", {
  x <- c(0.5, 1.3, -0.7, 2.1)
  expect_equal(.sinc_norm(x), sin(pi * x) / (pi * x))
})

# ─── .firwin_lowpass() ───────────────────────────────────────────────────────

test_that(".firwin_lowpass returns vector of correct length", {
  expect_length(.firwin_lowpass(101L, 0.4), 101L)
  expect_length(.firwin_lowpass(51L,  0.2), 51L)
})

test_that(".firwin_lowpass DC gain equals 1 (normalized)", {
  h <- .firwin_lowpass(101L, 0.4)
  expect_equal(sum(h), 1.0, tolerance = 1e-14)
})

test_that(".firwin_lowpass kernel is symmetric (linear phase FIR)", {
  h <- .firwin_lowpass(101L, 0.3)
  expect_equal(h, rev(h), tolerance = 1e-15)
})

test_that(".firwin_lowpass errors on unknown window", {
  expect_error(.firwin_lowpass(101L, 0.4, window = "unknown_win"))
})

test_that(".firwin_lowpass accepts hann and blackman windows", {
  expect_no_error(.firwin_lowpass(101L, 0.4, window = "hann"))
  expect_no_error(.firwin_lowpass(101L, 0.4, window = "blackman"))
})

# ─── .firwin_kernel() ────────────────────────────────────────────────────────

test_that(".firwin_kernel returns vector of length N", {
  h <- .firwin_kernel(N    = 201L,
                      freq = c(0, 0.1, 0.2, 0.9, 1.0),
                      gain = c(0,   0,   1,   1,   0))
  expect_length(h, 201L)
})

test_that(".firwin_kernel is symmetric (linear phase)", {
  h <- .firwin_kernel(201L,
                      freq = c(0, 0.1, 0.2, 0.9, 1.0),
                      gain = c(0,   0,   1,   1,   0))
  expect_equal(h, rev(h), tolerance = 1e-14)
})

test_that(".firwin_kernel errors when N is even", {
  expect_error(
    .firwin_kernel(200L, freq = c(0, 0.2, 0.8, 1), gain = c(0, 1, 1, 0))
  )
})

test_that(".firwin_kernel BP1 (1-40 Hz, hamming) matches MNE reference kernel", {
  l_tb  <- min(max(0.25 * 1.0,  2.0), 1.0)
  h_tb  <- min(max(0.25 * 40.0, 2.0), 128.0 - 40.0)
  N     <- .auto_filter_length(l_tb, h_tb, 256)
  h_r   <- .firwin_kernel(
    N    = N,
    freq = c(0,
             (1.0 - l_tb) / 128, 1.0 / 128,
             40.0 / 128,         (40.0 + h_tb) / 128,
             1),
    gain = c(0, 0, 1, 1, 0, 0)
  )
  h_mne <- load_ref("bp1_kernel.csv")
  expect_equal(length(h_r), length(h_mne))
  expect_lt(max(abs(h_r - h_mne)), 1e-14)
})

test_that(".firwin_kernel BP5 (hann window) matches MNE reference kernel", {
  l_tb  <- min(max(0.25 * 1.0,  2.0), 1.0)
  h_tb  <- min(max(0.25 * 40.0, 2.0), 88.0)
  N     <- .auto_filter_length(l_tb, h_tb, 256, window = "hann")
  h_r   <- .firwin_kernel(
    N      = N,
    freq   = c(0,
               (1.0 - l_tb) / 128, 1.0 / 128,
               40.0 / 128,         (40.0 + h_tb) / 128,
               1),
    gain   = c(0, 0, 1, 1, 0, 0),
    window = "hann"
  )
  h_mne <- load_ref("bp5_kernel.csv")
  expect_equal(length(h_r), length(h_mne))
  expect_lt(max(abs(h_r - h_mne)), 1e-14)
})

# ─── .reflect_limited_pad() ──────────────────────────────────────────────────

test_that(".reflect_limited_pad output length = n + n_left + n_right", {
  x <- as.double(1:10)
  expect_length(.reflect_limited_pad(x, 4L, 3L), 17L)
  expect_length(.reflect_limited_pad(x, 0L, 0L), 10L)
})

test_that(".reflect_limited_pad preserves original signal in middle", {
  x      <- as.double(1:10)
  result <- .reflect_limited_pad(x, 3L, 5L)
  expect_equal(result[4:13], x)
})

test_that(".reflect_limited_pad left reflection: 2*x[1] - x[reversed inner]", {
  x      <- as.double(c(10, 20, 30, 40, 50))
  result <- .reflect_limited_pad(x, 3L, 0L)
  # Python: 2*x[0] - x[3:0:-1] = 20 - c(40,30,20) = c(-20,-10,0)
  expect_equal(result[1:3], c(-20, -10, 0))
})

test_that(".reflect_limited_pad right reflection: 2*x[end] - x[reversed inner]", {
  x      <- as.double(c(10, 20, 30, 40, 50))
  result <- .reflect_limited_pad(x, 0L, 3L)
  # Python: 2*x[-1] - x[-2:-5:-1] = 100 - c(40,30,20) = c(60,70,80)
  expect_equal(result[6:8], c(60, 70, 80))
})

test_that(".reflect_limited_pad zero-fills when n_pad > len(x)-1", {
  x      <- as.double(c(1, 2, 3))
  result <- .reflect_limited_pad(x, 4L, 0L)
  expect_length(result, 7L)
  expect_equal(result[1], 0.0)  # zero fill at far left
})

test_that(".reflect_limited_pad zero padding works symmetrically", {
  x      <- as.double(c(5, 10))
  result <- .reflect_limited_pad(x, 3L, 3L)
  expect_length(result, 8L)
})

# ─── .overlap_add_filter() ───────────────────────────────────────────────────

test_that(".overlap_add_filter length-1 kernel scales x by h", {
  x <- rnorm(200)
  expect_equal(.overlap_add_filter(x, 1.0), x)
  expect_equal(.overlap_add_filter(x, 2.0), 2.0 * x)
})

test_that(".overlap_add_filter returns same-length vector as input", {
  set.seed(1)
  x <- rnorm(512)
  h <- .firwin_lowpass(101L, 0.4)
  expect_length(.overlap_add_filter(x, h), 512L)
})

test_that(".overlap_add_filter attenuates stop-band content", {
  sfreq <- 256
  t     <- seq(0, 4 - 1/sfreq, by = 1/sfreq)   # 1024 samples
  # 10 Hz (pass) + 80 Hz (stop for 40 Hz low-pass)
  x     <- sin(2 * pi * 10 * t) + sin(2 * pi * 80 * t)
  h     <- .firwin_lowpass(851L, 40 / 128)
  y     <- .overlap_add_filter(x, h)
  # Drop the first 10% (edge transient), check residual std vs a pure 10 Hz signal
  n_drop  <- ceiling(0.1 * length(y))
  y_mid   <- y[(n_drop + 1):length(y)]
  ref_mid <- sin(2 * pi * 10 * t[(n_drop + 1):length(t)])
  expect_lt(sd(y_mid - ref_mid), 0.01)
})

# ─── eeg_bandpass() — input validation ───────────────────────────────────────

test_that("eeg_bandpass rejects non-eeg input", {
  expect_error(eeg_bandpass(list()), "class 'eeg'")
  expect_error(eeg_bandpass(data.frame()), "class 'eeg'")
})

test_that("eeg_bandpass requires at least one of l_freq / h_freq", {
  expect_error(eeg_bandpass(eeg1ch), "At least one of 'l_freq' or 'h_freq'")
})

test_that("eeg_bandpass rejects non-positive l_freq", {
  expect_error(eeg_bandpass(eeg1ch, l_freq = 0),  "'l_freq' must be a positive")
  expect_error(eeg_bandpass(eeg1ch, l_freq = -5), "'l_freq' must be a positive")
})

test_that("eeg_bandpass rejects h_freq at or above Nyquist", {
  expect_error(eeg_bandpass(eeg1ch, h_freq = 128), "Nyquist")
  expect_error(eeg_bandpass(eeg1ch, h_freq = 200), "Nyquist")
})

test_that("eeg_bandpass rejects l_freq >= h_freq", {
  expect_error(eeg_bandpass(eeg1ch, l_freq = 40, h_freq = 1),
               "must be less than 'h_freq'")
  expect_error(eeg_bandpass(eeg1ch, l_freq = 40, h_freq = 40),
               "must be less than 'h_freq'")
})

test_that("eeg_bandpass rejects unknown channel name", {
  expect_error(
    eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40,
                 channels = "BOGUS", verbose = FALSE),
    "Channel\\(s\\) not found"
  )
})

test_that("eeg_bandpass rejects out-of-range channel index", {
  expect_error(
    eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40,
                 channels = 99L, verbose = FALSE),
    "Channel indices out of range"
  )
})

# ─── eeg_bandpass() — return structure ───────────────────────────────────────

test_that("eeg_bandpass returns class eeg", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40, verbose = FALSE)
  expect_s3_class(out, "eeg")
})

test_that("eeg_bandpass output data has identical dimensions to input", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40, verbose = FALSE)
  expect_equal(dim(out$data), dim(eeg1ch$data))
})

test_that("eeg_bandpass appends exactly one history entry labelled 'FIR bandpass'", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40, verbose = FALSE)
  expect_length(out$preprocessing_history, 1L)
  expect_true(grepl("FIR bandpass", out$preprocessing_history[[1]]))
})

test_that("eeg_bandpass preserves channels and sampling_rate", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40, verbose = FALSE)
  expect_equal(out$channels,      eeg1ch$channels)
  expect_equal(out$sampling_rate, eeg1ch$sampling_rate)
})

test_that("eeg_bandpass channel selection by name leaves others unchanged", {
  eeg3 <- new_eeg(
    data          = matrix(rnorm(512 * 3, sd = 10), nrow = 3),
    channels      = c("Fp1", "Fp2", "Fz"),
    sampling_rate = 256L
  )
  out <- eeg_bandpass(eeg3, l_freq = 1, h_freq = 40,
                      channels = "Fp2", verbose = FALSE)
  expect_equal(out$data[1, ], eeg3$data[1, ])
  expect_equal(out$data[3, ], eeg3$data[3, ])
})

test_that("eeg_bandpass channel selection by integer index leaves others unchanged", {
  eeg3 <- new_eeg(
    data          = matrix(rnorm(512 * 3, sd = 10), nrow = 3),
    channels      = c("Fp1", "Fp2", "Fz"),
    sampling_rate = 256L
  )
  out <- eeg_bandpass(eeg3, l_freq = 1, h_freq = 40,
                      channels = 1L, verbose = FALSE)
  expect_equal(out$data[2, ], eeg3$data[2, ])
  expect_equal(out$data[3, ], eeg3$data[3, ])
})

# ─── eeg_bandpass() — MNE numerical accuracy ─────────────────────────────────

test_that("eeg_bandpass BP1 (1-40 Hz, hamming) matches MNE to 1e-10", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1.0, h_freq = 40.0,
                      fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("bp1_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

test_that("eeg_bandpass BP2 (highpass 1 Hz, hamming) matches MNE to 1e-10", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1.0, h_freq = NULL,
                      fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("bp2_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

test_that("eeg_bandpass BP3 (lowpass 40 Hz, hamming) matches MNE to 1e-10", {
  out <- eeg_bandpass(eeg1ch, l_freq = NULL, h_freq = 40.0,
                      fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("bp3_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

test_that("eeg_bandpass BP4 (mixed trans bandwidth) matches MNE to 1e-10", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1.0, h_freq = 40.0,
                      l_trans_bandwidth = "auto", h_trans_bandwidth = 5.0,
                      fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("bp4_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

test_that("eeg_bandpass BP5 (hann window) matches MNE to 1e-10", {
  out <- eeg_bandpass(eeg1ch, l_freq = 1.0, h_freq = 40.0,
                      fir_window = "hann", verbose = FALSE)
  mne <- load_ref("bp5_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

# ─── eeg_notch() — input validation ──────────────────────────────────────────

test_that("eeg_notch rejects non-eeg input", {
  expect_error(eeg_notch(list()), "class 'eeg'")
})

test_that("eeg_notch rejects non-positive freqs", {
  expect_error(eeg_notch(eeg1ch, freqs = 0),   "positive values")
  expect_error(eeg_notch(eeg1ch, freqs = -50), "positive values")
})

test_that("eeg_notch rejects empty freqs vector", {
  expect_error(eeg_notch(eeg1ch, freqs = numeric(0)), "non-empty")
})

test_that("eeg_notch rejects freqs at or above Nyquist", {
  expect_error(eeg_notch(eeg1ch, freqs = 128), "Nyquist")
  expect_error(eeg_notch(eeg1ch, freqs = 200), "Nyquist")
})

test_that("eeg_notch rejects non-positive trans_bandwidth", {
  expect_error(eeg_notch(eeg1ch, freqs = 50, trans_bandwidth = 0),  "positive numeric")
  expect_error(eeg_notch(eeg1ch, freqs = 50, trans_bandwidth = -1), "positive numeric")
})

test_that("eeg_notch rejects unknown channel name", {
  expect_error(
    eeg_notch(eeg1ch, freqs = 50, channels = "BOGUS", verbose = FALSE),
    "Channel\\(s\\) not found"
  )
})

# ─── eeg_notch() — return structure ──────────────────────────────────────────

test_that("eeg_notch returns class eeg", {
  out <- eeg_notch(eeg1ch, freqs = 50, verbose = FALSE)
  expect_s3_class(out, "eeg")
})

test_that("eeg_notch output data has identical dimensions to input", {
  out <- eeg_notch(eeg1ch, freqs = 50, verbose = FALSE)
  expect_equal(dim(out$data), dim(eeg1ch$data))
})

test_that("eeg_notch appends exactly one history entry labelled 'FIR notch'", {
  out <- eeg_notch(eeg1ch, freqs = 50, verbose = FALSE)
  expect_length(out$preprocessing_history, 1L)
  expect_true(grepl("FIR notch", out$preprocessing_history[[1]]))
})

test_that("eeg_notch preserves channels and sampling_rate", {
  out <- eeg_notch(eeg1ch, freqs = 50, verbose = FALSE)
  expect_equal(out$channels,      eeg1ch$channels)
  expect_equal(out$sampling_rate, eeg1ch$sampling_rate)
})

test_that("eeg_notch works with multiple harmonic frequencies", {
  out <- eeg_notch(eeg1ch, freqs = c(20, 50, 80), verbose = FALSE)
  expect_s3_class(out, "eeg")
  expect_equal(dim(out$data), dim(eeg1ch$data))
})

test_that("eeg_notch NULL notch_widths uses freq/200 default without error", {
  expect_no_error(eeg_notch(eeg1ch, freqs = 50, notch_widths = NULL, verbose = FALSE))
})

test_that("eeg_notch custom notch_widths scalar is accepted", {
  out <- eeg_notch(eeg1ch, freqs = 50, notch_widths = 2.0, verbose = FALSE)
  expect_s3_class(out, "eeg")
})

test_that("eeg_notch custom notch_widths vector is accepted", {
  out <- eeg_notch(eeg1ch, freqs = c(20, 50), notch_widths = c(0.5, 1.0),
                   verbose = FALSE)
  expect_s3_class(out, "eeg")
})

# ─── eeg_notch() — MNE numerical accuracy ────────────────────────────────────

test_that("eeg_notch N1 (50 Hz, hamming) matches MNE to 1e-10", {
  out <- eeg_notch(eeg1ch, freqs = 50.0, trans_bandwidth = 1.0,
                   fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("n1_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

test_that("eeg_notch N2 ([20,50,80] Hz, hamming) matches MNE to 1e-10", {
  out <- eeg_notch(eeg1ch, freqs = c(20.0, 50.0, 80.0),
                   trans_bandwidth = 1.0, fir_window = "hamming", verbose = FALSE)
  mne <- load_ref("n2_filtered.csv")
  expect_equal(length(out$data[1, ]), length(mne))
  expect_lt(max(abs(out$data[1, ] - mne)), 1e-10)
})

# ─── filter chaining ─────────────────────────────────────────────────────────

test_that("bandpass then notch accumulates two history entries", {
  bp  <- eeg_bandpass(eeg1ch, l_freq = 1, h_freq = 40, verbose = FALSE)
  out <- eeg_notch(bp, freqs = 20, verbose = FALSE)
  expect_s3_class(out, "eeg")
  expect_equal(dim(out$data), dim(eeg1ch$data))
  expect_length(out$preprocessing_history, 2L)
})
