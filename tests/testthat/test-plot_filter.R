# ============================================================================
#                      Test File for plot_filter.R
# ============================================================================
#
# plot_filter() shows a filter of eeg_bandpass() or eeg_notch() in three
# pictures: the impulse response (what comes out when one single spike goes
# in), the magnitude response (how much of every frequency is let through, in
# dB) and the group delay (how late every frequency comes out). It measures
# the filter by sending a spike through the real filter function, so the
# pictures must show exactly the filter that is applied to data.
#
# Functions tested:
#   1. plot_filter()                 - the three pictures and the numbers behind
#   2. .measure_impulse_response()   - sends one spike through a filter function
#   3. .frequency_response()         - magnitude and group delay of a response
#
# Test suites:
#   1. Input validation
#   2. What comes back
#   3. The picture is the filter that is applied (sine waves through the real
#      filter functions)
#   4. Reference numbers
#   5. Notch filters and several notches
#   6. Drawing
#   7. .measure_impulse_response()
#   8. .frequency_response()
#
# The reference numbers in suite 4 come from the Python reference
# implementation of the same three plots (impulse response, magnitude in dB,
# group delay), run on the same filters. They sit in this file, so the tests
# need neither Python nor any reference folder. Where the reference draws the
# impulse response of a zero-phase IIR filter from a window that is cut at its
# own ringing estimate, the tests compare with a long run instead, so no number
# here depends on that cut.
#
# Author: Christos Dalamarinis
# ============================================================================

library(testthat)
library(eeganalysis)

# Runs `expr` against a throwaway PNG device so that no plot window is shown
# during the test run, and always cleans the device up afterwards.
with_null_device <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  force(expr)
}

# Draws `expr` on a throwaway PNG device and returns the bytes of the image,
# so that two drawings can be compared.
png_bytes <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit(unlink(tmp), add = TRUE)
  tryCatch(force(expr), finally = grDevices::dev.off())
  readBin(tmp, "raw", file.size(tmp))
}

# The gain, in dB, that the real filter gives to a sine wave of frequency f:
# the amplitude of the middle half of the output, found by projecting it onto
# a sine and a cosine of that frequency.
applied_gain_db <- function(filter_fun, f, fs, seconds = 80) {
  n  <- round(seconds * fs)
  tt <- (0:(n - 1)) / fs
  s  <- new_eeg(data = matrix(sin(2 * pi * f * tt), nrow = 1L),
                channels = "Cz", sampling_rate = fs)
  y   <- filter_fun(s)$data[1, ]
  mid <- (n %/% 4):(3 * n %/% 4)
  a <- 2 * sum(y[mid] * sin(2 * pi * f * tt[mid])) / length(mid)
  b <- 2 * sum(y[mid] * cos(2 * pi * f * tt[mid])) / length(mid)
  20 * log10(sqrt(a^2 + b^2))
}

# The magnitude (dB) that plot_filter() draws for frequency f, taken from the
# impulse response it measured.
plotted_gain_db <- function(res, f) {
  .frequency_response(res$impulse, round(res$time * res$sfreq), res$sfreq,
                      f)$magnitude_db
}

# ============================================================================
#                 TEST SUITE 1: Input Validation
# ============================================================================

test_that("plot_filter rejects an x that is neither an eeg object nor a sampling rate", {
  expect_error(plot_filter("512", h_freq = 40), "class 'eeg' or a sampling rate")
  expect_error(plot_filter(list(), h_freq = 40), "class 'eeg' or a sampling rate")
  expect_error(plot_filter(c(256, 512), h_freq = 40), "class 'eeg' or a sampling rate")
  expect_error(plot_filter(NA_real_, h_freq = 40), "class 'eeg' or a sampling rate")
  expect_error(plot_filter(0, h_freq = 40), "single positive number")
  expect_error(plot_filter(-256, h_freq = 40), "single positive number")
})

test_that("plot_filter needs exactly one kind of filter", {
  expect_error(plot_filter(512), "Give 'l_freq' and/or 'h_freq'")
  expect_error(plot_filter(512, l_freq = 1, freqs = 50), "not both")
  expect_error(plot_filter(512, h_freq = 40, freqs = 50), "not both")
})

test_that("plot_filter checks method, phase, panels and fscale", {
  expect_error(plot_filter(512, h_freq = 40, method = "butter"), "should be one of")
  expect_error(plot_filter(512, h_freq = 40, phase = "forward"),
               "'phase' must be one of \"zero\", \"minimum\" for an FIR filter")
  expect_error(plot_filter(512, h_freq = 40, method = "iir", phase = "minimum"),
               "'phase' must be one of \"zero\", \"forward\" for an IIR filter")
  expect_error(plot_filter(512, h_freq = 40, panels = "phase"), "'panels' must be")
  expect_error(plot_filter(512, h_freq = 40, panels = character(0)), "'panels' must be")
  expect_error(plot_filter(512, h_freq = 40, fscale = "db"), "should be one of")
})

test_that("plot_filter checks flim, alim, tlim and dlim", {
  expect_error(plot_filter(512, h_freq = 40, flim = c(10, 1)), "two increasing frequencies")
  expect_error(plot_filter(512, h_freq = 40, flim = c(5, 5)), "two increasing frequencies")
  expect_error(plot_filter(512, h_freq = 40, flim = 5), "two increasing frequencies")
  expect_error(plot_filter(512, h_freq = 40, flim = c(0, 100)), "above 0 Hz for a log")
  expect_error(plot_filter(512, h_freq = 40, flim = c(1, 300)), "Nyquist")
  expect_error(plot_filter(512, h_freq = 40, alim = c(0, -80)), "'alim' must be")
  expect_error(plot_filter(512, h_freq = 40, alim = c(5, 5)), "'alim' must be")
  expect_error(plot_filter(512, h_freq = 40, tlim = c(1, 0)), "'tlim' and 'dlim'")
  expect_error(plot_filter(512, h_freq = 40, dlim = 1), "'tlim' and 'dlim'")
})

test_that("the checks of the real filter apply: the arguments are passed on", {
  expect_error(plot_filter(512, l_freq = 40, h_freq = 1), "must be less than 'h_freq'")
  expect_error(plot_filter(512, h_freq = 300), "Nyquist")
  expect_error(plot_filter(512, freqs = 50, method = "iir", notch_widths = -1),
               "non-negative")
  expect_error(plot_filter(512, l_freq = 1, h_freq = 40, method = "iir",
                           iir_params = list(order = 0)), "whole number")
})

test_that("a warning of the real filter is shown once, not once per measurement", {
  # the IIR notch needs a window of 64 s, so the spike is sent in twice
  n_warn <- 0L
  with_null_device(withCallingHandlers(
    plot_filter(512, freqs = 50, method = "iir", fir_window = "hann"),
    warning = function(w) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    }))
  expect_equal(n_warn, 1L)
})

# ============================================================================
#                 TEST SUITE 2: What Comes Back
# ============================================================================

test_that("plot_filter returns the curves, invisibly", {
  vis <- with_null_device(withVisible(plot_filter(512, l_freq = 1, h_freq = 40)))
  expect_false(vis$visible)
  # it draws and returns, and prints nothing (the real filter functions are
  # told to be quiet)
  expect_silent(with_null_device(plot_filter(512, l_freq = 1, h_freq = 40)))
  expect_silent(with_null_device(plot_filter(512, freqs = 50, method = "iir")))
  res <- vis$value
  expect_named(res, c("time", "impulse", "freq", "magnitude_db", "delay",
                      "sfreq", "window"))
  expect_length(res$freq, 1000)
  expect_length(res$magnitude_db, 1000)
  expect_length(res$delay, 1000)
  expect_equal(length(res$time), length(res$impulse))
  expect_equal(res$sfreq, 512)
  expect_identical(res$freq[c(1, 1000)], c(0.1, 256))
  # log-spaced: every step is the same ratio
  expect_equal(diff(log(res$freq)), rep(diff(log(res$freq))[1], 999),
               tolerance = 1e-9)
})

test_that("a zero-phase filter is symmetric around the spike; a causal one starts at it", {
  zero <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40))
  expect_equal(zero$time, -rev(zero$time))
  expect_equal(zero$impulse, rev(zero$impulse), tolerance = 1e-12)
  expect_equal(zero$time[which.max(zero$impulse)], 0)
  expect_length(zero$impulse, 1691L)                 # the 1691-tap FIR kernel
  expect_true(all(abs(zero$delay) < 1e-9, na.rm = TRUE))   # no delay at all

  for (cs in list(list(method = "fir", phase = "minimum"),
                  list(method = "iir", phase = "forward"))) {
    r <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40,
                                      method = cs$method, phase = cs$phase))
    expect_equal(r$time[1], 0)
    expect_true(all(r$time >= 0))
    # a causal filter delays the passband by a few milliseconds at least
    in_band <- r$freq > 5 & r$freq < 30
    expect_true(all(r$delay[in_band] > 0))
  }
})

test_that("the frequency axis follows fscale and flim", {
  lin <- with_null_device(plot_filter(512, h_freq = 40, fscale = "linear"))
  expect_equal(lin$freq, seq(0, 256, length.out = 1000))
  sub <- with_null_device(plot_filter(512, h_freq = 40, flim = c(10, 100)))
  expect_equal(range(sub$freq), c(10, 100))
  sub2 <- with_null_device(plot_filter(512, h_freq = 40, fscale = "linear",
                                       flim = c(30, 50)))
  expect_equal(range(sub2$freq), c(30, 50))
})

test_that("the ends of the frequency grid are exactly the limits, also where 10^log10(x) is not x", {
  # on a log grid the last point of 0.01 ... 1024 Hz comes out as 1024.0000000000002
  # (above the Nyquist frequency of 2048 Hz data) unless the ends are set exactly
  res <- with_null_device(plot_filter(2048, h_freq = 40, flim = c(0.01, 1024)))
  expect_identical(res$freq[c(1, 1000)], c(0.01, 1024))
  expect_false(anyNA(res$magnitude_db))
  lin <- with_null_device(plot_filter(512, h_freq = 40, fscale = "linear", flim = c(0.3, 45)))
  expect_identical(lin$freq[c(1, 1000)], c(0.3, 45))
})

test_that("an eeg object and its sampling rate give the same numbers", {
  eeg <- new_eeg(data = matrix(0, nrow = 1, ncol = 10), channels = "Cz",
                 sampling_rate = 256)
  a <- with_null_device(plot_filter(eeg, h_freq = 40))
  b <- with_null_device(plot_filter(256, h_freq = 40))
  expect_identical(a, b)
})

test_that("a longer window is used until the response fits, and it is reported", {
  # a narrow IIR notch rings for tens of seconds: 32 s is not enough
  r <- with_null_device(plot_filter(512, freqs = 50, method = "iir"))
  expect_gt(r$window, 32)
  # the response is cut where it falls below a ten-billionth of its peak
  expect_lt(abs(r$impulse[1]) / max(abs(r$impulse)), 1e-9)
})

test_that("a response longer than the longest window gives a warning", {
  # 0.0001 Hz at 2 Hz: the FIR kernel is about 66000 samples = 9 hours
  expect_warning(
    with_null_device(plot_filter(2, l_freq = 1e-4, h_freq = 0.9)),
    "longer than the longest window tried"
  )
})

# ============================================================================
#     TEST SUITE 3: The Picture Is the Filter That Is Applied
# ============================================================================
#
# Each filter below is applied by the real filter function to sine waves of
# several frequencies; the gain it gives them must be the gain the picture
# shows. This is what keeps plot_filter() honest if the filter code changes.

test_that("the magnitude drawn is the gain of the filter that is applied (band-pass)", {
  fs <- 256
  settings <- list(
    list(label = "FIR zero",
         plot = function() plot_filter(fs, l_freq = 1, h_freq = 40),
         filt = function(e) eeg_bandpass(e, 1, 40, verbose = FALSE),
         f = c(0.7, 1.5, 12, 38, 42, 44)),
    list(label = "IIR zero",
         plot = function() plot_filter(fs, l_freq = 1, h_freq = 40, method = "iir"),
         filt = function(e) eeg_bandpass(e, 1, 40, method = "iir", verbose = FALSE),
         f = c(0.7, 1, 2, 12, 40, 60)),
    list(label = "IIR forward",
         plot = function() plot_filter(fs, l_freq = 1, h_freq = 40, method = "iir",
                                       phase = "forward"),
         filt = function(e) eeg_bandpass(e, 1, 40, method = "iir",
                                         phase = "forward", verbose = FALSE),
         f = c(0.7, 1, 2, 12, 40, 60)),
    list(label = "FIR minimum",
         plot = function() plot_filter(fs, l_freq = 1, h_freq = 40, phase = "minimum"),
         filt = function(e) eeg_bandpass(e, 1, 40, phase = "minimum", verbose = FALSE),
         f = c(0.7, 1.5, 12, 38, 42, 44))
  )
  for (st in settings) {
    res <- with_null_device(st$plot())
    for (f in st$f) {
      plotted <- plotted_gain_db(res, f)
      applied <- applied_gain_db(st$filt, f, fs)
      expect_lt(abs(plotted - applied), 0.05,
                label = paste0(st$label, " at ", f, " Hz: drawn ", round(plotted, 3),
                               " dB, applied ", round(applied, 3), " dB"))
    }
  }
})

test_that("the magnitude drawn is the gain of the filter that is applied (notch)", {
  fs <- 512
  settings <- list(
    list(label = "FIR notch",
         plot = function() plot_filter(fs, freqs = 50),
         filt = function(e) eeg_notch(e, 50, verbose = FALSE),
         f = c(48, 49.55, 49.8, 50.3, 51)),
    list(label = "IIR notch zero",
         plot = function() plot_filter(fs, freqs = 50, method = "iir"),
         filt = function(e) eeg_notch(e, 50, method = "iir", verbose = FALSE),
         f = c(48, 49.2, 49.5, 50.6, 51)),
    list(label = "IIR notch forward",
         plot = function() plot_filter(fs, freqs = 50, method = "iir", phase = "forward"),
         filt = function(e) eeg_notch(e, 50, method = "iir", phase = "forward",
                                      verbose = FALSE),
         f = c(48, 49.2, 49.5, 50.6, 51))
  )
  for (st in settings) {
    res <- with_null_device(st$plot())
    for (f in st$f) {
      plotted <- plotted_gain_db(res, f)
      applied <- applied_gain_db(st$filt, f, fs)
      expect_lt(abs(plotted - applied), 0.05,
                label = paste0(st$label, " at ", f, " Hz: drawn ", round(plotted, 3),
                               " dB, applied ", round(applied, 3), " dB"))
    }
  }
})

test_that("the delay drawn is the delay of a narrow burst through the real causal filter", {
  # A burst of 8 Hz with a smooth bell-shaped envelope arrives late by about the
  # group delay at 8 Hz: the centre of the envelope moves by that much.
  fs <- 256
  tt <- (0:(60 * fs - 1)) / fs
  burst <- exp(-((tt - 30)^2) / (2 * 0.4^2)) * sin(2 * pi * 8 * tt)
  centre <- function(y) {
    env <- Mod(.analytic_signal(y))^2
    sum(tt * env) / sum(env)
  }
  for (cs in list(list(method = "iir", phase = "forward"),
                  list(method = "fir", phase = "minimum"))) {
    eeg <- new_eeg(data = matrix(burst, nrow = 1L), channels = "Cz", sampling_rate = fs)
    y   <- eeg_bandpass(eeg, 1, 40, method = cs$method, phase = cs$phase,
                        verbose = FALSE)$data[1, ]
    res <- with_null_device(plot_filter(fs, l_freq = 1, h_freq = 40,
                                        method = cs$method, phase = cs$phase))
    drawn <- stats::approx(log(res$freq), res$delay, xout = log(8))$y
    expect_equal(centre(y) - centre(burst), drawn, tolerance = 0.1,
                 label = paste(cs$method, cs$phase))
  }
})

# ============================================================================
#                 TEST SUITE 4: Reference Numbers
# ============================================================================
#
# Magnitude (dB) and group delay (s) at nine points of the default log
# frequency grid (0.1 Hz to 256 Hz, 1000 points; the positions are 2, 100, 250,
# 400, 500, 600, 750, 900 and 1000), and a few samples of the impulse response.
# Where the reference magnitude is below -90 dB the number is dominated by the
# last digits of both programs and is not compared.

pf_pos  <- c(2, 100, 250, 400, 500, 600, 750, 900, 1000)
pf_freq <- c(0.100788655448, 0.217647413332, 0.70713324023, 2.29746548228,
             5.03980991788, 11.0555236647, 35.9192335522, 116.701060764, 256)

expect_matches_reference <- function(res, pos, freq, mag, delay,
                                     tol_mag = 1e-3, tol_delay = 1e-3) {
  expect_equal(res$freq[pos], freq, tolerance = 1e-9)
  ok <- mag > -90
  expect_lt(max(abs(res$magnitude_db[pos][ok] - mag[ok])), tol_mag)
  okd <- ok & mag > -60 & !is.na(delay)
  expect_lt(max(abs(res$delay[pos][okd] - delay[okd])), tol_delay)
}

test_that("plot_filter matches the reference: FIR band-pass 1-40 Hz, zero phase, 512 Hz", {
  res <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40))
  expect_matches_reference(res, pf_pos, pf_freq,
    mag   = c(-34.8931875666, -21.0063950082, -1.64184220765, -0.029567534739,
              0.00311160513865, -3.36586693196e-06, 0.00249912336944,
              -73.8679239841, -79.0362061658),
    delay = rep(0, 9),
    tol_mag = 1e-6, tol_delay = 1e-9)
  # the kernel itself, at lags 0, 1, 2, 5, 10, 50, 200 and 800
  lag <- c(0, 1, 2, 5, 10, 50, 200, 800)
  ref <- c(0.173650202098, 0.164805821613, 0.139879624223, 0.0214645984548,
           -0.0231912251423, -0.000334219530168, -0.00132296850275,
           3.38896656966e-05)
  expect_lt(max(abs(res$impulse[match(lag, round(res$time * 512))] - ref)), 1e-9)
})

test_that("plot_filter matches the reference: IIR band-pass 1-40 Hz, order 4, zero phase", {
  # both passes, so the magnitude is twice the one-pass dB value and there is no delay
  res <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40, method = "iir"))
  expect_matches_reference(res, pf_pos, pf_freq,
    mag   = c(-161.161077619, -107.602999311, -25.4252329121, -0.00450386186809,
              -1.04525238231e-08, -1.21568454903e-05, -2.85652108046,
              -88.1108644331, -200),
    delay = rep(0, 9))
})

test_that("plot_filter matches the reference: IIR band-pass 1-40 Hz, order 4, causal", {
  res <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40, method = "iir",
                                      phase = "forward"))
  expect_matches_reference(res, pf_pos, pf_freq,
    mag   = c(-80.5805388096, -53.8014996555, -12.7126164561, -0.00225193093405,
              -5.22626094409e-09, -6.07842274459e-06, -1.42826054023,
              -44.0554322165, -200),
    delay = c(0.40764446365, 0.415114013888, 0.564330862493, 0.0980480599176,
              0.0273109972726, 0.0141918946462, 0.0169107981157,
              0.00150816279408, NA),
    tol_mag = 1e-3, tol_delay = 2e-4)
  # the causal response itself, at lags 0, 1, 3, 6, 10, 20, 50, 100 and 300
  lag <- c(0, 1, 3, 6, 10, 20, 50, 100, 300)
  ref <- c(0.00188238492775, 0.0126758655597, 0.0831642408409, 0.167550384533,
           0.0220654628378, -0.0142997767943, -0.0135033225582,
           -0.00165407415924, 0.000701675988842)
  expect_lt(max(abs(res$impulse[match(lag, round(res$time * 512))] - ref)), 1e-9)
})

test_that("plot_filter matches the reference: FIR minimum-phase band-pass, 512 Hz", {
  # the rebuilt minimum-phase kernel is rounding-sensitive at 0 Hz (see the help of
  # eeg_bandpass()): two correct implementations agree to about 1e-3 dB
  res <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40, phase = "minimum"))
  expect_equal(res$freq[pf_pos], pf_freq, tolerance = 1e-9)
  ref <- c(-34.9139803496, -21.0113022135, -1.64229215347, -0.0296095710595,
           0.00310296983785, -5.47881878634e-06, 0.00249899945215,
           -73.8678445726, -79.0363655241)
  expect_lt(max(abs(res$magnitude_db[pf_pos] - ref)), 5e-3)
})

# A zero-phase IIR filter is the same filter run forward and then backward, so
# its impulse response has to equal one long forward-backward run of the filter
# (here the zero-phase band-pass of suite 4, 512 Hz): the values below come from
# such a run on 64 seconds of zeros around the spike.
test_that("the zero-phase IIR impulse response is a long forward-backward run", {
  res <- with_null_device(plot_filter(512, l_freq = 1, h_freq = 40, method = "iir"))
  lag <- round(res$time * 512)
  at  <- function(m) res$impulse[match(m, lag)]
  expect_equal(at(0:3), at(0:-3))                   # symmetric around the spike
  expect_gt(max(abs(res$impulse)), 0.1)
  expect_gt(length(res$impulse), 8000L)             # rings for well over 8000 samples
  expect_lt(abs(res$impulse[1]) / max(abs(res$impulse)), 1e-8)
})

# ============================================================================
#         TEST SUITE 5: Notch Filters and Several Notches
# ============================================================================

# nine-point grid near 50 Hz: positions 780, 786, 789, 791, 792, 793, 795, 800
pn_pos  <- c(780, 786, 789, 791, 792, 793, 795, 800)
pn_freq <- c(45.4649386706, 47.6591748936, 48.7956911884, 49.568385924,
             49.9593097, 50.3533165177, 51.1506767281, 53.1997561853)

test_that("plot_filter matches the reference: FIR notch at 50 Hz, 512 Hz", {
  res <- with_null_device(plot_filter(512, freqs = 50))
  expect_matches_reference(res, pn_pos, pn_freq,
    mag   = c(0.00097037359861, 0.00864794239448, 0.0027978374378, -3.17518043216,
              -52.9266495663, -7.49130169983, 0.0166839737232, 0.00240849528839),
    delay = rep(0, 8), tol_mag = 1e-6, tol_delay = 1e-9)
})

test_that("plot_filter matches the reference: IIR notch at 50 Hz, zero phase, 512 Hz", {
  res <- with_null_device(plot_filter(512, freqs = 50, method = "iir"))
  expect_matches_reference(res, pn_pos, pn_freq,
    mag   = c(-8.0569080174e-07, -0.000191443534674, -0.0428144947543,
              -26.4396343594, -195.870649401, -39.2792026899, -0.0692614753847,
              -2.25248786714e-05),
    delay = rep(0, 8))
})

test_that("plot_filter matches the reference: IIR notch at 50 Hz, causal, 512 Hz", {
  res <- with_null_device(plot_filter(512, freqs = 50, method = "iir",
                                      phase = "forward"))
  expect_matches_reference(res, pn_pos, pn_freq,
    mag   = c(-4.02845401232e-07, -9.57217673359e-05, -0.0214072473772,
              -13.2198171797, -97.9353247006, -19.6396013449, -0.0346307376924,
              -1.12624393357e-05),
    delay = c(0.0127864067757, 0.0490688662801, 0.207654144576, 0.890292715253,
              0.666832695659, 0.792455236257, 0.229997240451, 0.0258058068934),
    tol_mag = 1e-3, tol_delay = 2e-4)
})

test_that("an IIR notch is 6 dB down at its cutoffs (3 dB for one pass) and gone at its centre", {
  # cutoffs: 50 -/+ (width / 2 + transition / 2) = 49.375 and 50.625 Hz
  edge <- function(phase) with_null_device(
    plot_filter(512, freqs = 50, method = "iir", phase = phase,
                fscale = "linear", flim = c(49.375, 50.625)))
  zero <- edge("zero")
  fwd  <- edge("forward")
  expect_equal(zero$magnitude_db[c(1, 1000)], rep(-6.0206, 2), tolerance = 1e-4)
  expect_equal(fwd$magnitude_db[c(1, 1000)],  rep(-3.0103, 2), tolerance = 1e-4)
  centre <- with_null_device(plot_filter(512, freqs = 50, method = "iir",
                                         fscale = "linear", flim = c(49.99, 50.01)))
  expect_lt(max(centre$magnitude_db), -80)
})

test_that("an IIR band-pass is 6 dB down at its cutoffs (3 dB for one pass)", {
  edge <- function(phase) with_null_device(
    plot_filter(512, l_freq = 1, h_freq = 40, method = "iir", phase = phase,
                flim = c(1, 40)))
  expect_equal(edge("zero")$magnitude_db[c(1, 1000)],    rep(-6.0206, 2), tolerance = 1e-4)
  expect_equal(edge("forward")$magnitude_db[c(1, 1000)], rep(-3.0103, 2), tolerance = 1e-4)
})

test_that("several IIR notches are the single notches one after the other", {
  # dB add up, and so do the delays of a causal filter
  one <- function(f) with_null_device(plot_filter(512, freqs = f, method = "iir",
                                                  phase = "forward"))
  both <- with_null_device(plot_filter(512, freqs = c(50, 100), method = "iir",
                                       phase = "forward"))
  a <- one(50)
  b <- one(100)
  ok <- a$magnitude_db > -60 & b$magnitude_db > -60
  expect_lt(max(abs(both$magnitude_db[ok] - (a$magnitude_db + b$magnitude_db)[ok])), 1e-4)
  okd <- ok & !is.na(a$delay) & !is.na(b$delay) & !is.na(both$delay)
  expect_lt(max(abs(both$delay[okd] - (a$delay + b$delay)[okd])), 1e-5)
})

test_that("plot_filter shows the notch widths and transition width it is given", {
  # cutoffs 50 -/+ (1 + 0.5) = 48.5 and 51.5 Hz: 6 dB down there
  res <- with_null_device(
    plot_filter(512, freqs = 50, method = "iir", notch_widths = 2,
                trans_bandwidth = 1, fscale = "linear", flim = c(48.5, 51.5)))
  expect_equal(res$magnitude_db[c(1, 1000)], rep(-6.0206, 2), tolerance = 1e-4)
  # an argument of the FIR filter is passed on: the Blackman window needs a
  # longer kernel than the default Hamming window for the same transition
  blackman <- with_null_device(plot_filter(512, freqs = 50, fir_window = "blackman"))
  hamming  <- with_null_device(plot_filter(512, freqs = 50))
  expect_gt(length(blackman$impulse), length(hamming$impulse) + 1000L)
})

# ============================================================================
#                      TEST SUITE 6: Drawing
# ============================================================================

test_that("plot_filter draws every combination of panels on both frequency scales", {
  for (pn in list("time", "magnitude", "delay", c("magnitude", "delay"),
                  c("delay", "time", "magnitude"))) {
    for (sc in c("log", "linear")) {
      expect_no_error(with_null_device(
        plot_filter(256, l_freq = 1, h_freq = 40, panels = pn, fscale = sc)))
    }
  }
  expect_no_error(with_null_device(
    plot_filter(256, freqs = 50, method = "iir", phase = "forward",
                dlim = c(-0.1, 1), tlim = c(-0.2, 1), alim = c(-40, 5))))
})

test_that("plot_filter leaves the graphics settings as it found them", {
  with_null_device({
    graphics::par(mfrow = c(2, 2))
    before <- graphics::par("mfrow", "mar")
    plot_filter(256, h_freq = 40)
    expect_equal(graphics::par("mfrow", "mar"), before)
  })
})

test_that("the drawing is repeatable and changes with the filter and the title", {
  a1 <- png_bytes(plot_filter(256, l_freq = 1, h_freq = 40))
  a2 <- png_bytes(plot_filter(256, l_freq = 1, h_freq = 40))
  b  <- png_bytes(plot_filter(256, l_freq = 1, h_freq = 30))
  m  <- png_bytes(plot_filter(256, l_freq = 1, h_freq = 40, main = "my own title"))
  expect_identical(a1, a2)
  expect_false(identical(a1, b))
  expect_false(identical(a1, m))
})

# ============================================================================
#              TEST SUITE 7: .measure_impulse_response()
# ============================================================================

# a causal moving average of 60 samples, as a "filter function" for an eeg object
# (59 zeros are put in front, so that even a short recording is longer than the
# filter; they are cut off again)
moving_average_60 <- function(e) {
  y <- as.numeric(stats::filter(c(numeric(59), e$data[1, ]), rep(1 / 60, 60),
                                sides = 1))[-(1:59)]
  e$data[1, ] <- y
  e
}

test_that(".measure_impulse_response finds the taps of a known filter", {
  m <- .measure_impulse_response(20, moving_average_60, causal = TRUE)
  expect_equal(m$h, rep(1 / 60, 60), tolerance = 1e-12)
  expect_equal(m$lags, 0:59)
  expect_true(m$contained)
  expect_equal(m$window, 32)
})

test_that(".measure_impulse_response doubles the window until the response fits", {
  m <- .measure_impulse_response(20, moving_average_60, causal = TRUE, window = 2)
  expect_equal(m$h, rep(1 / 60, 60), tolerance = 1e-12)
  expect_equal(m$window, 16)                  # 2 s -> 4 s -> 8 s -> 16 s
  expect_true(m$contained)
})

test_that(".measure_impulse_response says so when the window cannot be made long enough", {
  m <- .measure_impulse_response(20, moving_average_60, causal = TRUE, window = 2,
                                 max_window = 2)
  expect_false(m$contained)
  expect_equal(m$window, 2)
})

test_that(".measure_impulse_response of a symmetric filter is symmetric around the spike", {
  centered <- function(e) {                    # 5-sample moving average, centred
    y <- as.numeric(stats::filter(e$data[1, ], rep(1 / 5, 5), sides = 2))
    y[is.na(y)] <- 0
    e$data[1, ] <- y
    e
  }
  m <- .measure_impulse_response(20, centered, causal = FALSE)
  expect_equal(m$lags, -2:2)
  expect_equal(m$h, rep(0.2, 5), tolerance = 1e-12)
})

test_that(".measure_impulse_response stops when the filter returns nothing", {
  kill <- function(e) {
    e$data[1, ] <- 0
    e
  }
  expect_error(.measure_impulse_response(20, kill, causal = TRUE),
               "no response to a single spike")
})

# ============================================================================
#                TEST SUITE 8: .frequency_response()
# ============================================================================

test_that(".frequency_response of a pure delay: gain 1, delay of 2 samples", {
  fs <- 100
  f  <- c(1, 10, 25, 40)
  r  <- .frequency_response(c(0, 0, 1), 0:2, fs, f)
  expect_equal(r$magnitude_db, rep(0, 4), tolerance = 1e-9)
  expect_equal(r$delay, rep(2 / fs, 4), tolerance = 1e-9)
})

test_that(".frequency_response of a 4-sample moving average matches the formula", {
  fs <- 100
  f  <- c(1, 5, 12, 20, 33, 45)
  w  <- 2 * pi * f / fs
  r  <- .frequency_response(rep(1 / 4, 4), 0:3, fs, f)
  expect_equal(10 * log10(abs(sin(2 * w) / (4 * sin(w / 2)))^2), r$magnitude_db,
               tolerance = 1e-9)
  expect_equal(r$delay, rep(1.5 / fs, length(f)), tolerance = 1e-9)
})

test_that(".frequency_response of a symmetric kernel around lag 0 has no delay", {
  fs <- 100
  f  <- c(2, 10, 30, 45)
  w  <- 2 * pi * f / fs
  r  <- .frequency_response(c(1, 2, 1) / 4, -1:1, fs, f)
  expect_equal(10 * log10(((2 + 2 * cos(w)) / 4)^2), r$magnitude_db, tolerance = 1e-9)
  expect_equal(r$delay, rep(0, 4), tolerance = 1e-9)
})

test_that(".frequency_response gives NA for the delay where almost nothing is left", {
  # a 4-sample moving average has a zero at 25 Hz when the sampling rate is 100 Hz
  r <- .frequency_response(rep(1 / 4, 4), 0:3, 100, c(10, 25, 30))
  expect_true(is.na(r$delay[2]))
  expect_false(anyNA(r$delay[c(1, 3)]))
  expect_lt(r$magnitude_db[2], -100)
})

test_that(".frequency_response: the FFT route agrees with the exact sums", {
  fs <- 1000
  f  <- 10^seq(log10(0.5), log10(499), length.out = 200)
  h  <- c(rep(1 / 200, 200))                     # moving average of 200 samples
  ex <- .frequency_response(h, 0:199, fs, f, exact = TRUE)
  ff <- .frequency_response(h, 0:199, fs, f, exact = FALSE)
  # compared where the signal is not tiny: right next to a deep dip of the
  # response even a tiny difference in the numbers shows up in dB
  ok <- ex$magnitude_db > -50
  expect_lt(max(abs(ex$magnitude_db[ok] - ff$magnitude_db[ok])), 1e-3)
  # (this moving average has a deep dip every 5 Hz, a hard case: the delay of
  # 99.5 ms is still right to 0.1 ms)
  ok_d <- ex$magnitude_db > -30
  expect_lt(max(abs(ex$delay[ok_d] - ff$delay[ok_d])), 1e-4)
  # and for a response centred on lag 0 (negative lags wrap around in the FFT)
  hs <- c(1, 4, 6, 4, 1) / 16
  ex2 <- .frequency_response(hs, -2:2, fs, f, exact = TRUE)
  ff2 <- .frequency_response(hs, -2:2, fs, f, exact = FALSE)
  expect_lt(max(abs(ex2$magnitude_db - ff2$magnitude_db)), 1e-3)
})

test_that(".frequency_response: the FFT route keeps a frequency a hair above the Nyquist frequency", {
  fs <- 1000
  r  <- .frequency_response(rep(1 / 4, 4), 0:3, fs, c(10, fs / 2 * (1 + 1e-15)),
                            exact = FALSE)
  expect_false(anyNA(r$magnitude_db))
  expect_false(anyNA(r$H))
})

test_that(".frequency_response picks the FFT route only for very long responses", {
  # a 30000-sample response takes the FFT route by default, and gives the same as exact
  h  <- exp(-(0:29999) / 4000) / 4000
  f  <- c(0.5, 2, 10, 50)
  d  <- .frequency_response(h, 0:29999, 1000, f)
  ex <- .frequency_response(h, 0:29999, 1000, f, exact = TRUE)
  expect_equal(d$magnitude_db, ex$magnitude_db, tolerance = 1e-3)
})

test_that(".freq_ticks gives 1-2-5 ticks inside the limits", {
  expect_equal(.freq_ticks(c(0.1, 256), "log"),
               c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200))
  expect_equal(.freq_ticks(c(1, 100), "log"), c(1, 2, 5, 10, 20, 50, 100))
  expect_equal(.freq_ticks(c(0, 100), "linear"), pretty(c(0, 100)))
})

test_that(".filter_title describes the filter", {
  expect_equal(.filter_title("bandpass", "fir", "zero", 1, 40, NULL, 4),
               "FIR band-pass 1-40 Hz, zero phase")
  expect_equal(.filter_title("bandpass", "iir", "forward", 0.1, NULL, NULL, 8),
               "IIR (Butterworth, order 8) high-pass 0.1 Hz, causal (one pass)")
  expect_equal(.filter_title("bandpass", "fir", "minimum", NULL, 40, NULL, 4),
               "FIR low-pass 40 Hz, causal (minimum phase)")
  expect_equal(.filter_title("notch", "iir", "zero", NULL, NULL, c(50, 100), 4),
               "IIR (Butterworth, order 4) notch 50, 100 Hz, zero phase (forward and backward)")
})
