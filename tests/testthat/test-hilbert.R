# ============================================================================
#                       Test File for hilbert.R
# ============================================================================
#
# The Hilbert transform gives, for every sample of a rhythm, two numbers: how
# big the rhythm is right now (the envelope) and where in its up-and-down cycle
# it is right now (the phase). hilbert.R replaces the chosen channels of an
# 'eeg' or 'eeg_epochs' object by one of the two.
#
# Functions tested:
#   1. apply_hilbert()      - Envelope or phase of each channel
#   2. .analytic_signal()   - Internal helper: the analytic signal through the
#                             FFT (the voltage plus a second number, per sample)
#   3. .hilbert_envelope(), .hilbert_phase()
#                           - Older internal helpers, now one-liners that call
#                             .analytic_signal(); they must give exactly what
#                             they gave before
#
# Test suites:
#   1. apply_hilbert() input validation
#   2. Known answers (waves whose envelope and phase are known)
#   3. The object that comes back
#   4. Reference numbers
#   5. Epoched data, and filter -> transform -> epoch
#   6. .analytic_signal()
#   7. Messages, history and the warning
#   8. The two older helpers
#
# The reference numbers in suite 4 were measured with an independent Python
# implementation of the Hilbert transform on the same made-up recordings that the
# fixtures below build, so no test needs Python. They are listed in
# notes/hilbert-build-spec.md, section 8.1.
#
# Author: Christos Dalamarinis
# Date: Oct 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

# The internal helper, reached through the package namespace.
analytic_signal <- eeganalysis:::.analytic_signal

# A small recording from a channels x time matrix. By default A and B are EEG
# channels and EOG is an eye channel.
make_hilbert_eeg <- function(X, types = c("eeg", "eeg", "eog"), sr = 256, bads = NULL) {
  new_eeg(data = X, channels = rownames(X), sampling_rate = sr,
          channel_types = types, bads = bads)
}

# 4 seconds at 256 Hz: 1024 samples are exactly 40 whole cycles of 10 Hz, so a
# 10 Hz wave looks the same at both ends and has no edge effect.
tt    <- (0:1023) / 256
x_sin <- 5 * sin(2 * pi * 10 * tt)
X0    <- rbind(A = x_sin, B = 2 * x_sin, EOG = x_sin)
eeg0  <- make_hilbert_eeg(X0)

# A recording with every field filled in, to check that nothing but the picked
# channels and the history changes.
eeg_rich <- new_eeg(
  data = rbind(Cz = x_sin, Fz = 2 * x_sin, VEOG = x_sin),
  channels = c("Cz", "Fz", "VEOG"), sampling_rate = 256,
  channel_types = c("eeg", "eeg", "eog"), bads = "Fz",
  events = data.frame(onset = c(300L, 700L), onset_time = (c(300, 700) - 1) / 256,
                      type = 1L, description = "Trigger: 1"),
  annotations = data.frame(onset = 1, duration = 0.5, description = "BAD_test",
                           channel = NA_character_, stringsAsFactors = FALSE),
  montage = create_montage(c("Cz", "Fz")),
  metadata = list(subject = "s01"), reference = "average",
  preprocessing_history = list("step one", "step two"))

# A small 'eeg_epochs' object: 4 channels (the 4th is an eye channel), 104
# samples per epoch, 3 epochs.
make_hilbert_epochs <- function() {
  set.seed(11)
  n  <- 1500
  ev <- data.frame(onset = c(300L, 700L, 1100L),
                   onset_time = (c(300, 700, 1100) - 1) / 256,
                   type = 1L, description = "Trigger: 1")
  cont <- new_eeg(data = matrix(rnorm(4 * n) * 10, 4),
                  channels = c("Cz", "Fz", "Pz", "Oz"), sampling_rate = 256,
                  events = ev, channel_types = c("eeg", "eeg", "eeg", "eog"))
  epoch_eeg(cont, events = "all", tmin = -0.1, tmax = 0.3,
            baseline = NULL, baseline_method = "none", verbose = FALSE)
}
ep0 <- make_hilbert_epochs()

# Distance between two phase curves, measured around the circle. A sample on the
# +pi/-pi seam can come out as +pi in one curve and -pi in the other, so the
# plain difference can be a whole turn (6.2832) for the same angle.
circ <- function(a, b) abs(Arg(exp(1i * (a - b))))

# The newest entry in the history of an object.
last_step <- function(eeg) {
  eeg$preprocessing_history[[length(eeg$preprocessing_history)]]
}

# Passes when no element differs from the reference by more than `tol`.
# (expect_equal() with a tolerance is relative, which is the wrong yardstick
# for reference values printed to a fixed number of decimals.)
expect_abs_close <- function(object, expected, tol, label = NULL) {
  stopifnot(length(object) == length(expected))
  worst <- max(abs(object - expected))
  first <- function(x) format(as.vector(x)[seq_len(min(6, length(x)))])
  expect(isTRUE(worst < tol),
         paste0(if (!is.null(label)) paste0(label, ": ") else "",
                "largest difference ", format(worst, digits = 3),
                " is not below ", format(tol, digits = 3),
                "\n  got:      ", paste(first(object), collapse = " "),
                "\n  expected: ", paste(first(expected), collapse = " ")))
}

# ============================================================================
# TEST SUITE 1: apply_hilbert() input validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Requires an 'eeg' or 'eeg_epochs' object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies apply_hilbert() refuses anything that is not an
# 'eeg' or 'eeg_epochs' object.
test_that("apply_hilbert errors when eeg_obj is not an 'eeg' or 'eeg_epochs' object", {
  msg <- "'eeg_obj' must be an object of class 'eeg' or 'eeg_epochs'."
  for (bad in list(list(), NULL, 1:3, matrix(1, 2, 2), "eeg", data.frame(a = 1))) {
    expect_error(apply_hilbert(bad), msg, fixed = TRUE,
                 info = paste(class(bad), collapse = "/"))
  }
})

# ----------------------------------------------------------------------------
# Test 1.2: Continuous data must be a matrix
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies apply_hilbert() errors when eeg_obj$data of a
# continuous object is missing or is not a channels x time matrix.
test_that("apply_hilbert errors when the data of an 'eeg' object are not a matrix", {
  msg <- "eeg_obj$data must be a numeric matrix (channels x timepoints)."

  no_data <- eeg0
  no_data$data <- NULL
  expect_error(apply_hilbert(no_data), msg, fixed = TRUE)

  flat <- eeg0
  flat$data <- as.vector(eeg0$data)
  expect_error(apply_hilbert(flat), msg, fixed = TRUE)

  cube <- eeg0
  cube$data <- array(1, c(2, 3, 4))
  expect_error(apply_hilbert(cube), msg, fixed = TRUE)

  frame <- eeg0
  frame$data <- as.data.frame(t(eeg0$data))
  expect_error(apply_hilbert(frame), msg, fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 1.3: Epoched data must be loaded and 3-dimensional
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies apply_hilbert() errors for epochs without data
# (preload = FALSE) and for epochs whose data are not channels x time x epochs.
test_that("apply_hilbert errors for epochs without data or with data that are not 3D", {
  no_data <- ep0
  no_data$data <- NULL
  expect_error(apply_hilbert(no_data),
               "epoch data not loaded - re-run epoch_eeg() with preload = TRUE.",
               fixed = TRUE)

  msg <- "eeg_obj$data must be a 3D array (channels x timepoints x epochs)."
  flat <- ep0
  flat$data <- ep0$data[, , 1]
  expect_error(apply_hilbert(flat), msg, fixed = TRUE)

  four <- ep0
  four$data <- array(1, c(2, 3, 4, 5))
  expect_error(apply_hilbert(four), msg, fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 1.4: 'output' must be exactly "envelope" or "phase"
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies output rejects other words, partial words, other
# cases, two values (unless it is the untouched default), numbers and NULL, and
# accepts both valid words and the default.
test_that("apply_hilbert validates 'output'", {
  msg <- "'output' must be \"envelope\" or \"phase\"."
  for (bad in list("power", "env", "Envelope", c("phase", "envelope"),
                   c("envelope", "phase", "power"), 1, NA, NULL, character(0))) {
    expect_error(apply_hilbert(eeg0, output = bad), msg, fixed = TRUE,
                 info = paste("output =", paste(bad, collapse = ", ")))
  }

  expect_no_error(apply_hilbert(eeg0, output = "envelope", verbose = FALSE))
  expect_no_error(apply_hilbert(eeg0, output = "phase", verbose = FALSE))
  expect_no_error(apply_hilbert(eeg0, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 1.5: 'n_fft' must be "auto", NULL or one whole number
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies n_fft rejects other words, vectors, fractions, zero,
# negative numbers, NA, NaN, Inf, logicals and lists.
test_that("apply_hilbert validates 'n_fft'", {
  msg <- "'n_fft' must be \"auto\", NULL or a single whole number."
  for (bad in list("fast", "AUTO", c(1024, 2048), 1.5, 0, -5, NA, NA_real_, NaN,
                   Inf, TRUE, list(1024))) {
    expect_error(apply_hilbert(eeg0, n_fft = bad), msg, fixed = TRUE,
                 info = paste("n_fft =", paste(unlist(bad), collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.6: 'n_fft' must be at least the number of time points
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a too-small n_fft is refused with both numbers in
# the message (for continuous and epoched data), and that exactly the number of
# time points, a bigger number and a whole number given as a double are fine.
test_that("apply_hilbert requires n_fft to be at least the number of time points", {
  expect_error(apply_hilbert(eeg0, n_fft = 500),
               "'n_fft' (500) must be at least the number of time points (1024).",
               fixed = TRUE)
  expect_error(apply_hilbert(eeg0, n_fft = 1023),
               "'n_fft' (1023) must be at least the number of time points (1024).",
               fixed = TRUE)
  expect_error(apply_hilbert(ep0, n_fft = 100),
               "'n_fft' (100) must be at least the number of time points (104).",
               fixed = TRUE)

  expect_no_error(apply_hilbert(eeg0, n_fft = 1024, verbose = FALSE))
  expect_no_error(apply_hilbert(eeg0, n_fft = 1024L, verbose = FALSE))
  expect_no_error(apply_hilbert(eeg0, n_fft = 2000, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 1.7: 'verbose' must be TRUE or FALSE
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies verbose rejects text, NA, two values, numbers and
# NULL.
test_that("apply_hilbert validates 'verbose'", {
  msg <- "'verbose' must be TRUE or FALSE."
  for (bad in list("yes", NA, c(TRUE, FALSE), 1, NULL)) {
    expect_error(apply_hilbert(eeg0, verbose = bad), msg, fixed = TRUE,
                 info = paste("verbose =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.8: 'picks' names must exist
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies unknown channel names are refused and every unknown
# name is listed (names are case-sensitive).
test_that("apply_hilbert errors for channel names that are not in the object", {
  expect_error(apply_hilbert(eeg0, picks = c("A", "Nope")),
               "These channels in 'picks' are not in eeg_obj: Nope.", fixed = TRUE)
  expect_error(apply_hilbert(eeg0, picks = c("X", "A", "Nope")),
               "not in eeg_obj: X, Nope.", fixed = TRUE)
  expect_error(apply_hilbert(eeg0, picks = "a"),
               "not in eeg_obj: a.", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 1.9: 'picks' numbers must be channel numbers
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies channel numbers outside 1..number of channels are
# refused, with the number of channels in the message.
test_that("apply_hilbert errors for channel numbers outside the recording", {
  msg <- "'picks' numbers must be between 1 and 3."
  for (bad in list(0, 4, -1, c(1, 10), 1e6)) {
    expect_error(apply_hilbert(eeg0, picks = bad), msg, fixed = TRUE,
                 info = paste("picks =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.10: 'picks' must be NULL, names or whole numbers
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies picks rejects fractions, logicals, lists, NA, NaN,
# vectors containing NA and factors.
test_that("apply_hilbert errors when picks is not NULL, names or channel numbers", {
  msg <- "'picks' must be NULL, channel names or channel numbers."
  for (bad in list(1.5, TRUE, list(), NA, NA_real_, NaN, c(1, NA), factor("A"))) {
    expect_error(apply_hilbert(eeg0, picks = bad), msg, fixed = TRUE,
                 info = paste("picks =", paste(unlist(bad), collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.11: There must be at least one channel to transform
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies empty picks, and a recording without any EEG
# channel (no picks given), are refused with the same message; and that giving
# picks explicitly is the way out. An old object without a 'channel_types'
# field behaves like one without EEG channels.
test_that("apply_hilbert errors when there is no channel to transform", {
  msg <- "No channels to transform (no channel is typed \"eeg\", or 'picks' is empty)."

  expect_error(apply_hilbert(eeg0, picks = character(0)), msg, fixed = TRUE)
  expect_error(apply_hilbert(eeg0, picks = integer(0)), msg, fixed = TRUE)
  expect_error(apply_hilbert(eeg0, picks = numeric(0)), msg, fixed = TRUE)

  no_eeg <- make_hilbert_eeg(X0, types = c("misc", "misc", "eog"))
  expect_error(apply_hilbert(no_eeg), msg, fixed = TRUE)
  expect_no_error(apply_hilbert(no_eeg, picks = "A", verbose = FALSE))

  old <- eeg0
  old$channel_types <- NULL
  expect_error(apply_hilbert(old), msg, fixed = TRUE)
  expect_no_error(apply_hilbert(old, picks = c("A", "B"), verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 1.12: Picked channels must be finite
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies NA, NaN and Inf in a picked channel are refused, that
# every such channel is named (in channel order), that the message points to
# annotate_nan(), that a non-finite value in a channel that is not picked does
# no harm, and that epoched data are checked in every epoch.
test_that("apply_hilbert errors when a picked channel has non-finite values", {
  nan_eeg <- make_hilbert_eeg(X0)
  nan_eeg$data["A", 5] <- NaN
  nan_eeg$data["B", 7] <- Inf
  expect_error(apply_hilbert(nan_eeg),
               "Channel(s) with NA or non-finite values cannot be transformed: A, B. Clean the data first (see annotate_nan()).",
               fixed = TRUE)
  expect_error(apply_hilbert(nan_eeg, picks = "A"),
               "cannot be transformed: A. Clean", fixed = TRUE)
  expect_error(apply_hilbert(nan_eeg, picks = c("B", "A")),
               "cannot be transformed: A, B. Clean", fixed = TRUE)

  na_eeg <- make_hilbert_eeg(X0)
  na_eeg$data["B", 3] <- NA
  expect_error(apply_hilbert(na_eeg), "cannot be transformed: B. Clean", fixed = TRUE)

  # a non-finite value in a channel that is not picked is fine
  expect_no_error(apply_hilbert(nan_eeg, picks = "EOG", verbose = FALSE))
  nan_eog <- make_hilbert_eeg(X0)
  nan_eog$data["EOG", 9] <- NA
  expect_no_error(apply_hilbert(nan_eog, verbose = FALSE))

  # epochs: one NA in one epoch of one channel is enough
  nan_ep <- ep0
  nan_ep$data[2, 10, 2] <- NA
  expect_error(apply_hilbert(nan_ep), "cannot be transformed: Fz.", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 1.13: The checks run in a fixed order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when two things are wrong the first check in
# the order class, data shape, output, n_fft, verbose, picks, finite values is
# the one reported.
test_that("apply_hilbert reports the first problem in a fixed order of checks", {
  expect_error(apply_hilbert(list(), output = "power"), "must be an object of class",
               fixed = TRUE)

  no_data <- eeg0
  no_data$data <- NULL
  expect_error(apply_hilbert(no_data, output = "power"), "numeric matrix", fixed = TRUE)

  expect_error(apply_hilbert(eeg0, output = "power", n_fft = "fast"),
               "'output' must be", fixed = TRUE)
  expect_error(apply_hilbert(eeg0, n_fft = "fast", verbose = "yes"),
               "'n_fft' must be", fixed = TRUE)
  expect_error(apply_hilbert(eeg0, n_fft = 10, verbose = "yes"),
               "'n_fft' (10) must be at least", fixed = TRUE)
  expect_error(apply_hilbert(eeg0, verbose = "yes", picks = "Nope"),
               "'verbose' must be", fixed = TRUE)

  nan_eeg <- make_hilbert_eeg(X0)
  nan_eeg$data["A", 5] <- NaN
  expect_error(apply_hilbert(nan_eeg, picks = c("A", "Nope")),
               "not in eeg_obj", fixed = TRUE)
})

# ============================================================================
# TEST SUITE 2: Known answers
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: A steady sine wave has a constant envelope
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a 5 microvolt, 10 Hz sine wave (40 whole cycles)
# gives an envelope of 5 at every sample (10 for the channel twice as big), and
# that the EOG channel is left exactly as it was.
test_that("the envelope of a steady sine wave is its amplitude at every sample", {
  env <- apply_hilbert(eeg0, "envelope", verbose = FALSE)
  expect_abs_close(env$data["A", ], rep(5, 1024), 1e-9, label = "channel A")
  expect_abs_close(env$data["B", ], rep(10, 1024), 1e-9, label = "channel B")
  expect_identical(env$data["EOG", ], X0["EOG", ])
})

# ----------------------------------------------------------------------------
# Test 2.2: The phase of a sine wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the phase of a sine wave starts at -pi/2 (going up
# through the middle) and follows the angle of (sin, -cos) all the way, compared
# around the circle.
test_that("the phase of a sine wave starts at -pi/2 and follows the angle of (sin, -cos)", {
  ph <- apply_hilbert(eeg0, "phase", verbose = FALSE)$data["A", ]
  expected <- Arg(complex(real = x_sin, imaginary = -5 * cos(2 * pi * 10 * tt)))
  expect_lt(abs(ph[1] + pi / 2), 1e-9)
  expect_lt(max(circ(ph, expected)), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.3: The phase of a cosine wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the phase of a cosine wave is 0 at the first sample
# (the top of a peak) and then equals 2 * pi * f * t, and its envelope is 5.
test_that("a cosine wave starts at phase 0 and its phase grows as 2 * pi * f * t", {
  cos_eeg <- make_hilbert_eeg(rbind(A = 5 * cos(2 * pi * 10 * tt),
                                    B = 5 * cos(2 * pi * 10 * tt), EOG = tt))
  ph  <- apply_hilbert(cos_eeg, "phase", verbose = FALSE)$data["A", ]
  env <- apply_hilbert(cos_eeg, "envelope", verbose = FALSE)$data["A", ]
  expect_lt(abs(ph[1]), 1e-9)
  expect_lt(max(circ(ph, 2 * pi * 10 * tt)), 1e-9)
  expect_abs_close(env, rep(5, 1024), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.4: The phase moves forward in time
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the phase of a 10 Hz wave sampled at 256 Hz grows by
# 2 * pi * 10 / 256 from one sample to the next (counting the jump at +pi -> -pi
# as part of the same turn). This fixes the direction of the angle: a wave that
# is later in its cycle has a larger phase.
test_that("the phase grows by 2 * pi * f / sampling rate per sample", {
  ph <- apply_hilbert(eeg0, "phase", verbose = FALSE)$data["A", ]
  step <- diff(ph)
  step[step < 0] <- step[step < 0] + 2 * pi
  expect_abs_close(step, rep(2 * pi * 10 / 256, 1023), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.5: A rhythm that swells and shrinks
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the envelope of a 10 Hz rhythm whose size swings
# between 2.5 and 7.5 microvolts (slowly, at 0.5 Hz) is exactly that swinging
# size.
test_that("the envelope of a swelling and shrinking rhythm is its swelling size", {
  t3   <- (0:2559) / 256
  size <- 5 + 2.5 * sin(2 * pi * 0.5 * t3)
  wave <- size * sin(2 * pi * 10 * t3)
  am   <- make_hilbert_eeg(rbind(A = wave, B = wave, EOG = t3))
  env  <- apply_hilbert(am, "envelope", verbose = FALSE)$data["A", ]
  expect_abs_close(env, size, 1e-9)
  expect_equal(range(env), c(2.5, 7.5), tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 2.6: Random data stay inside the allowed ranges
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the envelope is never negative and the phase never
# leaves -pi..pi, for noise.
test_that("the envelope is never negative and the phase stays between -pi and pi", {
  set.seed(3)
  rnd <- make_hilbert_eeg(rbind(A = rnorm(3000) * 10, B = rnorm(3000) * 10, EOG = rnorm(3000)))
  env <- apply_hilbert(rnd, "envelope", verbose = FALSE)
  ph  <- apply_hilbert(rnd, "phase", verbose = FALSE)
  expect_gte(min(env$data[1:2, ]), 0)
  expect_gte(min(ph$data[1:2, ]), -pi)
  expect_lte(max(ph$data[1:2, ]), pi)
})

# ----------------------------------------------------------------------------
# Test 2.7: Constants
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a constant +3 has envelope 3 and phase 0, and a
# constant -3 has envelope 3 and phase pi.
test_that("a constant has the envelope of its size and phase 0 (or pi when negative)", {
  const <- make_hilbert_eeg(rbind(A = rep(3, 64), B = rep(-3, 64), EOG = rep(0, 64)), sr = 64)
  env <- apply_hilbert(const, "envelope", verbose = FALSE)$data
  ph  <- apply_hilbert(const, "phase", verbose = FALSE)$data
  expect_abs_close(env[1:2, ], rep(3, 128), 1e-12)
  expect_abs_close(ph[1, ], rep(0, 64), 1e-12)
  expect_abs_close(ph[2, ], rep(pi, 64), 1e-12)
})

# ----------------------------------------------------------------------------
# Test 2.8: Scaling and flipping a wave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that multiplying a signal by a number multiplies its
# envelope by the size of that number and leaves the phase alone, and that
# flipping the signal (times -1) keeps the envelope and turns the phase by half a
# lap.
test_that("scaling a signal scales its envelope, and flipping it turns its phase by pi", {
  set.seed(8)
  x <- rnorm(2000) * 10
  trio <- make_hilbert_eeg(rbind(A = x, B = 3 * x, C = -2 * x), types = rep("eeg", 3))
  env <- apply_hilbert(trio, "envelope", verbose = FALSE)$data
  ph  <- apply_hilbert(trio, "phase", verbose = FALSE)$data

  expect_abs_close(env["B", ], 3 * env["A", ], 1e-9)
  expect_abs_close(env["C", ], 2 * env["A", ], 1e-9)
  expect_lt(max(circ(ph["B", ], ph["A", ])), 1e-9)
  expect_lt(max(circ(ph["C", ], ph["A", ] + pi)), 1e-9)
})

# ============================================================================
# TEST SUITE 3: The object that comes back
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Everything except the picked channels and the history is unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the result has the class, size and dimension names of
# the input and that every other field (bads, montage, annotations, events,
# metadata, reference, channel types, times, sampling rate) is identical, for
# both outputs.
test_that("apply_hilbert changes only the picked channels and the history", {
  keep <- setdiff(names(eeg_rich), c("data", "preprocessing_history"))
  for (out in c("envelope", "phase")) {
    r <- apply_hilbert(eeg_rich, out, verbose = FALSE)
    expect_identical(class(r), class(eeg_rich))
    expect_identical(dim(r$data), dim(eeg_rich$data))
    expect_identical(dimnames(r$data), dimnames(eeg_rich$data))
    expect_identical(unclass(r)[keep], unclass(eeg_rich)[keep], info = out)
  }
})

# ----------------------------------------------------------------------------
# Test 3.2: One step is added to the history
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies exactly one entry is appended at the end of the
# history, older entries are kept in place, and the entries have the documented
# text for both outputs.
test_that("apply_hilbert appends one step to the history and keeps the older steps", {
  r <- apply_hilbert(eeg0, "envelope", verbose = FALSE)
  expect_length(r$preprocessing_history, length(eeg0$preprocessing_history) + 1)
  expect_identical(last_step(r),
                   "apply_hilbert(): replaced 2 channel(s) by their amplitude envelope (no zero padding).")
  expect_identical(last_step(apply_hilbert(eeg0, "phase", verbose = FALSE)),
                   "apply_hilbert(): replaced 2 channel(s) by their instantaneous phase (radians) (no zero padding).")

  r2 <- apply_hilbert(eeg_rich, verbose = FALSE)
  expect_length(r2$preprocessing_history, 3)
  expect_identical(r2$preprocessing_history[1:2], list("step one", "step two"))
  expect_identical(last_step(r2),
                   "apply_hilbert(): replaced 2 channel(s) by their amplitude envelope (no zero padding).")
})

# ----------------------------------------------------------------------------
# Test 3.3: The object that was passed in is not changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies apply_hilbert() works on a copy.
test_that("apply_hilbert leaves the object that was passed in untouched", {
  before <- eeg0
  apply_hilbert(eeg0, "phase", verbose = FALSE)
  expect_identical(eeg0, before)

  before_ep <- ep0
  apply_hilbert(ep0, "envelope", verbose = FALSE)
  expect_identical(ep0, before_ep)
})

# ----------------------------------------------------------------------------
# Test 3.4: By default every EEG channel is transformed, bad ones included
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the default picks are the channels typed "eeg"
# (bads included, and 'bads' itself unchanged), and that channels of other types
# stay exactly as they were.
test_that("by default every channel typed 'eeg' is transformed, bad channels included", {
  env <- apply_hilbert(eeg0, "envelope", verbose = FALSE)
  bad_eeg <- make_hilbert_eeg(X0, bads = "B")
  r <- apply_hilbert(bad_eeg, verbose = FALSE)
  expect_identical(r$data, env$data)
  expect_identical(r$bads, "B")
  expect_abs_close(r$data["B", ], rep(10, 1024), 1e-9)

  mixed <- make_hilbert_eeg(X0, types = c("eog", "eeg", "misc"))
  rm <- apply_hilbert(mixed, verbose = FALSE)
  expect_identical(rm$data["A", ], X0["A", ])
  expect_identical(rm$data["EOG", ], X0["EOG", ])
  expect_abs_close(rm$data["B", ], rep(10, 1024), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 3.5: picks by names and by numbers
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies picks given as names (in any order) or as numbers
# (with repeats) give the same result as the default, that a repeated channel is
# transformed and counted once, and that mixing up the order does not matter.
test_that("picks by names or numbers, in any order and with repeats, work the same", {
  env <- apply_hilbert(eeg0, "envelope", verbose = FALSE)

  by_names   <- apply_hilbert(eeg0, picks = c("B", "A"), verbose = FALSE)
  by_numbers <- apply_hilbert(eeg0, picks = c(2, 1, 1), verbose = FALSE)
  by_ints    <- apply_hilbert(eeg0, picks = 1:2, verbose = FALSE)
  expect_identical(by_names$data, env$data)
  expect_identical(by_numbers$data, env$data)
  expect_identical(by_ints$data, env$data)
  expect_identical(last_step(by_numbers),
                   "apply_hilbert(): replaced 2 channel(s) by their amplitude envelope (no zero padding).")

  twice <- apply_hilbert(eeg0, picks = c("A", "A"), verbose = FALSE)
  expect_match(last_step(twice), "replaced 1 channel(s)", fixed = TRUE)
  expect_identical(twice$data["B", ], X0["B", ])
})

# ----------------------------------------------------------------------------
# Test 3.6: picks can choose channels that are not EEG
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies picks = "EOG" (or its number) transforms only that
# channel, and picks = "A" leaves every other channel alone.
test_that("picks chooses exactly the channels that are transformed", {
  r_eog <- apply_hilbert(eeg0, picks = "EOG", verbose = FALSE)
  expect_identical(r_eog$data[1:2, ], X0[1:2, ])
  expect_abs_close(r_eog$data["EOG", ], rep(5, 1024), 1e-9)
  expect_identical(apply_hilbert(eeg0, picks = 3, verbose = FALSE)$data, r_eog$data)

  r_a <- apply_hilbert(eeg0, picks = "A", verbose = FALSE)
  expect_identical(r_a$data["B", ], X0["B", ])
  expect_identical(r_a$data["EOG", ], X0["EOG", ])
  expect_abs_close(r_a$data["A", ], rep(5, 1024), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 3.7: Only real numbers come back
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the data of the result are ordinary real numbers for
# both outputs and for continuous and epoched data; the complex analytic signal
# is never stored in the object.
test_that("the data that come back are real numbers, never complex", {
  for (out in c("envelope", "phase")) {
    r <- apply_hilbert(eeg0, out, verbose = FALSE)
    expect_false(is.complex(r$data), info = out)
    expect_type(r$data, "double")

    re <- apply_hilbert(ep0, out, verbose = FALSE)
    expect_false(is.complex(re$data), info = out)
    expect_type(re$data, "double")
  }
})

# ----------------------------------------------------------------------------
# Test 3.8: The result is a valid 'eeg' object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the result can go straight into other functions of
# the package: printing it and the diagnostic summary run without an error.
test_that("the result can be used by other functions of the package", {
  r <- apply_hilbert(eeg0, "envelope", verbose = FALSE)
  expect_no_error(capture.output(print(r)))
  expect_no_error(suppressMessages(capture.output(eeg_summary(r))))
})

# ============================================================================
# TEST SUITE 4: Reference numbers
# ============================================================================
# Recordings of three channels made from one random vector (256 Hz):
#   Cz   = x            an EEG channel
#   Fz   = 0.5 * x      an EEG channel: exactly half the envelope, same phase
#   VEOG = -x           an eye channel: not picked by default, comes back as it went in
# The numbers are for channel Cz, with and without zero padding.

x1001  <- { set.seed(1001);  rnorm(1001)  * 10 }
x15361 <- { set.seed(15361); rnorm(15361) * 10 }
make_ref_eeg <- function(x) make_hilbert_eeg(rbind(Cz = x, Fz = 0.5 * x, VEOG = -x))

samples_1001  <- c(1, 2, 3, 501, 1000, 1001)
ref_1001 <- list(
  input      = c(21.886481, -1.775473, -1.852753, 8.461368, 0.923366, -0.708226),
  env_none   = c(24.779109, 17.324429, 18.656686, 8.940806, 19.885355, 11.738733),
  phase_none = c(0.488019, 1.67346, 1.670268, 0.328967, 1.524345, -1.631165),
  env_auto   = c(22.972588, 17.125577, 15.808027, 8.908653, 16.920783, 1.3661),
  phase_auto = c(0.308725, 1.674657, 1.68827, 0.318226, 1.516199, 2.115809)
)

samples_15361 <- c(1, 2, 7681, 15360, 15361)
ref_15361 <- list(
  env_none   = c(6.339527, 12.489048, 9.45765, 42.189886, 18.968029),
  phase_none = c(1.221838, -2.264664, -2.487432, 1.300124, 3.020815),
  env_auto   = c(8.386223, 13.333536, 9.477924, 38.768868, 19.62755),
  phase_auto = c(1.309356, -2.213062, -2.484647, 1.275552, 2.855515)
)

# ----------------------------------------------------------------------------
# Test 4.1: The input recipe
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the random recordings above are still the ones the
# reference numbers were measured on.
test_that("the random recordings of the reference numbers are unchanged", {
  expect_abs_close(x1001[samples_1001], ref_1001$input, 5e-6)
})

# ----------------------------------------------------------------------------
# Test 4.2: 1001 samples
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the envelope and the phase of channel Cz of a
# 1001-sample recording against the reference numbers, without padding
# (n_fft = NULL) and with the default padding (1001 -> 1024), and that Fz has half
# the envelope and the same phase.
test_that("1001 samples: envelope and phase match the reference numbers", {
  e <- make_ref_eeg(x1001)
  for (mode in c("none", "auto")) {
    nf <- if (mode == "none") NULL else "auto"
    env <- apply_hilbert(e, "envelope", n_fft = nf, verbose = FALSE)$data
    ph  <- apply_hilbert(e, "phase", n_fft = nf, verbose = FALSE)$data
    expect_abs_close(env["Cz", samples_1001], ref_1001[[paste0("env_", mode)]], 5e-6,
                     label = paste("envelope,", mode))
    expect_abs_close(ph["Cz", samples_1001], ref_1001[[paste0("phase_", mode)]], 5e-6,
                     label = paste("phase,", mode))
    expect_abs_close(env["Fz", ], 0.5 * env["Cz", ], 1e-9)
    expect_lt(max(circ(ph["Fz", ], ph["Cz", ])), 1e-9)
    expect_identical(env["VEOG", ], -x1001)
  }
})

# ----------------------------------------------------------------------------
# Test 4.3: 15361 samples
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The same for a 15361-sample recording (padding 15361 ->
# 15552).
test_that("15361 samples: envelope and phase match the reference numbers", {
  e <- make_ref_eeg(x15361)
  for (mode in c("none", "auto")) {
    nf <- if (mode == "none") NULL else "auto"
    env <- apply_hilbert(e, "envelope", n_fft = nf, verbose = FALSE)$data
    ph  <- apply_hilbert(e, "phase", n_fft = nf, verbose = FALSE)$data
    expect_abs_close(env["Cz", samples_15361], ref_15361[[paste0("env_", mode)]], 5e-6,
                     label = paste("envelope,", mode))
    expect_abs_close(ph["Cz", samples_15361], ref_15361[[paste0("phase_", mode)]], 5e-6,
                     label = paste("phase,", mode))
  }
})

# ----------------------------------------------------------------------------
# Test 4.4: What "auto" means
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies "auto" pads 1001 samples to 1024, an explicit whole
# number sets the padded length, and padding really changes the numbers.
test_that("n_fft = 'auto' pads to the next fast length and a number sets the length", {
  e <- make_ref_eeg(x1001)
  auto <- apply_hilbert(e, "envelope", n_fft = "auto", verbose = FALSE)$data["Cz", ]
  none <- apply_hilbert(e, "envelope", n_fft = NULL, verbose = FALSE)$data["Cz", ]
  big  <- apply_hilbert(e, "envelope", n_fft = 2048, verbose = FALSE)$data["Cz", ]

  expect_abs_close(auto, Mod(analytic_signal(x1001, 1024)), 1e-12)
  expect_abs_close(big, Mod(analytic_signal(x1001, 2048)), 1e-12)
  expect_abs_close(none, Mod(analytic_signal(x1001)), 1e-12)
  expect_gt(max(abs(auto - none)), 1)     # padding moves the numbers, mostly near the ends
})

# ----------------------------------------------------------------------------
# Test 4.5: A length that is already fast
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies "auto" and NULL give exactly the same numbers when
# the number of samples (15360 = 2^10 * 3 * 5) needs no padding, and that the
# history then says so.
test_that("when no padding is needed, 'auto' and NULL are exactly the same", {
  x <- x15361[1:15360]
  e <- make_hilbert_eeg(rbind(A = x, B = x, EOG = x))
  auto <- apply_hilbert(e, n_fft = "auto", verbose = FALSE)
  none <- apply_hilbert(e, n_fft = NULL, verbose = FALSE)
  expect_identical(auto$data, none$data)
  expect_match(last_step(auto), "no zero padding", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 4.6: The history says how the signal was padded
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the history entry reports the padded length for
# several recording lengths (1001 -> 1024, 4097 -> 4320, 1000 -> none) and for an
# explicit n_fft.
test_that("the history entry says how the signal was padded", {
  step_for <- function(n, ...) {
    x <- rnorm(n)
    last_step(apply_hilbert(make_hilbert_eeg(rbind(A = x, B = x, EOG = x)), verbose = FALSE, ...))
  }
  set.seed(2)
  expect_identical(step_for(1001),
                   "apply_hilbert(): replaced 2 channel(s) by their amplitude envelope (zero-padded from 1001 to 1024 samples).")
  expect_match(step_for(4097), "zero-padded from 4097 to 4320 samples", fixed = TRUE)
  expect_match(step_for(1000), "(no zero padding)", fixed = TRUE)
  expect_match(step_for(1001, n_fft = 1500), "zero-padded from 1001 to 1500 samples", fixed = TRUE)
  expect_match(step_for(1001, n_fft = NULL), "(no zero padding)", fixed = TRUE)
})

# ============================================================================
# TEST SUITE 5: Epoched data, and filter -> transform -> epoch
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Every epoch is transformed on its own
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies each epoch of each EEG channel is transformed by
# itself (104 samples padded to 108), the eye channel is left alone, the class
# and every other field are unchanged and the history says so.
test_that("apply_hilbert transforms every epoch on its own", {
  re <- apply_hilbert(ep0, "envelope", verbose = FALSE)
  manual <- ep0$data
  for (i in 1:3) for (k in 1:3) manual[i, , k] <- Mod(analytic_signal(ep0$data[i, , k], 108))

  expect_s3_class(re, "eeg_epochs")
  expect_identical(dim(re$data), c(4L, 104L, 3L))
  expect_abs_close(re$data, manual, 1e-12)
  expect_identical(re$data[4, , ], ep0$data[4, , ])

  keep <- setdiff(names(ep0), c("data", "preprocessing_history"))
  expect_identical(unclass(re)[keep], unclass(ep0)[keep])
  expect_identical(last_step(re),
                   "apply_hilbert(): replaced 3 channel(s) by their amplitude envelope (zero-padded from 104 to 108 samples).")
})

# ----------------------------------------------------------------------------
# Test 5.2: Phase of one picked channel in epochs
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies output = "phase" with picks = "Cz" changes only Cz in
# every epoch.
test_that("apply_hilbert on epochs with picks = 'Cz' changes only that channel", {
  rp <- apply_hilbert(ep0, "phase", picks = "Cz", verbose = FALSE)
  manual <- ep0$data[1, , ]
  for (k in 1:3) manual[, k] <- Arg(analytic_signal(ep0$data[1, , k], 108))
  expect_abs_close(rp$data[1, , ], manual, 1e-12)
  expect_identical(rp$data[2:4, , ], ep0$data[2:4, , ])
})

# ----------------------------------------------------------------------------
# Test 5.3: A single epoch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an object with one epoch (where R drops the epoch
# dimension when one channel is taken out) works.
test_that("apply_hilbert works for an object with a single epoch", {
  ep1 <- ep0
  ep1$data <- ep0$data[, , 1, drop = FALSE]
  r <- apply_hilbert(ep1, "envelope", verbose = FALSE)
  expect_identical(dim(r$data), c(4L, 104L, 1L))
  expect_abs_close(r$data[1, , 1], Mod(analytic_signal(ep0$data[1, , 1], 108)), 1e-12)
  expect_identical(r$data[4, , 1], ep0$data[4, , 1])
})

# ----------------------------------------------------------------------------
# Test 5.4: n_fft for epochs
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an explicit n_fft is used for every epoch (and
# reported), and NULL means no padding.
test_that("apply_hilbert on epochs uses the n_fft that was asked for", {
  r128 <- apply_hilbert(ep0, "envelope", n_fft = 128, verbose = FALSE)
  expect_abs_close(r128$data[2, , 3], Mod(analytic_signal(ep0$data[2, , 3], 128)), 1e-12)
  expect_match(last_step(r128), "zero-padded from 104 to 128 samples", fixed = TRUE)

  rnone <- apply_hilbert(ep0, "envelope", n_fft = NULL, verbose = FALSE)
  expect_abs_close(rnone$data[2, , 3], Mod(analytic_signal(ep0$data[2, , 3])), 1e-12)
  expect_match(last_step(rnone), "no zero padding", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 5.5: The example of the help page: filter, transform, epoch
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Runs the whole route of the help page on a made-up recording:
# band-pass, phase, then epoch_eeg() without a baseline. The phases stay inside
# -pi..pi; with epoch_eeg()'s default baseline (a mean is subtracted) they would
# not, which is why the help page says to switch it off for a phase.
test_that("filter, then phase, then epochs works and keeps the phase in -pi..pi", {
  set.seed(21)
  sr <- 256
  n  <- sr * 20
  t  <- (0:(n - 1)) / sr
  alpha <- (6 + 3 * sin(2 * pi * 0.2 * t)) * sin(2 * pi * 10 * t)
  X <- rbind(Cz = alpha + rnorm(n) * 2, Pz = alpha * 0.5 + rnorm(n) * 2, Oz = rnorm(n) * 2)
  ev <- data.frame(onset = c(1024L, 2048L, 3072L, 4096L),
                   onset_time = (c(1024, 2048, 3072, 4096) - 1) / sr,
                   type = 1L, description = "Trigger: 1")
  eeg <- new_eeg(data = X, channels = rownames(X), sampling_rate = sr, events = ev)

  theta_phase <- apply_hilbert(eeg_bandpass(eeg, l_freq = 4, h_freq = 8, verbose = FALSE),
                               output = "phase", verbose = FALSE)
  epochs <- epoch_eeg(theta_phase, events = "all", tmin = -0.2, tmax = 0.8,
                      baseline = NULL, baseline_method = "none", verbose = FALSE)
  expect_s3_class(epochs, "eeg_epochs")
  expect_identical(dim(epochs$data), c(3L, 257L, 4L))
  expect_lte(max(abs(epochs$data)), pi)

  with_baseline <- epoch_eeg(theta_phase, events = "all", tmin = -0.2, tmax = 0.8,
                             verbose = FALSE)
  expect_gt(max(abs(with_baseline$data)), pi)
})

# ----------------------------------------------------------------------------
# Test 5.6: Band power over time
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the envelope of an 8-12 Hz band-passed recording
# follows an alpha rhythm that grows from 2 to 6 microvolts at 10 s (noise of 3
# microvolts around it): the mean envelope is about 2.1 before and 5.9 after.
test_that("the envelope of band-passed data follows the size of the rhythm over time", {
  set.seed(10)
  sr <- 256
  t  <- (0:(sr * 20 - 1)) / sr
  size <- 2 + 4 * plogis((t - 10) * 3)
  x <- size * sin(2 * pi * 10 * t) + rnorm(length(t)) * 3
  eeg <- new_eeg(data = rbind(Oz = x), channels = "Oz", sampling_rate = sr)
  alpha <- eeg_bandpass(eeg, l_freq = 8, h_freq = 12, verbose = FALSE)
  env <- apply_hilbert(alpha, "envelope", verbose = FALSE)$data[1, ]

  early <- mean(env[t >= 1 & t <= 8])
  late  <- mean(env[t >= 12 & t <= 19])
  expect_abs_close(c(early, late), c(2.11, 5.91), 0.05)
  expect_gt(late / early, 2.5)
})

# ============================================================================
# TEST SUITE 6: .analytic_signal()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: A matrix is the same as its columns one by one
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a time x signals matrix gives, column by column,
# exactly what each column gives alone (also when padded).
test_that(".analytic_signal of a matrix equals the analytic signal of each column", {
  set.seed(5)
  M <- matrix(rnorm(300 * 4), 300)
  expect_identical(analytic_signal(M, 324),
                   sapply(1:4, function(j) analytic_signal(M[, j], 324)))
  expect_identical(analytic_signal(M),
                   sapply(1:4, function(j) analytic_signal(M[, j])))
})

# ----------------------------------------------------------------------------
# Test 6.2: Shapes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a vector gives a complex vector, a matrix a complex
# matrix of the same shape and a one-column matrix stays a matrix.
test_that(".analytic_signal keeps the shape of its input", {
  set.seed(5)
  M <- matrix(rnorm(300 * 4), 300)
  zv <- analytic_signal(M[, 1])
  expect_null(dim(zv))
  expect_length(zv, 300)
  expect_true(is.complex(zv))

  zm <- analytic_signal(M, 324)
  expect_identical(dim(zm), dim(M))
  expect_true(is.complex(zm))

  z1 <- analytic_signal(M[, 1, drop = FALSE])
  expect_identical(dim(z1), c(300L, 1L))
})

# ----------------------------------------------------------------------------
# Test 6.3: n_fft = n is the same as no padding
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies NULL and the signal's own length are the same thing.
test_that(".analytic_signal with n_fft equal to the length equals the default", {
  set.seed(5)
  x <- rnorm(300)
  expect_identical(analytic_signal(x, 300), analytic_signal(x))
})

# ----------------------------------------------------------------------------
# Test 6.4: The real part is the signal itself
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the real part of the analytic signal is the input,
# for even and odd lengths, padded or not.
test_that("the real part of the analytic signal is the signal itself", {
  set.seed(5)
  for (n in c(300, 301)) {
    x <- rnorm(n)
    expect_abs_close(Re(analytic_signal(x)), x, 1e-12, label = paste("n =", n))
    expect_abs_close(Re(analytic_signal(x, n + 24)), x, 1e-12, label = paste("padded, n =", n))
  }
})

# ----------------------------------------------------------------------------
# Test 6.5: The weights: zero frequency, Nyquist, positive and negative
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies, in the frequency domain, that the analytic signal
# keeps the zero frequency once, keeps the Nyquist frequency (even lengths) once,
# doubles every positive frequency and removes every negative one. This is the
# one thing a wrong Nyquist or zero-frequency weight would break, and it is
# nearly invisible on ordinary data.
test_that(".analytic_signal keeps DC and Nyquist once, doubles the positive frequencies and drops the negative", {
  set.seed(6)
  # even length 16: bins 1 (DC), 2:8 (positive), 9 (Nyquist), 10:16 (negative)
  x  <- rnorm(16)
  X  <- fft(x)
  Z  <- fft(analytic_signal(x))
  expect_abs_close(Z[1], X[1], 1e-12, label = "DC (even)")
  expect_abs_close(Z[2:8], 2 * X[2:8], 1e-12, label = "positive (even)")
  expect_abs_close(Z[9], X[9], 1e-12, label = "Nyquist (even)")
  expect_abs_close(Z[10:16], rep(0, 7), 1e-12, label = "negative (even)")

  # odd length 15: bins 1 (DC), 2:8 (positive), 9:15 (negative)
  x  <- rnorm(15)
  X  <- fft(x)
  Z  <- fft(analytic_signal(x))
  expect_abs_close(Z[1], X[1], 1e-12, label = "DC (odd)")
  expect_abs_close(Z[2:8], 2 * X[2:8], 1e-12, label = "positive (odd)")
  expect_abs_close(Z[9:15], rep(0, 7), 1e-12, label = "negative (odd)")
})

# ----------------------------------------------------------------------------
# Test 6.6: Very short signals
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies signals of 1, 2 and 3 samples do not fail, keep their
# length, and that one and two samples (which have no frequency to double) come
# back as they went in.
test_that(".analytic_signal copes with signals of one, two and three samples", {
  x <- c(1, -2, 4)
  for (k in 1:3) {
    z <- analytic_signal(x[seq_len(k)])
    expect_length(z, k)
    expect_abs_close(Re(z), x[seq_len(k)], 1e-12, label = paste("length", k))
  }
  expect_equal(analytic_signal(5), 5 + 0i)
  expect_equal(analytic_signal(c(3, -1)), complex(real = c(3, -1), imaginary = c(0, 0)))
})

# ----------------------------------------------------------------------------
# Test 6.7: Padding is added at the end and cut off again
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a padded analytic signal equals the first n samples
# of the analytic signal of the signal followed by zeros, has the length of the
# unpadded signal, and is not the same as the unpadded one.
test_that(".analytic_signal pads at the end and cuts the padding off", {
  set.seed(7)
  x <- rnorm(100)
  padded <- analytic_signal(x, 128)
  full   <- analytic_signal(c(x, rep(0, 28)))
  expect_length(padded, 100)
  expect_identical(padded, full[1:100])
  expect_gt(max(abs(padded - analytic_signal(x))), 1e-3)
})

# ============================================================================
# TEST SUITE 7: Messages, history and the warning
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: The one-line message
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies verbose = TRUE prints exactly one line with the
# documented text (for both outputs, with and without padding), and verbose =
# FALSE prints nothing.
test_that("apply_hilbert prints one summary line, and nothing when verbose = FALSE", {
  msgs <- capture_messages(apply_hilbert(eeg0, "envelope", verbose = TRUE))
  expect_identical(msgs, "apply_hilbert(): 2 channel(s) replaced by their amplitude envelope (no zero padding).\n")

  e1001 <- make_ref_eeg(x1001)
  msgs <- capture_messages(apply_hilbert(e1001, "phase", verbose = TRUE))
  expect_identical(msgs, "apply_hilbert(): 2 channel(s) replaced by their instantaneous phase (radians) (zero-padded from 1001 to 1024 samples).\n")

  expect_silent(apply_hilbert(eeg0, verbose = FALSE))
  expect_silent(apply_hilbert(ep0, "phase", verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 7.2: The default is verbose
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the summary is printed when verbose is not given.
test_that("apply_hilbert is verbose by default", {
  expect_message(apply_hilbert(eeg0), "apply_hilbert(): 2 channel(s) replaced", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 7.3: No warning the first time
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a fresh object, one with an empty history, and one
# whose history holds other steps do not trigger the warning for a second
# application.
test_that("apply_hilbert does not warn on an object that has not been through it", {
  expect_length(eeg0$preprocessing_history, 0)
  expect_no_warning(apply_hilbert(eeg0, verbose = FALSE))
  expect_no_warning(apply_hilbert(eeg_rich, verbose = FALSE))

  other <- eeg0
  other$preprocessing_history <- list("eeg_bandpass(): 8-12 Hz", "note: apply_hilbert(): not applied yet")
  expect_no_warning(apply_hilbert(other, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 7.4: A warning the second time
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an object whose history already has an
# apply_hilbert() step gets the documented warning (for continuous and epoched
# data), and that the second application is still carried out.
test_that("apply_hilbert warns when the object already went through apply_hilbert", {
  warn <- "apply_hilbert(): this object already went through apply_hilbert() (see preprocessing_history); transforming it again rarely makes sense."

  env <- apply_hilbert(eeg0, "envelope", verbose = FALSE)
  expect_warning(again <- apply_hilbert(env, "phase", verbose = FALSE), warn, fixed = TRUE)
  expect_length(again$preprocessing_history, 2)
  expect_identical(again$preprocessing_history[[1]], last_step(env))
  # the envelope of a steady wave is a constant, whose phase is 0
  expect_abs_close(again$data["A", ], rep(0, 1024), 1e-9)

  env_ep <- apply_hilbert(ep0, "envelope", verbose = FALSE)
  expect_warning(apply_hilbert(env_ep, "phase", verbose = FALSE), warn, fixed = TRUE)
})

# ============================================================================
# TEST SUITE 8: The two older helpers
# ============================================================================
# .hilbert_envelope() (annotations.R) and .hilbert_phase() (ica_detect.R) used to
# carry one copy of the recipe each. They are now one-liners that call
# .analytic_signal(). This is the recipe as they carried it:

old_analytic <- function(x) {
  n  <- length(x)
  Xf <- fft(x)
  h  <- numeric(n)
  if (n %% 2 == 0) {
    h[1] <- 1
    h[n / 2 + 1] <- 1
    if (n > 2) h[2:(n / 2)] <- 2
  } else {
    h[1] <- 1
    if (n > 1) h[2:((n + 1) / 2)] <- 2
  }
  fft(Xf * h, inverse = TRUE) / n
}

# ----------------------------------------------------------------------------
# Test 8.1: The envelope helper gives what it always gave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies .hilbert_envelope() gives exactly (bit for bit) the
# result of the old code, for even and odd lengths from one sample to 4097.
test_that(".hilbert_envelope gives exactly what the old code gave", {
  set.seed(99)
  for (n in c(1, 2, 3, 4, 5, 8, 9, 100, 101, 1000, 1001, 4097)) {
    x <- rnorm(n) * 10
    expect_identical(eeganalysis:::.hilbert_envelope(x), Mod(old_analytic(x)),
                     info = paste("n =", n))
  }
})

# ----------------------------------------------------------------------------
# Test 8.2: The phase helper gives what it always gave
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies .hilbert_phase() gives exactly the old result (the
# phase scaled to 0..1, one full cycle) and stays inside [0, 1].
test_that(".hilbert_phase gives exactly what the old code gave", {
  set.seed(99)
  for (n in c(1, 2, 3, 4, 5, 8, 9, 100, 101, 1000, 1001, 4097)) {
    x <- rnorm(n) * 10
    z <- old_analytic(x)
    ph <- eeganalysis:::.hilbert_phase(x)
    expect_identical(ph, (Arg(z) + pi) / (2 * pi), info = paste("n =", n))
    expect_true(all(ph >= 0 & ph <= 1), info = paste("n =", n))
  }
})

# ----------------------------------------------------------------------------
# Test 8.3: annotate_muscle() finds the same stretch as before
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: annotate_muscle() is the only user of .hilbert_envelope().
# Verifies it still finds the one stretch of a 125 Hz burst (between 15 s and
# 17 s on all 8 channels of a 40 s recording at 1024 Hz) with the same start and
# length it found before the helper became a one-liner.
test_that("annotate_muscle still finds the same stretch as before", {
  set.seed(5)
  sr <- 1024
  n  <- sr * 40
  X  <- matrix(rnorm(8 * n) * 5, 8)
  b  <- (15 * sr):(17 * sr)
  X[, b] <- X[, b] + matrix(rep(30 * sin(2 * pi * 125 * (b / sr)), each = 8), 8)
  eeg <- new_eeg(data = X, channels = paste0("Ch", 1:8), sampling_rate = sr)
  am  <- suppressMessages(annotate_muscle(eeg))

  expect_equal(nrow(am$annotations), 1)
  expect_identical(am$annotations$description, "BAD_muscle")
  expect_abs_close(am$annotations$onset, 14.98535, 1e-5)
  expect_abs_close(am$annotations$duration, 2.029297, 1e-6)
})
