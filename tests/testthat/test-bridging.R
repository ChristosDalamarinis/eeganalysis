# ============================================================================
#                       Test File for bridging.R
# ============================================================================
#
# Gel caps can short two neighbouring electrodes together through a puddle of
# gel ("bridging"). The two channels then record almost the same signal.
# bridging.R finds such pairs and rebuilds them.
#
# Functions tested:
#   1. find_bridged_electrodes()        - Finds bridged electrode pairs
#   2. print.eeg_bridges()              - Print method of the detection result
#   3. interpolate_bridged_electrodes() - Rebuilds bridged electrodes
#   4. .bridge_pairs_index(), .bridge_ed(), .bridge_keep_windows(),
#      .bridge_cutoff(), .bridge_groups(), .bridge_centroid()
#                                        - Internal helpers
#
# Test suites:
#   1. find_bridged_electrodes() input validation
#   2. find_bridged_electrodes() results (and print.eeg_bridges())
#   3. Internal helpers
#   4. interpolate_bridged_electrodes() repair
#   5. interpolate_bridged_electrodes() input validation
#
# The reference numbers in suites 2 to 4 were measured with MNE-Python 1.12.1
# (compute_bridged_electrodes() and interpolate_bridged_electrodes()) on the
# same made-up recordings that the generator below builds, so no test needs
# Python. They are listed in notes/bridging-build-spec.md, sections 9.1 to 9.4.
#
# Author: Christos Dalamarinis
# Date: Oct 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

# 64 EEG-like channels, 60 s at 256 Hz, in microvolts: 10 smooth "brain
# sources", each seen more strongly by nearby electrodes, plus 3 uV of
# independent noise per channel (neighbouring channels are correlated, as in
# real EEG). independent = TRUE gives plain independent noise (10 uV) instead.
# set.seed() and rnorm() give the same numbers on every platform, so these
# recordings - and the reference numbers measured on them - are the same
# everywhere.
make_bridge_sim <- function(seed, n_sec = 60, sr = 256, independent = FALSE) {
  set.seed(seed)
  mont <- create_montage()
  chs  <- as.character(mont$positions$channel)
  pos  <- as.matrix(mont$positions[, c("x", "y", "z")])
  unit <- pos / sqrt(rowSums(pos^2))
  n    <- sr * n_sec
  K    <- 10
  src  <- matrix(rnorm(K * 3), K)
  src  <- src / sqrt(rowSums(src^2))
  S    <- t(apply(matrix(rnorm(K * n), K), 1,
                  function(x) as.numeric(stats::filter(x, rep(1 / 8, 8), circular = TRUE))))
  S    <- S / apply(S, 1, sd) * 8
  W    <- exp(-(acos(pmin(pmax(unit %*% t(src), -1), 1)) / 0.9)^2)
  X    <- W %*% S + matrix(rnorm(length(chs) * n), length(chs)) * 3
  if (independent) X <- matrix(rnorm(length(chs) * n), length(chs)) * 10
  rownames(X) <- chs
  X
}

# A gel bridge: every channel in members[-1] becomes a copy of members[1] plus
# a little amplifier noise (0.3 uV), over the samples in `cols`.
add_bridge <- function(X, members, seed, noise = 0.3, cols = seq_len(ncol(X))) {
  set.seed(seed)
  for (m in members[-1]) X[m, cols] <- X[members[1], cols] + rnorm(length(cols)) * noise
  X
}

# An 'eeg' object (all 64 channels typed "eeg", montage attached) built from a
# channels x time matrix.
make_bridge_eeg <- function(X, annotations = NULL, sr = 256) {
  mont <- create_montage()
  new_eeg(data = X, channels = rownames(X), sampling_rate = sr, montage = mont,
          annotations = annotations)
}

# The recordings the reference numbers were measured on:
#   S1  one bridge, Cz and CPz
#   S2  no bridge at all
#   S3  two separate bridges, Cz-CPz and P3-P5
#   S3c a chain of three: Cz, CPz and Pz all read the same gel
#   S4  independent noise: neighbours are not similar either, and no bridge
#   S5  a Cz-CPz bridge during the first 40 % of the recording only
bridge_scenarios <- function() {
  X0 <- make_bridge_sim(42)
  n  <- ncol(X0)
  list(
    S1  = add_bridge(X0, c("Cz", "CPz"), seed = 1),
    S2  = X0,
    S3  = add_bridge(add_bridge(X0, c("Cz", "CPz"), seed = 1), c("P3", "P5"), seed = 2),
    S3c = add_bridge(X0, c("Cz", "CPz", "Pz"), seed = 3),
    S4  = make_bridge_sim(43, independent = TRUE),
    S5  = add_bridge(X0, c("Cz", "CPz"), seed = 4, cols = seq_len(round(0.4 * n)))
  )
}

# Annotation of scenario S6 (used together with S1): a bad stretch from 10 s
# to 13 s, which overlaps two of the thirty 2-second windows.
ann6 <- data.frame(onset = 10, duration = 3, description = "BAD_test",
                   channel = NA_character_, stringsAsFactors = FALSE)

sc <- bridge_scenarios()

# A small 'eeg_epochs' object, to check that epoched data are refused.
make_bridge_epochs <- function() {
  ev <- data.frame(onset = c(300L, 500L, 700L),
                   onset_time = (c(300, 500, 700) - 1) / 256,
                   type = 1L,
                   description = "Trigger: 1")
  eeg <- new_eeg(data = matrix(rnorm(3 * 1000), 3),
                 channels = c("Cz", "Fz", "Pz"),
                 sampling_rate = 256, events = ev)
  epoch_eeg(eeg, events = "all", tmin = -0.1, tmax = 0.3,
            baseline = NULL, baseline_method = "none", verbose = FALSE)
}

# Detection takes about 0.3 s, so a result is computed the first time a test
# asks for it and is reused after that.
.fx <- new.env()
fx_get <- function(key, expr) {
  if (is.null(.fx[[key]])) .fx[[key]] <- expr
  .fx[[key]]
}

# Detection result of one scenario, with the default settings.
bridge_res <- function(name) {
  fx_get(paste0("res_", name),
         find_bridged_electrodes(make_bridge_eeg(sc[[name]]), verbose = FALSE))
}

# S3 after the repair, from the detection result.
repaired_s3 <- function() {
  fx_get("rep_S3",
         interpolate_bridged_electrodes(make_bridge_eeg(sc$S3), bridge_res("S3")))
}

# The pairs of S3, written by hand.
s3_pairs <- rbind(c("P3", "P5"), c("CPz", "Cz"))
s3_rebuilt <- c("P3", "P5", "CPz", "Cz")

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

# Samples 1, 1001 and 15001 of every rebuilt channel, before and after the
# repair, measured with MNE-Python forced to the same sphere centre as this
# package, c(0, 0, 0) (build spec, section 9.4). "before" also checks that the
# generator above still builds the same recordings.
repair_samples <- c(1, 1001, 15001)

ref_s3 <- list(
  P3  = list(before = c(-2.38064, 4.15947, -0.55229),
             after  = c(-5.33763, 3.90452, 1.79855)),
  P5  = list(before = c(-1.99966, 4.78109, -0.56482),
             after  = c(-3.35391, 1.88593, 1.80842)),
  CPz = list(before = c(-12.9552, 1.02996, -1.23843),
             after  = c(-13.88668, 3.9013, -2.03375)),
  Cz  = list(before = c(-12.76727, 0.68947, -1.17426),
             after  = c(-14.57901, 3.14282, -3.12197))
)

ref_s3c <- list(
  Pz  = list(before = c(-12.94245, 0.86229, -1.46741),
             after  = c(-10.06764, 4.20896, 1.76202)),
  CPz = list(before = c(-13.05585, 0.82201, -1.31387),
             after  = c(-13.81252, 3.05631, -2.3223)),
  Cz  = list(before = c(-12.76727, 0.68947, -1.17426),
             after  = c(-15.23336, 3.65446, -3.68334))
)

# The newest entry in the history of an object.
last_step <- function(eeg) {
  eeg$preprocessing_history[[length(eeg$preprocessing_history)]]
}

expect_repair_matches <- function(before, after, ref) {
  for (ch in names(ref)) {
    expect_abs_close(before$data[ch, repair_samples], ref[[ch]]$before, 5e-5,
                     label = paste(ch, "before"))
    expect_abs_close(after$data[ch, repair_samples], ref[[ch]]$after, 5e-5,
                     label = paste(ch, "after"))
  }
}

# ============================================================================
# TEST SUITE 1: find_bridged_electrodes() input validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Requires a continuous 'eeg' object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies find_bridged_electrodes() refuses anything that is
# not an 'eeg' object, and says clearly that epoched data are not supported.
test_that("find_bridged_electrodes errors when eeg_obj is not a continuous 'eeg' object", {
  expect_error(find_bridged_electrodes(list()), "class 'eeg'")
  expect_error(find_bridged_electrodes(NULL), "class 'eeg'")
  expect_error(find_bridged_electrodes(matrix(1, 2, 2)), "class 'eeg'")

  expect_error(find_bridged_electrodes(make_bridge_epochs()),
               "Epoched data are not supported")
})

# ----------------------------------------------------------------------------
# Test 1.2: The data must be a matrix
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies find_bridged_electrodes() errors when
# eeg_obj$data is missing or is not a channels x time matrix.
test_that("find_bridged_electrodes errors when eeg_obj$data is not a matrix", {
  eeg <- make_bridge_eeg(sc$S1)

  no_data <- eeg
  no_data$data <- NULL
  expect_error(find_bridged_electrodes(no_data), "numeric matrix")

  flat <- eeg
  flat$data <- as.vector(eeg$data)
  expect_error(find_bridged_electrodes(flat), "numeric matrix")
})

# ----------------------------------------------------------------------------
# Test 1.3: 'lm_cutoff' must be one positive number
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies lm_cutoff rejects zero, negative, missing, infinite,
# non-numeric and vector values.
test_that("find_bridged_electrodes validates lm_cutoff", {
  eeg <- make_bridge_eeg(sc$S1)

  for (bad in list(-1, 0, NA_real_, Inf, c(16, 5), "a", TRUE, NULL)) {
    expect_error(find_bridged_electrodes(eeg, lm_cutoff = bad),
                 "'lm_cutoff' must be a single positive number",
                 info = paste("lm_cutoff =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.4: 'epoch_threshold' must be from 0 up to (not including) 1
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies epoch_threshold rejects 1 and above, negative,
# missing, non-numeric and vector values, and accepts 0.
test_that("find_bridged_electrodes validates epoch_threshold", {
  eeg <- make_bridge_eeg(sc$S1)

  for (bad in list(1, 1.5, -0.1, NA_real_, c(0.1, 0.2), "a", NULL)) {
    expect_error(find_bridged_electrodes(eeg, epoch_threshold = bad),
                 "'epoch_threshold' must be a single number from 0 up to",
                 info = paste("epoch_threshold =", paste(bad, collapse = ", ")))
  }

  # 0 is the lowest allowed value: every pair that dips below the cutoff in
  # any window is reported
  expect_no_error(find_bridged_electrodes(eeg, epoch_threshold = 0, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 1.5: 'l_freq' and 'h_freq' must make a band
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the band needs two single finite numbers with
# 0 < l_freq < h_freq.
test_that("find_bridged_electrodes validates l_freq and h_freq", {
  eeg <- make_bridge_eeg(sc$S1)

  # lower edge at or above the upper edge, or not positive
  expect_error(find_bridged_electrodes(eeg, l_freq = 40), "0 < l_freq < h_freq")
  expect_error(find_bridged_electrodes(eeg, l_freq = 30), "0 < l_freq < h_freq")
  expect_error(find_bridged_electrodes(eeg, l_freq = 0), "0 < l_freq < h_freq")
  expect_error(find_bridged_electrodes(eeg, l_freq = -1), "0 < l_freq < h_freq")
  expect_error(find_bridged_electrodes(eeg, h_freq = 0.2), "0 < l_freq < h_freq")

  # not single finite numbers
  expect_error(find_bridged_electrodes(eeg, l_freq = NA_real_), "'l_freq' and 'h_freq'")
  expect_error(find_bridged_electrodes(eeg, l_freq = c(0.5, 1)), "'l_freq' and 'h_freq'")
  expect_error(find_bridged_electrodes(eeg, h_freq = "a"), "'l_freq' and 'h_freq'")
  expect_error(find_bridged_electrodes(eeg, h_freq = Inf), "'l_freq' and 'h_freq'")
})

# ----------------------------------------------------------------------------
# Test 1.6: 'h_freq' must be below the Nyquist frequency
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies h_freq at or above half the sampling rate is
# refused, and the message names that limit.
test_that("find_bridged_electrodes requires h_freq below the Nyquist frequency", {
  eeg <- make_bridge_eeg(sc$S1)   # 256 Hz

  # (the message is matched in full: the filter has a similar check of its own)
  expect_error(find_bridged_electrodes(eeg, h_freq = 200),
               "'h_freq' must be below the Nyquist frequency (128 Hz)", fixed = TRUE)
  expect_error(find_bridged_electrodes(eeg, h_freq = 128),
               "'h_freq' must be below the Nyquist frequency (128 Hz)", fixed = TRUE)

  # an unusable band is reported before the Nyquist limit
  expect_error(find_bridged_electrodes(eeg, l_freq = 400, h_freq = 300),
               "0 < l_freq < h_freq")
})

# ----------------------------------------------------------------------------
# Test 1.7: 'epoch_duration' must be one positive number
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies epoch_duration rejects zero, negative, missing,
# infinite, non-numeric and vector values.
test_that("find_bridged_electrodes validates epoch_duration", {
  eeg <- make_bridge_eeg(sc$S1)

  for (bad in list(0, -2, NA_real_, Inf, c(2, 4), "a", NULL)) {
    expect_error(find_bridged_electrodes(eeg, epoch_duration = bad),
                 "'epoch_duration' must be a single positive number",
                 info = paste("epoch_duration =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.8: 'verbose' must be TRUE or FALSE
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies verbose rejects anything but a single TRUE/FALSE.
test_that("find_bridged_electrodes validates verbose", {
  eeg <- make_bridge_eeg(sc$S1)

  for (bad in list("yes", NA, c(TRUE, FALSE), 1, NULL)) {
    expect_error(find_bridged_electrodes(eeg, verbose = bad),
                 "'verbose' must be TRUE or FALSE",
                 info = paste("verbose =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 1.9: Needs at least two good EEG channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the function errors when fewer than 2 channels are
# typed "eeg" and not bad (a pair needs two channels), and works with exactly 2.
test_that("find_bridged_electrodes needs at least 2 good EEG channels", {
  eeg <- make_bridge_eeg(sc$S1)

  one_left <- eeg
  one_left$bads <- eeg$channels[-1]
  expect_error(find_bridged_electrodes(one_left), "(found 1)", fixed = TRUE)

  none_left <- eeg
  none_left$bads <- eeg$channels
  expect_error(find_bridged_electrodes(none_left), "(found 0)", fixed = TRUE)

  one_typed <- eeg
  one_typed$channel_types[-1] <- "misc"
  expect_error(find_bridged_electrodes(one_typed), "(found 1)", fixed = TRUE)

  two_left <- eeg
  two_left$bads <- eeg$channels[-(1:2)]
  res <- find_bridged_electrodes(two_left, verbose = FALSE)
  expect_equal(res$channels, eeg$channels[1:2])
  expect_equal(nrow(res$ed), 1)
})

# ----------------------------------------------------------------------------
# Test 1.10: Non-finite data in the examined channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies NA, NaN and Inf in an examined channel are refused
# (the filter would spread them), while the same values in a channel that is
# not examined (bad, or not EEG) are fine.
test_that("find_bridged_electrodes rejects non-finite values only in examined channels", {
  eeg <- make_bridge_eeg(sc$S1)

  for (bad in c(NaN, NA, Inf, -Inf)) {
    broken <- eeg
    broken$data[5, 100] <- bad
    expect_error(find_bridged_electrodes(broken), "NA or non-finite",
                 info = paste("value =", bad))
  }

  in_bad <- eeg
  in_bad$bads <- "Fz"
  in_bad$data["Fz", 100] <- NaN
  res <- find_bridged_electrodes(in_bad, verbose = FALSE)
  expect_false("Fz" %in% res$channels)

  in_eog <- set_channel_types(eeg, c(Fp1 = "eog"))
  in_eog$data["Fp1", 100] <- NaN
  res <- find_bridged_electrodes(in_eog, verbose = FALSE)
  expect_false("Fp1" %in% res$channels)
})

# ----------------------------------------------------------------------------
# Test 1.11: The recording must hold at least one whole window
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a recording shorter than one window is refused,
# and that exactly one whole window is enough.
test_that("find_bridged_electrodes errors when not even one window fits", {
  short <- make_bridge_eeg(sc$S1[, 1:256])   # 1 second
  expect_error(find_bridged_electrodes(short), "(it lasts 1 s)", fixed = TRUE)
  expect_error(find_bridged_electrodes(short, epoch_duration = 1.5),
               "window of 1.5 s", fixed = TRUE)

  # one sample short of a whole 2-second window
  expect_error(find_bridged_electrodes(make_bridge_eeg(sc$S1[, 1:511])),
               "too short for even one window")

  # exactly one window
  res <- find_bridged_electrodes(make_bridge_eeg(sc$S1[, 1:512]), verbose = FALSE)
  expect_equal(res$n_windows, 1)
})

# ----------------------------------------------------------------------------
# Test 1.12: Not every window may overlap a BAD annotation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the function errors when BAD annotations cover the
# whole recording (any letter case), but a single window left is enough, and
# an annotation that does not start with "BAD" is ignored.
test_that("find_bridged_electrodes errors when every window is BAD", {
  note <- function(description, onset = 0, duration = 100) {
    data.frame(onset = onset, duration = duration, description = description,
               channel = NA_character_, stringsAsFactors = FALSE)
  }

  expect_error(find_bridged_electrodes(make_bridge_eeg(sc$S1, note("BAD_all"))),
               "Every window overlaps a BAD annotation")
  expect_error(find_bridged_electrodes(make_bridge_eeg(sc$S1, note("bad_all"))),
               "Every window overlaps a BAD annotation")

  # the last window (58 s to 60 s) is the only one left
  res <- find_bridged_electrodes(
    make_bridge_eeg(sc$S1, note("BAD_most", duration = 58)), verbose = FALSE)
  expect_equal(res$n_windows, 1)
  expect_equal(res$n_windows_dropped, 29)

  # not a BAD description: nothing is dropped
  res <- find_bridged_electrodes(
    make_bridge_eeg(sc$S1, note("EDGE all")), verbose = FALSE)
  expect_equal(res$n_windows, 30)
})

# ----------------------------------------------------------------------------
# Test 1.13: Too few small distances to estimate the cutoff
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the error raised when the first check lets the
# search through (epoch_threshold = 0) but no distance is below lm_cutoff.
test_that("find_bridged_electrodes errors when there is nothing to estimate the cutoff from", {
  eeg <- make_bridge_eeg(sc$S1)
  expect_error(find_bridged_electrodes(eeg, lm_cutoff = 0.001, epoch_threshold = 0),
               "Too few small electrical distances")
})

# ----------------------------------------------------------------------------
# Test 1.14: The checks run in the documented order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when two things are wrong at once, the first
# check in the documented order is the one reported.
test_that("find_bridged_electrodes reports the first failing check", {
  eeg <- make_bridge_eeg(sc$S1)

  # class before data before the settings
  expect_error(find_bridged_electrodes(list(), lm_cutoff = -1), "class 'eeg'")
  no_data <- eeg
  no_data$data <- NULL
  expect_error(find_bridged_electrodes(no_data, lm_cutoff = -1), "numeric matrix")

  # settings before the channels and the data
  one_left <- eeg
  one_left$bads <- eeg$channels[-1]
  expect_error(find_bridged_electrodes(one_left, verbose = "yes"),
               "'verbose' must be TRUE or FALSE")

  # channels before non-finite values
  one_left$data[1, 100] <- NaN
  expect_error(find_bridged_electrodes(one_left), "At least 2 good EEG channels")

  # non-finite values before the length of the recording
  short <- make_bridge_eeg(sc$S1[, 1:256])
  short$data[5, 100] <- NaN
  expect_error(find_bridged_electrodes(short), "NA or non-finite")

  # length of the recording before the BAD annotations
  short_bad <- make_bridge_eeg(
    sc$S1[, 1:256],
    data.frame(onset = 0, duration = 100, description = "BAD_all",
               channel = NA_character_, stringsAsFactors = FALSE))
  expect_error(find_bridged_electrodes(short_bad), "too short")
})

# ============================================================================
# TEST SUITE 2: find_bridged_electrodes() results
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Scenario S1, the planted bridge is found
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the Cz-CPz bridge of S1 is the only pair found,
# reported in channel order, with the values measured with MNE-Python.
test_that("find_bridged_electrodes finds the planted Cz-CPz bridge (S1)", {
  res <- bridge_res("S1")

  expect_s3_class(res, "eeg_bridges")
  expect_equal(nrow(res$bridged), 1)

  # pairs follow the channel order of the object: CPz comes before Cz
  expect_equal(res$bridged$channel_1, "CPz")
  expect_equal(res$bridged$channel_2, "Cz")
  expect_equal(res$bridged$fraction_below, 1)
  expect_abs_close(res$bridged$median_ed, 0.02309, 1e-5)

  # the density estimate is almost flat at the bottom of the valley, so the
  # cutoff is only compared to 3 decimals
  expect_abs_close(res$local_minimum, 1.915696, 1e-3)

  expect_equal(res$n_windows, 30)
  expect_equal(res$n_windows_dropped, 0)
  expect_equal(length(res$channels), 64)
})

# ----------------------------------------------------------------------------
# Test 2.2: The structure of the result
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the fields, their order, their types and the
# settings stored with the result.
test_that("find_bridged_electrodes returns the documented structure", {
  res <- bridge_res("S1")

  expect_named(res, c("bridged", "ed", "ed_pairs", "local_minimum", "n_windows",
                      "n_windows_dropped", "channels", "lm_cutoff",
                      "epoch_threshold", "l_freq", "h_freq", "epoch_duration"))
  expect_named(res$bridged,
               c("channel_1", "channel_2", "fraction_below", "median_ed"))
  expect_named(res$ed_pairs, c("channel_1", "channel_2"))
  expect_true(is.matrix(res$ed))
  expect_equal(nrow(res$ed_pairs), nrow(res$ed))
  expect_equal(res$channels, make_bridge_eeg(sc$S1)$channels)

  # the default settings are stored
  expect_equal(res$lm_cutoff, 16)
  expect_equal(res$epoch_threshold, 0.5)
  expect_equal(res$l_freq, 0.5)
  expect_equal(res$h_freq, 30)
  expect_equal(res$epoch_duration, 2)
})

# ----------------------------------------------------------------------------
# Test 2.3: The table of electrical distances (S1)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the size of the distance table, which row belongs
# to which pair, and numbers measured with MNE-Python: the bridged pair's
# distances, how many distances are below 16, and the lowest medians.
test_that("find_bridged_electrodes gives the electrical distances of S1", {
  res <- bridge_res("S1")

  # 64 channels give 2016 pairs; 60 s of data give thirty 2-second windows
  expect_equal(dim(res$ed), c(2016, 30))

  row <- which(res$ed_pairs$channel_1 == "CPz" & res$ed_pairs$channel_2 == "Cz")
  expect_equal(row, 1504)
  expect_abs_close(res$ed[1504, 1:3], c(0.02657, 0.02190, 0.02563), 1e-5)

  expect_equal(sum(res$ed < 16), 7027)   # of 60480 values, 234.2 per window

  med <- apply(res$ed, 1, median)
  lowest <- order(med)[1:4]
  expect_equal(paste(res$ed_pairs$channel_1[lowest],
                     res$ed_pairs$channel_2[lowest], sep = " - "),
               c("CPz - Cz", "C3 - CP3", "C5 - CP5", "TP7 - CP5"))
  expect_abs_close(med[lowest], c(0.0231, 4.9301, 4.9308, 5.0753), 1e-4)
  expect_abs_close(median(med), 64.49, 0.01)
  expect_abs_close(max(res$ed), 285.76, 0.01)
  expect_gte(min(res$ed), 0)
})

# ----------------------------------------------------------------------------
# Test 2.4: Row order of the pairs
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the rows of 'ed' follow the pairs in order of the
# first channel and then the second: (1,2) (1,3) ... (1,64) (2,3) ... (63,64).
test_that("find_bridged_electrodes lists the pairs in channel order", {
  res <- bridge_res("S1")
  ch <- res$channels
  pairs <- res$ed_pairs

  expect_equal(c(pairs$channel_1[1], pairs$channel_2[1]), ch[c(1, 2)])
  expect_equal(c(pairs$channel_1[63], pairs$channel_2[63]), ch[c(1, 64)])
  expect_equal(c(pairs$channel_1[64], pairs$channel_2[64]), ch[c(2, 3)])
  expect_equal(c(pairs$channel_1[2016], pairs$channel_2[2016]), ch[c(63, 64)])

  # every pair appears once
  expect_false(anyDuplicated(paste(pairs$channel_1, pairs$channel_2)) > 0)
})

# ----------------------------------------------------------------------------
# Test 2.5: The distances equal an independent calculation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Recomputes three distances by hand from the filtered data
# (variance of the difference of two channels in one window, in uV^2) and
# compares them with the table: the bridged pair, a neighbouring pair and a
# distant pair, in three different windows.
test_that("find_bridged_electrodes distances match a calculation by hand", {
  eeg <- make_bridge_eeg(sc$S1)
  res <- bridge_res("S1")
  filtered <- eeg_bandpass(eeg, l_freq = 0.5, h_freq = 30, verbose = FALSE)$data

  # variance with n in the denominator (R's var() uses n - 1)
  pop_var <- function(v) mean((v - mean(v))^2)
  by_hand <- function(a, b, window) {
    cols <- (window - 1) * 512 + 1:512      # 2 s at 256 Hz
    pop_var(filtered[a, cols] - filtered[b, cols])
  }
  # the pair as it is listed: the first channel comes first in the channel order
  table_value <- function(a, b, window) {
    row <- which(res$ed_pairs$channel_1 == a & res$ed_pairs$channel_2 == b)
    stopifnot(length(row) == 1)
    res$ed[row, window]
  }

  expect_equal(table_value("CPz", "Cz", 1), by_hand("CPz", "Cz", 1), tolerance = 1e-9)
  expect_equal(table_value("C3", "CP3", 3), by_hand("C3", "CP3", 3), tolerance = 1e-9)
  expect_equal(table_value("Oz", "Fz", 30), by_hand("Oz", "Fz", 30), tolerance = 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.6: S2, similar neighbours are not bridges
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a recording with no bridge reports nothing, even
# though 233 of its 2016 pairs (mostly neighbouring electrodes, which are
# naturally similar) have a typical distance below 16 uV^2 - which is why the
# cutoff is found from the data and not simply set to 16.
test_that("find_bridged_electrodes finds nothing in a recording without a bridge (S2)", {
  res <- bridge_res("S2")

  expect_equal(nrow(res$bridged), 0)
  expect_s3_class(res$bridged, "data.frame")
  expect_named(res$bridged,
               c("channel_1", "channel_2", "fraction_below", "median_ed"))
  expect_type(res$bridged$median_ed, "double")

  # no pile of distances near zero: the valley is at the very bottom
  expect_lt(res$local_minimum, 1e-3)

  med <- apply(res$ed, 1, median)
  expect_equal(sum(med < 16), 233)
  expect_abs_close(min(med), 4.9301, 1e-4)
  expect_abs_close(sum(res$ed < 16) / ncol(res$ed), 235.9, 0.05)
})

# ----------------------------------------------------------------------------
# Test 2.7: S4, independent noise stops at the first check
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when hardly any distance is below lm_cutoff,
# the search stops: no cutoff (NA), no pairs.
test_that("find_bridged_electrodes stops at the first check when nothing is small (S4)", {
  res <- bridge_res("S4")

  expect_true(is.na(res$local_minimum))
  expect_equal(nrow(res$bridged), 0)
  expect_equal(sum(res$ed < 16), 0)
  expect_abs_close(min(res$ed), 28.59, 0.01)
})

# ----------------------------------------------------------------------------
# Test 2.8: S3 and S3c, several bridges
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies two separate bridges are both found, and a chain of
# three electrodes gives all three pairs. The order of the rows follows the
# channel order: "Pz - CPz", not the alphabetical "CPz - Pz".
test_that("find_bridged_electrodes finds two separate bridges (S3) and a chain (S3c)", {
  pair_names <- function(res) {
    paste(res$bridged$channel_1, res$bridged$channel_2, sep = " - ")
  }

  r3 <- bridge_res("S3")
  expect_equal(pair_names(r3), c("P3 - P5", "CPz - Cz"))
  expect_equal(r3$bridged$fraction_below, c(1, 1))
  expect_abs_close(r3$bridged$median_ed, c(0.022009, 0.023091), 1e-5)
  expect_abs_close(r3$local_minimum, 1.986222, 1e-3)

  r3c <- bridge_res("S3c")
  expect_equal(pair_names(r3c), c("Pz - CPz", "Pz - Cz", "CPz - Cz"))
  expect_equal(r3c$bridged$fraction_below, c(1, 1, 1))
  expect_abs_close(r3c$bridged$median_ed, c(0.045764, 0.023284, 0.023421), 1e-5)
  expect_abs_close(r3c$local_minimum, 2.038316, 1e-3)
})

# ----------------------------------------------------------------------------
# Test 2.9: S5, a bridge that is only there part of the time
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The Cz-CPz bridge of S5 exists for the first 40 % of the
# recording (12 of 30 windows). With the default epoch_threshold of 0.5 it is
# missed; with 0.3 it is found, with fraction_below 0.4. A pair needs MORE than
# epoch_threshold of the windows, so a threshold of exactly 0.4 still misses it.
test_that("find_bridged_electrodes needs a bridge in more than epoch_threshold of the windows (S5)", {
  eeg <- make_bridge_eeg(sc$S5)

  res <- bridge_res("S5")
  expect_equal(nrow(res$bridged), 0)
  expect_abs_close(res$local_minimum, 1.827988, 1e-3)

  res <- find_bridged_electrodes(eeg, epoch_threshold = 0.3, verbose = FALSE)
  expect_equal(nrow(res$bridged), 1)
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))
  expect_equal(res$bridged$fraction_below, 0.4)

  res <- find_bridged_electrodes(eeg, epoch_threshold = 0.39, verbose = FALSE)
  expect_equal(nrow(res$bridged), 1)

  res <- find_bridged_electrodes(eeg, epoch_threshold = 0.4, verbose = FALSE)
  expect_equal(nrow(res$bridged), 0)
})

# ----------------------------------------------------------------------------
# Test 2.10: S6, windows that overlap a BAD annotation are left out
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a BAD annotation from 10 s to 13 s drops the two
# windows it touches (10-12 s and 12-14 s, the 6th and 7th), and nothing else:
# the other 28 columns of the table are the same as without the annotation.
test_that("find_bridged_electrodes leaves out windows that overlap BAD annotations (S6)", {
  res <- find_bridged_electrodes(make_bridge_eeg(sc$S1, ann6), verbose = FALSE)

  expect_equal(res$n_windows, 28)
  expect_equal(res$n_windows_dropped, 2)
  expect_equal(ncol(res$ed), 28)
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))
  expect_abs_close(res$local_minimum, 1.913432, 1e-3)

  # the annotation only decides which windows are used; the filter still runs
  # over the whole recording
  expect_equal(res$ed, bridge_res("S1")$ed[, -c(6, 7)])
})

# ----------------------------------------------------------------------------
# Test 2.11: Which annotations count
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a description starting with "BAD" counts in any
# letter case, other descriptions are ignored, and an object with no
# annotations at all works.
test_that("find_bridged_electrodes only drops windows for annotations starting with BAD", {
  note <- function(description) {
    data.frame(onset = 10, duration = 3, description = description,
               channel = NA_character_, stringsAsFactors = FALSE)
  }

  res <- find_bridged_electrodes(make_bridge_eeg(sc$S1, note("bad_lowercase")),
                                 verbose = FALSE)
  expect_equal(res$n_windows_dropped, 2)

  res <- find_bridged_electrodes(make_bridge_eeg(sc$S1, note("EDGE x")),
                                 verbose = FALSE)
  expect_equal(res$n_windows_dropped, 0)

  no_notes <- make_bridge_eeg(sc$S1)
  no_notes$annotations <- NULL
  res <- find_bridged_electrodes(no_notes, verbose = FALSE)
  expect_equal(res$n_windows, 30)
})

# ----------------------------------------------------------------------------
# Test 2.12: Annotations written by annotate_amplitude()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the windows dropped are the ones that overlap the
# rows annotate_amplitude() really writes (onset and duration in seconds from
# the start), so the two functions agree on the annotation table.
test_that("find_bridged_electrodes respects annotations written by annotate_amplitude()", {
  eeg <- make_bridge_eeg(sc$S1)
  samples <- 5000:5200                       # about 19.5 s to 20.3 s
  eeg$data["Fz", samples] <- rep(c(300, -300), length.out = length(samples))
  eeg <- annotate_amplitude(eeg, peak = 200)
  expect_gt(nrow(eeg$annotations), 0)

  res <- find_bridged_electrodes(eeg, verbose = FALSE)

  # which of the thirty windows [t0, t0 + 2) overlap a row of the table
  note <- eeg$annotations
  t0 <- (0:29) * 2
  overlaps <- vapply(t0, function(s) {
    any(note$onset < s + 2 & note$onset + note$duration > s)
  }, logical(1))

  expect_equal(sum(overlaps), 2)             # the artifact straddles 20 s
  expect_equal(res$n_windows_dropped, sum(overlaps))
  expect_equal(res$n_windows, 30 - sum(overlaps))
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))
})

# ----------------------------------------------------------------------------
# Test 2.13: Which channels are examined
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies bad channels and channels not typed "eeg" are left
# out, and that leaving channels out does not change the distances of the pairs
# that remain (each channel is filtered on its own).
test_that("find_bridged_electrodes examines only good EEG channels", {
  eeg <- make_bridge_eeg(sc$S1)
  eeg$bads <- "Fz"
  eeg <- set_channel_types(eeg, c(Fp1 = "eog"))

  res <- find_bridged_electrodes(eeg, verbose = FALSE)
  expect_equal(length(res$channels), 62)
  expect_equal(nrow(res$ed), 1891)           # 62 * 61 / 2
  expect_false(any(c("Fz", "Fp1") %in% res$channels))
  expect_false(any(c("Fz", "Fp1") %in%
                     c(res$ed_pairs$channel_1, res$ed_pairs$channel_2)))
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))

  full <- bridge_res("S1")
  key <- function(pairs) paste(pairs$channel_1, pairs$channel_2)
  expect_equal(res$ed, full$ed[match(key(res$ed_pairs), key(full$ed_pairs)), ])
})

# ----------------------------------------------------------------------------
# Test 2.14: The object is not changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies detection leaves the object alone: nothing is
# written to bads (bridged electrodes still carry brain signal) or to the
# history.
test_that("find_bridged_electrodes does not modify the object", {
  eeg <- make_bridge_eeg(sc$S1)
  before <- eeg

  res <- find_bridged_electrodes(eeg, verbose = FALSE)

  expect_equal(nrow(res$bridged), 1)
  expect_identical(eeg, before)
  expect_equal(eeg$bads, character(0))
  expect_length(eeg$preprocessing_history, 0)
})

# ----------------------------------------------------------------------------
# Test 2.15: The reference does not matter
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the claim in the help: a re-reference subtracts the
# same signal from both channels of a pair, so the distances are the same
# before and after (here: common average, and everything re-referenced to Fz).
test_that("find_bridged_electrodes gives the same distances after re-referencing", {
  eeg <- make_bridge_eeg(sc$S1)
  ref <- bridge_res("S1")

  average <- eeg
  average$data <- sweep(eeg$data, 2, colMeans(eeg$data))
  res <- find_bridged_electrodes(average, verbose = FALSE)
  expect_abs_close(res$ed, ref$ed, 1e-8)
  expect_equal(res$bridged[, 1:2], ref$bridged[, 1:2])

  to_fz <- eeg
  to_fz$data <- sweep(eeg$data, 2, eeg$data["Fz", ])
  res <- find_bridged_electrodes(to_fz, verbose = FALSE)
  expect_abs_close(res$ed, ref$ed, 1e-8)
  expect_equal(res$bridged[, 1:2], ref$bridged[, 1:2])
})

# ----------------------------------------------------------------------------
# Test 2.16: Window length and other settings
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies epoch_duration sets the number of windows (the
# partial window at the end is dropped), the bridge is still found, and every
# setting is stored with the result.
test_that("find_bridged_electrodes honours epoch_duration and the other settings", {
  eeg <- make_bridge_eeg(sc$S1)

  res <- find_bridged_electrodes(eeg, epoch_duration = 4, verbose = FALSE)
  expect_equal(res$n_windows, 15)
  expect_equal(ncol(res$ed), 15)
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))

  # 60 s / 7 s = 8.57 windows: eight whole ones, the rest is dropped
  res <- find_bridged_electrodes(eeg, epoch_duration = 7, verbose = FALSE)
  expect_equal(res$n_windows, 8)
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))

  # a window length that is not a whole number of samples is rounded to the
  # nearest sample: 0.2 s is 51.2 samples (51), 0.3 s is 76.8 samples (77)
  res <- find_bridged_electrodes(eeg, epoch_duration = 0.2, verbose = FALSE)
  expect_equal(res$n_windows, floor(15360 / 51))
  res <- find_bridged_electrodes(eeg, epoch_duration = 0.3, verbose = FALSE)
  expect_equal(res$n_windows, floor(15360 / 77))

  res <- find_bridged_electrodes(eeg, lm_cutoff = 10, epoch_threshold = 0.4,
                                 l_freq = 1, h_freq = 40, epoch_duration = 4,
                                 verbose = FALSE)
  expect_equal(res$lm_cutoff, 10)
  expect_equal(res$epoch_threshold, 0.4)
  expect_equal(res$l_freq, 1)
  expect_equal(res$h_freq, 40)
  expect_equal(res$epoch_duration, 4)
  expect_equal(c(res$bridged$channel_1, res$bridged$channel_2), c("CPz", "Cz"))
})

# ----------------------------------------------------------------------------
# Test 2.17: 'lm_cutoff' below every distance
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that a tiny lm_cutoff (no distance below it) stops
# at the first check: no data-driven cutoff, no pairs.
test_that("find_bridged_electrodes stops at the first check for a tiny lm_cutoff", {
  eeg <- make_bridge_eeg(sc$S1)
  res <- find_bridged_electrodes(eeg, lm_cutoff = 0.001, verbose = FALSE)

  expect_true(is.na(res$local_minimum))
  expect_equal(nrow(res$bridged), 0)
  expect_equal(res$lm_cutoff, 0.001)
})

# ----------------------------------------------------------------------------
# Test 2.18: The console message
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the exact one-line summary for a found bridge, a
# chain, a recording with dropped windows and one that stops at the first
# check, and that verbose = FALSE prints nothing.
test_that("find_bridged_electrodes prints a one-line summary when verbose", {
  eeg <- make_bridge_eeg(sc$S1)

  expect_equal(
    capture_messages(find_bridged_electrodes(eeg)),
    paste0("find_bridged_electrodes(): 64 EEG channel(s), 2016 pair(s), ",
           "30 window(s) of 2 s. Data-driven cutoff: 1.916 uV^2. ",
           "1 bridged pair(s): CPz - Cz.\n"))

  expect_equal(
    capture_messages(find_bridged_electrodes(make_bridge_eeg(sc$S3c))),
    paste0("find_bridged_electrodes(): 64 EEG channel(s), 2016 pair(s), ",
           "30 window(s) of 2 s. Data-driven cutoff: 2.038 uV^2. ",
           "3 bridged pair(s): Pz - CPz, Pz - Cz, CPz - Cz.\n"))

  expect_equal(
    capture_messages(find_bridged_electrodes(make_bridge_eeg(sc$S1, ann6))),
    paste0("find_bridged_electrodes(): 64 EEG channel(s), 2016 pair(s), ",
           "28 window(s) of 2 s (2 dropped because of BAD annotations). ",
           "Data-driven cutoff: 1.913 uV^2. 1 bridged pair(s): CPz - Cz.\n"))

  expect_equal(
    capture_messages(find_bridged_electrodes(make_bridge_eeg(sc$S4))),
    paste0("find_bridged_electrodes(): 64 EEG channel(s), 2016 pair(s), ",
           "30 window(s) of 2 s. Too few small distances to suspect bridging. ",
           "No bridged electrodes found.\n"))

  expect_silent(find_bridged_electrodes(eeg, verbose = FALSE))
})

# ----------------------------------------------------------------------------
# Test 2.19: print.eeg_bridges() layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the printed report line by line for the chain of
# three (S3c), and that print() returns the object invisibly.
test_that("print.eeg_bridges shows the data, the cutoffs and the pairs", {
  res <- bridge_res("S3c")

  out <- capture.output(print(res))
  rule <- paste0(strrep("=", 70), " ")
  expect_equal(out, c(
    "",
    rule,
    "Electrode Bridging Check",
    rule,
    "",
    "DATA:",
    "  EEG channels examined: 64 (2016 pairs)",
    "  Windows used:          30 of 2 s",
    "  Filter:                0.5 - 30 Hz",
    "",
    "CUTOFFS:",
    "  lm_cutoff:             16 uV^2",
    "  Data-driven cutoff:    2.038 uV^2",
    "  epoch_threshold:       0.5",
    "",
    "BRIDGED PAIRS:",
    "  Pz       - CPz      below the cutoff in 100% of windows (median distance 0.046 uV^2)",
    "  Pz       - Cz       below the cutoff in 100% of windows (median distance 0.023 uV^2)",
    "  CPz      - Cz       below the cutoff in 100% of windows (median distance 0.023 uV^2)",
    ""))

  shown <- NULL
  capture.output(shown <- withVisible(print(res)))
  expect_false(shown$visible)
  expect_identical(shown$value, res)
})

# ----------------------------------------------------------------------------
# Test 2.20: print.eeg_bridges() special cases
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the lines that change: no pairs found, no cutoff
# needed, dropped windows, and a pair found in only part of the windows.
test_that("print.eeg_bridges handles no pairs, no cutoff, dropped windows and partial pairs", {
  out <- capture.output(print(bridge_res("S2")))
  expect_true("  none found" %in% out)

  out <- capture.output(print(bridge_res("S4")))
  expect_true("  Data-driven cutoff:    not needed (too few small distances)" %in% out)
  expect_true("  none found" %in% out)

  res <- find_bridged_electrodes(make_bridge_eeg(sc$S1, ann6), verbose = FALSE)
  out <- capture.output(print(res))
  expect_true("  Windows used:          28 of 2 s (2 dropped: BAD annotation)" %in% out)

  res <- find_bridged_electrodes(make_bridge_eeg(sc$S5), epoch_threshold = 0.3,
                                 verbose = FALSE)
  out <- capture.output(print(res))
  expect_true(paste0("  CPz      - Cz       below the cutoff in  40% of windows ",
                     "(median distance 5.728 uV^2)") %in% out)
})

# ============================================================================
# TEST SUITE 3: Internal helpers
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: .bridge_pairs_index() lists every pair once, in the right order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies all pairs i < j are listed, ordered by i and then
# j (the loop order of MNE-Python) - not the column-by-column order that
# upper.tri() would give.
test_that(".bridge_pairs_index lists all pairs i < j, ordered by i then j", {
  idx <- eeganalysis:::.bridge_pairs_index(4)
  expect_equal(colnames(idx), c("i", "j"))
  expect_equal(unname(idx),
               cbind(c(1, 1, 1, 2, 2, 3), c(2, 3, 4, 3, 4, 4)))

  expect_equal(dim(eeganalysis:::.bridge_pairs_index(2)), c(1, 2))

  idx64 <- eeganalysis:::.bridge_pairs_index(64)
  expect_equal(nrow(idx64), 2016)
  expect_true(all(idx64[, "i"] < idx64[, "j"]))
  expect_false(anyDuplicated(idx64) > 0)
  expect_equal(unname(idx64),
               unname(idx64[order(idx64[, "i"], idx64[, "j"]), ]))

  # same pairs as upper.tri(), but not in its order
  upper <- which(upper.tri(diag(64)), arr.ind = TRUE)
  expect_setequal(paste(idx64[, "i"], idx64[, "j"]),
                  paste(upper[, "row"], upper[, "col"]))
  expect_false(identical(unname(as.matrix(idx64)), unname(upper)))
})

# ----------------------------------------------------------------------------
# Test 3.2: .bridge_groups() merges pairs that share a channel
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies pairs sharing an electrode end up in one group,
# separate pairs stay separate, a pair joining two groups merges them, and the
# direction a pair is written in does not matter.
test_that(".bridge_groups merges pairs that share a channel", {
  # each group as a sorted "A+B+C" string, the groups in sorted order, so the
  # order they were met in does not matter
  as_sets <- function(groups) {
    sort(vapply(groups, function(g) paste(sort(g), collapse = "+"), character(1)))
  }

  g <- eeganalysis:::.bridge_groups(
    rbind(c("Cz", "CPz"), c("P3", "P5"), c("CPz", "Pz")))
  expect_equal(as_sets(g), sort(c("CPz+Cz+Pz", "P3+P5")))

  # two pairs that touch nothing stay two groups
  g <- eeganalysis:::.bridge_groups(rbind(c("A", "B"), c("C", "D")))
  expect_equal(as_sets(g), c("A+B", "C+D"))

  # a pair that joins two existing groups merges them
  g <- eeganalysis:::.bridge_groups(rbind(c("A", "B"), c("C", "D"), c("B", "C")))
  expect_equal(as_sets(g), "A+B+C+D")

  # A-B and B-A are the same pair
  g <- eeganalysis:::.bridge_groups(rbind(c("A", "B"), c("B", "A")))
  expect_equal(as_sets(g), "A+B")

  # no pairs, no groups
  expect_length(eeganalysis:::.bridge_groups(matrix(character(0), 0, 2)), 0)
})

# ----------------------------------------------------------------------------
# Test 3.3: .bridge_centroid() places the virtual electrode
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the position of the virtual electrode against
# values measured with MNE-Python's _find_centroid_sphere: the mean position
# pushed back onto the sphere at the group's mean radius. Two electrodes on
# opposite sides of the head have no usable centre (NA).
test_that(".bridge_centroid puts the virtual electrode on the sphere between the group", {
  mont <- create_montage()
  pos <- function(chs) {
    as.matrix(mont$positions[match(chs, mont$positions$channel), c("x", "y", "z")])
  }
  norm <- function(v) sqrt(sum(v^2))

  two <- unname(eeganalysis:::.bridge_centroid(pos(c("Cz", "CPz"))))
  expect_abs_close(two, c(0, -17.34124, 86.19615), 1e-4)
  expect_abs_close(norm(two), 87.92323, 1e-4)
  # ... which is the mean radius of the two electrodes
  expect_abs_close(norm(two), mean(sqrt(rowSums(pos(c("Cz", "CPz"))^2))), 1e-9)

  three <- unname(eeganalysis:::.bridge_centroid(pos(c("Cz", "CPz", "Pz"))))
  expect_abs_close(three, c(0, -34.13659, 80.94243), 1e-4)

  # one electrode: the centre is the electrode itself
  expect_abs_close(unname(eeganalysis:::.bridge_centroid(pos("Cz"))),
                   unname(pos("Cz")[1, ]), 1e-9)

  # only the shape of the positions matters: metres instead of millimetres
  small <- unname(eeganalysis:::.bridge_centroid(pos(c("Cz", "CPz")) / 1000))
  expect_abs_close(small * 1000, two, 1e-9)

  # opposite sides of the head
  opposite <- eeganalysis:::.bridge_centroid(rbind(c(0, 0, 87), c(0, 0, -87)))
  expect_equal(unname(opposite), c(NA_real_, NA_real_, NA_real_))
})

# ----------------------------------------------------------------------------
# Test 3.4: .bridge_keep_windows() at the edges of an annotation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies, on seven boundary cases measured with MNE-Python,
# which of thirty 2-second windows a BAD annotation removes. A window spans
# [start, start + 2) seconds, an annotation [onset, onset + duration), and they
# overlap only if each starts before the other ends.
test_that(".bridge_keep_windows drops exactly the windows an annotation overlaps", {
  starts <- (0:29) * 512 + 1                  # first sample of each window
  kept <- function(onset, duration, description = "BAD_x") {
    note <- data.frame(onset = onset, duration = duration,
                       description = description, stringsAsFactors = FALSE)
    eeganalysis:::.bridge_keep_windows(note, starts, 512, 256)
  }

  # onset, duration, number of windows kept
  cases <- list(
    c(10, 3, 28),          # 10-13 s: touches the windows 10-12 s and 12-14 s
    c(12, 1, 29),          # 12-13 s: inside one window
    c(8, 2, 29),           # 8-10 s: ends exactly where the next window starts
    c(11.998, 0.001, 29),  # in the 4 ms after the last sample of 10-12 s
    c(10, 0.001, 29),      # right at the start of a window
    c(1, 0.5, 29),         # inside the first window
    c(-1, 3, 29))          # starts before the recording
  for (case in cases) {
    expect_equal(sum(kept(case[1], case[2])), case[3],
                 info = paste("onset", case[1], "duration", case[2]))
  }

  # which windows go: the 6th and 7th for the first case
  expect_equal(which(!kept(10, 3)), c(6, 7))

  # "BAD" in any letter case; other descriptions are ignored
  expect_equal(sum(kept(10, 3, "bad_lowercase")), 28)
  expect_equal(sum(kept(10, 3, "EDGE x")), 30)

  # several rows: two BAD ones count, the other one does not
  several <- data.frame(onset = c(2, 20, 40), duration = c(1, 1, 1),
                        description = c("BAD_a", "EDGE", "BAD_b"),
                        stringsAsFactors = FALSE)
  expect_equal(which(!eeganalysis:::.bridge_keep_windows(several, starts, 512, 256)),
               c(2, 21))

  # no annotations at all
  expect_true(all(eeganalysis:::.bridge_keep_windows(NULL, starts, 512, 256)))
  empty <- data.frame(onset = numeric(0), duration = numeric(0),
                      description = character(0), stringsAsFactors = FALSE)
  expect_true(all(eeganalysis:::.bridge_keep_windows(empty, starts, 512, 256)))
})

# ----------------------------------------------------------------------------
# Test 3.5: .bridge_ed() computes the variance of the difference
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the distances against a pair-by-pair calculation
# (each window's own mean removed, variance with n in the denominator), that a
# constant offset between two channels changes nothing, that two copies of a
# signal are at distance 0, and that scaling the data by k scales the distances
# by k^2.
test_that(".bridge_ed equals the variance of the difference, pair by pair", {
  set.seed(5)
  x <- matrix(rnorm(5 * 1024), 5) + 3         # 5 channels, two windows of 512
  ed <- eeganalysis:::.bridge_ed(x, c(1, 513), 512)

  expect_equal(dim(ed), c(10, 2))             # 10 pairs, 2 windows

  idx <- eeganalysis:::.bridge_pairs_index(5)
  pop_var <- function(v) mean((v - mean(v))^2)
  by_hand <- sapply(1:2, function(w) {
    cols <- (w - 1) * 512 + 1:512
    apply(idx, 1, function(p) pop_var(x[p[1], cols] - x[p[2], cols]))
  })
  expect_abs_close(ed, by_hand, 1e-12)

  # a constant added to one channel is removed with the window mean
  shifted <- x
  shifted[2, ] <- shifted[2, ] + 7
  expect_abs_close(eeganalysis:::.bridge_ed(shifted, c(1, 513), 512), ed, 1e-9)

  # two copies of one signal
  copies <- rbind(x[1, 1:512], x[1, 1:512])
  expect_abs_close(eeganalysis:::.bridge_ed(copies, 1, 512), 0, 1e-12)

  # scaling the data by 3 scales the distances by 9
  expect_abs_close(eeganalysis:::.bridge_ed(3 * x, c(1, 513), 512), 9 * ed, 1e-9)
})

# ----------------------------------------------------------------------------
# Test 3.6: .bridge_ed() never returns a negative distance
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: var(a) + var(b) - 2 cov(a, b) can round to a hair below zero
# when two signals are almost exact copies. Here each of 20 large signals gets a
# copy scaled by 1 + 1e-9: the true distance is about 1e-10, far below the
# rounding noise of the formula (about 1e-7), so roughly half of these pairs
# come out slightly negative without the clamp. The helper must return 0 or more.
test_that(".bridge_ed clamps rounding noise at zero", {
  set.seed(11)
  signals <- matrix(rnorm(20 * 512, sd = 1e4), 20)       # one window of 512
  data <- rbind(signals, signals * (1 + 1e-9))           # rows 21-40 are the copies
  ed <- eeganalysis:::.bridge_ed(data, 1, 512)

  expect_gte(min(ed), 0)

  # the 20 pairs (signal, its copy) are at a distance of about 0
  idx <- eeganalysis:::.bridge_pairs_index(40)
  twin <- idx[, "i"] <= 20 & idx[, "j"] == idx[, "i"] + 20
  expect_equal(sum(twin), 20)
  expect_lt(max(ed[twin, ]), 1e-3)
})

# ----------------------------------------------------------------------------
# Test 3.7: .bridge_cutoff() first check and the cutoff
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the first check (NA when, on average, fewer than
# epoch_threshold distances per window are below lm_cutoff), that the cutoff
# found on S1 matches MNE-Python, and that it does not depend on the order of
# the windows.
test_that(".bridge_cutoff returns NA when too few distances are small, else the valley", {
  ed <- bridge_res("S1")$ed             # 234.2 distances below 16 per window

  # nothing below lm_cutoff
  expect_true(is.na(eeganalysis:::.bridge_cutoff(
    matrix(c(20, 30, 40), ncol = 1), 16, 0.5)))
  expect_true(is.na(eeganalysis:::.bridge_cutoff(ed, 0.001, 0.5)))

  # 234.2 per window is below 300 but not below 200
  expect_true(is.na(eeganalysis:::.bridge_cutoff(ed, 16, 300)))
  expect_false(is.na(eeganalysis:::.bridge_cutoff(ed, 16, 200)))

  cutoff <- eeganalysis:::.bridge_cutoff(ed, 16, 0.5)
  expect_abs_close(cutoff, 1.915696, 1e-3)

  # the distances are pooled, so the order of the windows is irrelevant
  set.seed(1)
  shuffled <- ed[, sample(ncol(ed))]
  expect_abs_close(eeganalysis:::.bridge_cutoff(shuffled, 16, 0.5), cutoff, 1e-6)
})

# ----------------------------------------------------------------------------
# Test 3.8: .bridge_cutoff() needs at least two different small distances
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the error for one small distance, and for two small
# distances that are equal (no spread to build a density from).
test_that(".bridge_cutoff errors when there are fewer than 2 different small distances", {
  expect_error(
    eeganalysis:::.bridge_cutoff(matrix(c(1, 20, 30), ncol = 1), 16, 0.5),
    "found 1 below lm_cutoff", fixed = TRUE)
  expect_error(
    eeganalysis:::.bridge_cutoff(matrix(c(5, 5, 30), ncol = 1), 16, 0.5),
    "Too few small electrical distances")
})

# ============================================================================
# TEST SUITE 4: interpolate_bridged_electrodes() repair
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Two separate bridges (S3)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the four bridged channels are rebuilt to the values
# measured with MNE-Python, and nothing else changes: not the other channels,
# not the other fields, not the input object. One entry is added to the
# history.
test_that("interpolate_bridged_electrodes rebuilds two separate bridges (S3)", {
  eeg <- make_bridge_eeg(sc$S3)
  before <- eeg
  out <- interpolate_bridged_electrodes(eeg, bridge_res("S3"))

  changed <- rownames(eeg$data)[rowSums(abs(out$data - eeg$data)) > 0]
  expect_setequal(changed, s3_rebuilt)
  expect_repair_matches(eeg, out, ref_s3)
  expect_abs_close(max(abs(out$data - eeg$data)), 11.0, 0.05)

  # every other field is the same
  for (field in setdiff(names(eeg), c("data", "preprocessing_history"))) {
    expect_identical(out[[field]], eeg[[field]], info = field)
  }
  expect_identical(dimnames(out$data), dimnames(eeg$data))
  expect_s3_class(out, "eeg")

  # the input object was not touched
  expect_identical(eeg, before)

  expect_length(out$preprocessing_history, 1)
  expect_equal(
    last_step(out),
    paste0("interpolate_bridged_electrodes(): rebuilt 4 bridged channel(s) in ",
           "2 group(s) (P3+P5, CPz+Cz) from 60 good channel(s) plus ",
           "2 virtual electrode(s)."))
})

# ----------------------------------------------------------------------------
# Test 4.2: Every way of naming the pairs gives the same result
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the detection result, a matrix, a data frame (text
# or factor columns), pairs written the other way round, pairs in another
# order, and a pair listed twice all give exactly the same repair.
test_that("interpolate_bridged_electrodes accepts the pairs in any form", {
  eeg <- make_bridge_eeg(sc$S3)
  base <- repaired_s3()

  as_matrix  <- interpolate_bridged_electrodes(eeg, s3_pairs)
  as_text_df <- interpolate_bridged_electrodes(
    eeg, data.frame(a = c("P3", "CPz"), b = c("P5", "Cz"), stringsAsFactors = FALSE))
  as_factors <- interpolate_bridged_electrodes(
    eeg, data.frame(a = factor(c("P3", "CPz")), b = factor(c("P5", "Cz"))))
  reversed   <- interpolate_bridged_electrodes(
    eeg, rbind(c("P5", "P3"), c("Cz", "CPz")))
  other_order <- interpolate_bridged_electrodes(
    eeg, rbind(c("CPz", "Cz"), c("P3", "P5")))
  twice      <- interpolate_bridged_electrodes(
    eeg, rbind(c("P3", "P5"), c("CPz", "Cz"), c("P5", "P3")))

  expect_identical(as_matrix, base)
  expect_identical(as_text_df, base)
  expect_identical(as_factors, base)
  expect_identical(reversed, base)
  expect_identical(other_order, base)
  expect_identical(twice, base)
})

# ----------------------------------------------------------------------------
# Test 4.3: A chain of three (S3c)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the three pairs of a chain become one group of
# three with one virtual electrode, and that the result matches MNE-Python.
# Two pairs of the chain are enough to give the same group.
test_that("interpolate_bridged_electrodes rebuilds a chain of three as one group (S3c)", {
  eeg <- make_bridge_eeg(sc$S3c)
  out <- interpolate_bridged_electrodes(eeg, bridge_res("S3c"))

  changed <- rownames(eeg$data)[rowSums(abs(out$data - eeg$data)) > 0]
  expect_setequal(changed, c("Pz", "CPz", "Cz"))
  expect_repair_matches(eeg, out, ref_s3c)
  expect_abs_close(max(abs(out$data - eeg$data)), 25.1, 0.05)

  expect_equal(
    last_step(out),
    paste0("interpolate_bridged_electrodes(): rebuilt 3 bridged channel(s) in ",
           "1 group(s) (Pz+CPz+Cz) from 61 good channel(s) plus ",
           "1 virtual electrode(s)."))

  # Cz-CPz and CPz-Pz already connect all three
  two_pairs <- interpolate_bridged_electrodes(
    eeg, rbind(c("Cz", "CPz"), c("CPz", "Pz")))
  expect_identical(two_pairs, out)
})

# ----------------------------------------------------------------------------
# Test 4.4: 'bad_limit' caps the size of a group
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a group larger than bad_limit is refused (the
# message names the channels in channel order and the limit), the default
# limit is 4, and raising the limit lets the repair through.
test_that("interpolate_bridged_electrodes refuses groups larger than bad_limit", {
  eeg <- make_bridge_eeg(sc$S3c)
  res <- bridge_res("S3c")

  expect_error(
    interpolate_bridged_electrodes(eeg, res, bad_limit = 2),
    "The channels Pz, CPz, Cz are bridged together and form a group of 3 electrodes (limit: 2)",
    fixed = TRUE)
  expect_no_error(interpolate_bridged_electrodes(eeg, res, bad_limit = 3))

  # the default limit is 4: a chain of five is too much
  chain5 <- rbind(c("Cz", "CPz"), c("CPz", "Pz"), c("Pz", "POz"), c("POz", "Oz"))
  expect_error(
    interpolate_bridged_electrodes(eeg, chain5),
    "The channels Oz, POz, Pz, CPz, Cz are bridged together and form a group of 5 electrodes (limit: 4)",
    fixed = TRUE)
  out <- interpolate_bridged_electrodes(eeg, chain5, bad_limit = 5)
  expect_match(last_step(out),
               "rebuilt 5 bridged channel(s) in 1 group(s) (Oz+POz+Pz+CPz+Cz)",
               fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 4.5: The repair is linear and keeps a constant
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The rebuilt channels are weighted averages of the donors and
# the virtual electrodes, with weights that add up to 1. So a recording where
# every channel holds the same constant stays that constant, and scaling the
# data scales the result.
test_that("interpolate_bridged_electrodes keeps a constant and is linear", {
  eeg <- make_bridge_eeg(sc$S3)

  flat <- eeg
  flat$data[] <- 7.5
  out <- interpolate_bridged_electrodes(flat, s3_pairs)
  expect_abs_close(out$data, flat$data, 1e-9)

  base <- repaired_s3()
  scaled <- eeg
  scaled$data <- 3 * eeg$data
  out <- interpolate_bridged_electrodes(scaled, s3_pairs)
  expect_abs_close(out$data, 3 * base$data, 1e-9)
})

# ----------------------------------------------------------------------------
# Test 4.6: Which channels are donors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the channels the bridged electrodes are rebuilt
# from are EEG channels with a montage position that are not bad and not in a
# group. A channel that is bad, not typed "eeg" or missing from the montage
# does not influence the result at all (here: its data are replaced by values
# that would ruin the repair if they were used), and is not counted in the
# history. Leaving a good channel out does change the result.
test_that("interpolate_bridged_electrodes only uses good EEG channels with a position as donors", {
  eeg <- make_bridge_eeg(sc$S3)
  base <- repaired_s3()
  wild <- function(x) x * 50 + 1000

  # a bad channel
  bad <- eeg
  bad$bads <- "Fz"
  out <- interpolate_bridged_electrodes(bad, s3_pairs)
  bad$data["Fz", ] <- wild(bad$data["Fz", ])
  out_wild <- interpolate_bridged_electrodes(bad, s3_pairs)
  expect_equal(out_wild$data[s3_rebuilt, ], out$data[s3_rebuilt, ], tolerance = 1e-12)
  expect_match(last_step(out), "from 59 good channel(s)", fixed = TRUE)
  expect_gt(max(abs(out$data[s3_rebuilt, ] - base$data[s3_rebuilt, ])), 0.01)
  # the bad channel itself is not touched
  expect_identical(out_wild$data["Fz", ], bad$data["Fz", ])

  # a channel that is not an EEG channel
  not_eeg <- set_channel_types(eeg, c(Fp1 = "eog"))
  out <- interpolate_bridged_electrodes(not_eeg, s3_pairs)
  not_eeg$data["Fp1", ] <- wild(not_eeg$data["Fp1", ])
  out_wild <- interpolate_bridged_electrodes(not_eeg, s3_pairs)
  expect_equal(out_wild$data[s3_rebuilt, ], out$data[s3_rebuilt, ], tolerance = 1e-12)
  expect_match(last_step(out), "from 59 good channel(s)", fixed = TRUE)

  # a channel with no position in the montage
  no_position <- eeg
  no_position$montage <- create_montage(setdiff(create_montage()$channels, "Fp1"))
  out <- interpolate_bridged_electrodes(no_position, s3_pairs)
  no_position$data["Fp1", ] <- wild(no_position$data["Fp1", ])
  out_wild <- interpolate_bridged_electrodes(no_position, s3_pairs)
  expect_equal(out_wild$data[s3_rebuilt, ], out$data[s3_rebuilt, ], tolerance = 1e-12)
  expect_match(last_step(out), "from 59 good channel(s)", fixed = TRUE)
})

# ----------------------------------------------------------------------------
# Test 4.7: Bad channels stay as they are
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the repaired channels are not marked bad and the
# existing list of bad channels is unchanged, in the usual order of work:
# detect on an object that has a bad channel, then repair that same object.
test_that("interpolate_bridged_electrodes leaves the list of bad channels alone", {
  eeg <- make_bridge_eeg(sc$S3)
  eeg$bads <- "Fz"

  found <- find_bridged_electrodes(eeg, verbose = FALSE)
  expect_equal(nrow(found$bridged), 2)
  out <- interpolate_bridged_electrodes(eeg, found)

  expect_equal(out$bads, "Fz")
  expect_setequal(rownames(eeg$data)[rowSums(abs(out$data - eeg$data)) > 0],
                  s3_rebuilt)
})

# ----------------------------------------------------------------------------
# Test 4.8: The sphere centre and the unit of the positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies 'origin' is the centre the positions are measured
# from: a montage moved by a fixed shift gives the same repair when origin is
# that shift, and a sphericity warning when it is not. Positions in metres
# instead of millimetres give the same repair.
test_that("interpolate_bridged_electrodes measures positions from 'origin' and ignores their unit", {
  eeg <- make_bridge_eeg(sc$S3)
  base <- repaired_s3()

  shift <- c(30, -20, 10)
  moved <- eeg
  moved_montage <- create_montage()
  moved_montage$positions$x <- moved_montage$positions$x + shift[1]
  moved_montage$positions$y <- moved_montage$positions$y + shift[2]
  moved_montage$positions$z <- moved_montage$positions$z + shift[3]
  moved$montage <- moved_montage

  expect_warning(interpolate_bridged_electrodes(moved, s3_pairs),
                 "not close to spherical")
  expect_no_warning(out <- interpolate_bridged_electrodes(moved, s3_pairs,
                                                          origin = shift))
  expect_abs_close(out$data, base$data, 1e-9)

  metres <- eeg
  metres_montage <- create_montage()
  metres_montage$positions[, c("x", "y", "z")] <-
    metres_montage$positions[, c("x", "y", "z")] / 1000
  metres$montage <- metres_montage
  out <- interpolate_bridged_electrodes(metres, s3_pairs)
  expect_abs_close(out$data, base$data, 1e-8)
})

# ----------------------------------------------------------------------------
# Test 4.9: The history and data without row names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the new history entry is added after the existing
# ones, and that the repair works when the data matrix has no row names (the
# channels are found through eeg_obj$channels).
test_that("interpolate_bridged_electrodes appends to the history and works without row names", {
  eeg <- make_bridge_eeg(sc$S3)
  base <- repaired_s3()

  eeg$preprocessing_history <- list("step one", "step two")
  out <- interpolate_bridged_electrodes(eeg, s3_pairs)
  expect_length(out$preprocessing_history, 3)
  expect_equal(out$preprocessing_history[1:2], list("step one", "step two"))
  expect_match(out$preprocessing_history[[3]], "^interpolate_bridged_electrodes\\(\\)")

  plain <- make_bridge_eeg(sc$S3)
  dimnames(plain$data) <- NULL
  out <- interpolate_bridged_electrodes(plain, s3_pairs)
  expect_null(rownames(out$data))
  expect_equal(unname(out$data), unname(base$data))
})

# ----------------------------------------------------------------------------
# Test 4.10: No pairs, nothing to do
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a result with no pairs (also as an empty table or
# an empty data frame) gives a message and the object back unchanged.
test_that("interpolate_bridged_electrodes returns the object unchanged when there are no pairs", {
  eeg <- make_bridge_eeg(sc$S3)

  no_pairs <- list(
    bridge_res("S2"),                                  # a detection result
    matrix(character(0), 0, 2),                        # an empty table
    data.frame(channel_1 = character(0), channel_2 = character(0),
               stringsAsFactors = FALSE))              # an empty data frame

  for (empty in no_pairs) {
    expect_message(out <- interpolate_bridged_electrodes(eeg, empty),
                   "no bridged pairs given")
    expect_identical(out, eeg)
  }
})

# ----------------------------------------------------------------------------
# Test 4.11: After the repair the detector finds nothing
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the repaired channels are no longer copies of each
# other: the detector finds no pair, and the typical (median) distance of each
# repaired pair now lies above the cutoff the original recording gave. (The
# median is used because a single window can still dip lower.)
test_that("a repaired recording shows no bridges", {
  before <- bridge_res("S3")
  after <- find_bridged_electrodes(repaired_s3(), verbose = FALSE)

  expect_equal(nrow(after$bridged), 0)

  median_after <- function(a, b) {
    row <- which(after$ed_pairs$channel_1 == a & after$ed_pairs$channel_2 == b)
    median(after$ed[row, ])
  }
  expect_gt(median_after("CPz", "Cz"), before$local_minimum)
  expect_gt(median_after("P3", "P5"), before$local_minimum)

  # the distances measured after the repair (0.023 and 0.022 before it)
  expect_abs_close(median_after("CPz", "Cz"), 2.724, 0.005)
  expect_abs_close(median_after("P3", "P5"), 2.189, 0.005)
})

# ============================================================================
# TEST SUITE 5: interpolate_bridged_electrodes() input validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Requires a continuous 'eeg' object with a data matrix
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies interpolate_bridged_electrodes() refuses an object
# that is not 'eeg' (epoched data included) or has no data matrix.
test_that("interpolate_bridged_electrodes errors when eeg_obj is not a usable 'eeg' object", {
  expect_error(interpolate_bridged_electrodes(list(), s3_pairs), "class 'eeg'")
  expect_error(interpolate_bridged_electrodes(NULL, s3_pairs), "class 'eeg'")
  expect_error(interpolate_bridged_electrodes(make_bridge_epochs(), s3_pairs),
               "class 'eeg'")

  eeg <- make_bridge_eeg(sc$S3)
  eeg$data <- NULL
  expect_error(interpolate_bridged_electrodes(eeg, s3_pairs), "numeric matrix")
})

# ----------------------------------------------------------------------------
# Test 5.2: 'bad_limit' and 'origin'
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies bad_limit must be a single whole number of at least
# 1, and origin a finite numeric vector of length 3.
test_that("interpolate_bridged_electrodes validates bad_limit and origin", {
  eeg <- make_bridge_eeg(sc$S3)

  for (bad in list(0, 1.5, -1, NA_real_, Inf, c(2, 3), "a", TRUE, NULL)) {
    expect_error(interpolate_bridged_electrodes(eeg, s3_pairs, bad_limit = bad),
                 "'bad_limit' must be a single whole number",
                 info = paste("bad_limit =", paste(bad, collapse = ", ")))
  }

  for (bad in list(c(0, 0), c(0, 0, 0, 0), c(0, 0, NA), c(0, 0, Inf), "a", NULL)) {
    expect_error(interpolate_bridged_electrodes(eeg, s3_pairs, origin = bad),
                 "'origin' must be a numeric vector of length 3",
                 info = paste("origin =", paste(bad, collapse = ", ")))
  }
})

# ----------------------------------------------------------------------------
# Test 5.3: 'bridged' must be a detection result or a two-column table
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies vectors, lists, tables with the wrong number of
# columns and tables that are not made of channel names are all refused, and a
# row naming the same channel twice is refused.
test_that("interpolate_bridged_electrodes validates 'bridged'", {
  eeg <- make_bridge_eeg(sc$S3)

  not_tables <- list(
    c("Cz", "CPz"),                           # a plain vector
    "Cz",
    NULL,
    list("Cz", "CPz"),
    matrix(1:4, 2),                           # numbers, not names
    matrix("Cz", 1, 3),                       # three columns
    data.frame(a = "Cz", stringsAsFactors = FALSE),
    data.frame(a = 1, b = 2))
  for (i in seq_along(not_tables)) {
    expect_error(interpolate_bridged_electrodes(eeg, not_tables[[i]]),
                 "two-column table of channel names",
                 info = paste("input", i))
  }

  expect_error(interpolate_bridged_electrodes(eeg, rbind(c("Cz", "Cz"))),
               "two different channels")
  expect_error(interpolate_bridged_electrodes(eeg, rbind(c("P3", "P5"), c("Cz", "Cz"))),
               "two different channels")
})

# ----------------------------------------------------------------------------
# Test 5.4: Every named channel must be an EEG channel of the object
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an unknown channel and a channel that is not typed
# "eeg" are refused, and the message names them.
test_that("interpolate_bridged_electrodes errors for channels that are not EEG channels", {
  eeg <- make_bridge_eeg(sc$S3)

  expect_error(interpolate_bridged_electrodes(eeg, rbind(c("Cz", "Nope"))),
               "are not EEG channels of eeg_obj: Nope")

  not_eeg <- set_channel_types(eeg, c(Fp1 = "eog"))
  expect_error(interpolate_bridged_electrodes(not_eeg, rbind(c("Fp1", "Fpz"))),
               "are not EEG channels of eeg_obj: Fp1")
})

# ----------------------------------------------------------------------------
# Test 5.5: A montage is needed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the repair errors without a montage, and with
# something that is not a montage.
test_that("interpolate_bridged_electrodes errors without a montage", {
  eeg <- make_bridge_eeg(sc$S3)

  no_montage <- eeg
  no_montage$montage <- NULL
  expect_error(interpolate_bridged_electrodes(no_montage, s3_pairs),
               "No montage attached")

  plain_list <- eeg
  plain_list$montage <- unclass(create_montage())
  expect_error(interpolate_bridged_electrodes(plain_list, s3_pairs),
               "No montage attached")
})

# ----------------------------------------------------------------------------
# Test 5.6: Every bridged channel needs a position
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a bridged channel missing from the montage is
# refused, and named.
test_that("interpolate_bridged_electrodes errors for a bridged channel with no position", {
  eeg <- make_bridge_eeg(sc$S3)
  eeg$montage <- create_montage(setdiff(create_montage()$channels, "P5"))

  expect_error(interpolate_bridged_electrodes(eeg, s3_pairs),
               "no position in the montage: P5")
})

# ----------------------------------------------------------------------------
# Test 5.7: Bridged channels must not already be bad
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a bridged channel that is already marked bad is
# refused, and the message points to interpolate_bads().
test_that("interpolate_bridged_electrodes errors for a bridged channel that is already bad", {
  eeg <- make_bridge_eeg(sc$S3)
  eeg$bads <- "Cz"

  expect_error(interpolate_bridged_electrodes(eeg, s3_pairs),
               "already marked bad: Cz")
  expect_error(interpolate_bridged_electrodes(eeg, s3_pairs), "interpolate_bads")
})

# ----------------------------------------------------------------------------
# Test 5.8: There must be donors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the repair errors when no good EEG channel with a
# position is left to rebuild from (all others typed non-EEG, or all bad).
test_that("interpolate_bridged_electrodes errors when there are no donor channels", {
  eeg <- make_bridge_eeg(sc$S3)
  pair <- rbind(c("CPz", "Cz"))

  only_pair <- eeg
  only_pair$channel_types[!(eeg$channels %in% c("CPz", "Cz"))] <- "misc"
  expect_error(interpolate_bridged_electrodes(only_pair, pair),
               "No good EEG channels with a montage position")

  all_bad <- eeg
  all_bad$bads <- setdiff(eeg$channels, c("CPz", "Cz"))
  expect_error(interpolate_bridged_electrodes(all_bad, pair),
               "No good EEG channels with a montage position")
})

# ----------------------------------------------------------------------------
# Test 5.9: Positions that are not close to a sphere
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a warning when an electrode lies far from the
# sphere the others are on (here: one position 1.6 times too far out), whether
# it is a channel the repair uses (Fz) or one it rebuilds (Cz), that the repair
# still runs, and that the normal montage gives no warning.
test_that("interpolate_bridged_electrodes warns when positions are not close to a sphere", {
  eeg <- make_bridge_eeg(sc$S3)

  stretch <- function(channel) {
    moved <- eeg
    montage <- create_montage()
    row <- montage$positions$channel == channel
    montage$positions[row, c("x", "y", "z")] <-
      1.6 * montage$positions[row, c("x", "y", "z")]
    moved$montage <- montage
    moved
  }

  expect_warning(out <- interpolate_bridged_electrodes(stretch("Fz"), s3_pairs),
                 "not close to spherical around 'origin'")
  expect_equal(dim(out$data), dim(eeg$data))

  expect_warning(interpolate_bridged_electrodes(stretch("Cz"), s3_pairs),
                 "not close to spherical around 'origin'")

  expect_no_warning(interpolate_bridged_electrodes(eeg, s3_pairs))
})

# ----------------------------------------------------------------------------
# Test 5.10: The virtual electrode needs a place
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the error when the electrodes of a group lie on
# opposite sides of the head, so no virtual electrode can be put between them
# (here: P5 made the mirror image of P3).
test_that("interpolate_bridged_electrodes errors when a group has no centre", {
  eeg <- make_bridge_eeg(sc$S3)
  montage <- create_montage()
  p3 <- montage$positions$channel == "P3"
  p5 <- montage$positions$channel == "P5"
  montage$positions[p5, c("x", "y", "z")] <- -montage$positions[p3, c("x", "y", "z")]
  eeg$montage <- montage

  expect_error(interpolate_bridged_electrodes(eeg, rbind(c("P3", "P5"))),
               "Cannot place a virtual electrode for P3, P5")
})

# ----------------------------------------------------------------------------
# Test 5.11: The checks run in the documented order
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when two things are wrong at once, the first
# check in the documented order is the one reported.
test_that("interpolate_bridged_electrodes reports the first failing check", {
  eeg <- make_bridge_eeg(sc$S3)

  # class, data, bad_limit, origin, then 'bridged'
  expect_error(interpolate_bridged_electrodes(list(), c("Cz", "CPz"), bad_limit = 0),
               "class 'eeg'")
  expect_error(interpolate_bridged_electrodes(eeg, c("Cz", "CPz"), bad_limit = 0),
               "'bad_limit' must be a single whole number")
  expect_error(interpolate_bridged_electrodes(eeg, c("Cz", "CPz"), origin = 1),
               "'origin' must be a numeric vector of length 3")

  # the same channel twice, before the channels are looked up
  expect_error(
    interpolate_bridged_electrodes(eeg, rbind(c("Cz", "Cz"), c("Cz", "Nope"))),
    "two different channels")

  # unknown channels, then the montage
  no_montage <- eeg
  no_montage$montage <- NULL
  expect_error(interpolate_bridged_electrodes(no_montage, rbind(c("Cz", "Nope"))),
               "are not EEG channels")

  # the montage, then positions, then bad channels
  no_montage$bads <- "Cz"
  expect_error(interpolate_bridged_electrodes(no_montage, rbind(c("Cz", "CPz"))),
               "No montage attached")

  no_p5 <- eeg
  no_p5$montage <- create_montage(setdiff(create_montage()$channels, "P5"))
  no_p5$bads <- "P3"
  expect_error(interpolate_bridged_electrodes(no_p5, s3_pairs),
               "no position in the montage: P5")

  # bad channels, then the size of the groups, then the donors
  bad_cz <- make_bridge_eeg(sc$S3c)
  bad_cz$bads <- "Cz"
  expect_error(interpolate_bridged_electrodes(bad_cz, bridge_res("S3c"), bad_limit = 2),
               "already marked bad")

  no_donors <- make_bridge_eeg(sc$S3c)
  no_donors$channel_types[!(no_donors$channels %in% c("Cz", "CPz", "Pz"))] <- "misc"
  expect_error(interpolate_bridged_electrodes(no_donors, bridge_res("S3c"), bad_limit = 2),
               "form a group of 3 electrodes")
})
