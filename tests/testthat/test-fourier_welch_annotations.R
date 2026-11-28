# ============================================================================
#      Test File for eeg_psd_welch() leaving out "bad" annotated stretches
# ============================================================================
#
# For continuous 'eeg' input, eeg_psd_welch() now leaves out every sample
# covered by an annotation whose description starts with "bad" (any case),
# unless reject_by_annotation = FALSE. The recording is split into the clean
# stretches between the bad ones; each clean stretch gets its own Welch
# estimate, and the stretches are averaged, weighted by the samples each one
# used. A clean stretch shorter than one segment is analysed as one shorter
# segment, zero-padded to the same FFT length.
#
# The main test rebuilds that by hand: cut the clean stretches out as
# separate recordings, run eeg_psd_welch() on each, and take the weighted
# average.
#
# Layout used by most tests: 12 s at 512 Hz = 6144 samples, 2 channels.
#   bad  4.0 - 4.5 s  = samples 2049 .. 2304
#   bad  5.0 - 5.5 s  = samples 2561 .. 2816
#   clean stretches: 1 .. 2048     (2048 samples, 1 segment, weight 2048)
#                    2305 .. 2560  ( 256 samples, short,     weight  256)
#                    2817 .. 6144  (3328 samples, 1 segment, weight 2048;
#                                   the last 1280 samples are not used)
#
# Test suites:
#   1. The estimate uses only the clean stretches
#   2. What does not count as "bad"
#   3. Edge cases (all bad, a one-sample stretch, epochs)
#   4. The stretch helper
#
# Author: Christos Dalamarinis
# Date: Nov 2026
# ============================================================================

library(testthat)
library(eeganalysis)

SR_A <- 512   # sampling rate (Hz)

ann_rows <- function(onset, duration, description = "BAD_test") {
  data.frame(onset = onset, duration = duration, description = description,
             channel = NA_character_, stringsAsFactors = FALSE)
}

ann_eeg <- function(X, ann = NULL) {
  new_eeg(X, paste0("Ch", seq_len(nrow(X))), SR_A, annotations = ann)
}

two_bad <- function(description = "BAD_test") {
  ann_rows(c(4.0, 5.0), c(0.5, 0.5), description)
}

set.seed(10L)
X_A <- rbind(20 * sin(2 * pi * 10 * (0:6143) / SR_A) + rnorm(6144, sd = 5),
             rnorm(6144, sd = 8))

psd_of <- function(X, ann = NULL, ...) {
  eeg_psd_welch(ann_eeg(X, ann), verbose = FALSE, ...)
}


# ============================================================================
# Suite 1: the estimate uses only the clean stretches
# ============================================================================

test_that("the result is the weighted average of the clean stretches", {
  full <- psd_of(X_A, two_bad())
  # each clean stretch on its own, same 2048-point FFT
  stretches <- list(1:2048, 2305:2560, 2817:6144)
  weights   <- c(2048, 256, 2048)
  parts <- Map(function(idx, w) {
    w * psd_of(X_A[, idx, drop = FALSE], n_fft = 2048)$power
  }, stretches, weights)
  expected <- Reduce(`+`, parts) / sum(weights)
  expect_equal(full$power, expected, tolerance = 1e-12)
  expect_equal(full$n_segments, 3)
})

test_that("what is inside a bad stretch does not matter", {
  X_big <- X_A
  X_big[, 2049:2304] <- 5000                     # huge artefact
  X_nan <- X_A
  X_nan[1L, 2600:2700] <- NaN                    # missing samples
  ref <- psd_of(X_A, two_bad())$power
  expect_identical(psd_of(X_big, two_bad())$power, ref)
  expect_no_warning(nan_res <- psd_of(X_nan, two_bad())$power)
  expect_identical(nan_res, ref)
})

test_that("reject_by_annotation = FALSE uses every sample", {
  X_nan <- X_A
  X_nan[1L, 2600:2700] <- NaN
  expect_identical(psd_of(X_A, two_bad(), reject_by_annotation = FALSE)$power,
                   psd_of(X_A)$power)
  # the missing samples are used too, so channel 1 becomes NA
  expect_warning(res <- psd_of(X_nan, two_bad(), reject_by_annotation = FALSE),
                 "NaN/Inf")
  expect_true(all(is.na(res$power[1L, ])))
})


# ============================================================================
# Suite 2: what does not count as "bad"
# ============================================================================

test_that("only descriptions starting with 'bad' count, in any case", {
  no_ann <- psd_of(X_A)$power
  expect_identical(psd_of(X_A, two_bad("eyes_closed"))$power, no_ann)
  expect_identical(psd_of(X_A, two_bad("my_BAD_note"))$power, no_ann)
  lower <- psd_of(X_A, two_bad("bad_blink"))$power
  expect_false(isTRUE(all.equal(lower, no_ann)))   # it does leave something out
  expect_identical(lower, psd_of(X_A, two_bad("BAD_test"))$power)
})


# ============================================================================
# Suite 3: edge cases
# ============================================================================

test_that("a recording that is bad everywhere is an error", {
  expect_error(psd_of(X_A, ann_rows(0, 12)), "Every sample")
})

test_that("a one-sample clean stretch gives zeros, not NaN", {
  # samples 1 .. 2047 and 2049 .. 6144 are bad; sample 2048 alone is clean
  ann <- ann_rows(c(0, 2048 / SR_A), c(2047 / SR_A, 4096 / SR_A))
  res <- psd_of(X_A, ann)$power
  expect_true(all(is.finite(res)))
  expect_true(all(res == 0))                     # one sample minus its mean
  for (type in c("hann", "hamming", "blackman", "none")) {
    expect_equal(eeganalysis:::.make_window(1L, type), 1)
  }
})

test_that("epochs are not affected", {
  set.seed(11L)
  ep <- structure(list(data = array(rnorm(2 * 513 * 3), c(2L, 513L, 3L)),
                       channels = c("Ch1", "Ch2"), sampling_rate = SR_A),
                  class = "eeg_epochs")
  expect_identical(
    eeg_psd_welch(ep, reject_by_annotation = TRUE, verbose = FALSE)$power,
    eeg_psd_welch(ep, reject_by_annotation = FALSE, verbose = FALSE)$power)
})

test_that("reject_by_annotation must be TRUE or FALSE", {
  expect_error(psd_of(X_A, reject_by_annotation = NA), "reject_by_annotation")
  expect_error(psd_of(X_A, reject_by_annotation = "yes"), "reject_by_annotation")
})


# ============================================================================
# Suite 4: the stretch helper
# ============================================================================

test_that(".welch_spans() finds the stretches, segments and weights", {
  bad <- rep(FALSE, 5120)
  bad[1537:2048] <- TRUE
  sp <- eeganalysis:::.welch_spans(bad, win_samp = 1024, step = 1024)
  expect_equal(sp$start,    c(1, 2049))
  expect_equal(sp$end,      c(1536, 5120))
  expect_equal(sp$win_samp, c(1024, 1024))
  expect_equal(sp$n_seg,    c(1, 3))
  expect_equal(sp$weight,   c(1024, 3072))       # the last 512 of stretch 1 unused
  # a stretch shorter than one segment becomes one segment of its own length
  bad2 <- rep(FALSE, 3000)
  bad2[1001:2900] <- TRUE
  sp2 <- eeganalysis:::.welch_spans(bad2, win_samp = 2048, step = 2048)
  expect_equal(sp2$win_samp, c(1000, 100))
  expect_equal(sp2$weight,   c(1000, 100))
})
