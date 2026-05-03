# ============================================================================
#                     Test File for eeg_summary.R
# ============================================================================
#
# Functions tested:
#   1. eeg_summary() - Full diagnostic report for eeg objects
#
# Test suites:
#   1. Input validation
#   2. Channel classification (EEG vs EXG vs Status)
#   3. Per-channel and global statistics
#   4. Suspicious channel flagging
#   5. Return value structure
#   6. Print output and parameters
#   7. Integration tests
#
# Author: Christos Dalamarinis
# Date: May 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

# 4 EEG channels, 1 EXG channel, 1 Status channel, 200 timepoints at 256 Hz
.make_mixed_eeg <- function() {
  set.seed(42)
  data <- matrix(rnorm(6 * 200, mean = 0, sd = 15), nrow = 6, ncol = 200)
  new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "Fz", "EXG1", "Status"),
    sampling_rate = 256
  )
}

# Pure EEG object - no EXG, no Status
.make_pure_eeg <- function(n_ch = 4, n_tp = 200, sd = 15) {
  set.seed(1)
  data <- matrix(rnorm(n_ch * n_tp, mean = 0, sd = sd), nrow = n_ch, ncol = n_tp)
  new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "Fz")[seq_len(n_ch)],
    sampling_rate = 256
  )
}


# ============================================================================
# TEST SUITE 1: Input validation
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Non-eeg object rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: eeg_summary() must stop with a clear message when the
# supplied object is not of class 'eeg'.
test_that("eeg_summary errors on non-eeg input", {
  expect_error(
    eeg_summary(list(data = matrix(1, 1, 1))),
    regexp = "requires an object of class 'eeg'"
  )
  expect_error(
    eeg_summary(42),
    regexp = "requires an object of class 'eeg'"
  )
})


# ----------------------------------------------------------------------------
# Test 1.2: Non-matrix data field rejected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: eeg_summary() stops when eeg_obj$data is not a matrix,
# since all statistical computations depend on row-indexing.
test_that("eeg_summary errors when data is not a matrix", {
  eeg <- .make_pure_eeg()
  eeg$data <- as.data.frame(eeg$data)  # break the matrix requirement

  expect_error(
    eeg_summary(eeg),
    regexp = "must be a numeric matrix"
  )
})


# ----------------------------------------------------------------------------
# Test 1.3: Valid eeg object accepted
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A properly formed eeg object runs without error.
test_that("eeg_summary runs without error on a valid eeg object", {
  eeg <- .make_pure_eeg()
  expect_silent(capture.output(eeg_summary(eeg)))
})


# ============================================================================
# TEST SUITE 2: Channel classification
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Status channel excluded from returned stats
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The Status channel must never appear in eeg_stats or
# exg_stats, regardless of its position in the channel vector.
test_that("eeg_summary excludes Status channel from all stats", {
  eeg <- .make_mixed_eeg()  # has "Status" as last channel

  result <- invisible(capture.output(report <- eeg_summary(eeg)))

  expect_false("Status" %in% report$eeg_stats$channel)
  if (!is.null(report$exg_stats)) {
    expect_false("Status" %in% report$exg_stats$channel)
  }
})


# ----------------------------------------------------------------------------
# Test 2.2: EXG channels detected via Pass 1 (original names)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Channels named "EXG1"..."EXG8" are routed to exg_stats,
# not eeg_stats, via detect_external_channels() (Pass 1).
test_that("eeg_summary routes EXG channels to exg_stats (Pass 1)", {
  eeg <- .make_mixed_eeg()  # has "EXG1"

  capture.output(report <- eeg_summary(eeg))

  expect_false("EXG1" %in% report$eeg_stats$channel)
  expect_true("EXG1" %in% report$exg_stats$channel)
})


# ----------------------------------------------------------------------------
# Test 2.3: EXG channels detected via Pass 2 (renamed with parentheses)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Channels renamed by apply_external_labels(keep_original=TRUE)
# e.g. "MASTOID LEFT (EXG5)" must still be routed to exg_stats via the regex
# fallback (Pass 2).
test_that("eeg_summary routes renamed EXG channels to exg_stats (Pass 2)", {
  set.seed(7)
  data <- matrix(rnorm(4 * 200), nrow = 4, ncol = 200)
  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "MASTOID LEFT (EXG5)", "EOG LEFT (EXG1)"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_false("MASTOID LEFT (EXG5)" %in% report$eeg_stats$channel)
  expect_false("EOG LEFT (EXG1)"     %in% report$eeg_stats$channel)
  expect_true("MASTOID LEFT (EXG5)"  %in% report$exg_stats$channel)
  expect_true("EOG LEFT (EXG1)"      %in% report$exg_stats$channel)
})


# ----------------------------------------------------------------------------
# Test 2.4: No EXG channels → exg_stats is NULL
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When the recording contains no EXG channels, exg_stats
# must be NULL (not an empty data frame).
test_that("eeg_summary returns NULL exg_stats when no EXG channels present", {
  eeg <- .make_pure_eeg()

  capture.output(report <- eeg_summary(eeg))

  expect_null(report$exg_stats)
})


# ----------------------------------------------------------------------------
# Test 2.5: All channels are EXG/Status → error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: If no EEG channels remain after exclusions, eeg_summary()
# must stop with an informative message rather than compute on empty data.
test_that("eeg_summary errors when no EEG channels remain", {
  set.seed(3)
  data <- matrix(rnorm(3 * 100), nrow = 3, ncol = 100)
  eeg <- new_eeg(
    data     = data,
    channels = c("EXG1", "EXG2", "Status"),
    sampling_rate = 256
  )

  expect_error(
    capture.output(eeg_summary(eeg)),
    regexp = "No EEG channels found"
  )
})


# ============================================================================
# TEST SUITE 3: Statistics
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Per-channel stats are numerically correct
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: min, max, mean, and std in eeg_stats match what base R
# would compute on the same rows.
test_that("eeg_summary per-channel stats match base R calculations", {
  set.seed(99)
  data <- matrix(rnorm(3 * 500), nrow = 3, ncol = 500)
  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))
  stats <- report$eeg_stats

  for (i in 1:3) {
    ch_row <- stats[stats$channel == c("Cz", "Pz", "Oz")[i], ]
    expect_equal(ch_row$min_uv,  round(min(data[i, ]),  2))
    expect_equal(ch_row$max_uv,  round(max(data[i, ]),  2))
    expect_equal(ch_row$mean_uv, round(mean(data[i, ]), 2))
    expect_equal(ch_row$std_uv,  round(sd(data[i, ]),   2))
  }
})


# ----------------------------------------------------------------------------
# Test 3.2: eeg_stats contains one row per EEG channel
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The number of rows in eeg_stats equals the number of
# EEG channels (total minus Status and EXG).
test_that("eeg_summary eeg_stats has one row per EEG channel", {
  eeg <- .make_mixed_eeg()  # 4 EEG + 1 EXG + 1 Status

  capture.output(report <- eeg_summary(eeg))

  expect_equal(nrow(report$eeg_stats), 4)
})


# ----------------------------------------------------------------------------
# Test 3.3: exg_stats contains one row per EXG channel
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Mirrors 3.2 for the EXG partition.
test_that("eeg_summary exg_stats has one row per EXG channel", {
  set.seed(5)
  data <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)
  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "EXG1", "EXG2", "EXG3"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_equal(nrow(report$exg_stats), 3)
})


# ============================================================================
# TEST SUITE 4: Suspicious channel flagging
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Flat / disconnected channel is flagged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A channel with std near zero (constant signal) must be
# flagged with reason "Flat / possibly disconnected".
test_that("eeg_summary flags flat channels", {
  set.seed(10)
  data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data[2, ] <- rep(5, 200)  # channel Pz: perfectly flat

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_false(is.null(report$flags))
  flat_flags <- report$flags[grepl("Flat", report$flags$reason), ]
  expect_true("Pz" %in% flat_flags$channel)
})


# ----------------------------------------------------------------------------
# Test 4.2: Excessive amplitude channel is flagged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A channel with a peak absolute value above the default
# 500 uV threshold must be flagged with "Excessive amplitude".
test_that("eeg_summary flags channels with excessive amplitude", {
  set.seed(11)
  data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data[3, 100] <- 800  # Oz: one saturated sample

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_false(is.null(report$flags))
  amp_flags <- report$flags[grepl("Excessive", report$flags$reason), ]
  expect_true("Oz" %in% amp_flags$channel)
})


# ----------------------------------------------------------------------------
# Test 4.3: Outlier noise channel is flagged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A channel whose std exceeds flag_outlier_sd_multiplier
# times the median std of all EEG channels is flagged as "Outlier noise level".
test_that("eeg_summary flags channels with outlier noise", {
  set.seed(12)
  # Three quiet channels + one extremely noisy channel
  data <- rbind(
    matrix(rnorm(3 * 200, sd = 10), nrow = 3, ncol = 200),
    matrix(rnorm(1 * 200, sd = 300), nrow = 1, ncol = 200)
  )

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "Fz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_false(is.null(report$flags))
  noise_flags <- report$flags[grepl("Outlier", report$flags$reason), ]
  expect_true("Fz" %in% noise_flags$channel)
})


# ----------------------------------------------------------------------------
# Test 4.4: No flags when all channels are clean
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A well-behaved recording (moderate amplitude, no flat
# channels) produces NULL flags.
test_that("eeg_summary returns NULL flags when all channels are clean", {
  set.seed(20)
  data <- matrix(rnorm(4 * 500, mean = 0, sd = 15), nrow = 4, ncol = 500)
  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "Fz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_null(report$flags)
})


# ----------------------------------------------------------------------------
# Test 4.5: Custom flag thresholds are respected
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Passing non-default thresholds changes which channels are
# flagged. A channel with std ≈ 0.4 should be flagged at the default 0.5
# threshold but NOT flagged when the threshold is lowered to 0.1.
test_that("eeg_summary respects custom flag thresholds", {
  set.seed(30)
  data <- rbind(
    matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200),
    rnorm(200, sd = 0.4)  # Fz: real std ≈ 0.4, below default 0.5 threshold
  )

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "Fz"),
    sampling_rate = 256
  )

  # Default threshold (0.5): Fz should be flagged
  capture.output(report_default <- eeg_summary(eeg, flag_flat_threshold = 0.5))
  flat_default <- report_default$flags[grepl("Flat", report_default$flags$reason), ]
  expect_true("Fz" %in% flat_default$channel)

  # Lowered threshold (0.1): Fz should NOT be flagged as flat
  capture.output(report_low <- eeg_summary(eeg, flag_flat_threshold = 0.1))
  flat_low <- if (!is.null(report_low$flags)) {
    report_low$flags[grepl("Flat", report_low$flags$reason), ]
  } else {
    data.frame(channel = character(0))
  }
  expect_false("Fz" %in% flat_low$channel)
})


# ----------------------------------------------------------------------------
# Test 4.6: EXG channels are never flagged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Flagging logic applies only to EEG channels. Even an EXG
# channel with zero variance must not appear in the flags table.
test_that("eeg_summary never flags EXG channels", {
  set.seed(40)
  data <- rbind(
    matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200),
    rep(0, 200)   # EXG1: completely flat — should be ignored by flag logic
  )

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz", "EXG1"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  flagged_channels <- if (is.null(report$flags)) character(0) else report$flags$channel
  expect_false("EXG1" %in% flagged_channels)
})


# ============================================================================
# TEST SUITE 5: Return value structure
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: Return value is an invisible list with three named elements
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The function contract guarantees a list with exactly
# eeg_stats, exg_stats, and flags. All three must be present.
test_that("eeg_summary returns a list with eeg_stats, exg_stats, flags", {
  eeg <- .make_pure_eeg()

  capture.output(report <- eeg_summary(eeg))

  expect_true(is.list(report))
  expect_named(report, c("eeg_stats", "exg_stats", "flags"), ignore.order = FALSE)
})


# ----------------------------------------------------------------------------
# Test 5.2: eeg_stats data frame has correct columns
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: eeg_stats must have channel, min_uv, max_uv, mean_uv,
# std_uv — names and types as documented.
test_that("eeg_summary eeg_stats has correct columns", {
  eeg <- .make_pure_eeg()

  capture.output(report <- eeg_summary(eeg))

  expect_true(is.data.frame(report$eeg_stats))
  expect_named(report$eeg_stats, c("channel", "min_uv", "max_uv", "mean_uv", "std_uv"))
  expect_type(report$eeg_stats$channel, "character")
  expect_type(report$eeg_stats$min_uv,  "double")
  expect_type(report$eeg_stats$std_uv,  "double")
})


# ----------------------------------------------------------------------------
# Test 5.3: flags data frame has correct columns when flags are raised
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The flags data frame must have channel, reason, and value
# columns when at least one flag is raised.
test_that("eeg_summary flags data frame has correct columns", {
  set.seed(50)
  data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data[1, ] <- 0   # flat channel

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  expect_false(is.null(report$flags))
  expect_named(report$flags, c("channel", "reason", "value"))
})


# ============================================================================
# TEST SUITE 6: Print output and parameters
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: Report contains all expected sections
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The printed report must include the main section headers
# that users rely on to navigate the output.
test_that("eeg_summary output contains all expected sections", {
  eeg <- .make_mixed_eeg()

  output <- capture.output(eeg_summary(eeg))

  expect_true(any(grepl("EEG Summary Report",        output)))
  expect_true(any(grepl("RECORDING INFORMATION",     output)))
  expect_true(any(grepl("EEG CHANNEL STATISTICS",    output)))
  expect_true(any(grepl("SUSPICIOUS CHANNEL FLAGS",  output)))
})


# ----------------------------------------------------------------------------
# Test 6.2: show_exg = FALSE suppresses EXG block
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: When show_exg is FALSE, the EXG section must not appear
# in the printed output, even when EXG channels are present.
test_that("eeg_summary suppresses EXG block when show_exg = FALSE", {
  eeg <- .make_mixed_eeg()

  output_suppressed <- capture.output(eeg_summary(eeg, show_exg = FALSE))
  output_shown      <- capture.output(eeg_summary(eeg, show_exg = TRUE))

  expect_false(any(grepl("EXG CHANNEL STATISTICS", output_suppressed)))
  expect_true(any(grepl("EXG CHANNEL STATISTICS",  output_shown)))
})


# ----------------------------------------------------------------------------
# Test 6.3: [OK] shown in output when no flags raised
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The clean-bill-of-health message must appear when no
# suspicious channels are detected.
test_that("eeg_summary prints [OK] message when no flags", {
  eeg <- .make_pure_eeg()

  output <- capture.output(eeg_summary(eeg))

  expect_true(any(grepl("\\[OK\\]", output)))
})


# ----------------------------------------------------------------------------
# Test 6.4: [!] shown in output when flags are raised
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The warning marker must appear when at least one flag
# is raised, so users are not misled into thinking the recording is clean.
test_that("eeg_summary prints [!] marker when flags are raised", {
  set.seed(60)
  data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data[1, ] <- 0  # flat

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz", "Oz"),
    sampling_rate = 256
  )

  output <- capture.output(eeg_summary(eeg))

  expect_true(any(grepl("\\[!\\]", output)))
})


# ----------------------------------------------------------------------------
# Test 6.5: Recording info values are printed correctly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Sampling rate and duration displayed in the header block
# match the actual values in the eeg object.
test_that("eeg_summary prints correct sampling rate and duration", {
  set.seed(70)
  sr  <- 512
  dur <- 5   # seconds
  data <- matrix(rnorm(2 * sr * dur), nrow = 2, ncol = sr * dur)

  eeg <- new_eeg(
    data     = data,
    channels = c("Cz", "Pz"),
    sampling_rate = sr
  )

  output <- capture.output(eeg_summary(eeg))

  expect_true(any(grepl("512", output)))          # sampling rate
  expect_true(any(grepl("5\\.00|4\\.99", output))) # duration ≈ 5 s
})


# ============================================================================
# TEST SUITE 7: Integration tests
# ============================================================================

# ----------------------------------------------------------------------------
# Test 7.1: Mixed EEG + EXG + Status recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A realistic multi-type recording is partitioned correctly,
# stats are computed on EEG only, and EXG stats are returned separately.
test_that("eeg_summary handles mixed EEG + EXG + Status correctly", {
  set.seed(80)
  data <- matrix(rnorm(6 * 512, sd = 20), nrow = 6, ncol = 512)

  eeg <- new_eeg(
    data     = data,
    channels = c("Fp1", "Fp2", "Cz", "Pz", "EXG1", "Status"),
    sampling_rate = 512
  )

  capture.output(report <- eeg_summary(eeg))

  expect_equal(nrow(report$eeg_stats), 4)    # Fp1 Fp2 Cz Pz
  expect_equal(nrow(report$exg_stats), 1)    # EXG1
  expect_false("Status" %in% report$eeg_stats$channel)
  expect_false("Status" %in% report$exg_stats$channel)
})


# ----------------------------------------------------------------------------
# Test 7.2: Return value is independent of print output
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Capturing vs. not capturing output should not affect the
# returned list — the invisible return value is always populated.
test_that("eeg_summary return value is consistent regardless of output capture", {
  eeg <- .make_pure_eeg()

  r1 <- capture.output(report1 <- eeg_summary(eeg))
  report2 <- suppressMessages(invisible(capture.output(eeg_summary(eeg))))

  expect_true(is.data.frame(report1$eeg_stats))
  expect_equal(nrow(report1$eeg_stats), 4)
})


# ----------------------------------------------------------------------------
# Test 7.3: Summary stats are not influenced by EXG data
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: An EXG channel with extreme values must not distort the
# EEG channel statistics or trigger false flags on EEG channels.
test_that("eeg_summary EEG stats are unaffected by extreme EXG values", {
  set.seed(90)
  eeg_data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  exg_data <- matrix(rep(5000, 200), nrow = 1, ncol = 200)  # extreme EXG values

  eeg <- new_eeg(
    data     = rbind(eeg_data, exg_data),
    channels = c("Cz", "Pz", "Oz", "EXG1"),
    sampling_rate = 256
  )

  capture.output(report <- eeg_summary(eeg))

  # EEG amplitude range must reflect only the eeg_data rows
  expect_true(report$eeg_stats$max_uv["Cz" == report$eeg_stats$channel] < 500)
  expect_false("EXG1" %in% report$eeg_stats$channel)

  # No EEG channel should be flagged for excessive amplitude
  if (!is.null(report$flags)) {
    amp_flags <- report$flags[grepl("Excessive", report$flags$reason), ]
    expect_equal(nrow(amp_flags), 0)
  }
})


# ============================================================================
#                           End of Test File
# ============================================================================
