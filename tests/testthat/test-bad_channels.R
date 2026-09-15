# ============================================================================
#                       Test File for bad_channels.R
# ============================================================================
#
# Functions tested:
#   1. find_bad_channels()             - Exported orchestrator
#   2. flag_channel_quality()          - Internal: flat/amplitude/outlier
#   3. find_channel_neighbors()        - Internal: spatial neighbor geometry
#   4. flag_low_neighbor_correlation() - Internal: spatial-neighbor correlation
#   5. flag_lof_outlier_channels()     - Internal: global Local Outlier Factor
#
# Test suites:
#   1. flag_channel_quality()
#   2. find_channel_neighbors()
#   3. flag_low_neighbor_correlation()
#   4. flag_lof_outlier_channels()
#   5. find_bad_channels() input validation
#   6. find_bad_channels() with montage attached
#   7. find_bad_channels() without montage
#   8. Additive bads / preprocessing_history
#   9. return_details
#
# Author: Christos Dalamarinis
# Date: Sep 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
# Shared test fixtures
# ============================================================================

# 13 real 10-20 channels, spread across the scalp - enough for create_montage()
# to resolve and for k = 4 spatial neighbors to be meaningful.
.standard_channels <- c("Fp1", "Fp2", "F3", "F4", "Fz",
                         "C3", "C4", "Cz", "P3", "P4", "Pz", "O1", "O2")

# Plain, well-behaved eeg object: every channel is the same shared pattern
# plus independent small noise, so channels are mutually correlated and no
# check should flag anything.
.make_clean_eeg <- function(n_tp = 300, seed = 1) {
  set.seed(seed)
  shared <- rnorm(n_tp) * 15
  data <- t(vapply(seq_along(.standard_channels), function(i) {
    shared + rnorm(n_tp, sd = 2)
  }, numeric(n_tp)))
  new_eeg(data = data, channels = .standard_channels, sampling_rate = 256)
}

# Same as above, but with a montage already attached.
.make_clean_eeg_with_montage <- function(n_tp = 300, seed = 1) {
  eeg <- .make_clean_eeg(n_tp = n_tp, seed = seed)
  set_montage(eeg, create_montage(.standard_channels))
}

# A "rich" fixture with one problem per check:
#   O1 - flat (near-zero variance)
#   F3 - spatially decorrelated (independent noise, normal amplitude/variance)
#   Pz - global outlier (much larger variance than everything else)
# The remaining 10 channels are all the same shared pattern + small noise.
.make_problem_eeg <- function(n_tp = 400, seed = 42) {
  set.seed(seed)
  shared <- rnorm(n_tp) * 15

  data <- t(vapply(.standard_channels, function(ch) {
    if (ch == "O1") {
      rep(0.01, n_tp)                       # flat
    } else if (ch == "F3") {
      rnorm(n_tp, sd = 15)                  # independent, normal amplitude
    } else if (ch == "Pz") {
      rnorm(n_tp, sd = 200)                 # globally huge variance
    } else {
      shared + rnorm(n_tp, sd = 2)          # correlated with the group
    }
  }, numeric(n_tp)))

  new_eeg(data = data, channels = .standard_channels, sampling_rate = 256)
}

.make_problem_eeg_with_montage <- function(n_tp = 400, seed = 42) {
  eeg <- .make_problem_eeg(n_tp = n_tp, seed = seed)
  set_montage(eeg, create_montage(.standard_channels))
}


# ============================================================================
# TEST SUITE 1: flag_channel_quality()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1-1.3: flags flat, excessive amplitude, and outlier channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The three ported checks fire on the same kind of fixtures
# eeg_summary()'s own tests use (see test-eeg_summary.R suite 4).
test_that("flag_channel_quality flags flat, excessive-amplitude, and outlier channels", {
  set.seed(10)
  data <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data[2, ] <- rep(5, 200)  # Pz: perfectly flat
  eeg_flat <- new_eeg(data = data, channels = c("Cz", "Pz", "Oz"), sampling_rate = 256)
  flags_flat <- eeganalysis:::flag_channel_quality(eeg_flat)
  expect_true("Pz" %in% flags_flat$channel[grepl("Flat", flags_flat$reason)])

  set.seed(11)
  data2 <- matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200)
  data2[3, 100] <- 800  # Oz: one saturated sample
  eeg_amp <- new_eeg(data = data2, channels = c("Cz", "Pz", "Oz"), sampling_rate = 256)
  flags_amp <- eeganalysis:::flag_channel_quality(eeg_amp)
  expect_true("Oz" %in% flags_amp$channel[grepl("Excessive", flags_amp$reason)])

  set.seed(12)
  data3 <- rbind(matrix(rnorm(3 * 200, sd = 10), nrow = 3, ncol = 200),
                 matrix(rnorm(1 * 200, sd = 300), nrow = 1, ncol = 200))
  eeg_out <- new_eeg(data = data3, channels = c("Cz", "Pz", "Oz", "Fz"), sampling_rate = 256)
  flags_out <- eeganalysis:::flag_channel_quality(eeg_out)
  expect_true("Fz" %in% flags_out$channel[grepl("Outlier", flags_out$reason)])
})

# ----------------------------------------------------------------------------
# Test 1.4: Returns a correctly-shaped 0-row frame, never NULL, when clean
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Unlike eeg_summary()'s $flags (which becomes NULL when
# clean), this internal helper's contract is a 0-row data frame with the
# right columns - find_bad_channels() relies on never seeing NULL here.
test_that("flag_channel_quality returns a 0-row channel/reason/value frame when clean", {
  eeg <- .make_clean_eeg()
  flags <- eeganalysis:::flag_channel_quality(eeg)

  expect_s3_class(flags, "data.frame")
  expect_equal(names(flags), c("channel", "reason", "value"))
  expect_equal(nrow(flags), 0)
})

# ----------------------------------------------------------------------------
# Test 1.5: Excludes channels already in eeg_obj$bads, and non-EEG channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A flat channel already marked bad is not re-flagged, and
# EXG/Status channels are never flagged regardless of their signal.
test_that("flag_channel_quality excludes existing bads and non-EEG channels", {
  set.seed(13)
  data <- matrix(rnorm(4 * 200, sd = 15), nrow = 4, ncol = 200)
  data[2, ] <- rep(0, 200)  # Pz: flat, but will be pre-marked bad
  eeg <- new_eeg(data = data, channels = c("Cz", "Pz", "Oz", "EXG1"),
                 sampling_rate = 256, bads = "Pz")

  flags <- eeganalysis:::flag_channel_quality(eeg)

  expect_false("Pz" %in% flags$channel)
  expect_false("EXG1" %in% flags$channel)
})

# ----------------------------------------------------------------------------
# Test 1.6: Respects custom thresholds
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Mirrors test-eeg_summary.R 4.5 - a channel with std ~ 0.4
# is flagged at the default 0.5 threshold but not at a lowered 0.1 threshold.
test_that("flag_channel_quality respects custom flag_flat_threshold", {
  set.seed(30)
  data <- rbind(matrix(rnorm(3 * 200, sd = 15), nrow = 3, ncol = 200),
                rnorm(200, sd = 0.4))
  eeg <- new_eeg(data = data, channels = c("Cz", "Pz", "Oz", "Fz"), sampling_rate = 256)

  flags_default <- eeganalysis:::flag_channel_quality(eeg, flag_flat_threshold = 0.5)
  expect_true("Fz" %in% flags_default$channel[grepl("Flat", flags_default$reason)])

  flags_strict <- eeganalysis:::flag_channel_quality(eeg, flag_flat_threshold = 0.1)
  expect_false("Fz" %in% flags_strict$channel[grepl("Flat", flags_strict$reason)])
})

# ----------------------------------------------------------------------------
# Test 1.7: Parity with eeg_summary()'s current flagging
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Locks in behavioral equivalence between this helper and
# eeg_summary()'s own (still separate, print-only) flagging logic, so a
# future refactor of eeg_summary() to call this helper doesn't silently
# change behavior.
test_that("flag_channel_quality matches eeg_summary()'s flagged channel set", {
  fixtures <- list(
    .make_clean_eeg(),
    .make_problem_eeg()
  )

  for (eeg in fixtures) {
    helper_channels  <- sort(unique(eeganalysis:::flag_channel_quality(eeg)$channel))
    capture.output(report <- eeg_summary(eeg))
    summary_channels <- sort(unique(if (is.null(report$flags)) character(0) else report$flags$channel))

    expect_equal(helper_channels, summary_channels)
  }
})


# ============================================================================
# TEST SUITE 2: find_channel_neighbors()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1-2.3: correct neighbor count, ordering, never self
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With enough channels available, each gets exactly k
# neighbors, nearest-to-farthest, and never lists itself.
test_that("find_channel_neighbors returns k correctly-ordered neighbors, never self", {
  m <- create_montage(.standard_channels)
  result <- eeganalysis:::find_channel_neighbors(m$positions, k = 4)

  expect_equal(nrow(result), length(.standard_channels) * 4)

  for (ch in .standard_channels) {
    sub <- result[result$channel == ch, ]
    expect_equal(nrow(sub), 4)
    expect_equal(sub$rank, 1:4)
    expect_equal(sub$distance, sort(sub$distance))
    expect_false(ch %in% sub$neighbor)
  }
})

# ----------------------------------------------------------------------------
# Test 2.4: Domain plausibility on real BioSemi geometry
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: On the full 64-channel montage, Cz's nearest neighbors are
# central channels, never occipital ones on the opposite side of the head.
test_that("find_channel_neighbors gives spatially plausible neighbors on real geometry", {
  m <- create_montage()  # all 64 standard electrodes
  result <- eeganalysis:::find_channel_neighbors(m$positions, k = 4)

  cz_neighbors <- result$neighbor[result$channel == "Cz"]
  expect_false(any(cz_neighbors %in% c("O1", "O2", "Oz", "Iz")))
})

# ----------------------------------------------------------------------------
# Test 2.5: Caps k and warns when too few channels are available
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With only 3 channels and k = 4, each channel gets exactly
# 2 neighbors (the other 2 channels) and exactly one warning is raised.
test_that("find_channel_neighbors caps k and warns with too few channels", {
  m <- create_montage(c("Cz", "Fz", "Pz"))

  expect_warning(
    result <- eeganalysis:::find_channel_neighbors(m$positions, k = 4),
    "only 2 other channel"
  )
  expect_true(all(table(result$channel) == 2))
})

# ----------------------------------------------------------------------------
# Test 2.6: Empty-but-correctly-shaped result for fewer than 2 channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A single-channel input returns a 0-row frame with the
# right columns rather than erroring.
test_that("find_channel_neighbors returns 0 rows for fewer than 2 channels", {
  m <- create_montage("Cz")
  result <- eeganalysis:::find_channel_neighbors(m$positions, k = 4)

  expect_equal(names(result), c("channel", "neighbor", "rank", "distance"))
  expect_equal(nrow(result), 0)
})

# ----------------------------------------------------------------------------
# Test 2.7: Regression - ranking follows Cartesian x/y/z, not spherical cols
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: theta/phi/radius are present but deliberately contradict
# the x/y/z-based ranking; the function must follow x/y/z only (radius is a
# fixed constant in the real electrode database and carries no information).
test_that("find_channel_neighbors uses Cartesian x/y/z, ignoring theta/phi/radius", {
  positions <- data.frame(
    channel = c("A", "B", "C"),
    x = c(0, 1, 10), y = c(0, 0, 0), z = c(0, 0, 0),
    # theta/phi deliberately contradict the x/y/z-based ranking above
    theta = c(0, 180, 90), phi = c(0, 90, 45), radius = c(87.54, 87.54, 87.54),
    stringsAsFactors = FALSE
  )

  result <- eeganalysis:::find_channel_neighbors(positions, k = 1)

  # By x/y/z, A's nearest is B (distance 1), not C (distance 10)
  expect_equal(result$neighbor[result$channel == "A"], "B")
})


# ============================================================================
# TEST SUITE 3: flag_low_neighbor_correlation()
# ============================================================================

# Shared neighbor map: A's neighbors are B, C, D, E (hand-built, so this
# suite is isolated from find_channel_neighbors()).
.corr_neighbor_map <- data.frame(
  channel = "A", neighbor = c("B", "C", "D", "E"),
  rank = 1:4, distance = 1:4, stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# Test 3.1-3.2: flags a decorrelated channel, not a well-correlated one
# ----------------------------------------------------------------------------
test_that("flag_low_neighbor_correlation flags decorrelated, not correlated, channels", {
  set.seed(50)
  n_tp <- 300
  shared <- rnorm(n_tp) * 10
  neighbor_data <- rbind(
    B = shared + rnorm(n_tp, sd = 1),
    C = shared + rnorm(n_tp, sd = 1),
    D = shared + rnorm(n_tp, sd = 1),
    E = shared + rnorm(n_tp, sd = 1)
  )

  a_correlated   <- shared + rnorm(n_tp, sd = 1)
  a_decorrelated <- rnorm(n_tp, sd = 10)  # unrelated to shared

  data_good <- rbind(A = a_correlated, neighbor_data)
  flags_good <- eeganalysis:::flag_low_neighbor_correlation(
    data = data_good, channels = rownames(data_good), neighbor_map = .corr_neighbor_map)
  expect_false("A" %in% flags_good$channel)

  data_bad <- rbind(A = a_decorrelated, neighbor_data)
  flags_bad <- eeganalysis:::flag_low_neighbor_correlation(
    data = data_bad, channels = rownames(data_bad), neighbor_map = .corr_neighbor_map)
  expect_true("A" %in% flags_bad$channel)
  expect_true(grepl("Low correlation", flags_bad$reason[flags_bad$channel == "A"]))
})

# ----------------------------------------------------------------------------
# Test 3.3-3.4: robust (median) vs non-robust (mean) give different results
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: One neighbor (E) is a huge-amplitude, uncorrelated
# channel. The median reference (4 neighbors, so the outlier lands at an
# extreme rank) stays close to the shared pattern; the mean reference is
# dragged toward E's noise. A should therefore correlate better with the
# robust (median) reference than the non-robust (mean) one.
test_that("robust = TRUE (median) is less affected by one extreme neighbor than robust = FALSE (mean)", {
  set.seed(51)
  n_tp <- 400
  shared <- rnorm(n_tp) * 10

  neighbor_data <- rbind(
    B = shared + rnorm(n_tp, sd = 1),
    C = shared + rnorm(n_tp, sd = 1),
    D = shared + rnorm(n_tp, sd = 1),
    E = rnorm(n_tp, sd = 200)  # huge, unrelated to shared
  )
  a <- shared + rnorm(n_tp, sd = 1)
  data <- rbind(A = a, neighbor_data)
  channels <- rownames(data)

  ref_median <- apply(neighbor_data, 2, median)
  ref_mean   <- colMeans(neighbor_data)

  cor_median <- cor(a, ref_median)
  cor_mean   <- cor(a, ref_mean)

  expect_gt(cor_median, cor_mean)

  # And the helper itself reflects this: at a threshold between the two
  # correlations, robust = TRUE does not flag A while robust = FALSE does.
  threshold <- mean(c(cor_median, cor_mean))
  flags_robust <- eeganalysis:::flag_low_neighbor_correlation(
    data, channels, .corr_neighbor_map, correlation_threshold = threshold, robust = TRUE)
  flags_plain <- eeganalysis:::flag_low_neighbor_correlation(
    data, channels, .corr_neighbor_map, correlation_threshold = threshold, robust = FALSE)

  expect_false("A" %in% flags_robust$channel)
  expect_true("A" %in% flags_plain$channel)
})

# ----------------------------------------------------------------------------
# Test 3.5: Respects a custom correlation_threshold
# ----------------------------------------------------------------------------
test_that("flag_low_neighbor_correlation respects a custom correlation_threshold", {
  set.seed(52)
  n_tp <- 300
  shared <- rnorm(n_tp) * 10
  neighbor_data <- rbind(
    B = shared + rnorm(n_tp, sd = 1), C = shared + rnorm(n_tp, sd = 1),
    D = shared + rnorm(n_tp, sd = 1), E = shared + rnorm(n_tp, sd = 1)
  )
  a <- 0.5 * shared + rnorm(n_tp, sd = 5)  # moderate correlation
  data <- rbind(A = a, neighbor_data)
  channels <- rownames(data)

  actual_cor <- cor(a, apply(neighbor_data, 2, median))

  flags_lenient <- eeganalysis:::flag_low_neighbor_correlation(
    data, channels, .corr_neighbor_map, correlation_threshold = actual_cor - 0.2)
  flags_strict <- eeganalysis:::flag_low_neighbor_correlation(
    data, channels, .corr_neighbor_map, correlation_threshold = actual_cor + 0.2)

  expect_false("A" %in% flags_lenient$channel)
  expect_true("A" %in% flags_strict$channel)
})

# ----------------------------------------------------------------------------
# Test 3.6: Does not crash on a zero-variance channel or neighbor
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: cor() returns NA for a constant vector; the helper must
# treat that as flag-worthy rather than crashing on an unguarded comparison.
test_that("flag_low_neighbor_correlation does not crash on zero-variance signals", {
  n_tp <- 100
  data <- rbind(
    A = rep(0, n_tp),               # constant channel
    B = rnorm(n_tp), C = rnorm(n_tp), D = rnorm(n_tp), E = rnorm(n_tp)
  )
  channels <- rownames(data)

  expect_no_error(
    flags <- eeganalysis:::flag_low_neighbor_correlation(data, channels, .corr_neighbor_map)
  )
  expect_true("A" %in% flags$channel)
})

# ----------------------------------------------------------------------------
# Test 3.7: A channel with no available neighbors is skipped, not an error
# ----------------------------------------------------------------------------
test_that("flag_low_neighbor_correlation skips a channel whose neighbors aren't in the data", {
  n_tp <- 100
  data <- rbind(A = rnorm(n_tp), B = rnorm(n_tp))
  channels <- rownames(data)
  # .corr_neighbor_map's neighbors for "A" (B, C, D, E) are only partially
  # present in `channels` here (C, D, E are missing) - B alone remains, so
  # this should just run using B, not error.
  expect_no_error(
    eeganalysis:::flag_low_neighbor_correlation(data, channels, .corr_neighbor_map)
  )
})


# ============================================================================
# TEST SUITE 4: flag_lof_outlier_channels()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1-4.2: flags an obvious global outlier, clean on homogeneous data
# ----------------------------------------------------------------------------
test_that("flag_lof_outlier_channels flags an obvious global outlier and nothing on clean data", {
  set.seed(60)
  n_tp <- 500
  normal_data <- matrix(rnorm(7 * n_tp, sd = 10), nrow = 7, ncol = n_tp)
  channels_clean <- paste0("Ch", 1:7)
  flags_clean <- eeganalysis:::flag_lof_outlier_channels(normal_data, channels_clean, n_neighbors = 4)
  expect_equal(nrow(flags_clean), 0)

  outlier_data <- rbind(normal_data, matrix(rnorm(1 * n_tp, sd = 300), nrow = 1, ncol = n_tp))
  channels_outlier <- paste0("Ch", 1:8)
  flags_outlier <- eeganalysis:::flag_lof_outlier_channels(outlier_data, channels_outlier, n_neighbors = 4)
  expect_true("Ch8" %in% flags_outlier$channel)
})

# ----------------------------------------------------------------------------
# Test 4.3-4.4: caps and warns with too few channels; no warning otherwise
# ----------------------------------------------------------------------------
test_that("flag_lof_outlier_channels caps n_neighbors only when necessary", {
  set.seed(61)
  small_data <- matrix(rnorm(4 * 100), nrow = 4, ncol = 100)
  expect_warning(
    eeganalysis:::flag_lof_outlier_channels(small_data, paste0("Ch", 1:4), n_neighbors = 20),
    "only 3 other channel"
  )

  enough_data <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
  expect_no_warning(
    eeganalysis:::flag_lof_outlier_channels(enough_data, paste0("Ch", 1:10), n_neighbors = 4)
  )
})

# ----------------------------------------------------------------------------
# Test 4.5: Orientation regression - flags a channel (row), not a timepoint
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A wide matrix (far more timepoints than channels) with
# one outlier row must flag exactly one entry from the channel vocabulary -
# guards against an accidental transpose in flag_lof_outlier_channels().
test_that("flag_lof_outlier_channels flags channels, not timepoints (orientation guard)", {
  set.seed(62)
  n_tp <- 500
  data <- rbind(
    matrix(rnorm(5 * n_tp, sd = 10), nrow = 5, ncol = n_tp),
    matrix(rnorm(1 * n_tp, sd = 300), nrow = 1, ncol = n_tp)
  )
  channels <- c("Cz", "Pz", "Oz", "Fz", "Fp1", "OddOne")

  flags <- eeganalysis:::flag_lof_outlier_channels(data, channels, n_neighbors = 3)

  expect_equal(flags$channel, "OddOne")
  expect_true(all(flags$channel %in% channels))
})


# ============================================================================
# TEST SUITE 5: find_bad_channels() input validation
# ============================================================================

test_that("find_bad_channels errors on non-eeg input", {
  expect_error(find_bad_channels(list()), "class 'eeg'")
})

test_that("find_bad_channels errors when eeg_obj$data is not a matrix", {
  eeg <- .make_clean_eeg()
  eeg$data <- as.data.frame(eeg$data)
  expect_error(find_bad_channels(eeg), "must be a numeric matrix")
})

test_that("find_bad_channels errors when no EEG channels remain", {
  set.seed(2)
  data <- matrix(rnorm(2 * 100), nrow = 2, ncol = 100)
  eeg_all_exg <- new_eeg(data = data, channels = c("EXG1", "EXG2"), sampling_rate = 256)
  expect_error(find_bad_channels(eeg_all_exg), "No EEG channels")

  eeg_all_bad <- new_eeg(data = data, channels = c("Cz", "Pz"), sampling_rate = 256,
                          bads = c("Cz", "Pz"))
  expect_error(find_bad_channels(eeg_all_bad), "No EEG channels")
})


# ============================================================================
# TEST SUITE 6: find_bad_channels() with montage attached
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1-6.4: all four checks run and catch their respective problems
# ----------------------------------------------------------------------------
test_that("find_bad_channels catches flat, decorrelated, and globally-outlying channels", {
  eeg <- .make_problem_eeg_with_montage()

  result <- find_bad_channels(eeg, lof_n_neighbors = 6)

  expect_s3_class(result, "eeg")
  expect_true("O1" %in% result$bads)  # flat
  expect_true("F3" %in% result$bads)  # spatially decorrelated
  expect_true("Pz" %in% result$bads)  # global outlier
})

# ----------------------------------------------------------------------------
# Test 6.5: Warns but completes when some channels lack a montage position
# ----------------------------------------------------------------------------
test_that("find_bad_channels warns but completes when some channels lack a montage position", {
  eeg <- .make_clean_eeg()
  partial_montage <- create_montage(.standard_channels[1:6])
  eeg <- suppressWarnings(set_montage(eeg, partial_montage))
  # set_montage() restricts $montage to the overlap, so eeg's remaining
  # channels (7-13) now have no montage position at all.

  expect_warning(
    result <- find_bad_channels(eeg, lof_n_neighbors = 6),
    "no position in the attached montage"
  )
  expect_s3_class(result, "eeg")
})


# ============================================================================
# TEST SUITE 7: find_bad_channels() without montage
# ============================================================================

test_that("find_bad_channels messages (not errors) and still runs the other checks with no montage", {
  eeg <- .make_problem_eeg()  # no set_montage() call

  expect_message(
    result <- find_bad_channels(eeg, lof_n_neighbors = 6),
    "no montage attached"
  )

  expect_true("O1" %in% result$bads)  # flat - doesn't need montage
  expect_true("Pz" %in% result$bads)  # global LOF - doesn't need montage
})


# ============================================================================
# TEST SUITE 8: Additive bads / preprocessing_history
# ============================================================================

test_that("find_bad_channels adds to existing bads without removing them, and logs history", {
  eeg <- .make_clean_eeg_with_montage()
  eeg$bads <- "Fp2"  # manually marked, would not itself be re-flagged (clean data)

  n_history_before <- length(eeg$preprocessing_history)
  result <- find_bad_channels(eeg, lof_n_neighbors = 6)

  expect_true("Fp2" %in% result$bads)  # survived
  expect_equal(length(result$preprocessing_history), n_history_before + 1)
  expect_true(grepl("find_bad_channels", result$preprocessing_history[[length(result$preprocessing_history)]]))
})

test_that("find_bad_channels appends exactly one history entry even when nothing new is flagged", {
  eeg <- .make_clean_eeg_with_montage()
  result <- find_bad_channels(eeg, lof_n_neighbors = 6)

  last_entry <- result$preprocessing_history[[length(result$preprocessing_history)]]
  expect_true(grepl("0 channel", last_entry))
  expect_equal(length(result$bads), 0)
})


# ============================================================================
# TEST SUITE 9: return_details
# ============================================================================

test_that("return_details controls the return shape", {
  eeg <- .make_problem_eeg_with_montage()

  default_result <- find_bad_channels(eeg, lof_n_neighbors = 6)
  expect_s3_class(default_result, "eeg")

  detailed_result <- find_bad_channels(eeg, lof_n_neighbors = 6, return_details = TRUE)
  expect_type(detailed_result, "list")
  expect_s3_class(detailed_result$eeg_obj, "eeg")
  expect_s3_class(detailed_result$flags, "data.frame")

  newly_flagged <- setdiff(detailed_result$eeg_obj$bads, eeg$bads)
  expect_setequal(unique(detailed_result$flags$channel), newly_flagged)
})


# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - flag_channel_quality(): flat/amplitude/outlier flagging, 0-row contract,
#   bads/non-EEG exclusion, custom thresholds, parity with eeg_summary()
# - find_channel_neighbors(): correct k/ordering/no-self, spatial plausibility,
#   capping + warning, 0-row contract, Cartesian-not-spherical regression
# - flag_low_neighbor_correlation(): decorrelated flagging, robust vs plain
#   aggregation, custom threshold, zero-variance crash guard, missing-neighbor
#   skip
# - flag_lof_outlier_channels(): global outlier flagging, capping + warning,
#   orientation (channel- not timepoint-) regression
# - find_bad_channels(): input validation, all four checks end-to-end with
#   and without a montage, additive bads, preprocessing_history, return_details
# ============================================================================
