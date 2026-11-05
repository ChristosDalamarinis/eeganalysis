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
#   6. .column_medians()               - Internal: middle value of every column at once
#   7. .centre_by_stretch()            - Internal: subtract each channel's own level
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
#  10. .column_medians()
#  11. .centre_by_stretch()
#  12. The amplitude check measures from each channel's own level
#  13. Raw BioSemi-style data, end to end
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
  eeg_all_exg <- new_eeg(data = data, channels = c("EXG1", "EXG2"), sampling_rate = 256,
                         channel_types = c("eog", "eog"))
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
# Fixtures for suites 10-13: raw BioSemi-style data
# ============================================================================
#
# Raw BioSemi data keeps a constant offset on every electrode (tens of mV =
# thousands of uV) that also drifts slowly; the brain signal is only a wiggle of
# a few tens of uV around it. The two fixtures below build that:
#   .make_offset_eeg()  - n_ch channels "Ch1".."ChN", each = its own offset
#                         (spread from -28000 to +28000 uV) + an optional slow
#                         drift (0.02 Hz sine, own phase per channel) + a small
#                         normal wiggle (default sd 20 uV)
#   .make_biosemi_eeg() - the 64 BioSemi electrodes with the real layout and
#                         montage, neighbours that move together, one broken
#                         channel (wiggle sd 300 uV), offsets and optional drift

.make_offset_eeg <- function(n_ch = 10, seconds = 10, fs = 256, drift = 0,
                             seed = 7, wiggle_sd = rep(20, n_ch)) {
  set.seed(seed)
  n  <- seconds * fs
  tt <- (0:(n - 1)) / fs
  data <- t(vapply(seq_len(n_ch), function(i) {
    w <- as.numeric(stats::filter(rnorm(n), 0.9, "recursive"))
    w <- w / sd(w) * wiggle_sd[i]
    offset <- -28000 + (i - 1) * 56000 / (n_ch - 1)
    offset + drift * sin(2 * pi * 0.02 * tt + i) + w
  }, numeric(n)))
  new_eeg(data = data, channels = paste0("Ch", seq_len(n_ch)),
          sampling_rate = fs, channel_types = rep("eeg", n_ch))
}

.make_biosemi_eeg <- function(seconds = 5, drift = 0, offsets = TRUE,
                              bad = "T8", seed = 2) {
  fs  <- 512
  lay <- eeganalysis:::.plot_layout_table("biosemi64")
  nm  <- lay$electrode
  n   <- seconds * fs
  tt  <- (0:(n - 1)) / fs
  pos <- as.matrix(lay[, c("x", "y", "z")])
  pos <- pos / max(abs(pos))
  set.seed(seed)
  src  <- scale(vapply(1:6, function(j) {
    as.numeric(stats::filter(rnorm(n), 0.97, "recursive"))
  }, numeric(n))) * 12
  load <- vapply(1:6, function(j) {
    as.numeric(cos(pos %*% rnorm(3, sd = 2) + runif(1, 0, 2 * pi)))
  }, numeric(64))
  X <- load %*% t(src) + matrix(rnorm(64 * n, sd = 1.5), 64, n)
  w <- as.numeric(stats::filter(rnorm(n), 0.95, "recursive"))
  X[match(bad, nm), ] <- w / sd(w) * 300
  if (drift > 0) {
    X <- X + t(vapply(1:64, function(i) drift * sin(2 * pi * 0.02 * tt + i), numeric(n)))
  }
  if (offsets) {
    X <- X + runif(64, -30000, 30000)
  }
  eeg <- new_eeg(data = X, channels = nm, sampling_rate = fs,
                 channel_types = rep("eeg", 64))
  set_montage(eeg, "biosemi64")
}

.excessive_channels <- function(flags) {
  flags$channel[grepl("Excessive", flags$reason)]
}


# ============================================================================
# TEST SUITE 10: .column_medians()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 10.1: hand-worked middle values
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The middle value of 4 neighbours at one moment, with one
# broken neighbour (500): 10, 11, 12, 500 -> (11 + 12) / 2 = 11.5, while the
# mean would be 133.25. Odd counts take the middle one.
test_that(".column_medians takes the middle value of every column", {
  m <- cbind(c(10, 12, 11, 500), c(4, 3, 2, 1), c(5, 5, 5, 5))
  expect_equal(eeganalysis:::.column_medians(m), c(11.5, 2.5, 5))

  odd <- cbind(c(3, 1, 2), c(-1, -5, 9))
  expect_equal(eeganalysis:::.column_medians(odd), c(2, -1))

  expect_equal(eeganalysis:::.column_medians(matrix(c(7, 8, 9), nrow = 1)), c(7, 8, 9))
  expect_equal(eeganalysis:::.column_medians(rbind(c(1, 2), c(3, 6))), c(2, 4))
})

# ----------------------------------------------------------------------------
# Test 10.2: same answer as apply(m, 2, median)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: It is a faster way to the same numbers, for 1 to 8 rows,
# unsorted values, negative values, ties, whole numbers, a single column.
test_that(".column_medians gives the same numbers as apply(m, 2, median)", {
  set.seed(3)
  for (k in 1:8) {
    m <- matrix(rnorm(k * 60, sd = 50), nrow = k)
    expect_equal(eeganalysis:::.column_medians(m), apply(m, 2, median),
                 tolerance = 1e-14, info = paste("k =", k))
  }
  ties <- matrix(sample(1:3, 4 * 40, replace = TRUE), nrow = 4)
  expect_equal(eeganalysis:::.column_medians(ties), apply(ties, 2, median))

  ints <- matrix(1:12, nrow = 3)
  expect_equal(eeganalysis:::.column_medians(ints), apply(ints, 2, median))

  one_col <- matrix(c(9, 1, 5, 3), ncol = 1)
  expect_equal(eeganalysis:::.column_medians(one_col), 4)

  m <- matrix(rnorm(5 * 30), nrow = 5)
  expect_equal(eeganalysis:::.column_medians(m[sample(5), ]),
               eeganalysis:::.column_medians(m))            # row order does not matter
})

# ----------------------------------------------------------------------------
# Test 10.3: missing values, infinite values, empty matrices
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: A column with a missing value gives NA (like median() with
# its default), only for that column - the neighbor check relies on it, because
# an NA reference makes the channel get flagged. Infinite values and empty
# matrices behave like median().
test_that(".column_medians handles NA, NaN, Inf and empty matrices like median()", {
  m <- cbind(c(1, 2, 3, 4), c(1, NA, 3, 4), c(5, 6, NaN, 8), c(2, 4, 6, 8))
  res <- eeganalysis:::.column_medians(m)
  expect_equal(res[c(1, 4)], c(2.5, 5))
  expect_true(is.na(res[2]) && is.na(res[3]))
  expect_equal(res, apply(m, 2, median))

  odd <- cbind(c(1, NA, 3), c(1, 2, 3))
  expect_equal(eeganalysis:::.column_medians(odd), c(NA, 2))

  inf <- cbind(c(-Inf, Inf, 1, 2), c(Inf, Inf, 5, 6), c(-Inf, Inf, 0, 0)[1:4])
  expect_equal(eeganalysis:::.column_medians(inf), apply(inf, 2, median))
  expect_equal(eeganalysis:::.column_medians(rbind(c(-Inf, Inf), c(Inf, Inf))),
               apply(rbind(c(-Inf, Inf), c(Inf, Inf)), 2, median))

  expect_equal(eeganalysis:::.column_medians(matrix(0, nrow = 3, ncol = 0)), numeric(0))
  expect_equal(eeganalysis:::.column_medians(matrix(0, nrow = 0, ncol = 4)), rep(NA_real_, 4))
})


# ============================================================================
# TEST SUITE 11: .centre_by_stretch()
# ============================================================================

# ----------------------------------------------------------------------------
# Test 11.1: the worked examples (hand numbers)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: One stretch: 8 numbers standing near 12630 have the average
# 12628.125 subtracted. Three stretches of 4 numbers (the level creeps up):
# each stretch loses its OWN average (12001.75, 12684.75, 13382.5), so what is
# left is only the wiggle - one average for everything (12689.67) would leave
# -690 .. +720.
test_that(".centre_by_stretch subtracts the channel's own average in every stretch", {
  a <- c(12609, 12665, 12599, 12632, 12673, 12615, 12618, 12614)
  out <- eeganalysis:::.centre_by_stretch(matrix(a, nrow = 1), window_samples = 8)
  expect_equal(as.numeric(out), a - 12628.125)
  expect_equal(as.numeric(out)[1:3], c(-19.125, 36.875, -29.125))

  s1 <- c(11993, 12003, 12031, 11980)
  s2 <- c(12673, 12696, 12673, 12697)
  s3 <- c(13385, 13345, 13406, 13394)
  out3 <- eeganalysis:::.centre_by_stretch(matrix(c(s1, s2, s3), nrow = 1), window_samples = 4)
  expect_equal(as.numeric(out3),
               c(s1 - 12001.75, s2 - 12684.75, s3 - 13382.5))
  expect_lt(max(abs(out3)), 40)

  one_average <- c(s1, s2, s3) - mean(c(s1, s2, s3))
  expect_gt(max(abs(one_average)), 700)                      # why one average is not enough
})

# ----------------------------------------------------------------------------
# Test 11.2: where the stretches start and end
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: As many whole stretches as fit; the leftover joins the last
# stretch (never a tiny stretch of its own); a recording shorter than one
# stretch is one stretch; a window of one sample leaves nothing.
test_that(".centre_by_stretch: leftover joins the last stretch, short recordings are one stretch", {
  x <- matrix(1:11, nrow = 1)
  out <- eeganalysis:::.centre_by_stretch(x, window_samples = 4)         # 2 stretches: 1:4 and 5:11
  expect_equal(as.numeric(out), c(-1.5, -0.5, 0.5, 1.5, -3, -2, -1, 0, 1, 2, 3))

  expect_equal(as.numeric(eeganalysis:::.centre_by_stretch(matrix(c(1, 2, 6), nrow = 1), 10)),
               c(-2, -1, 3))                                             # shorter than the window
  expect_equal(as.numeric(eeganalysis:::.centre_by_stretch(matrix(1:8, nrow = 1), 4)),
               c(-1.5, -0.5, 0.5, 1.5, -1.5, -0.5, 0.5, 1.5))            # exactly 2 stretches
  expect_equal(as.numeric(eeganalysis:::.centre_by_stretch(matrix(c(5, 9, 7), nrow = 1), 3)),
               c(-2, 2, 0))                                              # window = recording
  expect_equal(as.numeric(eeganalysis:::.centre_by_stretch(matrix(c(5, 9, 7), nrow = 1), 1)),
               c(0, 0, 0))                                               # one sample: nothing left
})

# ----------------------------------------------------------------------------
# Test 11.3: shape, rows, properties
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Same size and names; every row is handled on its own; a
# constant channel becomes all zero; in every stretch of every row the average
# is 0; the wiggle itself is unchanged.
test_that(".centre_by_stretch keeps the shape and handles every row on its own", {
  set.seed(9)
  x <- matrix(rnorm(3 * 1000, sd = 10), nrow = 3,
              dimnames = list(c("a", "b", "c"), NULL)) + c(5000, -20000, 1)
  x["c", ] <- 42                                              # a constant channel
  out <- eeganalysis:::.centre_by_stretch(x, window_samples = 128)

  expect_equal(dim(out), dim(x))
  expect_equal(rownames(out), c("a", "b", "c"))
  expect_equal(as.numeric(out["c", ]), rep(0, 1000))

  starts <- c(seq(1, 7 * 128, by = 128))                      # 7 stretches; the last one is 232 long
  ends   <- c(starts[-1] - 1, 1000)
  for (s in seq_along(starts)) {
    idx <- starts[s]:ends[s]
    expect_equal(unname(rowMeans(out[, idx])), c(0, 0, 0), tolerance = 1e-9, info = paste("stretch", s))
  }
  # only the level is removed: the wiggle inside a stretch is untouched
  expect_equal(out["a", 1:128] - out["a", 1], x["a", 1:128] - x["a", 1])

  # the other rows do not change a row's result
  alone <- eeganalysis:::.centre_by_stretch(x["a", , drop = FALSE], 128)
  expect_equal(as.numeric(alone), as.numeric(out["a", ]))

  expect_equal(dim(eeganalysis:::.centre_by_stretch(matrix(0, nrow = 2, ncol = 0), 4)), c(2L, 0L))
})


# ============================================================================
# TEST SUITE 12: the amplitude check measures from each channel's own level
# ============================================================================

# ----------------------------------------------------------------------------
# Test 12.1: offsets alone must not flag anything
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Ten good channels standing at -28000 .. +28000 uV with a
# normal 20 uV wiggle: before this check measured from the channel's own level,
# all ten were flagged "Excessive amplitude" because their raw numbers are far
# above 500 uV.
test_that("channels that only stand on a large offset are not flagged", {
  eeg <- .make_offset_eeg()
  expect_gt(max(abs(eeg$data)), 20000)                         # the offsets are really there
  flags <- eeganalysis:::flag_channel_quality(eeg)
  expect_equal(nrow(flags), 0)
  expect_equal(nrow(eeganalysis:::flag_channel_quality(eeg, amplitude_window = Inf)), 0)
})

# ----------------------------------------------------------------------------
# Test 12.2: a slow drift of the level, and why it is re-measured
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The level drifts (2000 uV up and down, 0.02 Hz). Measured
# every 2 s, the channel is fine; measured once for the whole 10 s recording,
# the drift swings several good channels past 500 uV.
test_that("a slowly drifting level is followed by the 2 s stretches, but not by one average", {
  eeg <- .make_offset_eeg(drift = 2000)
  expect_equal(length(.excessive_channels(eeganalysis:::flag_channel_quality(eeg))), 0)
  expect_equal(length(.excessive_channels(eeganalysis:::flag_channel_quality(eeg, amplitude_window = 2))), 0)
  one_average <- .excessive_channels(eeganalysis:::flag_channel_quality(eeg, amplitude_window = Inf))
  expect_gte(length(one_average), 3)
  # a window longer than the recording is one stretch, the same as Inf (same
  # channels, same reasons, same numbers; only the words of the flag differ)
  long_window <- eeganalysis:::flag_channel_quality(eeg, amplitude_window = 100)
  one_stretch <- eeganalysis:::flag_channel_quality(eeg, amplitude_window = Inf)
  expect_equal(long_window[, c("channel", "reason")], one_stretch[, c("channel", "reason")])
  expect_equal(sub(".* = ", "", long_window$value), sub(".* = ", "", one_stretch$value))
})

# ----------------------------------------------------------------------------
# Test 12.3: what must still be flagged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With offsets everywhere, a really broken channel (wiggle of
# several hundred uV) and a single 800 uV spike are still "Excessive
# amplitude", and nothing else is.
test_that("a broken channel and a big spike are still flagged when every channel has an offset", {
  sds <- rep(20, 10)
  sds[4] <- 300
  broken <- .make_offset_eeg(wiggle_sd = sds)
  expect_equal(.excessive_channels(eeganalysis:::flag_channel_quality(broken)), "Ch4")

  spike <- .make_offset_eeg()
  spike$data[7, 1000] <- spike$data[7, 1000] + 800
  flags <- eeganalysis:::flag_channel_quality(spike)
  expect_equal(.excessive_channels(flags), "Ch7")
  peak <- as.numeric(sub(" uV.*", "", sub(".* = ", "", flags$value[grepl("Excessive", flags$reason)])))
  expect_gt(peak, 700)                                         # about the 800 uV of the spike,
  expect_lt(peak, 900)                                         # not the 12000+ of the offset
})

# ----------------------------------------------------------------------------
# Test 12.4: the other two checks do not depend on the offset
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Flat and outlier-noise checks use the std, which a constant
# level does not change: a flat channel standing at 12000 uV is flagged flat and
# only flat; the flags for the same data with and without offsets are the same.
test_that("the flat and noise checks give the same answer with or without offsets", {
  sds <- rep(20, 10)
  sds[2] <- 0.01
  sds[9] <- 200
  eeg <- .make_offset_eeg(wiggle_sd = sds)
  flags <- eeganalysis:::flag_channel_quality(eeg)
  expect_equal(flags$channel[grepl("Flat", flags$reason)], "Ch2")
  expect_false("Ch2" %in% .excessive_channels(flags))

  centred <- eeg
  centred$data <- eeg$data - rowMeans(eeg$data)
  flags_centred <- eeganalysis:::flag_channel_quality(centred)
  keep <- function(f) f[!grepl("Excessive", f$reason), c("channel", "reason")]
  expect_equal(keep(flags), keep(flags_centred))
  expect_true("Ch9" %in% flags$channel[grepl("Outlier", flags$reason)])
})

# ----------------------------------------------------------------------------
# Test 12.4b: the noise check still looks at the whole signal
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Only the amplitude check measures from the channel's own
# level; the std (flat and noise checks) is taken from the signal as it is, as
# before. One channel gets a slow drift of 400 uV: its std (about 280 uV) stands
# out from the others (about 20 uV) and it is flagged "Outlier noise level",
# while its peak amplitude from its own level (about 100 uV) is not excessive.
test_that("the noise check still judges the whole signal, a strong slow drift included", {
  eeg <- .make_offset_eeg()
  tt  <- (0:(ncol(eeg$data) - 1)) / 256
  eeg$data[5, ] <- eeg$data[5, ] + 400 * sin(2 * pi * 0.02 * tt)

  flags <- eeganalysis:::flag_channel_quality(eeg)
  expect_equal(flags$channel[grepl("Outlier", flags$reason)], "Ch5")
  expect_equal(length(.excessive_channels(flags)), 0)
})

# ----------------------------------------------------------------------------
# Test 12.5: the wording of the flag
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The flag text says what is measured: from the channel's own
# level, re-measured every window; with Inf, only "its own level".
test_that("the amplitude flag says it is measured from the channel's own level", {
  sds <- rep(20, 10)
  sds[4] <- 300
  eeg <- .make_offset_eeg(wiggle_sd = sds)
  value <- function(...) {
    f <- eeganalysis:::flag_channel_quality(eeg, ...)
    f$value[grepl("Excessive", f$reason)]
  }
  expect_match(value(), "^peak \\|amplitude\\| from its own level, re-measured every 2 s = ")
  expect_match(value(), "uV  \\(threshold: > 500 uV\\)$")
  expect_match(value(amplitude_window = 0.5), "re-measured every 0.5 s")
  expect_match(value(amplitude_window = Inf), "^peak \\|amplitude\\| from its own level = ")
  expect_false(grepl("re-measured", value(amplitude_window = Inf)))
})

# ----------------------------------------------------------------------------
# Test 12.6: amplitude_window: default and checks
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: 2 is the default; wrong values are refused with a clear
# message, also when no channel is eligible; the shortest window is 2 samples.
test_that("amplitude_window defaults to 2 s and is validated", {
  eeg <- .make_offset_eeg(drift = 2000, wiggle_sd = c(rep(20, 3), 300, rep(20, 6)))
  expect_equal(eeganalysis:::flag_channel_quality(eeg),
               eeganalysis:::flag_channel_quality(eeg, amplitude_window = 2))

  bad_values <- list(0, -1, NA, NA_real_, NaN, "a", c(2, 3), numeric(0), NULL, TRUE)
  for (v in bad_values) {
    expect_error(eeganalysis:::flag_channel_quality(eeg, amplitude_window = v),
                 "'amplitude_window' must be a single number above 0", info = deparse(v))
  }
  expect_error(eeganalysis:::flag_channel_quality(eeg, amplitude_window = 0.001),
               "shorter than 2 samples at 256 Hz")
  expect_error(eeganalysis:::flag_channel_quality(eeg, amplitude_window = 1 / 256),
               "shorter than 2 samples")
  expect_error(eeganalysis:::flag_channel_quality(eeg, amplitude_window = 2 / 256), NA)   # 2 samples: allowed

  all_bad <- eeg
  all_bad$bads <- all_bad$channels
  expect_error(eeganalysis:::flag_channel_quality(all_bad, amplitude_window = -1),
               "'amplitude_window' must be a single number above 0")
  expect_equal(nrow(eeganalysis:::flag_channel_quality(all_bad)), 0)
})


# ============================================================================
# TEST SUITE 13: raw BioSemi-style data, end to end
# ============================================================================

# ----------------------------------------------------------------------------
# Test 13.1: find_bad_channels() on 64 channels with offsets
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: 64 BioSemi electrodes, every one on an offset of up to
# +-30000 uV, one broken channel (T8, wiggle sd 300 uV). Before the fix the
# amplitude check flagged all 64; now it flags T8 only, and T8 reaches bads.
test_that("find_bad_channels flags only the broken channel for amplitude on raw offset data", {
  eeg <- .make_biosemi_eeg()
  expect_gt(max(abs(eeg$data)), 20000)

  result <- suppressMessages(find_bad_channels(eeg, return_details = TRUE))
  expect_equal(.excessive_channels(result$flags), "T8")
  expect_true("T8" %in% result$eeg_obj$bads)
})

# ----------------------------------------------------------------------------
# Test 13.2: amplitude_window is passed on by find_bad_channels()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: With a strong slow drift, 2 s stretches still give T8 only;
# one average for the whole recording flags many more channels. So the argument
# really reaches the check.
test_that("find_bad_channels passes amplitude_window on to the amplitude check", {
  eeg <- .make_biosemi_eeg(drift = 2000)
  two  <- suppressMessages(find_bad_channels(eeg, return_details = TRUE))
  inf  <- suppressMessages(find_bad_channels(eeg, amplitude_window = Inf, return_details = TRUE))
  expect_equal(.excessive_channels(two$flags), "T8")
  expect_gt(length(.excessive_channels(inf$flags)), 5)
  expect_error(suppressMessages(find_bad_channels(eeg, amplitude_window = -1)),
               "'amplitude_window' must be a single number above 0")
})

# ----------------------------------------------------------------------------
# Test 13.3: the neighbor check gives the numbers it gave before (median in one go)
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Taking all the middle values at once instead of one moment at
# a time must not change any result. The correlations below were recorded from
# the previous code on the existing problem fixture (robust = TRUE uses the
# median, FALSE the mean). Also: a missing value in a neighbour still gives a
# flagged "correlation = NA".
test_that("the neighbor-correlation check gives the same correlations as before", {
  eeg <- .make_problem_eeg_with_montage()

  med <- suppressWarnings(find_bad_channels(eeg, robust = TRUE, return_details = TRUE))$flags
  med <- med[grepl("correlation", med$reason), ]
  expect_equal(med$channel, c("F3", "Pz", "O1"))
  expect_match(med$value[1], "^correlation = -0.031 ")
  expect_match(med$value[2], "^correlation = 0.074 ")
  expect_match(med$value[3], "^correlation = NA ")

  avg <- suppressWarnings(find_bad_channels(eeg, robust = FALSE, return_details = TRUE))$flags
  avg <- avg[grepl("correlation", avg$reason), ]
  expect_equal(avg$channel, c("F3", "Cz", "P3", "P4", "Pz", "O1", "O2"))
  expect_equal(sub(" .*", "", sub("^correlation = ", "", avg$value)),
               c("-0.031", "0.282", "0.213", "0.276", "0.076", "NA", "0.218"))
})

test_that("a missing value in the neighbors still makes the check report NA", {
  lay <- c("Cz", "Fz", "Pz", "C3", "C4")
  set.seed(5)
  data <- matrix(rnorm(5 * 200, sd = 10), nrow = 5)
  neighbor_map <- data.frame(channel = "Cz", neighbor = c("Fz", "Pz", "C3", "C4"),
                             rank = 1:4, distance = 1:4)
  data[2, 50] <- NA                                           # one neighbor has a missing value
  flags <- eeganalysis:::flag_low_neighbor_correlation(data, lay, neighbor_map, robust = TRUE)
  expect_equal(flags$channel, "Cz")
  expect_match(flags$value, "^correlation = NA ")
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
# - .column_medians(): hand-worked middle values, same numbers as
#   apply(m, 2, median), NA/NaN/Inf/empty handling
# - .centre_by_stretch(): worked examples, stretch edges (leftover joins the
#   last stretch), every stretch averages 0, shape and row independence
# - amplitude check on raw offset data: offsets alone not flagged, drift
#   followed by 2 s stretches, broken channel and spike still flagged, flat
#   and noise checks unchanged, flag wording, amplitude_window default and
#   validation, end to end on 64 BioSemi electrodes, neighbor-check numbers
#   identical to the previous implementation
# ============================================================================
