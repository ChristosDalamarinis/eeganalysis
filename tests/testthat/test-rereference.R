# ============================================================================
#                    Test File for rereference.R
# ============================================================================
#
# This test file provides comprehensive testing for the rereference.R script,
# which re-references EEG data to a new reference scheme.
#
# Function tested:
#   eeg_rereference() - Change reference scheme (average / custom channel(s))
#
# HOW THE FUNCTION WORKS (and therefore how each suite tests it):
# -----------------------------------------------------------------------
# eeg_rereference(eeg, ref, exclude, copy) operates in 8 logical steps:
#
#  1. Copies the eeg object (R's copy-on-modify semantics protect the
#     original when copy = TRUE; the code always does eeg_out <- eeg).
#  2. Extracts eeg$data (channels x time matrix) and eeg$channels (names).
#  3. VALIDATES that nrow(data) == length(channels) → Suite 1.
#  4. Resolves 'exclude' to channel indices (by name via match(), or directly
#     if numeric); removes NAs; builds excl_idx.
#     contrib_idx = all_channels − excl_idx → Suite 7.
#     Stops if contrib_idx is empty → Suite 1.
#  5. Resolves 'ref':
#       - "average"  → ref_idx = contrib_idx
#       - character  → match() to chan_names; stops on NA → Suite 1.
#       - integer    → used directly
#     ref_idx = intersect(ref_idx, contrib_idx); stops if empty → Suite 1.
#  6. Computes ref_signal = colMeans(data[ref_idx,], na.rm=TRUE)
#     → vector of length ncol (one value per time point) → Suite 9.
#  7. Subtracts ref_signal from all channels in apply_idx (= contrib_idx)
#     via sweep(). Excluded channels are untouched → Suites 4–7.
#  8. If eeg$meta != NULL, updates eeg$meta$reference → Suite 8.
#
# NOTE: new_eeg() stores metadata as $metadata (not $meta). The function
#   updates $meta; for freshly constructed objects $meta is NULL, so the
#   update is a no-op unless the caller adds $meta manually.
#
# Author: Christos Dalamarinis
# Date: February 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                        SHARED HELPER UTILITIES
# ============================================================================

# Helper: create a minimal valid eeg object with deterministic data.
# Data layout: rows = channels, columns = time points.
make_eeg <- function(n_channels = 4,
                     n_timepoints = 50,
                     channel_names = NULL,
                     seed = 42) {
  set.seed(seed)
  if (is.null(channel_names)) {
    channel_names <- paste0("Ch", seq_len(n_channels))
  } else {
    # Update n_channels to match the provided channel_names
    n_channels <- length(channel_names)
  }
  data_mat <- matrix(rnorm(n_channels * n_timepoints, mean = 50, sd = 10),
                     nrow = n_channels,
                     ncol = n_timepoints)
  new_eeg(
    data = data_mat,
    channels = channel_names,
    sampling_rate = 256
  )
}

# Helper: create a small eeg object with hand-crafted values so mathematical
# assertions can be stated exactly without relying on rnorm().
#
#   Data (3 channels × 4 time points):
#         t1   t2   t3   t4
#   Ch1:  10   20   30   40
#   Ch2:   0   10   20   30
#   Ch3:  -5    5   15   25
make_known_eeg <- function() {
  data_mat <- matrix(
    c( 10, 20, 30, 40,   # Ch1
       0, 10, 20, 30,   # Ch2
       -5,  5, 15, 25),  # Ch3
    nrow = 3, ncol = 4, byrow = TRUE
  )
  new_eeg(
    data          = data_mat,
    channels      = c("Ch1", "Ch2", "Ch3"),
    sampling_rate = 256
  )
}


# ============================================================================
#                       TEST SUITE 1 – Input Validation
# ============================================================================
# WHAT THESE TESTS DO:
#   Each test deliberately triggers a different error branch inside
#   eeg_rereference(). The tests call the function with invalid arguments
#   and assert that the correct error message is emitted, confirming that
#   the validation guards execute before any computation begins.

test_that("eeg_rereference stops when nrow(data) != length(channels)", {
  # Build a broken eeg-like list: 3 data rows but 5 channel names.
  bad_eeg <- list(
    data          = matrix(1:300, nrow = 3, ncol = 100),
    channels      = paste0("Ch", 1:5),   # 5 names but only 3 rows
    sampling_rate = 256,
    times         = seq(0, 99/256, length.out = 100),
    events        = data.frame(),
    metadata      = list(),
    reference     = "original",
    preprocessing_history = list()
  )
  class(bad_eeg) <- "eeg"
  
  expect_error(
    eeg_rereference(bad_eeg, ref = "average"),
    "Number of columns in signals does not match number of channel names"
  )
})

test_that("eeg_rereference stops when all channels are excluded", {
  # Excluding all channels makes contrib_idx empty; the function must stop
  # before attempting colMeans on a 0-row matrix.
  eeg <- make_eeg(n_channels = 3, channel_names = c("Ch1", "Ch2", "Ch3"))
  
  expect_error(
    eeg_rereference(eeg,
                    ref     = "average",
                    exclude = c("Ch1", "Ch2", "Ch3")),
    "No channels left to compute reference after applying 'exclude'"
  )
})

test_that("eeg_rereference stops when a named ref channel is not found", {
  # match() returns NA for unrecognised names; any(is.na(ref_idx)) triggers
  # the error guard.
  eeg <- make_eeg()
  
  expect_error(
    eeg_rereference(eeg, ref = "NonExistentElectrode"),
    "Some reference channels specified in 'ref' were not found"
  )
})

test_that("eeg_rereference stops when all ref channels are excluded", {
  # After intersect(ref_idx, contrib_idx), if nothing remains, the function
  # must stop with a dedicated message.
  eeg <- make_eeg(n_channels = 4,
                  channel_names = c("Cz", "Pz", "M1", "M2"))
  
  expect_error(
    eeg_rereference(eeg,
                    ref     = c("M1", "M2"),
                    exclude = c("M1", "M2")),
    "Reference channels are all excluded by 'exclude'"
  )
})


# ============================================================================
#               TEST SUITE 2 – Return Type and Object Integrity
# ============================================================================
# WHAT THESE TESTS DO:
#   Call the function with valid input and inspect the returned object.
#   Tests confirm the class, matrix dimensions, channel order, and that
#   non-data fields (sampling_rate, times, events) were not accidentally
#   mutated during the re-referencing computation.

test_that("eeg_rereference returns an object of class 'eeg'", {
  eeg    <- make_eeg()
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_s3_class(result, "eeg")
})

test_that("eeg_rereference preserves data matrix dimensions", {
  eeg    <- make_eeg(n_channels = 5, n_timepoints = 80)
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_equal(dim(result$data), c(5, 80))
})

test_that("eeg_rereference preserves channel names in order", {
  eeg    <- make_eeg(channel_names = c("Cz", "Pz", "Fz", "Oz"))
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_equal(result$channels, c("Cz", "Pz", "Fz", "Oz"))
})

test_that("eeg_rereference preserves sampling_rate, times, and events", {
  eeg    <- make_eeg(n_channels = 3, n_timepoints = 50)
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_equal(result$sampling_rate, eeg$sampling_rate)
  expect_equal(result$times,         eeg$times)
  expect_equal(result$events,        eeg$events)
})


# ============================================================================
#                 TEST SUITE 3 – copy Parameter Behaviour
# ============================================================================
# WHAT THESE TESTS DO:
#   With copy = TRUE (default), the function does eeg_out <- eeg before
#   modifying anything. R's copy-on-modify semantics mean the caller's
#   object is never mutated. These tests snapshot eeg$data before the call
#   and compare it afterwards. They also verify the returned data truly
#   differs (i.e., re-referencing happened).

test_that("copy = TRUE does not modify the original eeg object", {
  eeg           <- make_eeg()
  original_data <- eeg$data   # snapshot
  
  eeg_rereference(eeg, ref = "average", copy = TRUE)
  
  expect_equal(eeg$data, original_data,
               label = "Original eeg$data unchanged after copy = TRUE call")
})

test_that("copy = TRUE returns modified (rereferenced) data", {
  eeg    <- make_eeg()
  result <- eeg_rereference(eeg, ref = "average", copy = TRUE)
  
  # Data must actually change — rereferencing is not a no-op.
  expect_false(identical(result$data, eeg$data))
})

test_that("copy = FALSE still returns a valid modified eeg object", {
  # The current implementation always copies internally (eeg_out <- eeg),
  # so copy = FALSE behaves the same as copy = TRUE in practice.
  eeg    <- make_eeg()
  result <- eeg_rereference(eeg, ref = "average", copy = FALSE)
  
  expect_s3_class(result, "eeg")
  expect_equal(dim(result$data), dim(eeg$data))
})


# ============================================================================
#         TEST SUITE 4 – Average Reference: Mathematical Correctness
# ============================================================================
# WHAT THESE TESTS DO:
#   When ref = "average", the reference signal at time t is
#       ref_signal[t] = mean(data[contrib_idx, t])
#   and every contributing channel i is updated:
#       data_new[i, t] = data[i, t] − ref_signal[t]
#
#   Defining property: colMeans of all contributing channels = 0 after
#   the operation. Tests verify this property, verify exact numeric values
#   on the known fixture, and verify idempotency (applying twice = once).

test_that("average reference makes column means of all channels ~0", {
  eeg    <- make_eeg(n_channels = 6, n_timepoints = 100, seed = 7)
  result <- eeg_rereference(eeg, ref = "average")
  
  col_means_after <- colMeans(result$data)
  expect_true(all(abs(col_means_after) < 1e-10),
              label = "Column means are ~0 at every time point")
})

test_that("average reference values match hand-computed expectations", {
  # Input (3 channels × 4 time points):
  #   Ch1: 10  20  30  40
  #   Ch2:  0  10  20  30
  #   Ch3: -5   5  15  25
  #
  # Column means (reference signal):
  #   t1: 5/3,  t2: 35/3,  t3: 65/3,  t4: 95/3
  #
  # Expected output (each row = original row − ref_signal):
  #   Ch1: 25/3   25/3   25/3   25/3
  #   Ch2: -5/3   -5/3   -5/3   -5/3
  #   Ch3: -20/3  -20/3  -20/3  -20/3
  eeg        <- make_known_eeg()
  result     <- eeg_rereference(eeg, ref = "average")
  ref_signal <- c(5/3, 35/3, 65/3, 95/3)
  
  expect_equal(result$data[1, ], c(10, 20, 30, 40) - ref_signal, tolerance = 1e-10)
  expect_equal(result$data[2, ], c( 0, 10, 20, 30) - ref_signal, tolerance = 1e-10)
  expect_equal(result$data[3, ], c(-5,  5, 15, 25) - ref_signal, tolerance = 1e-10)
})

test_that("average reference is idempotent: applying it twice leaves data unchanged", {
  # After the first pass, column means are 0. A second pass subtracts those
  # zero means, changing nothing.
  eeg   <- make_eeg(n_channels = 5, n_timepoints = 60)
  once  <- eeg_rereference(eeg, ref = "average")
  twice <- eeg_rereference(once, ref = "average")
  
  expect_equal(once$data, twice$data, tolerance = 1e-10)
})


# ============================================================================
#               TEST SUITE 5 – Single Named Channel Reference
# ============================================================================
# WHAT THESE TESTS DO:
#   When ref = "Ch1" (one channel name), ref_signal = data[ch1_idx, ].
#   All contrib channels including Ch1 itself get ref_signal subtracted.
#   Tests verify: (a) the reference channel row becomes all zeros;
#   (b) exact numeric values for other channels on the known fixture;
#   (c) numeric index ref gives the same result as named ref.

test_that("single named ref channel becomes all zeros after rereferencing", {
  eeg    <- make_eeg(channel_names = c("Cz", "Pz", "Fz", "Oz"))
  result <- eeg_rereference(eeg, ref = "Cz")
  
  cz_row <- which(result$channels == "Cz")
  expect_true(all(abs(result$data[cz_row, ]) < 1e-10),
              label = "Reference channel Cz becomes all zeros")
})

test_that("single named ref: other channels equal original minus ref channel", {
  # Input:
  #   Ch1: 10  20  30  40   ← reference
  #   Ch2:  0  10  20  30
  #   Ch3: -5   5  15  25
  #
  # After ref = "Ch1":
  #   Ch1:  0   0   0   0
  #   Ch2: -10 -10 -10 -10
  #   Ch3: -15 -15 -15 -15
  eeg    <- make_known_eeg()
  result <- eeg_rereference(eeg, ref = "Ch1")
  
  expect_equal(result$data[1, ], c(  0,   0,   0,   0), tolerance = 1e-10)
  expect_equal(result$data[2, ], c(-10, -10, -10, -10), tolerance = 1e-10)
  expect_equal(result$data[3, ], c(-15, -15, -15, -15), tolerance = 1e-10)
})

test_that("numeric ref index gives identical result to named ref", {
  # ref = 1L (first channel) must equal ref = "Ch1".
  eeg            <- make_known_eeg()
  result_named   <- eeg_rereference(eeg, ref = "Ch1")
  result_numeric <- eeg_rereference(eeg, ref = 1L)
  
  expect_equal(result_named$data, result_numeric$data, tolerance = 1e-10)
})


# ============================================================================
#             TEST SUITE 6 – Multi-Channel (Linked) Reference
# ============================================================================
# WHAT THESE TESTS DO:
#   When ref is a vector of two or more channels, ref_signal = rowwise mean
#   of those channels at each time point (i.e., colMeans of a 2-row subset).
#   Tests verify: (a) exact numeric values; (b) the "mean of ref channels = 0"
#   property; (c) numeric vector ref equals named vector ref.

test_that("linked two-channel ref produces correct values on known fixture", {
  # Linked ref = mean(Ch1, Ch2):
  #   t1: 5,  t2: 15,  t3: 25,  t4: 35
  #
  # After ref = c("Ch1", "Ch2"):
  #   Ch1:  5   5   5   5
  #   Ch2: -5  -5  -5  -5
  #   Ch3: -10 -10 -10 -10
  eeg    <- make_known_eeg()
  result <- eeg_rereference(eeg, ref = c("Ch1", "Ch2"))
  
  expect_equal(result$data[1, ], c(  5,  5,  5,  5), tolerance = 1e-10)
  expect_equal(result$data[2, ], c( -5, -5, -5, -5), tolerance = 1e-10)
  expect_equal(result$data[3, ], c(-10,-10,-10,-10),  tolerance = 1e-10)
})

test_that("linked ref: mean of the two ref channels is ~0 at every time point", {
  eeg    <- make_eeg(channel_names = c("Cz", "Pz", "M1", "M2", "Fz", "Oz"),
                     seed = 3)
  result <- eeg_rereference(eeg, ref = c("M1", "M2"))
  
  m1_idx      <- which(result$channels == "M1")
  m2_idx      <- which(result$channels == "M2")
  mean_of_ref <- (result$data[m1_idx, ] + result$data[m2_idx, ]) / 2
  
  expect_true(all(abs(mean_of_ref) < 1e-10),
              label = "Mean of linked ref channels is ~0 at every time point")
})

test_that("numeric vector ref gives same result as named vector ref", {
  eeg            <- make_eeg(channel_names = c("Fp1", "Fp2", "M1", "M2", "Cz"),
                             seed = 9)
  result_named   <- eeg_rereference(eeg, ref = c("M1", "M2"))
  result_numeric <- eeg_rereference(eeg, ref = c(3L, 4L))
  
  expect_equal(result_named$data, result_numeric$data, tolerance = 1e-10)
})


# ============================================================================
#                     TEST SUITE 7 – Exclude Parameter
# ============================================================================
# WHAT THESE TESTS DO:
#   Channels in 'exclude' are removed from contrib_idx, which means:
#     (a) They do NOT contribute to the reference signal computation.
#     (b) They are NOT present in apply_idx, so their rows in $data are
#         left untouched.
#   Each test takes a snapshot of the excluded channel's row before the
#   call and compares it to the result, confirming exact preservation.
#   One test on the known fixture verifies that excluding Ch3 causes the
#   reference to be computed from Ch1 and Ch2 only.

test_that("excluded channel by name is not modified in output", {
  eeg          <- make_eeg(channel_names = c("Cz", "Pz", "Fz", "EOG"))
  original_eog <- eeg$data[4, ]  # snapshot before call
  
  result <- eeg_rereference(eeg, ref = "average", exclude = "EOG")
  
  expect_equal(result$data[4, ], original_eog,
               label = "Excluded EOG channel row is unchanged")
})

test_that("excluded channel by numeric index is not modified", {
  eeg          <- make_eeg(channel_names = c("Cz", "Pz", "Fz", "EOG"))
  original_eog <- eeg$data[4, ]
  
  result <- eeg_rereference(eeg, ref = "average", exclude = 4L)
  
  expect_equal(result$data[4, ], original_eog)
})

test_that("excluded channel is not included in average reference computation", {
  # With Ch3 excluded, the reference = mean(Ch1, Ch2) only.
  #
  # Input:
  #   Ch1: 10  20  30  40
  #   Ch2:  0  10  20  30
  #   Ch3: -5   5  15  25  ← excluded
  #
  # Ref = mean(Ch1, Ch2):  t1:5  t2:15  t3:25  t4:35
  #
  # After rereferencing:
  #   Ch1:  5   5   5   5
  #   Ch2: -5  -5  -5  -5
  #   Ch3: -5   5  15  25   (unchanged – excluded)
  eeg    <- make_known_eeg()
  result <- eeg_rereference(eeg, ref = "average", exclude = "Ch3")
  
  expect_equal(result$data[1, ], c( 5,  5,  5,  5), tolerance = 1e-10,
               label = "Ch1 rereferenced to mean(Ch1, Ch2)")
  expect_equal(result$data[2, ], c(-5, -5, -5, -5), tolerance = 1e-10,
               label = "Ch2 rereferenced to mean(Ch1, Ch2)")
  expect_equal(result$data[3, ], c(-5,  5, 15, 25), tolerance = 1e-10,
               label = "Ch3 (excluded) retains original values")
})

test_that("multiple channels can be excluded simultaneously", {
  eeg           <- make_eeg(channel_names = c("Cz", "Pz", "Fz", "EOG1", "EOG2"),
                            seed = 15)
  original_eog1 <- eeg$data[4, ]
  original_eog2 <- eeg$data[5, ]
  
  result <- eeg_rereference(eeg,
                            ref     = "average",
                            exclude = c("EOG1", "EOG2"))
  
  expect_equal(result$data[4, ], original_eog1, label = "EOG1 unchanged")
  expect_equal(result$data[5, ], original_eog2, label = "EOG2 unchanged")
})

test_that("exclude accepts numeric indices for multiple channels", {
  eeg           <- make_eeg(channel_names = c("Cz", "Pz", "EOG1", "EOG2"))
  original_eog1 <- eeg$data[3, ]
  original_eog2 <- eeg$data[4, ]
  
  result <- eeg_rereference(eeg, ref = "average", exclude = c(3L, 4L))
  
  expect_equal(result$data[3, ], original_eog1)
  expect_equal(result$data[4, ], original_eog2)
})

test_that("unrecognised channel name in exclude is silently ignored", {
  # match() returns NA for unknown names; excl_idx[!is.na(excl_idx)] drops them.
  eeg <- make_eeg(n_channels = 3)
  
  expect_no_error(
    eeg_rereference(eeg, ref = "average", exclude = "DOES_NOT_EXIST")
  )
})


# ============================================================================
#                       TEST SUITE 8 – Metadata Update
# ============================================================================
# WHAT THESE TESTS DO:
#   The function runs: if (!is.null(eeg_out$meta)) { eeg_out$meta$reference <- ... }
#   new_eeg() creates $metadata (not $meta), so $meta is NULL on fresh objects.
#   To test the update branch, we manually add $meta to the eeg object.
#   We also verify that the function runs without error on standard new_eeg
#   objects (where $meta is NULL and the update is skipped).

test_that("average ref sets eeg$meta$reference to 'average' when $meta exists", {
  eeg      <- make_eeg()
  eeg$meta <- list()   # inject the field the function will update
  
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_equal(result$meta$reference, "average")
})

test_that("named channel ref sets eeg$meta$reference to the channel name", {
  eeg      <- make_eeg(channel_names = c("Cz", "Pz", "M1", "M2"))
  eeg$meta <- list()
  
  result <- eeg_rereference(eeg, ref = "M1")
  
  expect_equal(result$meta$reference, "M1")
})

test_that("linked ref sets eeg$meta$reference to 'M1+M2'", {
  eeg      <- make_eeg(channel_names = c("Cz", "Pz", "M1", "M2"))
  eeg$meta <- list()
  
  result <- eeg_rereference(eeg, ref = c("M1", "M2"))
  
  expect_equal(result$meta$reference, "M1+M2")
})

test_that("function runs without error when eeg$meta is NULL (standard object)", {
  # The !is.null guard prevents the update from failing on standard objects.
  eeg <- make_eeg()
  
  expect_no_error(eeg_rereference(eeg, ref = "average"))
})


# ============================================================================
#                   TEST SUITE 9 – Edge Cases and Robustness
# ============================================================================

test_that("single-channel eeg referenced to average becomes all zeros", {
  # Degenerate but valid: only one channel. Average ref = itself → all zeros.
  single_eeg <- new_eeg(
    data          = matrix(c(5, 10, 15, 20), nrow = 1),
    channels      = "Cz",
    sampling_rate = 256
  )
  result <- eeg_rereference(single_eeg, ref = "average")
  
  expect_equal(result$data[1, ], c(0, 0, 0, 0), tolerance = 1e-10)
})

test_that("channel differences are invariant to rereferencing", {
  # Re-referencing subtracts the same ref_signal[t] from all channels at t.
  # Therefore (Ch1[t] − ref[t]) − (Ch2[t] − ref[t]) = Ch1[t] − Ch2[t].
  # The pairwise difference between any two channels is preserved.
  eeg    <- make_eeg(n_channels = 4, n_timepoints = 100, seed = 101)
  result <- eeg_rereference(eeg, ref = "average")
  
  diff_before <- eeg$data[1, ] - eeg$data[2, ]
  diff_after  <- result$data[1, ] - result$data[2, ]
  
  expect_equal(diff_before, diff_after, tolerance = 1e-10,
               label = "Ch1 − Ch2 difference is invariant to rereferencing")
})

test_that("colMeans uses na.rm = TRUE: NA in one channel does not propagate", {
  # Construct a 3-channel × 4-timepoint matrix with NA at (Ch2, t2).
  # The reference at t2 should come from Ch1 and Ch3 only.
  #
  # t2 reference = mean(20, NA, 5, na.rm=TRUE) = 12.5
  # Ch1[t2] after: 20 − 12.5 = 7.5
  # Ch3[t2] after:  5 − 12.5 = -7.5
  data_with_na <- matrix(c(10,  0, -5,
                           20, NA,  5,
                           30, 20, 15,
                           40, 30, 25),
                         nrow = 3, ncol = 4)
  eeg_na <- new_eeg(
    data          = data_with_na,
    channels      = c("Ch1", "Ch2", "Ch3"),
    sampling_rate = 256
  )
  result <- eeg_rereference(eeg_na, ref = "average")
  
  expect_equal(result$data[1, 2],  7.5, tolerance = 1e-10)
  expect_equal(result$data[3, 2], -7.5, tolerance = 1e-10)
})

test_that("rereferencing produces no NaN or Inf values on well-formed data", {
  eeg    <- make_eeg(n_channels = 6, n_timepoints = 200, seed = 77)
  result <- eeg_rereference(eeg, ref = "average")
  
  expect_false(any(is.nan(result$data)),      label = "No NaN values")
  expect_false(any(is.infinite(result$data)), label = "No Inf values")
})

test_that("64-channel average reference satisfies the column-mean = 0 property", {
  # Verify the defining property at realistic BioSemi 64-channel scale.
  set.seed(55)
  n_chan  <- 64
  n_time  <- 512
  eeg_64  <- new_eeg(
    data          = matrix(rnorm(n_chan * n_time), nrow = n_chan),
    channels      = paste0("Ch", seq_len(n_chan)),
    sampling_rate = 2048
  )
  result <- eeg_rereference(eeg_64, ref = "average")
  
  col_means <- colMeans(result$data)
  expect_true(all(abs(col_means) < 1e-10),
              label = "Column means ~0 for 64-channel average reference")
})

test_that("single time-point matrix is handled without error", {
  # Edge case: ncol = 1. ref_signal is a scalar; should still work.
  eeg_1tp <- new_eeg(
    data          = matrix(c(10, 0, -5), nrow = 3, ncol = 1),
    channels      = c("Ch1", "Ch2", "Ch3"),
    sampling_rate = 256
  )
  result <- eeg_rereference(eeg_1tp, ref = "average")
  
  # ref = mean(10, 0, −5) = 5/3
  expect_equal(result$data[1, 1], 10 - 5/3, tolerance = 1e-10)
  expect_equal(result$data[2, 1],  0 - 5/3, tolerance = 1e-10)
  expect_equal(result$data[3, 1], -5 - 5/3, tolerance = 1e-10)
})