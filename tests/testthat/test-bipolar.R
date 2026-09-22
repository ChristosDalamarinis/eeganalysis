# ============================================================================
#                    Test File for bipolar.R
# ============================================================================
#
# Function tested:
#   set_bipolar_reference() - derive a new channel as anode minus cathode
#
# HOW THE FUNCTION WORKS (and therefore how each suite tests it):
# -----------------------------------------------------------------------
# set_bipolar_reference(eeg, anode, cathode, ch_name, drop = TRUE):
#
#  1. Validates 'eeg' is class 'eeg'; 'anode'/'cathode'/'ch_name' are all
#     supplied (no defaults) and are single, non-missing character strings;
#     'anode' != 'cathode'; both exist in eeg$channels.
#  2. Validates 'ch_name' does not collide with a channel that will still
#     exist afterward (eeg$channels minus {anode, cathode} if drop = TRUE,
#     else the full list - so reusing anode's own name works once drop
#     frees it up, but not otherwise).
#  3. new_row = eeg$data[anode, ] - eeg$data[cathode, ]; appended to $data,
#     $channels (as ch_name) and $channel_types (as "external", set
#     directly - not re-derived via classify_channels()).
#  4. If drop = TRUE: removes anode/cathode from data/channels/channel_types
#     and from $bads (cleans up now-stale bad-channel names).
#  5. If anode or cathode was in eeg$bads (checked against the ORIGINAL
#     object, before step 4's cleanup): warns, and ch_name is added to the
#     result's $bads - it was built from a known-noisy electrode.
#  6. Appends one entry to $preprocessing_history.
#  7. Everything else (times, sampling_rate, montage, annotations, events,
#     metadata, reference) passes through unchanged.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

# 5 channels x 4 time points, hand-crafted so anode - cathode is exact and
# easy to check by hand.
#
#          t1   t2   t3   t4
#   Fp1:   10   20   30   40    (regular EEG channel)
#   Cz:     5   15   25   35    (regular EEG channel)
#   EXG1:  100  110  120  130   (anode   - "above the eye")
#   EXG2:   40   45   50   55   (cathode - "below the eye")
#   Status:  1    1    1    1   (BioSemi trigger line)
#
# EXG1 - EXG2 = (60, 65, 70, 75) at every time point - the value every
# "derives the new channel correctly" assertion below checks against.
make_bipolar_fixture <- function() {
  data_mat <- matrix(
    c( 10,  20,  30,  40,     # Fp1
        5,  15,  25,  35,     # Cz
      100, 110, 120, 130,     # EXG1 (anode)
       40,  45,  50,  55,     # EXG2 (cathode)
        1,   1,   1,   1),    # Status
    nrow = 5, ncol = 4, byrow = TRUE
  )
  new_eeg(
    data          = data_mat,
    channels      = c("Fp1", "Cz", "EXG1", "EXG2", "Status"),
    sampling_rate = 256
  )
}

VEOG_EXPECTED <- c(60, 65, 70, 75)

# ============================================================================
# Suite 1: the new channel's values
# ============================================================================

test_that("the new channel equals anode minus cathode exactly", {
  eeg <- make_bipolar_fixture()
  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "VEOG")
  idx <- match("VEOG", out$channels)
  expect_equal(out$data[idx, ], VEOG_EXPECTED)
})

test_that("the new channel is classified as external", {
  eeg <- make_bipolar_fixture()
  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "VEOG")
  idx <- match("VEOG", out$channels)
  expect_equal(out$channel_types[idx], "external")
})

# ============================================================================
# Suite 2: drop behavior
# ============================================================================

test_that("drop = TRUE (default) removes the source channels", {
  eeg <- make_bipolar_fixture()
  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "VEOG")

  expect_false("EXG1" %in% out$channels)
  expect_false("EXG2" %in% out$channels)
  expect_true("VEOG" %in% out$channels)
  expect_equal(length(out$channels), length(eeg$channels) - 1)
  expect_equal(nrow(out$data), length(out$channels))
})

test_that("drop = FALSE keeps the source channels untouched", {
  eeg <- make_bipolar_fixture()
  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "VEOG", drop = FALSE)

  expect_true(all(c("EXG1", "EXG2", "VEOG") %in% out$channels))
  expect_equal(length(out$channels), length(eeg$channels) + 1)
  expect_equal(nrow(out$data), length(out$channels))

  expect_equal(out$data[match("EXG1", out$channels), ],
               eeg$data[match("EXG1", eeg$channels), ])
  expect_equal(out$data[match("EXG2", out$channels), ],
               eeg$data[match("EXG2", eeg$channels), ])
})

test_that("reusing anode's own name works once drop = TRUE frees it up", {
  eeg <- make_bipolar_fixture()
  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "EXG1")

  expect_equal(sum(out$channels == "EXG1"), 1)
  expect_false("EXG2" %in% out$channels)
  idx <- match("EXG1", out$channels)
  expect_equal(out$data[idx, ], VEOG_EXPECTED)   # it's the DERIVED channel
})

# ============================================================================
# Suite 3: validation errors
# ============================================================================

test_that("an unknown anode or cathode channel errors clearly", {
  eeg <- make_bipolar_fixture()
  expect_error(
    set_bipolar_reference(eeg, anode = "NoSuch", cathode = "EXG2", ch_name = "VEOG"),
    "anode.*NoSuch"
  )
  expect_error(
    set_bipolar_reference(eeg, anode = "EXG1", cathode = "NoSuch", ch_name = "VEOG"),
    "cathode.*NoSuch"
  )
})

test_that("anode and cathode must differ", {
  eeg <- make_bipolar_fixture()
  expect_error(
    set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG1", ch_name = "VEOG"),
    "cannot be derived from itself"
  )
})

test_that("anode, cathode and ch_name are required - no defaults", {
  eeg <- make_bipolar_fixture()
  expect_error(set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2"),
               "ch_name")
  expect_error(set_bipolar_reference(eeg, cathode = "EXG2", ch_name = "VEOG"),
               "anode")
})

test_that("ch_name colliding with a surviving channel errors", {
  eeg <- make_bipolar_fixture()
  expect_error(
    set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "Cz"),
    "already names a channel"
  )
  # EXG1 is still present when drop = FALSE, so reusing its name collides
  expect_error(
    set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "EXG1",
                          drop = FALSE),
    "already names a channel"
  )
})

test_that("a non-'eeg' input is rejected with a clear error", {
  expect_error(
    set_bipolar_reference(list(), anode = "A", cathode = "B", ch_name = "VEOG"),
    "class 'eeg'"
  )
})

# ============================================================================
# Suite 4: bad-channel handling
# ============================================================================

test_that("a bad source channel warns and the new channel is marked bad too", {
  eeg <- make_bipolar_fixture()
  eeg$bads <- "EXG1"

  expect_warning(
    out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "VEOG"),
    "EXG1"
  )
  expect_true("VEOG" %in% out$bads)
})

test_that("no warning and no bads entry when neither source is bad", {
  eeg <- make_bipolar_fixture()

  expect_no_warning(
    out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "VEOG")
  )
  expect_false("VEOG" %in% out$bads)
})

test_that("dropping a source channel removes its now-stale name from bads", {
  eeg <- make_bipolar_fixture()
  eeg$bads <- c("EXG2", "Cz")     # EXG2 bad (will be dropped); Cz unrelated

  out <- suppressWarnings(
    set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "VEOG")
  )

  expect_false("EXG2" %in% out$bads)   # stale - EXG2 no longer exists
  expect_true("Cz" %in% out$bads)      # unrelated bad channel untouched
  expect_true("VEOG" %in% out$bads)    # EXG2 was bad -> propagated
})

# ============================================================================
# Suite 5: preprocessing_history
# ============================================================================

test_that("exactly one preprocessing_history entry is appended, naming the derivation", {
  eeg <- make_bipolar_fixture()
  n0  <- length(eeg$preprocessing_history)

  out   <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "VEOG")
  entry <- out$preprocessing_history[[length(out$preprocessing_history)]]

  expect_equal(length(out$preprocessing_history), n0 + 1)
  expect_match(entry, "VEOG")
  expect_match(entry, "EXG1")
  expect_match(entry, "EXG2")
})

# ============================================================================
# Suite 6: reference invariance
# ============================================================================

test_that("the derived channel is the same whether built before or after re-referencing", {
  eeg <- make_bipolar_fixture()

  # Bipolar-derive first, from the original (unreferenced) data.
  before <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                                  ch_name = "VEOG")

  # Re-reference EVERY channel (including EXG1/EXG2) to the common average
  # first, THEN derive the bipolar channel from the shifted EXG1/EXG2.
  reref_first <- eeg_rereference(eeg, ref = "average")
  after <- set_bipolar_reference(reref_first, anode = "EXG1", cathode = "EXG2",
                                 ch_name = "VEOG")

  # Both EXG1 and EXG2 shifted by the same average-reference signal, so the
  # shift cancels out of their difference algebraically.
  expect_equal(before$data[match("VEOG", before$channels), ],
               after$data[match("VEOG", after$channels), ])
})

# ============================================================================
# Suite 7: integration with fit_eog_regression()
# ============================================================================

test_that("fit_eog_regression() finds the derived channel automatically by name", {
  eeg <- eeg_rereference(make_bipolar_fixture(), ref = "average",
                         exclude = c("EXG1", "EXG2", "Status"))
  eeg <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
                               ch_name = "VEOG")

  model <- fit_eog_regression(eeg)
  expect_equal(model$ch_names_artifact, "VEOG")
})

# ============================================================================
# Suite 8: everything else passes through unchanged
# ============================================================================

test_that("non-channel fields pass through unchanged", {
  eeg <- make_bipolar_fixture()
  eeg$metadata <- list(subject = "S01")

  out <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2", ch_name = "VEOG")

  expect_equal(out$times, eeg$times)
  expect_equal(out$sampling_rate, eeg$sampling_rate)
  expect_null(out$montage)
  expect_equal(out$annotations, eeg$annotations)
  expect_equal(out$events, eeg$events)
  expect_equal(out$metadata, eeg$metadata)
  expect_equal(out$reference, eeg$reference)
})
