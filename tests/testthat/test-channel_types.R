# ============================================================================
#                       Test File for channel_types.R
# ============================================================================
#
# Function tested:
#   set_channel_types() - change the type of the channels you name, leave
#                         every other channel exactly as it is
#
# HOW THE FUNCTION WORKS (and therefore how each suite tests it):
# -----------------------------------------------------------------------
# set_channel_types(eeg, types):
#
#  1. Validates 'eeg' is class 'eeg'; 'types' is a named character vector
#     with every name non-missing and non-empty.
#  2. Every name in 'types' must exist in eeg$channels; every value must be
#     one of .valid_channel_types() ("eeg", "eog", "ecg", "emg", "resp",
#     "gsr", "temp", "bio", "misc", "status").
#  3. If eeg$channel_types is NULL (an old or hand-built object), it is
#     worked out first via classify_channels(), the same default new_eeg()
#     would use, so the channels NOT named in 'types' still get a type.
#  4. The named channels' types are overwritten; everything else (data,
#     channels, bads, montage, ...) is copied through unchanged.
#  5. Appends one entry to $preprocessing_history naming every channel
#     changed and its new type.
#
# Author: Christos Dalamarinis
# Date: Oct - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPER
# ============================================================================

make_fixture <- function() {
  set.seed(1)
  new_eeg(
    data = matrix(rnorm(5 * 100), nrow = 5, ncol = 100),
    channels = c("Fz", "Cz", "VEOG", "EXG5", "Status"),
    sampling_rate = 256
  )
}

# ============================================================================
# Suite 1: input validation
# ============================================================================

test_that("errors on non-eeg input", {
  expect_error(set_channel_types(list(), c(VEOG = "eog")),
               "class 'eeg'")
  expect_error(set_channel_types(data.frame(), c(VEOG = "eog")),
               "class 'eeg'")
})

test_that("errors when types is not a named character vector", {
  eeg <- make_fixture()

  expect_error(set_channel_types(eeg, "eog"),
               "named character vector")
  expect_error(set_channel_types(eeg, character(0)),
               "named character vector")
  expect_error(set_channel_types(eeg, c(eog = 1)),
               "named character vector")

  unnamed <- c("eog", "ecg")
  names(unnamed) <- c("VEOG", "")
  expect_error(set_channel_types(eeg, unnamed),
               "named character vector")
})

test_that("errors when a named channel is not in eeg$channels", {
  eeg <- make_fixture()
  expect_error(set_channel_types(eeg, c(BOGUS = "eog")),
               "not found in eeg\\$channels: BOGUS")
})

test_that("errors when a type value is not one of the valid ones", {
  eeg <- make_fixture()
  expect_error(set_channel_types(eeg, c(VEOG = "external")),
               "must be one of")
  expect_error(set_channel_types(eeg, c(VEOG = "EOG")),
               "must be one of")  # case-sensitive: "EOG" is not "eog"
})

# ============================================================================
# Suite 2: core behaviour
# ============================================================================

test_that("changes exactly the named channel and leaves every other channel alone", {
  eeg <- make_fixture()
  out <- set_channel_types(eeg, c(VEOG = "eog"))

  expect_equal(unname(out$channel_types),
               replace(eeg$channel_types, eeg$channels == "VEOG", "eog"))
})

test_that("changes multiple channels in one call", {
  eeg <- make_fixture()
  out <- set_channel_types(eeg, c(VEOG = "eog", EXG5 = "gsr"))

  types <- setNames(out$channel_types, out$channels)
  expect_equal(unname(types["VEOG"]), "eog")
  expect_equal(unname(types["EXG5"]), "gsr")
  expect_equal(unname(types["Fz"]),   "eeg")
  expect_equal(unname(types["Cz"]),   "eeg")
  expect_equal(unname(types["Status"]), "status")
})

test_that("works regardless of the order channels are named in", {
  eeg <- make_fixture()
  a <- set_channel_types(eeg, c(VEOG = "eog", EXG5 = "gsr"))
  b <- set_channel_types(eeg, c(EXG5 = "gsr", VEOG = "eog"))
  expect_equal(a$channel_types, b$channel_types)
})

test_that("every other field is copied through unchanged", {
  eeg <- make_fixture()
  eeg$bads <- "Fz"
  eeg <- suppressWarnings(set_montage(eeg, create_montage(c("Fz", "Cz"))))

  out <- set_channel_types(eeg, c(VEOG = "eog"))

  expect_equal(out$data, eeg$data)
  expect_equal(out$channels, eeg$channels)
  expect_equal(out$bads, eeg$bads)
  expect_equal(out$montage, eeg$montage)
  expect_equal(out$sampling_rate, eeg$sampling_rate)
})

test_that("does not modify the input object", {
  eeg <- make_fixture()
  before <- eeg$channel_types
  invisible(set_channel_types(eeg, c(VEOG = "eog")))
  expect_equal(eeg$channel_types, before)
})

# ============================================================================
# Suite 3: preprocessing_history
# ============================================================================

test_that("appends exactly one entry to preprocessing_history", {
  eeg <- make_fixture()
  out <- set_channel_types(eeg, c(VEOG = "eog"))
  expect_equal(length(out$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
})

test_that("the history entry names the channel and its new type", {
  eeg <- make_fixture()
  out <- set_channel_types(eeg, c(VEOG = "eog"))
  entry <- out$preprocessing_history[[length(out$preprocessing_history)]]
  expect_true(grepl("VEOG", entry, fixed = TRUE))
  expect_true(grepl("eog", entry, fixed = TRUE))
})

test_that("the history entry names every channel when several are set at once", {
  eeg <- make_fixture()
  out <- set_channel_types(eeg, c(VEOG = "eog", EXG5 = "gsr"))
  entry <- out$preprocessing_history[[length(out$preprocessing_history)]]
  expect_true(grepl("VEOG", entry, fixed = TRUE) && grepl("eog", entry, fixed = TRUE))
  expect_true(grepl("EXG5", entry, fixed = TRUE) && grepl("gsr", entry, fixed = TRUE))
})

test_that("stacks on a pre-existing history instead of overwriting it", {
  eeg <- make_fixture()
  eeg$preprocessing_history <- list("Prior step: band-pass filter applied")

  out <- set_channel_types(eeg, c(VEOG = "eog"))

  expect_equal(length(out$preprocessing_history), 2)
  expect_equal(out$preprocessing_history[[1]], "Prior step: band-pass filter applied")
})

# ============================================================================
# Suite 4: objects with no channel_types field yet
# ============================================================================

test_that("an object with no channel_types still gets one for every channel", {
  eeg <- make_fixture()
  eeg$channel_types <- NULL

  out <- set_channel_types(eeg, c(VEOG = "eog"))

  types <- setNames(out$channel_types, out$channels)
  expect_equal(unname(types["VEOG"]),   "eog")     # the one explicitly set
  expect_equal(unname(types["Fz"]),     "eeg")     # worked out, like new_eeg() would
  expect_equal(unname(types["Status"]), "status")  # worked out, like new_eeg() would
})

# ============================================================================
# Suite 5: integration with the rest of the pipeline
# ============================================================================

test_that("a type set this way survives eeg_bandpass(), eeg_notch() and downsample()", {
  eeg <- make_fixture()
  eeg <- set_channel_types(eeg, c(VEOG = "eog", EXG5 = "gsr"))

  bp <- eeg_bandpass(eeg, l_freq = 1, h_freq = 40, verbose = FALSE)
  nt <- eeg_notch(eeg, freqs = 50, verbose = FALSE)
  ds <- downsample(eeg, target_rate = 128, verbose = FALSE)

  expect_equal(bp$channel_types, eeg$channel_types)
  expect_equal(nt$channel_types, eeg$channel_types)
  expect_equal(ds$channel_types, eeg$channel_types)
})

test_that("a channel typed 'eog' this way is then found automatically", {
  eeg <- make_fixture()
  eeg <- set_channel_types(eeg, c(VEOG = "eog"))

  expect_equal(.resolve_reference_channels(eeg, NULL, "EOG"), "VEOG")
})
