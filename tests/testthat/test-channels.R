# ============================================================================
#                    Test File for channels.R
# ============================================================================
#
# Functions tested:
#   drop_channels()    - remove channels by name
#   pick_channels()    - keep channels by name (in the order given)
#   pick_types()       - keep channels by type
#   rename_channels()  - change channel names
#   add_channels()     - append the channels of another eeg object
#   combine_channels() - merge groups of channels into one channel each
#   eeg_crop()         - cut time off the start and the end
#   eeg_concatenate()  - join recordings one after the other
#   (and the shared helper .select_channels() behind the first four)
#
# THE TEST RECORDING
#   The numbers in the expectations come from the reference implementation
#   (run once on this exact recording and written in here as literals, so
#   these tests need no other software). The recording: 7 channels
#   (Fz Cz Pz O1 = eeg, VEOG = eog, EXG7 = ecg, Status = status), 100 Hz,
#   1000 samples (10 s). The value of channel i at sample s (s = 1..1000) is
#   i * 100 + (s - 1), so every number is easy to check by eye: Fz starts at
#   100, Cz at 200, ... Pz and O1 are bad. Three annotations: a (2 s, 1 s long,
#   all channels), b (7.5 s, 1 s, aimed at Cz) and c (3.5 s, 2 s, all
#   channels). Three events at samples 101, 301 and 801.
#
# Test suites:
#   1. drop_channels()
#   2. pick_channels()
#   3. pick_types()
#   4. rename_channels()
#   5. add_channels()
#   6. combine_channels()
#   7. eeg_crop()
#   8. eeg_concatenate()
#   9. Things that must always hold (round trips, both object classes)
#  10. Chains with the EOG tools of gaps 3.3 / 3.4
#  11. Wrong input
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

toy_names <- c("Fz", "Cz", "Pz", "O1", "VEOG", "EXG7", "Status")
toy_types <- c("eeg", "eeg", "eeg", "eeg", "eog", "ecg", "status")

toy_annotations <- function() {
  data.frame(onset = c(2.0, 7.5, 3.5), duration = c(1, 1, 2),
             description = c("BAD_a", "BAD_b", "BAD_c"),
             channel = c(NA, "Cz", NA), stringsAsFactors = FALSE)
}

toy_eeg <- function(annotations = toy_annotations()) {
  d <- t(vapply(1:7, function(i) i * 100 + (0:999), numeric(1000)))
  eeganalysis::new_eeg(
    data = d, channels = toy_names, sampling_rate = 100,
    channel_types = toy_types, bads = c("Pz", "O1"),
    annotations = annotations,
    events = data.frame(onset = c(101L, 301L, 801L), onset_time = c(1, 3, 8),
                        type = c(1L, 2L, 3L),
                        description = c("t1", "t2", "t3"),
                        stringsAsFactors = FALSE),
    montage = create_montage(c("Fz", "Cz", "Pz", "O1")),
    metadata = list(subject = "S01"))
}

# random data, a reference that the EOG regression accepts
toy_noise <- function() {
  set.seed(1)
  eeganalysis::new_eeg(
    data = matrix(rnorm(7 * 2000), 7, 2000), channels = toy_names,
    sampling_rate = 100, channel_types = toy_types, bads = c("Pz", "O1"),
    reference = "Common Average")
}

# 7 channels x 50 samples x 4 epochs; value = i * 100 + sample + 1000 * epoch
toy_epochs <- function() {
  arr <- array(0, dim = c(7, 50, 4), dimnames = list(toy_names, NULL, NULL))
  for (i in 1:7) arr[i, , ] <- i * 100 + outer(0:49, 1:4, function(s, e) s + 1000 * e)
  structure(list(
    data = arr, channels = toy_names, times = (-10:39) / 100,
    events = data.frame(onset = 1:4, type = 1:4), sampling_rate = 100,
    tmin = -0.1, tmax = 0.39, baseline = NULL, baseline_method = "none",
    n_epochs = 4L, rejected = rep(FALSE, 4), rejection_log = data.frame(),
    metadata = list(), channel_types = toy_types, bads = c("Pz", "O1"),
    reference = "original", preprocessing_history = list(),
    montage = create_montage(c("Fz", "Cz", "Pz", "O1"))),
    class = "eeg_epochs")
}

# Every field that mentions channels must agree with the others.
check_consistent <- function(e) {
  n <- length(e$channels)
  if (inherits(e, "eeg_epochs")) {
    expect_equal(dim(e$data)[1], n)
    expect_equal(dimnames(e$data)[[1]], e$channels)
  } else {
    expect_equal(nrow(e$data), n)
  }
  expect_length(e$channel_types, n)
  expect_true(all(e$bads %in% e$channels))
  expect_false(anyDuplicated(e$channels) > 0)
  if (!is.null(e$montage)) {
    expect_true(all(e$montage$channels %in% e$channels))
    expect_equal(nrow(e$montage$positions), length(e$montage$channels))
    expect_equal(e$montage$positions$channel, e$montage$channels)
  }
  if (!is.null(e$annotations) && nrow(e$annotations) > 0) {
    ch <- e$annotations$channel[!is.na(e$annotations$channel)]
    expect_true(all(ch %in% e$channels))
  }
}

# The fields of two objects that must be equal, leaving out the history.
same_but_history <- function(a, b) {
  a$preprocessing_history <- NULL
  b$preprocessing_history <- NULL
  expect_equal(a, b)
}


# ============================================================================
# SUITE 1: drop_channels()
# ============================================================================

test_that("drop_channels(Pz) keeps every field in step (reference numbers)", {
  x <- toy_eeg()
  y <- drop_channels(x, "Pz")
  expect_equal(y$channels, c("Fz", "Cz", "O1", "VEOG", "EXG7", "Status"))
  expect_equal(y$channel_types, c("eeg", "eeg", "eeg", "eog", "ecg", "status"))
  expect_equal(y$bads, "O1")                     # bads were Pz and O1
  expect_equal(dim(y$data), c(6L, 1000L))
  expect_equal(y$data[3, 1], 400)                # O1 is the third row now
  expect_equal(y$data[1, 1000], 100 + 999)
  expect_equal(y$montage$channels, c("Fz", "Cz", "O1"))
  expect_equal(nrow(y$annotations), 3L)          # b is aimed at Cz: still there
  expect_identical(y$times, x$times)
  expect_identical(y$events, x$events)
  expect_identical(y$sampling_rate, x$sampling_rate)
  expect_identical(y$reference, x$reference)
  expect_identical(y$metadata, x$metadata)
  expect_length(y$preprocessing_history, 1L)
  expect_match(y$preprocessing_history[[1]], "^Channels dropped: Pz \\(1 of 7\\)$")
  check_consistent(y)
  expect_identical(x, toy_eeg())                 # the input is not changed
})

test_that("drop_channels() deletes the annotation rows aimed at a dropped channel", {
  y <- drop_channels(toy_eeg(), "Cz")
  expect_equal(y$annotations$description, c("BAD_a", "BAD_c"))   # b was aimed at Cz
  expect_equal(nrow(drop_channels(toy_eeg(), c("Fz", "VEOG"))$annotations), 3L)
  check_consistent(y)
})

test_that("drop_channels() can drop several channels and remove all bad ones", {
  y <- drop_channels(toy_eeg(), c("Pz", "O1", "VEOG"))
  expect_equal(y$channels, c("Fz", "Cz", "EXG7", "Status"))
  expect_equal(y$bads, character(0))
  expect_equal(y$montage$channels, c("Fz", "Cz"))
  expect_match(y$preprocessing_history[[1]], "^Channels dropped: Pz, O1, VEOG \\(3 of 7\\)$")
  check_consistent(y)
})

test_that("a montage with no channel left becomes NULL", {
  y <- drop_channels(toy_eeg(), c("Fz", "Cz", "Pz", "O1"))
  expect_null(y$montage)
  expect_true("montage" %in% names(y))           # the field itself stays
})

test_that("drop_channels() stops, with the reference messages, for unknown names", {
  x <- toy_eeg()
  expect_error(drop_channels(x, "Nope"), "Channel(s) Nope not found, nothing dropped.",
               fixed = TRUE)
  expect_error(drop_channels(x, c("Fz", "Nope")), "Nope not found", fixed = TRUE)
  z <- drop_channels(x, c("Nope", "Fz"), on_missing = "ignore")
  expect_equal(z$channels, setdiff(toy_names, "Fz"))
  expect_identical(drop_channels(x, "Nope", on_missing = "ignore"), x)   # nothing to do
  expect_error(drop_channels(x, toy_names), "All channels would be dropped.", fixed = TRUE)
  expect_error(drop_channels(x, c("Fz", "Fz")), "repeats: Fz")
  expect_error(drop_channels(x, character(0)), "non-empty character vector")
  expect_error(drop_channels(x, NA_character_), "non-empty character vector")
  expect_error(drop_channels(x, 1), "non-empty character vector")
  expect_error(drop_channels(x, "Fz", on_missing = "maybe"), "should be one of")
})

test_that("drop_channels() works on epochs (data array and its dimnames)", {
  e <- toy_epochs()
  y <- drop_channels(e, c("Pz", "VEOG"))
  expect_s3_class(y, "eeg_epochs")
  expect_equal(dim(y$data), c(5L, 50L, 4L))
  expect_equal(dimnames(y$data)[[1]], c("Fz", "Cz", "O1", "EXG7", "Status"))
  expect_equal(unname(y$data["O1", 3, 2]), 400 + 2 + 2000)       # O1 keeps its own numbers
  expect_equal(y$bads, "O1")
  expect_equal(y$channel_types, c("eeg", "eeg", "eeg", "ecg", "status"))
  expect_equal(y$n_epochs, 4L)
  check_consistent(y)
})

test_that("an object without channel_types stays without them", {
  x <- toy_eeg(); x$channel_types <- NULL
  y <- drop_channels(x, "Pz")
  expect_null(y$channel_types)
  expect_equal(y$channels, setdiff(toy_names, "Pz"))
})


# ============================================================================
# SUITE 2: pick_channels()
# ============================================================================

test_that("pick_channels(c('Cz','Fz')) follows the order asked for (reference numbers)", {
  x <- toy_eeg()
  y <- pick_channels(x, c("Cz", "Fz"))
  expect_equal(y$channels, c("Cz", "Fz"))
  expect_equal(y$data[1, 1], 200)                # the first row is Cz now
  expect_equal(y$data[2, 1], 100)
  expect_equal(y$channel_types, c("eeg", "eeg"))
  expect_equal(y$bads, character(0))
  expect_equal(y$montage$channels, c("Fz", "Cz"))   # the montage keeps its own order
  expect_match(y$preprocessing_history[[1]], "^Channels picked: Cz, Fz \\(2 of 7\\)$")
  check_consistent(y)
  expect_identical(x, toy_eeg())
})

test_that("pick_channels(ordered = FALSE) keeps the original order", {
  y <- pick_channels(toy_eeg(), c("Cz", "Fz"), ordered = FALSE)
  expect_equal(y$channels, c("Fz", "Cz"))
  expect_equal(y$data[1, 1], 100)
})

test_that("pick_channels() keeps the bad channels among the picked ones", {
  y <- pick_channels(toy_eeg(), c("Cz", "Pz"), ordered = FALSE)
  expect_equal(y$bads, "Pz")
  check_consistent(y)
})

test_that("pick_channels() removes the annotation rows aimed at a channel that is gone", {
  expect_equal(nrow(pick_channels(toy_eeg(), c("Fz", "O1"))$annotations), 2L)  # b was aimed at Cz
  expect_equal(nrow(pick_channels(toy_eeg(), c("Fz", "Cz"))$annotations), 3L)
})

test_that("pick_channels() is strict about names, and on_missing = 'ignore' skips them", {
  x <- toy_eeg()
  expect_error(pick_channels(x, c("Cz", "Nope")),
               "Missing channel(s) required by 'channels': Nope", fixed = TRUE)
  expect_error(pick_channels(x, c("Cz", "Nope"), ordered = FALSE),
               "Missing channel(s) required by 'channels': Nope", fixed = TRUE)
  expect_equal(pick_channels(x, c("Nope", "Cz"), on_missing = "ignore")$channels, "Cz")
  expect_error(pick_channels(x, "Nope", on_missing = "ignore"),
               "No channels match the selection.", fixed = TRUE)
  expect_error(pick_channels(x, c("Fz", "Fz")), "repeats: Fz")
  expect_error(pick_channels(x, "Fz", ordered = NA), "'ordered' must be TRUE or FALSE")
  expect_error(pick_channels(x, "Fz", ordered = "yes"), "'ordered' must be TRUE or FALSE")
})

test_that("pick_channels() works on epochs", {
  y <- pick_channels(toy_epochs(), c("O1", "Fz"))
  expect_equal(dim(y$data), c(2L, 50L, 4L))
  expect_equal(dimnames(y$data)[[1]], c("O1", "Fz"))
  expect_equal(unname(y$data["O1", 1, 1]), 400 + 1000)
  expect_equal(y$bads, "O1")
  check_consistent(y)
})


# ============================================================================
# SUITE 3: pick_types()
# ============================================================================

test_that("pick_types(eeg) leaves the bad channels out and empties bads (reference numbers)", {
  x <- toy_eeg()
  y <- pick_types(x, "eeg")
  expect_equal(y$channels, c("Fz", "Cz"))        # Pz and O1 were bad
  expect_equal(y$bads, character(0))
  expect_equal(y$channel_types, c("eeg", "eeg"))
  expect_equal(y$montage$channels, c("Fz", "Cz"))
  expect_match(y$preprocessing_history[[1]],
               "^Channels picked by type \\(eeg\\): 2 of 7 kept; left out \\(excluded\\): Pz, O1$")
  check_consistent(y)
  expect_identical(x, toy_eeg())
})

test_that("pick_types(exclude = character(0)) keeps the bad channels and bads", {
  y <- pick_types(toy_eeg(), "eeg", exclude = character(0))
  expect_equal(y$channels, c("Fz", "Cz", "Pz", "O1"))
  expect_equal(y$bads, c("Pz", "O1"))
  expect_equal(pick_types(toy_eeg(), "eeg", exclude = NULL)$channels, y$channels)
  expect_match(y$preprocessing_history[[1]], "^Channels picked by type \\(eeg\\): 4 of 7 kept$")
})

test_that("pick_types() with several kinds keeps the original order", {
  y <- pick_types(toy_eeg(), c("eeg", "eog", "ecg", "status"))
  expect_equal(y$channels, c("Fz", "Cz", "VEOG", "EXG7", "Status"))
  expect_equal(y$channel_types, c("eeg", "eeg", "eog", "ecg", "status"))
  expect_equal(pick_types(toy_eeg(), c("status", "eog"))$channels, c("VEOG", "Status"))
})

test_that("pick_types() 'exclude' takes names as well as the word bads", {
  expect_equal(pick_types(toy_eeg(), "eeg", exclude = "Fz")$channels, c("Cz", "Pz", "O1"))
  expect_equal(pick_types(toy_eeg(), "eeg", exclude = "Fz")$bads, c("Pz", "O1"))
  expect_equal(pick_types(toy_eeg(), "eeg", exclude = c("bads", "Fz"))$channels, "Cz")
  expect_equal(pick_types(toy_eeg(), "eog")$channels, "VEOG")     # a bad eeg channel does not matter
})

test_that("pick_types() uses the name-based default when channel_types is missing", {
  x <- toy_eeg(); x$channel_types <- NULL
  y <- pick_types(x, "status")
  expect_equal(y$channels, "Status")
  expect_null(y$channel_types)
})

test_that("pick_types() stops when nothing matches, and for wrong arguments", {
  x <- toy_eeg()
  expect_error(pick_types(x, "gsr"), "No channels match the selection.", fixed = TRUE)
  expect_error(pick_types(x, "eeg", exclude = toy_names), "No channels match", fixed = TRUE)
  expect_error(pick_types(x, "banana"), "'types' must be from: eeg, eog, ecg")
  expect_error(pick_types(x), "'types' must be a non-empty character vector")
  expect_error(pick_types(x, character(0)), "'types' must be a non-empty character vector")
  expect_error(pick_types(x, 1), "'types' must be a non-empty character vector")
  expect_error(pick_types(x, "eeg", exclude = 5), "'exclude' must be")
})

test_that("pick_types() works on epochs", {
  y <- pick_types(toy_epochs(), "eeg", exclude = character(0))
  expect_equal(dim(y$data), c(4L, 50L, 4L))
  expect_equal(dimnames(y$data)[[1]], c("Fz", "Cz", "Pz", "O1"))
  check_consistent(y)
})


# ============================================================================
# SUITE 4: rename_channels()
# ============================================================================

test_that("rename_channels() carries the new names everywhere (reference numbers)", {
  x <- toy_eeg()
  y <- rename_channels(x, c(Fz = "Fz_new", Pz = "Pz_new"))
  expect_equal(y$channels, c("Fz_new", "Cz", "Pz_new", "O1", "VEOG", "EXG7", "Status"))
  expect_equal(y$bads, c("Pz_new", "O1"))
  expect_equal(y$montage$channels, c("Fz_new", "Cz", "Pz_new", "O1"))
  expect_equal(y$montage$positions$channel, y$montage$channels)
  expect_equal(y$montage$positions[, c("x", "y", "z")], x$montage$positions[, c("x", "y", "z")])
  expect_equal(y$channel_types, x$channel_types)
  expect_equal(y$data, x$data)                    # the signals do not move
  expect_match(y$preprocessing_history[[1]], "^Channels renamed: Fz -> Fz_new, Pz -> Pz_new$")
  check_consistent(y)
  expect_identical(x, toy_eeg())
})

test_that("rename_channels() renames the channel in the annotations aimed at it", {
  y <- rename_channels(toy_eeg(), c(Cz = "Cz_new"))
  expect_equal(y$annotations$channel, c(NA, "Cz_new", NA))
  check_consistent(y)
})

test_that("rename_channels() accepts a function", {
  y <- rename_channels(toy_eeg(), toupper)
  expect_equal(y$channels, c("FZ", "CZ", "PZ", "O1", "VEOG", "EXG7", "STATUS"))
  expect_equal(y$bads, c("PZ", "O1"))
  expect_match(y$preprocessing_history[[1]], "^Channels renamed with a function: ")
  check_consistent(y)
  expect_error(rename_channels(toy_eeg(), function(x) x[-1]), "must return one non-empty")
  expect_error(rename_channels(toy_eeg(), function(x) rep(NA_character_, length(x))),
               "must return one non-empty")
})

test_that("two names can be swapped: each channel keeps its own signal", {
  x <- toy_eeg()
  y <- rename_channels(x, c(Fz = "Cz", Cz = "Fz"))
  expect_equal(y$channels, c("Cz", "Fz", "Pz", "O1", "VEOG", "EXG7", "Status"))
  expect_equal(y$data, x$data)
  expect_equal(y$data[2, 1], 200)                 # the channel now called Fz holds the old Cz signal
  expect_equal(y$annotations$channel, c(NA, "Fz", NA))   # b was aimed at the signal of old Cz
  check_consistent(y)
})

test_that("rename_channels() refuses duplicates and unknown names (reference messages)", {
  x <- toy_eeg()
  expect_error(rename_channels(x, c(Fz = "Cz")),
               "New channel names are not unique, renaming failed.", fixed = TRUE)
  expect_error(rename_channels(x, c(Fz = "A", Cz = "A")), "not unique")
  expect_error(rename_channels(x, c(Nope = "X")), "Invalid channel name(s) Nope", fixed = TRUE)
  expect_error(rename_channels(x, c("A", "B")), "must be a named character vector")
  expect_error(rename_channels(x, c(Fz = NA_character_)), "must be a named character vector")
  expect_error(rename_channels(x, c(Fz = "")), "must be a named character vector")
  expect_error(rename_channels(x, c(Fz = "A", Fz = "B")), "more than once")
  expect_error(rename_channels(x, character(0)), "must be a named character vector")
  expect_error(rename_channels(x, 5), "must be a named character vector")
  expect_identical(rename_channels(x, c(Fz = "Fz")), x)          # nothing changes: nothing to do
})

test_that("rename_channels() works on epochs (data dimnames included)", {
  y <- rename_channels(toy_epochs(), c(Cz = "CZ"))
  expect_equal(dimnames(y$data)[[1]][2], "CZ")
  expect_equal(y$channels[2], "CZ")
  expect_equal(unname(y$data["CZ", 1, 1]), 200 + 1000)
  check_consistent(y)
})


# ============================================================================
# SUITE 5: add_channels()
# ============================================================================

heart <- function(n = 1000, name = "NEW", ...) {
  eeganalysis::new_eeg(data = matrix(5, 1, n), channels = name,
                       sampling_rate = 100, ...)
}

test_that("add_channels() appends the channel at the end (reference numbers)", {
  x <- toy_eeg()
  y <- add_channels(x, heart())
  expect_equal(y$channels, c(toy_names, "NEW"))
  expect_equal(y$channel_types, c(toy_types, "eeg"))
  expect_equal(y$bads, c("Pz", "O1"))
  expect_equal(dim(y$data), c(8L, 1000L))
  expect_equal(y$data[8, ], rep(5, 1000))
  expect_equal(y$data[1:7, ], x$data)             # the old rows are unchanged
  expect_null(rownames(y$data))                   # no "new_row" style labels
  expect_identical(y$montage, x$montage)          # the montage is left as it is
  expect_identical(y$annotations, x$annotations)
  expect_identical(y$events, x$events)
  expect_match(y$preprocessing_history[[1]], "^Channels added: NEW \\(1\\)$")
  check_consistent(y)
  expect_identical(x, toy_eeg())
})

test_that("add_channels() takes the types and the bad channels of the added object", {
  new <- eeganalysis::new_eeg(data = matrix(5, 2, 1000), channels = c("ECG", "RESP"),
                              sampling_rate = 100, channel_types = c("ecg", "resp"),
                              bads = "RESP")
  y <- add_channels(toy_eeg(), new)
  expect_equal(tail(y$channel_types, 2), c("ecg", "resp"))
  expect_equal(y$bads, c("Pz", "O1", "RESP"))
  expect_match(y$preprocessing_history[[1]], "^Channels added: ECG, RESP \\(2\\)$")
})

test_that("add_channels() says when events or annotations of the added object are not copied", {
  new <- heart(events = data.frame(onset = 5L, onset_time = 0.04, type = 1L,
                                   description = "x", stringsAsFactors = FALSE),
               annotations = data.frame(onset = 1, duration = 1, description = "BAD_x",
                                        channel = NA_character_, stringsAsFactors = FALSE))
  y <- add_channels(toy_eeg(), new)
  expect_match(y$preprocessing_history[[1]], "the events and annotations of the added object were not copied")
  expect_equal(nrow(y$events), 3L)
})

test_that("add_channels() keeps row names when the data have them", {
  x <- toy_eeg(); rownames(x$data) <- toy_names
  y <- add_channels(x, heart())
  expect_equal(rownames(y$data), c(toy_names, "NEW"))
})

test_that("add_channels() stops for different length, name clash and sampling rate", {
  x <- toy_eeg()
  expect_error(add_channels(x, heart(n = 500)),
               "All data dimensions except channels must match, got [1000] != [500].",
               fixed = TRUE)
  expect_error(add_channels(x, heart(name = "Fz")), "channel(s) already in eeg$channels: Fz",
               fixed = TRUE)
  other <- eeganalysis::new_eeg(data = matrix(5, 1, 1000), channels = "NEW", sampling_rate = 200)
  expect_error(add_channels(x, other), "The sampling rates differ (100 Hz in 'eeg', 200 Hz in 'new')",
               fixed = TRUE)
  twice <- eeganalysis::new_eeg(data = matrix(5, 2, 1000), channels = c("A", "A"), sampling_rate = 100)
  expect_error(add_channels(x, twice), "'new' repeats channel name(s): A", fixed = TRUE)
  expect_error(add_channels(x, "NEW"), "'new' must be an object of class 'eeg'")
})


# ============================================================================
# SUITE 6: combine_channels()
# ============================================================================

test_that("combine_channels(mean) of Fz, Cz, Pz (reference numbers)", {
  x <- toy_eeg()
  y <- combine_channels(x, list(ROI = c("Fz", "Cz", "Pz")))
  expect_equal(y$channels, "ROI")
  expect_equal(y$channel_types, "eeg")
  expect_equal(dim(y$data), c(1L, 1000L))
  expect_equal(y$data[1, 1], 200)                 # mean of 100, 200, 300: Pz (bad) is included
  expect_equal(y$data[1, ], 200 + (0:999))        # and so for every sample
  expect_equal(y$bads, character(0))
  expect_null(y$montage)
  expect_match(y$preprocessing_history[[1]],
               "^Channels combined \\(mean\\): ROI = Fz \\+ Cz \\+ Pz; the original channels were not kept$")
  check_consistent(y)
  expect_identical(x, toy_eeg())
})

test_that("combine_channels(drop_bad = TRUE) leaves the bad channels out (150, not 200)", {
  y <- combine_channels(toy_eeg(), list(ROI = c("Fz", "Cz", "Pz")), drop_bad = TRUE)
  expect_equal(y$data[1, 1], 150)                 # mean of 100 and 200
  expect_match(y$preprocessing_history[[1]], "ROI = Fz \\+ Cz;")
})

test_that("combine_channels() median, population std and a function (reference numbers)", {
  x <- toy_eeg()
  expect_equal(combine_channels(x, list(ROI = 1:3), method = "median")$data[1, 1], 200)
  expect_equal(combine_channels(x, list(ROI = 1:3), method = "std")$data[1, 1],
               81.64965809277261, tolerance = 1e-12)     # divides by 3, not by 2 (that would give 100)
  expect_equal(combine_channels(x, list(ROI = 1:3), method = "std")$data[1, 500],
               81.64965809277261, tolerance = 1e-12)     # the spread does not change along the ramp
  mx <- combine_channels(x, list(ROI = 1:3), method = function(b) apply(b, 2, max))
  expect_equal(mx$data[1, 1], 300)
  expect_match(mx$preprocessing_history[[1]], "^Channels combined \\(custom function\\)")
  med4 <- combine_channels(x, list(ROI = 1:4), method = "median")
  expect_equal(med4$data[1, 1], 250)              # even number of channels: mean of the middle two
})

test_that("the median is not the mean: a group with an outlier", {
  d <- rbind(rep(1, 10), rep(2, 10), rep(10, 10), rep(7, 10))
  x <- eeganalysis::new_eeg(data = d, channels = c("A", "B", "C", "D"), sampling_rate = 10)
  expect_equal(combine_channels(x, list(G = 1:3), method = "median")$data[1, ], rep(2, 10))
  expect_equal(combine_channels(x, list(G = 1:3), method = "mean")$data[1, ], rep(13 / 3, 10))
  expect_equal(combine_channels(x, list(G = 1:4), method = "median")$data[1, 1], 4.5)   # (2 + 7) / 2
  # the population standard deviation of 1, 2, 10: sqrt(((1 - 13/3)^2 + (2 - 13/3)^2 + (10 - 13/3)^2) / 3)
  expect_equal(combine_channels(x, list(G = 1:3), method = "std")$data[1, 1],
               sqrt(((1 - 13 / 3)^2 + (2 - 13 / 3)^2 + (10 - 13 / 3)^2) / 3))
})

test_that("combine_channels() with several groups, channel numbers or names (reference numbers)", {
  x <- toy_eeg()
  y <- combine_channels(x, list(A = 1:2, B = 3:4))
  expect_equal(y$channels, c("A", "B"))
  expect_equal(y$data[, 1], c(150, 350))
  z <- combine_channels(x, list(A = c("Fz", "Cz"), B = c("Pz", "O1")))
  expect_equal(z$data, y$data)
  expect_equal(combine_channels(x, list(V = "VEOG"))$channel_types, "eog")
  expect_equal(combine_channels(x, list(V = "VEOG"))$data[1, ], x$data[5, ])   # a group of one
})

test_that("combine_channels(keep_status = TRUE) puts the status channel first (reference order)", {
  y <- combine_channels(toy_eeg(), list(A = 1:2), keep_status = TRUE)
  expect_equal(y$channels, c("Status", "A"))
  expect_equal(y$channel_types, c("status", "eeg"))
  expect_equal(y$data[1, 1], 700)
  expect_equal(y$data[2, 1], 150)
  expect_match(y$preprocessing_history[[1]], "\\(status kept\\)$")
  check_consistent(y)
})

test_that("combine_channels() deletes annotation rows aimed at channels that are gone", {
  y <- combine_channels(toy_eeg(), list(A = 1:2))
  expect_equal(y$annotations$description, c("BAD_a", "BAD_c"))
})

test_that("combine_channels() stops with the reference messages", {
  x <- toy_eeg()
  expect_error(combine_channels(x, list(X = c(1, 5))),
               "Cannot combine sensors of different types; \"X\" contains types eeg, eog.",
               fixed = TRUE)
  expect_error(combine_channels(x, list(X = c(1, 8))), "Some channel indices are out of bounds")
  expect_error(combine_channels(x, list(X = c(0, 1))), "Some channel indices are out of bounds")
  expect_error(combine_channels(x, list(X = c("Fz", "Nope"))),
               "Channel(s) Nope of group 'X' not found.", fixed = TRUE)
  expect_error(combine_channels(x, list(A = 1:2), method = "banana"),
               "got \"banana\"", fixed = TRUE)
  expect_error(combine_channels(x, list(A = 1:2), method = 5), "must be a function, or one of")
  expect_error(combine_channels(x, list(A = 1:2), method = function(b) 1:3),
               "must return one number per sample")
  expect_error(combine_channels(x, list(1:2)), "'groups' must be a named list")
  expect_error(combine_channels(x, "Fz"), "'groups' must be a named list")
  expect_error(combine_channels(x, list(A = 1:2, A = 3:4)), "has a name more than once")
  expect_error(combine_channels(x, list(A = character(0))), "non-empty vector of channel names")
  expect_error(combine_channels(x, list(A = c("Pz", "O1")), drop_bad = TRUE),
               "no channels left after leaving out the bad channels")
  expect_error(combine_channels(x, list(Status = 1:2), keep_status = TRUE), "would exist twice")
  expect_error(combine_channels(x, list(A = 1:2), keep_status = NA), "'keep_status' must be TRUE or FALSE")
  expect_error(combine_channels(x, list(A = 1:2), drop_bad = "yes"), "'drop_bad' must be TRUE or FALSE")
})


# ============================================================================
# SUITE 7: eeg_crop()
# ============================================================================

test_that("eeg_crop(2, 5) gives the reference numbers", {
  x <- toy_eeg()
  y <- eeg_crop(x, 2, 5)
  expect_equal(ncol(y$data), 301L)                 # 2.00 s to 5.00 s, both included
  expect_equal(y$data[1, 1], 300)                  # Fz at sample 201: 100 + 200
  expect_equal(y$data[1, 301], 600)
  expect_equal(y$data[7, 1], 700 + 200)
  expect_equal(nrow(y$data), 7L)
  expect_equal(y$times, (0:300) / 100)             # times start at 0 again; the last is 3.0
  expect_equal(y$times[301], 3)
  expect_identical(y$channels, x$channels)
  expect_match(y$preprocessing_history[[1]],
               "^Cropped: samples 201 to 501 of 1000 kept \\(2 s to 5 s of the original; 301 samples\\)$")
  expect_identical(x, toy_eeg())
})

test_that("eeg_crop() shifts and cuts the annotations like the reference does", {
  y <- eeg_crop(toy_eeg(), 2, 5)       # a: 2-3 s, b: 7.5-8.5 s, c: 3.5-5.5 s
  expect_equal(y$annotations$description, c("BAD_a", "BAD_c"))   # b is outside and gone
  expect_equal(y$annotations$onset, c(0, 1.5))
  expect_equal(y$annotations$duration, c(1, 1.51), tolerance = 1e-9)   # c is cut at the last sample + 1
  expect_true(all(is.na(y$annotations$channel)))
})

test_that("eeg_crop() handles annotations that only touch the kept part (reference behaviour)", {
  an <- data.frame(onset = c(1, 5, 4, 2.5), duration = c(1, 1, 0, 0),
                   description = c("ends_at_tmin", "starts_at_tmax", "zero_in", "zero_at_start"),
                   channel = NA_character_, stringsAsFactors = FALSE)
  y <- eeg_crop(toy_eeg(an), 2, 5)
  expect_equal(nrow(y$annotations), 4L)            # all four touch the kept part
  get <- function(d) y$annotations[y$annotations$description == d, ]
  expect_equal(get("ends_at_tmin")$onset, 0);      expect_equal(get("ends_at_tmin")$duration, 0)
  expect_equal(get("zero_at_start")$onset, 0.5);   expect_equal(get("zero_at_start")$duration, 0)
  expect_equal(get("zero_in")$onset, 2);           expect_equal(get("zero_in")$duration, 0)
  expect_equal(get("starts_at_tmax")$onset, 3);    expect_equal(get("starts_at_tmax")$duration, 0.01, tolerance = 1e-9)
  # completely outside: removed
  an2 <- data.frame(onset = c(0, 5.2), duration = c(1.5, 1), description = c("before", "after"),
                    channel = NA_character_, stringsAsFactors = FALSE)
  expect_equal(nrow(eeg_crop(toy_eeg(an2), 2, 5)$annotations), 0L)
})

test_that("eeg_crop() shifts the events, in samples and in seconds", {
  y <- eeg_crop(toy_eeg(), 2, 5)                   # events were at samples 101, 301, 801
  expect_equal(nrow(y$events), 1L)                 # only the one at sample 301 is inside
  expect_equal(y$events$onset, 101L)               # 301 - 200
  expect_equal(y$events$onset_time, 1)
  expect_equal(y$events$description, "t2")
  expect_true(is.integer(y$events$onset))
  z <- eeg_crop(toy_eeg(), 0, 4.99)
  expect_equal(z$events$onset, c(101L, 301L))      # unchanged: the crop starts at 0
})

test_that("eeg_crop() counts samples like the reference: nearest sample, both ends kept", {
  x <- toy_eeg()
  expect_equal(ncol(eeg_crop(x, 2, 5, include_tmax = FALSE)$data), 300L)
  expect_equal(ncol(eeg_crop(x, 2.004, 2.996)$data), 101L)       # 200.4 and 299.6 round to 200 and 300
  expect_equal(eeg_crop(x, 2.004, 2.996)$data[1, 1], 300)        # it starts at sample 201, as crop(2, ..)
  expect_equal(ncol(eeg_crop(x, 3, 3)$data), 1L)                 # tmin == tmax is allowed: one sample
  expect_equal(eeg_crop(x, 3, 3)$data[1, 1], 400)
  expect_equal(ncol(eeg_crop(x, tmin = 2)$data), 800L)
  expect_equal(ncol(eeg_crop(x, tmin = 2, include_tmax = FALSE)$data), 799L)   # also without tmax
  expect_equal(ncol(eeg_crop(x, 2, 2.01)$data), 2L)
  expect_equal(ncol(eeg_crop(x, 0, 9.99)$data), 1000L)
  expect_equal(ncol(eeg_crop(x)$data), 1000L)
})

test_that("cropping the whole recording changes nothing but the history", {
  x <- toy_eeg()
  y <- eeg_crop(x)
  expect_identical(y$data, x$data)
  expect_identical(y$times, x$times)
  expect_equal(y$annotations, x$annotations)
  expect_equal(y$events, x$events)
  same_but_history(y, x)
})

test_that("eeg_crop() stops with the reference messages", {
  x <- toy_eeg()
  expect_error(eeg_crop(x, -1), "tmin (-1) must be >= 0.", fixed = TRUE)
  expect_error(eeg_crop(x, 0, 20), "tmax (20) must be less than or equal to the max time (9.99 s).",
               fixed = TRUE)
  expect_error(eeg_crop(x, 5, 2), "tmin (5) must not be greater than tmax (2).", fixed = TRUE)
  expect_error(eeg_crop(x, 12), "tmin (12) must not be greater than tmax (9.99).", fixed = TRUE)
  expect_error(eeg_crop(x, 9.995), "must not be greater than tmax")          # half a sample after the end
  expect_error(eeg_crop(x, 3, 3, include_tmax = FALSE), "No samples would be left after cropping.",
               fixed = TRUE)
  expect_error(eeg_crop(x, "a"), "'tmin' must be one number")
  expect_error(eeg_crop(x, c(1, 2)), "'tmin' must be one number")
  expect_error(eeg_crop(x, NA_real_), "'tmin' must be one number")
  expect_error(eeg_crop(x, 0, "b"), "'tmax' must be one number")
  expect_error(eeg_crop(x, 0, 5, include_tmax = NA), "'include_tmax' must be TRUE or FALSE")
})

test_that("eeg_crop() works without events and annotations", {
  x <- toy_eeg(); x$events <- NULL; x$annotations <- NULL
  y <- eeg_crop(x, 2, 5)
  expect_equal(ncol(y$data), 301L)
  x2 <- eeganalysis::new_eeg(data = matrix(1:20, 2, 10), channels = c("A", "B"), sampling_rate = 10)
  expect_equal(ncol(eeg_crop(x2, 0.2, 0.5)$data), 4L)
  expect_equal(nrow(eeg_crop(x2, 0.2, 0.5)$events), 0L)
})


# ============================================================================
# SUITE 8: eeg_concatenate()
# ============================================================================

test_that("eeg_concatenate(x, x) gives the reference numbers", {
  x <- toy_eeg()
  y <- eeg_concatenate(x, x)
  expect_equal(dim(y$data), c(7L, 2000L))
  expect_equal(y$data[, 1:1000], x$data)
  expect_equal(y$data[, 1001:2000], x$data)
  expect_equal(y$times, (0:1999) / 100)
  expect_equal(y$times[2000], 19.99)
  expect_identical(y$channels, x$channels)
  expect_identical(y$bads, x$bads)
  expect_identical(y$montage, x$montage)
  expect_match(y$preprocessing_history[[1]], "^Concatenated 2 recordings \\(1000 \\+ 1000 samples\\)$")
  expect_identical(x, toy_eeg())
})

test_that("eeg_concatenate() shifts the annotations and adds one boundary row per join", {
  y <- eeg_concatenate(toy_eeg(), toy_eeg())
  # the reference list: a@2, c@3.5, b@7.5, boundary@10, a@12, c@13.5, b@17.5
  expect_equal(y$annotations$onset, c(2, 3.5, 7.5, 10, 12, 13.5, 17.5))
  expect_equal(y$annotations$duration, c(1, 2, 1, 0, 1, 2, 1))
  expect_equal(y$annotations$description,
               c("BAD_a", "BAD_c", "BAD_b", "BAD_boundary", "BAD_a", "BAD_c", "BAD_b"))
  expect_equal(y$annotations$channel, c(NA, NA, "Cz", NA, NA, NA, "Cz"))
  z <- eeg_concatenate(toy_eeg(), toy_eeg(), toy_eeg())
  expect_equal(z$annotations$onset[z$annotations$description == "BAD_boundary"], c(10, 20))
  expect_equal(ncol(z$data), 3000L)
})

test_that("a boundary annotation makes epochs that span a join get rejected", {
  y <- eeg_concatenate(toy_eeg(), toy_eeg())
  b <- y$annotations[y$annotations$description == "BAD_boundary", ]
  epoch_start <- 9.8; epoch_end <- 10.2
  # the rule the epoching uses (R/epoch2.R): onset < end and onset + duration > start
  expect_true(b$onset < epoch_end && b$onset + b$duration > epoch_start)
  expect_false(b$onset < 9.0 && b$onset + b$duration > 8.0)       # an epoch before the join is fine
})

test_that("epoch_eeg() on a joined recording rejects the epochs that stretch across the join", {
  rec <- function(first, last, ev) {                    # samples first..last of a ramp recording, events at 'ev'
    eeganalysis::new_eeg(
      data = t(vapply(1:3, function(i) i * 100 + (first:last) - 1, numeric(last - first + 1L))),
      channels = c("Fz", "Cz", "Pz"), sampling_rate = 100,
      events = data.frame(onset = ev, onset_time = (ev - 1) / 100, type = 1L,
                          description = "s", stringsAsFactors = FALSE))
  }
  window <- function(e, ...) epoch_eeg(e, events = "all", tmin = -0.2, tmax = 0.8,
                                       baseline = NULL, verbose = FALSE, ...)
  # two recordings of 600 samples (6 s each); events at 1.00 s and 5.49 s, then at 0.04 s and 1.00 s of the second
  joined <- eeg_concatenate(rec(1, 600, c(101L, 550L)), rec(601, 1200, c(5L, 101L)))
  expect_equal(joined$events$onset, c(101L, 550L, 605L, 701L))
  ep <- window(joined)
  # the epoch of 5.49 s ends at 6.29 s and the epoch of 6.04 s starts at 5.84 s: both cross the join at 6 s
  expect_equal(dim(ep$data)[3], 2L)
  expect_equal(ep$events$onset, c(101L, 701L))
  # with a list of descriptions, the join counts only when "BAD_boundary" is in the list
  expect_equal(dim(window(joined, reject_by_annotation = "BAD_boundary")$data)[3], 2L)
  expect_equal(dim(window(joined, reject_by_annotation = "BAD_muscle")$data)[3], 4L)
  # the same four events in one recording that was never joined: nothing is rejected
  whole <- window(rec(1, 1200, c(101L, 550L, 605L, 701L)))
  expect_equal(dim(whole$data)[3], 4L)
})

test_that("eeg_concatenate() shifts the events by the samples before each recording", {
  y <- eeg_concatenate(toy_eeg(), toy_eeg())
  expect_equal(y$events$onset, c(101L, 301L, 801L, 1101L, 1301L, 1801L))
  expect_equal(y$events$onset_time, c(1, 3, 8, 11, 13, 18))
  expect_equal(y$events$description, rep(c("t1", "t2", "t3"), 2))
  short <- eeg_crop(toy_eeg(), 0, 4.99)                           # 500 samples, events at 101 and 301
  z <- eeg_concatenate(short, toy_eeg())
  expect_equal(z$events$onset, c(101L, 301L, 601L, 801L, 1301L))
  expect_match(z$preprocessing_history[[length(z$preprocessing_history)]],
               "^Concatenated 2 recordings \\(500 \\+ 1000 samples\\)$")
})

test_that("eeg_concatenate() takes a list as well, and recordings without events", {
  x <- toy_eeg()
  expect_identical(eeg_concatenate(list(x, x)), eeg_concatenate(x, x))
  a <- toy_eeg(); a$events <- NULL; a$annotations <- NULL
  y <- eeg_concatenate(a, a)
  expect_equal(ncol(y$data), 2000L)
  expect_equal(y$annotations$description, "BAD_boundary")
})

test_that("eeg_concatenate() stops when the recordings do not match", {
  x <- toy_eeg()
  expect_error(eeg_concatenate(x), "Give at least two recordings to join.", fixed = TRUE)
  expect_error(eeg_concatenate(), "Give at least two recordings")
  expect_error(eeg_concatenate(x, 1), "Recording 2 is not an object of class 'eeg'")
  expect_error(eeg_concatenate(x, toy_epochs()), "Recording 2 is not an object of class 'eeg'")

  renamed <- rename_channels(x, c(Fz = "F1"))
  expect_error(eeg_concatenate(x, renamed), "the channel names or their order differ")
  expect_error(eeg_concatenate(x, renamed), "'F1' versus 'Fz'", fixed = TRUE)
  reordered <- pick_channels(x, rev(toy_names))
  expect_error(eeg_concatenate(x, reordered), "order differ")
  fewer <- drop_channels(x, "Fz")
  expect_error(eeg_concatenate(x, fewer), "it has 6 channels, recording 1 has 7", fixed = TRUE)

  other_rate <- x; other_rate$sampling_rate <- 200
  expect_error(eeg_concatenate(x, other_rate), "the sampling rates differ (200 Hz versus 100 Hz)",
               fixed = TRUE)
  other_types <- set_channel_types(x, c(EXG7 = "misc"))
  expect_error(eeg_concatenate(x, other_types), "the channel types differ")
  other_bads <- x; other_bads$bads <- "Fz"
  expect_error(eeg_concatenate(x, other_bads),
               "the bad channels differ (recording 1: Pz, O1; recording 2: Fz)", fixed = TRUE)
  expect_error(eeg_concatenate(x, other_bads), "Mark the same channels as bad", fixed = TRUE)
  no_bads <- x; no_bads$bads <- character(0)
  expect_error(eeg_concatenate(x, no_bads), "recording 2: none", fixed = TRUE)
  other_ref <- x; other_ref$reference <- "average"
  expect_error(eeg_concatenate(x, other_ref), "the reference differs ('average' versus 'original')",
               fixed = TRUE)
  same_bads_other_order <- x; same_bads_other_order$bads <- c("O1", "Pz")
  expect_equal(ncol(eeg_concatenate(x, same_bads_other_order)$data), 2000L)   # the same set is fine
})


# ============================================================================
# SUITE 9: Things that must always hold
# ============================================================================

test_that("picking all channels in the original order changes nothing", {
  x <- toy_eeg()
  y <- pick_channels(x, x$channels)
  expect_identical(y$data, x$data)
  expect_identical(y$channels, x$channels)
  expect_identical(y$channel_types, x$channel_types)
  expect_equal(y$bads, x$bads)
  expect_equal(y$montage, x$montage)
  expect_equal(y$annotations, x$annotations)
  same_but_history(y, x)
})

test_that("renaming a channel and renaming it back gives the original object", {
  x <- toy_eeg()
  y <- rename_channels(rename_channels(x, c(Fz = "A", Cz = "B")), c(A = "Fz", B = "Cz"))
  same_but_history(y, x)
  expect_length(y$preprocessing_history, 2L)
})

test_that("dropping a channel and adding its signal back restores it at the end", {
  x <- toy_eeg()
  veog <- eeganalysis::new_eeg(data = x$data[5, , drop = FALSE], channels = "VEOG",
                               sampling_rate = 100, channel_types = "eog")
  y <- add_channels(drop_channels(x, "VEOG"), veog)
  expect_equal(y$channels, c("Fz", "Cz", "Pz", "O1", "EXG7", "Status", "VEOG"))
  expect_equal(y$data[7, ], x$data[5, ])
  expect_equal(y$data[1:6, ], x$data[-5, ])
  expect_equal(y$channel_types[7], "eog")
  check_consistent(y)
})

test_that("cutting a recording in two and joining the parts gives the recording back", {
  x <- toy_eeg()
  a <- eeg_crop(x, 0, 4.99)
  b <- eeg_crop(x, 5, NULL)
  y <- eeg_concatenate(a, b)
  expect_equal(ncol(y$data), 1000L)
  expect_identical(y$data, x$data)
  expect_equal(y$times, x$times)
  expect_equal(y$events$onset, x$events$onset)                    # the events find their place again
  expect_equal(y$events$onset_time, x$events$onset_time)
  expect_true("BAD_boundary" %in% y$annotations$description)      # the join is marked, nothing else is added
  expect_equal(sum(y$annotations$description == "BAD_boundary"), 1L)
})

test_that("a channel function never changes the input object", {
  x <- toy_eeg(); before <- toy_eeg()
  invisible(drop_channels(x, "Fz")); invisible(pick_channels(x, "Fz"))
  invisible(pick_types(x, "eeg")); invisible(rename_channels(x, toupper))
  invisible(combine_channels(x, list(A = 1:2))); invisible(eeg_crop(x, 1, 2))
  invisible(add_channels(x, heart())); invisible(eeg_concatenate(x, x))
  expect_identical(x, before)
})

test_that("every channel function leaves fields it does not know about alone", {
  x <- toy_eeg(); x$extra <- list(a = 1)
  expect_equal(drop_channels(x, "Fz")$extra, list(a = 1))
  expect_equal(rename_channels(x, toupper)$extra, list(a = 1))
  expect_equal(eeg_crop(x, 1, 2)$extra, list(a = 1))
})

test_that("the result of a channel function is still an eeg object that other functions accept", {
  y <- drop_channels(toy_eeg(), "VEOG")
  expect_s3_class(y, "eeg")
  expect_output(print(y))
  z <- eeg_concatenate(eeg_crop(toy_eeg(), 0, 3), eeg_crop(toy_eeg(), 4, 7))
  expect_s3_class(z, "eeg")
  expect_output(print(z))
})


test_that("row names of the data, when there are any, follow the channels", {
  x <- toy_eeg(); rownames(x$data) <- toy_names
  y <- pick_channels(x, c("Cz", "Fz"))
  expect_equal(rownames(y$data), c("Cz", "Fz"))
  z <- rename_channels(x, c(Fz = "A"))
  expect_equal(rownames(z$data)[1], "A")
  expect_equal(rownames(drop_channels(x, "Pz")$data), setdiff(toy_names, "Pz"))
  expect_null(rownames(drop_channels(toy_eeg(), "Fz")$data))     # none before, none after
  cmb <- combine_channels(x, list(A = 1:2))
  expect_equal(rownames(cmb$data), "A")
  expect_equal(rownames(eeg_concatenate(x, x)$data), toy_names)
})

test_that("a long list of names in a history line is shortened", {
  expect_equal(eeganalysis:::.name_list(letters[1:10]), "a, b, c, d, e, f, g, h, i, j")
  expect_equal(eeganalysis:::.name_list(letters[1:11]), "a, b, c, d, e, f, g, h, i, j ... (1 more)")
  expect_equal(eeganalysis:::.name_list(letters[1:4], max = 2L), "a, b ... (2 more)")
  d <- matrix(rnorm(15 * 20), 15, 20)
  x <- eeganalysis::new_eeg(data = d, channels = paste0("c", 1:15), sampling_rate = 10)
  y <- drop_channels(x, paste0("c", 1:12))
  expect_equal(y$preprocessing_history[[1]],
               "Channels dropped: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 ... (2 more) (12 of 15)")
})

# ============================================================================
# SUITE 10: Chains with the EOG tools (gaps 3.3 / 3.4)
# ============================================================================
# The eye channel is found through channel_types == "eog". It must be there
# when the model is fitted AND when it is applied; after that it can go.

test_that("the EOG regression still finds the eye channel after picking channels", {
  x <- toy_noise()
  model <- fit_eog_regression(x)
  expect_equal(model$ch_names_artifact, "VEOG")
  expect_equal(model$ch_names, c("Fz", "Cz"))                     # the good scalp channels
  y <- pick_channels(x, c("Fz", "Cz", "VEOG"))
  expect_equal(fit_eog_regression(y)$ch_names_artifact, "VEOG")
  expect_equal(fit_eog_regression(y)$ch_names, c("Fz", "Cz"))
})

test_that("after pick_types(eeg) the EOG regression says there is no eye channel", {
  expect_error(fit_eog_regression(pick_types(toy_noise(), "eeg")), "no EOG channel found")
})

test_that("a model can be applied after the other channels are dropped, but not after the eye channel", {
  x <- toy_noise()
  model <- fit_eog_regression(x)
  ok <- apply_eog_regression(model, pick_types(x, c("eeg", "eog")))
  expect_s3_class(ok, "eeg")
  expect_equal(ok$channels, c("Fz", "Cz", "VEOG"))
  expect_error(apply_eog_regression(model, drop_channels(x, "VEOG")),
               "missing EOG channel(s) that 'model' was fit on: VEOG", fixed = TRUE)
  expect_error(apply_eog_regression(model, drop_channels(x, "Cz")),
               "missing target channel(s) that 'model' was fit on: Cz", fixed = TRUE)
})

test_that("renaming a channel after fitting makes the model complain, a new fit works", {
  x <- toy_noise()
  model <- fit_eog_regression(x)
  y <- rename_channels(x, c(VEOG = "VEOG2"))
  expect_error(apply_eog_regression(model, y), "missing EOG channel(s)", fixed = TRUE)
  expect_equal(fit_eog_regression(y)$ch_names_artifact, "VEOG2")  # the type came along with the channel
})

test_that("set_channel_types() then pick_types() is the way to drop eye and heart channels", {
  x <- eeganalysis::new_eeg(data = matrix(rnorm(5 * 100), 5, 100),
                            channels = c("Fz", "Cz", "EXG1", "EXG2", "EXG7"), sampling_rate = 100)
  x <- set_channel_types(x, c(EXG1 = "eog", EXG2 = "eog", EXG7 = "ecg"))
  expect_equal(pick_types(x, c("eog", "ecg"))$channels, c("EXG1", "EXG2", "EXG7"))
  expect_equal(pick_types(x, "eeg")$channels, c("Fz", "Cz"))
})


# ============================================================================
# SUITE 11: Wrong input
# ============================================================================

test_that("every function refuses an object that is not an eeg object", {
  expect_error(drop_channels(1:3, "a"), "'eeg' must be an object of class 'eeg' (see new_eeg()) or 'eeg_epochs'.",
               fixed = TRUE)
  expect_error(pick_channels("x", "a"), "'eeg' must be an object of class 'eeg'")
  expect_error(pick_types(list(), "eeg"), "'eeg' must be an object of class 'eeg'")
  expect_error(rename_channels(NULL, toupper), "'eeg' must be an object of class 'eeg'")
  expect_error(add_channels(1, toy_eeg()), "'eeg' must be an object of class 'eeg' (see new_eeg()).",
               fixed = TRUE)
  expect_error(combine_channels(1, list(A = 1)), "'eeg' must be an object of class 'eeg'")
  expect_error(eeg_crop(1), "'eeg' must be an object of class 'eeg'")
  structure_evoked <- structure(list(channels = "A"), class = "eeg_evoked")
  expect_error(drop_channels(structure_evoked, "A"), "'eeg' must be an object of class 'eeg'")
})

test_that("functions that need continuous data say so when given epochs", {
  e <- toy_epochs()
  expect_error(add_channels(e, toy_eeg()), "add_channels() works on continuous 'eeg' objects only",
               fixed = TRUE)
  expect_error(combine_channels(e, list(A = 1:2)), "combine_channels() works on continuous 'eeg' objects only",
               fixed = TRUE)
  expect_error(eeg_crop(e), "eeg_crop() works on continuous 'eeg' objects only", fixed = TRUE)
})
