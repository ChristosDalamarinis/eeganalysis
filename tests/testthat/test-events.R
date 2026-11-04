# ============================================================================
#                    Test File for events.R
# ============================================================================
#
# Functions tested:
#   add_events()   - write events (a time and a code each) into the events
#                    table of a continuous recording
#   (and the helper behind it: .events_existing())
#
# THE TEST RECORDING
#   ev_eeg(): 3 channels (Fz, Cz, Pz, all eeg), 100 Hz and 1000 samples (10 s)
#   unless a block says otherwise. The values are smooth curves; only the
#   events matter in most blocks. The expected numbers are worked out by hand
#   (sample = round(time * rate) + 1, onset_time = (sample - 1) / rate), so
#   they can be checked by eye.
#
# Test suites:
#   1. Time to sample
#   2. Codes and descriptions
#   3. The old events and the new rows
#   4. The empty events table (type stored as text)
#   5. Nothing else changes
#   6. Two events on one sample
#   7. Wrong input
#   8. Used by the other steps (epoching, stimulation artifact, crop, downsample)
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

ev_table <- function(onsets, types = rep(1L, length(onsets)), fs = 100) {
  data.frame(onset = as.integer(onsets), onset_time = (onsets - 1) / fs,
             type = as.integer(types), description = paste0("Trigger: ", types),
             stringsAsFactors = FALSE)
}

ev_data <- function(n = 1000L) {
  s <- seq_len(n)
  rbind(Fz = sin(0.1 * s), Cz = cos(0.07 * s), Pz = 0.001 * s)
}

ev_eeg <- function(events = NULL, fs = 100, n = 1000L, data = ev_data(n), ...) {
  eeganalysis::new_eeg(data = data, channels = rownames(data),
                       sampling_rate = fs, events = events, ...)
}


# ============================================================================
#   1. TIME TO SAMPLE
# ============================================================================

test_that("512 Hz: time 2 s is sample 1025, time 0 is sample 1, and the table is sorted", {
  x   <- ev_eeg(fs = 512, n = 2000L)
  out <- add_events(x, times = c(2, 0, 1.5), code = 9)
  ev  <- out$events

  expect_equal(ev$onset, c(1L, 769L, 1025L))            # 0 * 512 + 1, 1.5 * 512 + 1, 2 * 512 + 1
  expect_equal(ev$onset_time, c(0, 1.5, 2))
  expect_equal(ev$type, c(9L, 9L, 9L))
  expect_equal(ev$description, rep("Trigger: 9", 3))
  expect_type(ev$onset, "integer")
  expect_type(ev$type, "integer")
  expect_type(ev$onset_time, "double")
  expect_equal(rownames(ev), as.character(1:3))
})

test_that("a time is moved to the nearest sample and onset_time is the time of that sample", {
  x   <- ev_eeg()                                        # 100 Hz: one sample = 0.01 s
  out <- add_events(x, times = c(0.004, 0.006, 0.024, 0.0361), code = 1)
  # 0.4 -> 0, 0.6 -> 1, 2.4 -> 2, 3.61 -> 4 samples after the first one
  expect_equal(out$events$onset, c(1L, 2L, 3L, 5L))
  expect_equal(out$events$onset_time, c(0, 0.01, 0.02, 0.04))
})

test_that("the last sample is accepted, one past it is refused, and the message says so", {
  x <- ev_eeg()                                          # last sample 1000 = 9.99 s
  expect_equal(add_events(x, 9.99, 1)$events$onset, 1000L)
  expect_equal(add_events(x, 9.994, 1)$events$onset, 1000L)   # 999.4 rounds to 999
  expect_error(add_events(x, 9.996, 1), "outside the recording \\(0 to 9.99 s\\)")
  expect_error(add_events(x, 10, 1), "the time 10 s is outside the recording")
  expect_error(add_events(x, c(5, 10, 11), 1), "time 10 s .* 1 more time\\(s\\)")
})

test_that("times may be whole numbers, and need not be in order", {
  x   <- ev_eeg()
  out <- add_events(x, times = c(7L, 2L, 5L), code = 3L)
  expect_equal(out$events$onset, c(201L, 501L, 701L))
})


# ============================================================================
#   2. CODES AND DESCRIPTIONS
# ============================================================================

test_that("one code for all events, or one code per time; codes are stored as integers", {
  x <- ev_eeg()
  same <- add_events(x, times = c(1, 2, 3), code = 9)
  expect_equal(same$events$type, c(9L, 9L, 9L))
  each <- add_events(x, times = c(3, 1, 2), code = c(30, 10, 20))
  expect_equal(each$events$type, c(10L, 20L, 30L))           # sorted with their times
  expect_equal(each$events$onset, c(101L, 201L, 301L))
  expect_equal(add_events(x, 1, 65535)$events$type, 65535L)
  expect_equal(add_events(x, 1, 1)$events$type, 1L)
  expect_type(each$events$type, "integer")
})

test_that("the default description is the reader's text, own descriptions are used as given", {
  x <- ev_eeg()
  expect_equal(add_events(x, c(1, 2), code = c(1, 2))$events$description,
               c("Trigger: 1", "Trigger: 2"))
  expect_equal(add_events(x, c(1, 2), code = 9, description = "pulse")$events$description,
               c("pulse", "pulse"))
  expect_equal(add_events(x, c(1, 2), code = 9, description = c("first", "second"))$events$description,
               c("first", "second"))
  # descriptions are sorted along with their events
  expect_equal(add_events(x, c(2, 1), code = 9, description = c("late", "early"))$events$description,
               c("early", "late"))
})


# ============================================================================
#   3. THE OLD EVENTS AND THE NEW ROWS
# ============================================================================

test_that("the old events stay, the new ones are merged in, the result is sorted, row names reset", {
  x   <- ev_eeg(events = ev_table(c(501L, 101L), c(2L, 1L)))        # not in order on purpose
  out <- add_events(x, times = c(3, 7), code = 9)
  expect_equal(out$events$onset, c(101L, 301L, 501L, 701L))
  expect_equal(out$events$type, c(1L, 9L, 2L, 9L))
  expect_equal(out$events$description, c("Trigger: 1", "Trigger: 9", "Trigger: 2", "Trigger: 9"))
  expect_equal(rownames(out$events), as.character(1:4))
  expect_equal(out$events$onset_time, c(1, 3, 5, 7))
})

test_that("a column the new rows do not have (label) gets NA of the same kind; old values stay", {
  ev <- ev_table(c(101L, 501L), c(1L, 2L))
  ev$label <- c("start", "stop")
  ev$score <- c(0.5, 1.5)
  ev$group <- factor(c("a", "b"))
  x   <- ev_eeg(events = ev)
  out <- add_events(x, times = 3, code = 9)
  expect_equal(names(out$events), names(ev))                       # same columns, same order
  expect_equal(out$events$label, c("start", NA, "stop"))
  expect_equal(out$events$score, c(0.5, NA, 1.5))
  expect_s3_class(out$events$group, "factor")
  expect_equal(as.character(out$events$group), c("a", NA, "b"))
  expect_type(out$events$label, "character")
})

test_that("replace = TRUE throws the old table away, with all its extra columns", {
  ev <- ev_table(c(101L, 501L), c(1L, 2L))
  ev$label <- c("start", "stop")
  x   <- ev_eeg(events = ev)
  out <- add_events(x, times = c(1, 4, 6), code = 9, replace = TRUE)
  expect_equal(names(out$events), c("onset", "onset_time", "type", "description"))
  expect_equal(out$events$onset, c(101L, 401L, 601L))             # time 1 = the old event's sample: no clash now
  expect_equal(out$events$type, c(9L, 9L, 9L))
  expect_equal(out$preprocessing_history[[1]],
               "Events replaced: 2 old event(s) removed, 3 new event(s), code(s) 9")
})


# ============================================================================
#   4. THE EMPTY EVENTS TABLE (type stored as text)
# ============================================================================

test_that("an object made without events gets integer codes, and matching by code works", {
  x <- ev_eeg()
  expect_type(x$events$type, "character")                         # why this test exists
  expect_equal(nrow(x$events), 0L)

  out <- add_events(x, times = c(3, 6), code = 9)
  expect_type(out$events$type, "integer")
  ep <- epoch_eeg(out, events = 9, tmin = -0.2, tmax = 0.5, baseline = NULL, verbose = FALSE)
  expect_equal(dim(ep$data)[3], 2L)
})

test_that("a missing events table (NULL) counts as empty", {
  x <- ev_eeg()
  x$events <- NULL
  out <- add_events(x, times = 3, code = 9)
  expect_equal(out$events$onset, 301L)
  expect_equal(names(out$events), c("onset", "onset_time", "type", "description"))
  # and with replace = TRUE the history says 0 old events
  expect_match(add_events(x, times = 3, code = 9, replace = TRUE)$preprocessing_history[[1]],
               "0 old event\\(s\\) removed")
})

test_that("codes stored as text that are numbers are turned into numbers; other text is refused", {
  ev <- ev_table(c(101L, 301L), c(1L, 2L))
  ev$type <- c("1", "2")
  x   <- ev_eeg(events = ev)
  out <- add_events(x, times = 2, code = 9)
  expect_equal(out$events$type, c(1L, 9L, 2L))
  expect_type(out$events$type, "integer")

  ev$type <- c("go", "stop")
  expect_error(add_events(ev_eeg(events = ev), times = 2, code = 9), "not numbers \\(go, stop\\)")
})

test_that("columns of the old table that are stored in other kinds come out as integer, number and text", {
  ev <- data.frame(onset = c(101, 501), onset_time = c(1L, 5L), type = c(1, 2),
                   description = factor(c("a", "b")))
  out <- add_events(ev_eeg(events = ev), times = 3, code = 9)
  expect_type(out$events$onset, "integer")
  expect_type(out$events$onset_time, "double")
  expect_type(out$events$type, "integer")
  expect_type(out$events$description, "character")
  expect_equal(out$events$description, c("a", "Trigger: 9", "b"))
  expect_equal(out$events$onset, c(101L, 301L, 501L))
})

test_that("an events table without the four columns is refused", {
  x <- ev_eeg()
  x$events <- data.frame(onset = 101L, type = 1L)
  expect_error(add_events(x, times = 2, code = 9), "columns onset, onset_time, type and description")
  x$events <- "not a table"
  expect_error(add_events(x, times = 2, code = 9), "columns onset, onset_time, type and description")
  # with replace = TRUE the old table is thrown away, so its shape does not matter
  expect_equal(add_events(x, times = 2, code = 9, replace = TRUE)$events$onset, 201L)
})


# ============================================================================
#   5. NOTHING ELSE CHANGES
# ============================================================================

test_that("only events and the history change; the input is untouched", {
  ann <- data.frame(onset = 2, duration = 0.5, description = "BAD_x",
                    channel = NA_character_, stringsAsFactors = FALSE)
  x <- ev_eeg(events = ev_table(101L, 1L), bads = "Pz", annotations = ann,
              montage = create_montage(c("Fz", "Cz", "Pz")),
              metadata = list(subject = "S01"), reference = "Common Average")
  x$preprocessing_history <- list("An earlier step")
  copy <- x
  out  <- add_events(x, times = c(3, 4), code = 9)

  expect_identical(x, copy)
  for (f in setdiff(names(x), c("events", "preprocessing_history"))) {
    expect_identical(out[[f]], x[[f]], info = f)
  }
  expect_equal(length(out$preprocessing_history), 2L)
  expect_equal(out$preprocessing_history[[1]], "An earlier step")
  expect_equal(out$preprocessing_history[[2]], "Events added: 2 event(s), code(s) 9")
})

test_that("the history line names the codes", {
  x <- ev_eeg()
  expect_equal(add_events(x, c(1, 2, 3), code = c(2, 1, 2))$preprocessing_history[[1]],
               "Events added: 3 event(s), code(s) 1, 2")
  expect_equal(add_events(x, 1, code = 9)$preprocessing_history[[1]],
               "Events added: 1 event(s), code(s) 9")
})


# ============================================================================
#   6. TWO EVENTS ON ONE SAMPLE
# ============================================================================

test_that("two new events that round to the same sample are refused", {
  x <- ev_eeg()
  expect_error(add_events(x, c(1.001, 1.004), 1), "two of the new events fall on the same sample \\(sample 101\\)")
  expect_error(add_events(x, c(1, 1), c(1, 2)), "same sample")
  expect_equal(add_events(x, c(1.001, 1.02), 1)$events$onset, c(101L, 103L))   # neighbours are fine
})

test_that("a new event on a sample that already has an event is refused (unless replace = TRUE)", {
  x <- ev_eeg(events = ev_table(c(101L, 501L), c(1L, 2L)))
  expect_error(add_events(x, 1, 9), "an event already sits at sample 101 \\(time 1 s\\)")
  expect_error(add_events(x, c(3, 5), 9), "sample 501 \\(time 5 s\\)")
  expect_equal(add_events(x, 1.01, 9)$events$onset, c(101L, 102L, 501L))
  expect_equal(add_events(x, 1, 9, replace = TRUE)$events$type, 9L)
})


# ============================================================================
#   7. WRONG INPUT
# ============================================================================

test_that("epochs and other objects are refused", {
  rec <- add_events(ev_eeg(), times = c(3, 6), code = 9)
  ep  <- epoch_eeg(rec, events = 9, tmin = -0.2, tmax = 0.5, baseline = NULL, verbose = FALSE)
  expect_error(add_events(ep, times = 1, code = 1), "continuous 'eeg' objects only")
  expect_error(add_events(list(a = 1), times = 1, code = 1), "object of class 'eeg'")
})

test_that("'times' is checked", {
  x <- ev_eeg()
  expect_error(add_events(x, times = "5", code = 1), "non-empty numeric vector")
  expect_error(add_events(x, times = numeric(0), code = 1), "non-empty numeric vector")
  expect_error(add_events(x, times = c(1, NA), code = 1), "must not hold NA, NaN or Inf")
  expect_error(add_events(x, times = Inf, code = 1), "must not hold NA, NaN or Inf")
  expect_error(add_events(x, times = c(2, -0.5), code = 1), "must not be negative \\(first one: -0.5 s\\)")
})

test_that("'code' is checked", {
  x <- ev_eeg()
  expect_error(add_events(x, 1, "a"), "whole numbers from 1 to 65535")
  expect_error(add_events(x, 1, 0), "whole numbers from 1 to 65535")
  expect_error(add_events(x, 1, 65536), "whole numbers from 1 to 65535")
  expect_error(add_events(x, 1, 1.5), "whole numbers from 1 to 65535")
  expect_error(add_events(x, 1, NA_real_), "whole numbers from 1 to 65535")
  expect_error(add_events(x, 1, numeric(0)), "whole numbers from 1 to 65535")
  expect_error(add_events(x, c(1, 2), c(1, 2, 3)), "one number or one number per time \\(2 time\\(s\\), 3 code\\(s\\)\\)")
})

test_that("'description' and 'replace' are checked", {
  x <- ev_eeg()
  expect_error(add_events(x, 1, 1, description = 5), "'description' must be text")
  expect_error(add_events(x, 1, 1, description = NA_character_), "'description' must be text")
  expect_error(add_events(x, c(1, 2), 1, description = c("a", "b", "c")), "'description' must be text")
  expect_error(add_events(x, 1, 1, replace = NA), "'replace' must be TRUE or FALSE")
  expect_error(add_events(x, 1, 1, replace = "yes"), "'replace' must be TRUE or FALSE")
  expect_error(add_events(x, 1, 1, replace = c(TRUE, FALSE)), "'replace' must be TRUE or FALSE")
})


# ============================================================================
#   8. USED BY THE OTHER STEPS
# ============================================================================

test_that("resting state: one event every 2 s gives equal windows for epoch_eeg()", {
  x   <- ev_eeg()
  out <- add_events(x, times = seq(0, 8, by = 2), code = 1)
  expect_equal(out$events$onset, c(1L, 201L, 401L, 601L, 801L))
  ep <- epoch_eeg(out, events = 1, tmin = 0, tmax = 1.5, baseline = NULL, verbose = FALSE)
  expect_equal(dim(ep$data), c(3L, 151L, 5L))
  # the third window starts at 4 s = sample 401
  expect_equal(unname(ep$data[1, , 3]), unname(x$data["Fz", 401:551]))
})

test_that("added stimulation times work with fix_stim_artifact()", {
  d <- ev_data(2000L)
  d["Fz", 1001:1007] <- c(0.36, 5000, 4000, 3000, 2000, 1500, 0.96)
  x   <- ev_eeg(fs = 512, data = d)
  out <- add_events(x, times = 1000 / 512, code = 9)            # exactly sample 1001
  expect_equal(out$events$onset, 1001L)
  fixed <- fix_stim_artifact(out, events = 9)
  expect_equal(unname(fixed$data["Fz", 1001:1007]),
               c(0.36, 0.46, 0.56, 0.66, 0.76, 0.86, 0.96), tolerance = 1e-12)
})

test_that("added events move with eeg_crop() and downsample()", {
  x   <- add_events(ev_eeg(), times = c(2, 4, 6), code = c(1, 2, 3))
  cut <- eeg_crop(x, tmin = 3, tmax = 8)
  expect_equal(cut$events$onset, c(101L, 301L))                  # the 4 s and 6 s events, now at 1 s and 3 s
  expect_equal(cut$events$onset_time, c(1, 3))
  expect_equal(cut$events$type, c(2L, 3L))

  low <- downsample(x, target_rate = 50, verbose = FALSE)
  expect_equal(low$events$onset_time, c(2, 4, 6))
  expect_equal(low$events$onset, c(101L, 201L, 301L))
  expect_equal(low$events$type, c(1L, 2L, 3L))
})

test_that("the object with added events prints", {
  out <- add_events(ev_eeg(), times = c(2, 4), code = 9)
  expect_output(print(out), "Total events")
})
