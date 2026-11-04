# ============================================================================
#                    Test File for stim_artifact.R
# ============================================================================
#
# Functions tested:
#   fix_stim_artifact()   - replace the short spike a stimulator leaves
#   (and the helpers behind it: .stim_window_samples(), .stim_event_samples(),
#    .stim_block())
#
# THE TEST RECORDING
#   stim_data(): 4 channels (Fz and Cz = eeg, VEOG = eog, Status = status),
#   512 Hz, 2000 samples. Fz and Cz are smooth, non-straight curves
#   (50 * sin(0.37 s) + 0.01 s and 20 * cos(0.11 s) - 0.003 s, s = 1..2000),
#   so a straight-line fill can never be mistaken for the original data. Most
#   blocks write a few hand-picked numbers into the window so that the
#   expected result can be checked by eye. The "matches the reference"
#   blocks use the curves as they are; their expected numbers come from the
#   reference implementation (run once on this exact recording and written in
#   here as literals, so these tests need no other software).
#
# Test suites:
#   1. The three fills, with numbers that can be checked by eye
#   2. Matches the reference implementation
#   3. Nothing else changes
#   4. Which channels
#   5. Which events
#   6. Epochs
#   7. Wrong input
#   8. Why it has to come first (filter, epoch thresholds)
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

stim_data <- function(n = 2000L) {
  s <- seq_len(n)
  rbind(Fz     = 50 * sin(0.37 * s) + 0.01 * s,
        Cz     = 20 * cos(0.11 * s) - 0.003 * s,
        VEOG   = 5 * sin(0.05 * s),
        Status = rep(0, n))
}

stim_types <- c("eeg", "eeg", "eog", "status")

stim_events <- function(onsets, types = rep(1L, length(onsets)), fs = 512) {
  data.frame(onset = as.integer(onsets), onset_time = (onsets - 1) / fs,
             type = as.integer(types),
             description = paste0("Trigger: ", types),
             stringsAsFactors = FALSE)
}

stim_eeg <- function(onsets = 1001L, types = rep(1L, length(onsets)),
                     fs = 512, data = stim_data(), ...) {
  eeganalysis::new_eeg(data = data, channels = rownames(data),
                       sampling_rate = fs, channel_types = stim_types,
                       events = stim_events(onsets, types, fs), ...)
}

# the two-channel recording that the reference numbers belong to
ref_eeg <- function(onsets, types = rep(1L, length(onsets))) {
  eeganalysis::new_eeg(data = stim_data()[1:2, ], channels = c("Fz", "Cz"),
                       sampling_rate = 512,
                       events = stim_events(onsets, types))
}

# the sample numbers where one channel differs between two recordings
changed <- function(a, b, row = "Fz") unname(which(a$data[row, ] != b$data[row, ]))


# ============================================================================
#   1. THE THREE FILLS
# ============================================================================

test_that("linear: a straight line between the two ends of the window (512 Hz, 0 to 0.01 s)", {
  d <- stim_data()
  d["Fz", 1001:1008] <- c(0.36, 5000, 4000, 3000, 2000, 1500, 0.96, 777)
  d["Cz", 1001:1007] <- c(10, 9000, 9000, 9000, 9000, 9000, 40)
  x   <- stim_eeg(1001L, data = d)
  out <- fix_stim_artifact(x, events = 1)

  # 0.01 s at 512 Hz = 6 samples: the line runs from sample 1001 to sample 1007
  expect_equal(unname(out$data["Fz", 1001:1008]),
               c(0.36, 0.46, 0.56, 0.66, 0.76, 0.86, 0.96, 777), tolerance = 1e-12)
  expect_equal(unname(out$data["Cz", 1001:1007]),
               c(10, 15, 20, 25, 30, 35, 40), tolerance = 1e-12)
  # the first sample of the window is replaced by itself and the last is only
  # read: 5 of the 6 samples really change
  expect_equal(changed(out, x), 1002:1006)
})

test_that("linear: at 1000 Hz the default window is 10 samples wide and 9 of them change", {
  d <- stim_data(3000L)
  d["Fz", 1001:1012] <- c(0, rep(5000, 9), 100, 777)
  x   <- stim_eeg(1001L, fs = 1000, data = d)
  out <- fix_stim_artifact(x, events = 1)
  expect_equal(unname(out$data["Fz", 1001:1012]), c(seq(0, 100, by = 10), 777), tolerance = 1e-12)
  expect_equal(changed(out, x), 1002:1010)
})

test_that("window: the middle goes to zero, the two edge samples stay", {
  d <- stim_data()
  d["Fz", 1001:1008] <- 5000
  d["Cz", 1001:1008] <- -2000
  x   <- stim_eeg(1001L, data = d)
  out <- fix_stim_artifact(x, events = 1, mode = "window")

  # six samples are multiplied by 1, 0.25, 0, 0, 0.25, 1; sample 1007 is outside
  expect_equal(unname(out$data["Fz", 1001:1008]), c(5000, 1250, 0, 0, 1250, 5000, 5000, 5000))
  expect_equal(unname(out$data["Cz", 1001:1008]), c(-2000, -500, 0, 0, -500, -2000, -2000, -2000))
})

test_that("window: a window of exactly 4 samples works and one of 3 is refused", {
  d <- stim_data()
  d["Fz", 1001:1006] <- 5000
  x   <- stim_eeg(1001L, data = d)
  out <- fix_stim_artifact(x, events = 1, mode = "window", tmax = 4 / 512)
  expect_equal(unname(out$data["Fz", 1001:1006]), c(5000, 1250, 1250, 5000, 5000, 5000))
  expect_error(fix_stim_artifact(x, events = 1, mode = "window", tmax = 3 / 512),
               "needs at least 4")
})

test_that("constant: every sample of the window becomes the mean of the baseline", {
  d <- stim_data()
  d["Fz", 991:998]   <- 1:8                 # baseline -0.02 to -0.005 s = 10 to 3 samples before the event
  d["Cz", 991:998]   <- 10 * (1:8)
  d["Fz", 1001:1008] <- 5000
  d["Cz", 1001:1008] <- 5000
  x   <- stim_eeg(1001L, data = d)
  out <- fix_stim_artifact(x, events = 1, mode = "constant", baseline = c(-0.02, -0.005))

  # the end of the baseline is not included: the 8 samples 991 to 998
  expect_equal(unname(out$data["Fz", 1001:1008]), c(rep(4.5, 6), 5000, 5000))
  expect_equal(unname(out$data["Cz", 1001:1008]), c(rep(45, 6), 5000, 5000))
  # all six samples of the window change (the first one too)
  expect_equal(changed(out, x), 1001:1006)
})

test_that("a window that starts before the event: ceiling of rate * tmin", {
  d <- stim_data()
  d["Fz", 999:1007] <- c(1, rep(5000, 7), 9)
  x <- stim_eeg(1001L, data = d)

  # -0.005 s * 512 = -2.56 -> -2: the line runs from sample 999 to sample 1007
  out <- fix_stim_artifact(x, events = 1, tmin = -0.005, tmax = 0.01)
  expect_equal(unname(out$data["Fz", 999:1007]), 1:9, tolerance = 1e-12)
  expect_equal(changed(out, x), 1000:1006)

  # -0.002 s * 512 = -1.024 -> -1: the line starts one sample before the event
  out2 <- fix_stim_artifact(stim_eeg(1001L), events = 1, tmin = -0.002, tmax = 0.01)
  expect_equal(changed(out2, stim_eeg(1001L)), 1001:1006)
})

test_that("a spike that already sits on the first sample of the window is not repaired by the default window", {
  d <- stim_data()
  d["Fz", 1001:1007] <- c(5500, 5000, 4000, 3000, 2000, 1500, 0.96)
  x <- stim_eeg(1001L, data = d)

  # the line starts from the spike value: a ramp (5500 down to 0.96 in six
  # equal steps of -916.5067), not a repair
  out <- fix_stim_artifact(x, events = 1)
  expect_equal(unname(out$data["Fz", 1001:1007]),
               c(5500, 4583.4933, 3666.9867, 2750.48, 1833.9733, 917.4667, 0.96),
               tolerance = 1e-7)

  # starting the window one sample earlier gives a line between clean samples
  d["Fz", 1000]      <- 0.36
  out <- fix_stim_artifact(stim_eeg(1001L, data = d), events = 1, tmin = -0.002, tmax = 0.01)
  expect_lt(max(abs(out$data["Fz", 1000:1007])), 1)
})

test_that("two events closer together than the window are handled one after the other", {
  d <- stim_data()
  d["Fz", 1001:1010] <- c(0, 100, 100, 100, 100, 100, 600, 100, 100, 6)
  x   <- stim_eeg(c(1001L, 1004L), data = d)
  out <- fix_stim_artifact(x, events = 1)

  # first event: line 0 -> 600 over 1001..1007; second event: from the value now
  # at 1004 (300) to the value at 1010 (6), slope -49
  expect_equal(unname(out$data["Fz", 1001:1010]),
               c(0, 100, 200, 300, 251, 202, 153, 104, 55, 6), tolerance = 1e-9)
})


# ============================================================================
#   2. MATCHES THE REFERENCE IMPLEMENTATION
# ============================================================================

# The reference numbers: Fz and Cz, samples 998 to 1010, event at sample 1001.
ref_window <- function(x, from = 998L, to = 1010L) unname(x$data[, from:to])

test_that("matches the reference: linear, default window", {
  expected <- rbind(Fz = c(-39.643090323, -34.0590873126, -22.513246978, -6.5668911915, 4.1388007503, 14.8444926921, 25.5501846339, 36.2558765756, 46.9615685174, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942),
                    Cz = c(-22.6859385301, -22.9537943081, -22.9804162663, -22.7655188671, -21.7627951018, -20.7600713366, -19.7573475713, -18.754623806, -17.7519000408, -16.7491762755, -15.0725540435, -13.2502912507, -11.3044513161))
  out <- fix_stim_artifact(ref_eeg(1001L), events = 1)
  expect_equal(ref_window(out), unname(expected), tolerance = 1e-9)
})

test_that("matches the reference: linear, window -0.005 to 0.012 s", {
  expected <- rbind(Fz = c(-39.643090323, -34.0590873126, -24.83933773, -15.6195881474, -6.3998385648, 2.8199110178, 12.0396606005, 21.2594101831, 30.4791597657, 39.6989093483, 48.9186589309, 34.9134271167, 17.5484608942),
                    Cz = c(-22.6859385301, -22.9537943081, -22.0781009454, -21.2024075826, -20.3267142199, -19.4510208572, -18.5753274944, -17.6996341317, -16.823940769, -15.9482474062, -15.0725540435, -13.2502912507, -11.3044513161))
  out <- fix_stim_artifact(ref_eeg(1001L), events = 1, tmin = -0.005, tmax = 0.012)
  expect_equal(ref_window(out), unname(expected), tolerance = 1e-9)
})

test_that("matches the reference: window mode, 0 to 0.012 s", {
  expected <- rbind(Fz = c(-39.643090323, -34.0590873126, -22.513246978, -6.5668911915, 2.905767263, 0.0, 0.0, 0.0, 14.9934490183, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942),
                    Cz = c(-22.6859385301, -22.9537943081, -22.9804162663, -22.7655188671, -5.5779340029, 0.0, 0.0, 0.0, -4.5649637506, -16.7491762755, -15.0725540435, -13.2502912507, -11.3044513161))
  out <- fix_stim_artifact(ref_eeg(1001L), events = 1, mode = "window", tmax = 0.012)
  expect_equal(ref_window(out), unname(expected), tolerance = 1e-9)
})

test_that("matches the reference: constant mode, baseline -0.02 to -0.005 s", {
  expected <- rbind(Fz = c(-39.643090323, -34.0590873126, -22.513246978, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942),
                    Cz = c(-22.6859385301, -22.9537943081, -22.9804162663, -19.3877333637, -19.3877333637, -19.3877333637, -19.3877333637, -19.3877333637, -19.3877333637, -16.7491762755, -15.0725540435, -13.2502912507, -11.3044513161))
  out <- fix_stim_artifact(ref_eeg(1001L), events = 1, mode = "constant",
                           baseline = c(-0.02, -0.005))
  expect_equal(ref_window(out), unname(expected), tolerance = 1e-9)
})

test_that("matches the reference: two events 3 samples apart", {
  expected <- rbind(Fz = c(-39.643090323, -34.0590873126, -22.513246978, -6.5668911915, 4.1388007503, 14.8444926921, 25.5501846339, 24.2165640106, 22.8829433873, 21.549322764, 20.2157021408, 18.8820815175, 17.5484608942),
                    Cz = c(-22.6859385301, -22.9537943081, -22.9804162663, -22.7655188671, -21.7627951018, -20.7600713366, -19.7573475713, -18.3485315288, -16.9397154862, -15.5308994437, -14.1220834012, -12.7132673586, -11.3044513161))
  out <- fix_stim_artifact(ref_eeg(c(1001L, 1004L)), events = 1)
  expect_equal(ref_window(out), unname(expected), tolerance = 1e-9)
})

test_that("matches the reference: only the events with the chosen code", {
  expected <- rbind(Fz = c(63.6507130712, 64.5084679766, 58.664130543, 46.9100578781, 33.2776427007, 19.6452275233, 6.0128123458, -7.6196028316, -21.252018009, -34.8844331864, -32.265626414, -23.2388112149, -9.0243712264),
                    Cz = c(-1.4337824552, -3.6249902495, -5.8267387267, -8.0124498462, -9.9046276772, -11.7968055082, -13.6889833392, -15.5811611703, -17.4733390013, -19.3655168323, -20.7501562284, -21.9386570269, -22.9166891265))
  x   <- ref_eeg(c(1001L, 1501L), types = c(1L, 2L))
  out <- fix_stim_artifact(x, events = 2)
  expect_equal(ref_window(out, 1498L, 1510L), unname(expected), tolerance = 1e-9)
  # the event with code 1 was left alone
  expect_equal(ref_window(out, 998L, 1010L), ref_window(x, 998L, 1010L))
})


# ============================================================================
#   3. NOTHING ELSE CHANGES
# ============================================================================

test_that("only the data inside the windows of the chosen channels changes, and the input is untouched", {
  ann <- data.frame(onset = 1.5, duration = 0.2, description = "BAD_x",
                    channel = NA_character_, stringsAsFactors = FALSE)
  x <- stim_eeg(c(1001L, 1601L), annotations = ann, bads = "Cz",
                montage = create_montage(c("Fz", "Cz")),
                metadata = list(subject = "S01"), reference = "Common Average")
  copy <- x
  out  <- fix_stim_artifact(x, events = 1)

  expect_identical(x, copy)                                   # the input is untouched
  for (f in c("events", "annotations", "bads", "montage", "metadata", "reference",
              "channels", "channel_types", "times", "sampling_rate")) {
    expect_identical(out[[f]], x[[f]], info = f)
  }
  expect_equal(nrow(out$annotations), 1L)                     # no annotation is added
  inside <- c(1002:1006, 1602:1606)
  for (ch in c("Fz", "Cz")) {
    expect_equal(changed(out, x, ch), inside, info = ch)      # bad channel Cz is repaired too
  }
  expect_identical(out$data["VEOG", ], x$data["VEOG", ])      # eye and trigger channels
  expect_identical(out$data["Status", ], x$data["Status", ])  # are left alone
  expect_identical(dim(out$data), dim(x$data))
  expect_identical(dimnames(out$data), dimnames(x$data))
})

test_that("one line is added to the history and says what was done", {
  x <- stim_eeg(c(1001L, 1601L, 1801L))
  x$preprocessing_history <- list("An earlier step")
  out <- fix_stim_artifact(x)
  expect_equal(length(out$preprocessing_history), 2L)
  expect_equal(out$preprocessing_history[[1]], "An earlier step")
  expect_equal(out$preprocessing_history[[2]],
               "Stimulation artifact fixed: 3 event(s), window 0 to 0.01 s, linear, 2 channel(s)")

  h <- fix_stim_artifact(x, events = 1, tmin = -0.002, tmax = 0.012, mode = "constant",
                         baseline = c(-0.02, -0.005))$preprocessing_history[[2]]
  expect_equal(h, paste0("Stimulation artifact fixed: 3 event(s) (code 1), window -0.002 to 0.012 s, ",
                         "constant (baseline -0.02 to -0.005 s), 2 channel(s)"))
})


# ============================================================================
#   4. WHICH CHANNELS
# ============================================================================

test_that("by default only the channels typed eeg are fixed, bad ones included", {
  x   <- stim_eeg(1001L, bads = "Cz")
  out <- fix_stim_artifact(x, events = 1)
  expect_equal(changed(out, x, "Fz"), 1002:1006)
  expect_equal(changed(out, x, "Cz"), 1002:1006)
  expect_equal(changed(out, x, "VEOG"), integer(0))
  expect_equal(changed(out, x, "Status"), integer(0))
})

test_that("picks: names or numbers, and eye channels only when asked for", {
  x <- stim_eeg(1001L)
  by_name   <- fix_stim_artifact(x, events = 1, picks = c("Fz", "VEOG"))
  by_number <- fix_stim_artifact(x, events = 1, picks = c(1, 3))
  expect_equal(by_name$data, by_number$data)
  expect_equal(changed(by_name, x, "Fz"), 1002:1006)
  expect_equal(changed(by_name, x, "VEOG"), 1002:1006)
  expect_equal(changed(by_name, x, "Cz"), integer(0))        # not picked
})

test_that("picks: wrong channels and a recording without eeg channels are refused", {
  x <- stim_eeg(1001L)
  expect_error(fix_stim_artifact(x, events = 1, picks = "Nope"), "not found")
  expect_error(fix_stim_artifact(x, events = 1, picks = c("Fz", "Fz")), "more than once")
  expect_error(fix_stim_artifact(x, events = 1, picks = character(0)), "empty")
  expect_error(fix_stim_artifact(x, events = 1, picks = 9), "valid channel indices")

  no_eeg <- x
  no_eeg$channel_types <- rep("misc", 4)
  expect_error(fix_stim_artifact(no_eeg, events = 1), "no channel is typed")
  expect_equal(changed(fix_stim_artifact(no_eeg, events = 1, picks = "Fz"), x, "Fz"), 1002:1006)
})


# ============================================================================
#   5. WHICH EVENTS
# ============================================================================

test_that("\"all\" uses the experimental triggers only (lower 16 bits, not 0), the same rule as epoch_eeg()", {
  # codes: 1, 65536 (only a system bit: lower 16 bits are 0), 1048577 (system bit + code 1)
  x   <- stim_eeg(c(1001L, 1301L, 1601L), types = c(1L, 65536L, 1048577L))
  out <- fix_stim_artifact(x)
  expect_equal(changed(out, x), c(1002:1006, 1602:1606))
  expect_match(out$preprocessing_history[[1]], "2 event\\(s\\), window")
})

test_that("events = a code uses only the events with that code (also when the upper bits are set)", {
  x   <- stim_eeg(c(1001L, 1301L, 1601L), types = c(9L, 2L, 1048585L))   # 1048585 = system bit + 9
  out <- fix_stim_artifact(x, events = 9)
  expect_equal(changed(out, x), c(1002:1006, 1602:1606))
  out2 <- fix_stim_artifact(x, events = c(2, 9))
  expect_equal(changed(out2, x), c(1002:1006, 1302:1306, 1602:1606))
})

test_that("trigger codes above 255 are kept whole: all 16 bits count", {
  # 256 would look like 0 (a system trigger) if only 8 bits were used
  x <- stim_eeg(c(1001L, 1301L, 1601L), types = c(300L, 256L, 65535L))
  expect_equal(changed(fix_stim_artifact(x), x), c(1002:1006, 1302:1306, 1602:1606))
  expect_equal(changed(fix_stim_artifact(x, events = 300), x), 1002:1006)
  expect_equal(changed(fix_stim_artifact(x, events = 65535), x), 1602:1606)
})

test_that("events: codes that do not exist, and wrong values", {
  x <- stim_eeg(c(1001L, 1301L), types = c(9L, 2L))
  expect_error(fix_stim_artifact(x, events = 5), "no event has the code\\(s\\) 5")
  expect_warning(out <- fix_stim_artifact(x, events = c(9, 5)), "no event has the code\\(s\\) 5")
  expect_equal(changed(out, x), 1002:1006)
  expect_error(fix_stim_artifact(x, events = "nine"), "must be \"all\" or a vector of event codes")
  expect_error(fix_stim_artifact(x, events = c(9, NA)), "must be \"all\" or a vector of event codes")
})

test_that("events: none at all, or only system triggers", {
  x <- stim_eeg(1001L)
  none <- x
  none$events <- none$events[0, ]
  expect_error(fix_stim_artifact(none), "has no events")
  system_only <- stim_eeg(1001L, types = 65536L)
  expect_error(fix_stim_artifact(system_only), "no experimental trigger")
})

test_that("events given twice or out of order give the same result as the plain table", {
  plain <- fix_stim_artifact(stim_eeg(c(1001L, 1601L)), events = 1)
  messy <- fix_stim_artifact(stim_eeg(c(1601L, 1001L, 1001L)), events = 1)
  expect_equal(messy$data, plain$data)
  expect_match(messy$preprocessing_history[[1]], "2 event\\(s\\)")
})


# ============================================================================
#   6. EPOCHS
# ============================================================================

epoch_rec <- function() {
  stim_eeg(c(1001L, 1501L, 3001L), data = stim_data(4000L))
}

make_epochs <- function(rec) {
  epoch_eeg(rec, events = 1, tmin = -0.1, tmax = 0.3, baseline = NULL, verbose = FALSE)
}

test_that("epochs: the window is counted from time zero of every epoch", {
  ep  <- make_epochs(epoch_rec())
  out <- fix_stim_artifact(ep)
  expect_s3_class(out, "eeg_epochs")
  expect_equal(dim(out$data), c(4L, 206L, 3L))
  # -0.1 s at 512 Hz: 51 samples before the event, so time zero is column 52
  cols <- which(apply(out$data[1, , ] != ep$data[1, , ], 1, any))
  expect_equal(unname(cols), 53:57)
  expect_equal(out$data[3, , ], ep$data[3, , ])               # VEOG not picked
})

test_that("epochs: fixing the epochs equals fixing the continuous recording and then epoching", {
  rec <- epoch_rec()
  cases <- list(list(),
                list(tmin = -0.002, tmax = 0.012),
                list(mode = "window", tmax = 0.012),
                list(mode = "constant", baseline = c(-0.02, -0.005)))
  for (k in seq_along(cases)) {
    args <- cases[[k]]
    a <- make_epochs(do.call(fix_stim_artifact, c(list(rec, events = 1), args)))
    b <- do.call(fix_stim_artifact, c(list(make_epochs(rec)), args))
    expect_equal(a$data, b$data, tolerance = 1e-12, info = paste("case", k))
  }
})

test_that("epochs: matches the reference (epochs -0.1 to 0.3 s, Fz, columns 50 to 62)", {
  rec <- ref_eeg(c(1001L, 1501L))
  ep  <- epoch_eeg(rec, events = 1, tmin = -0.1, tmax = 0.3, baseline = NULL, verbose = FALSE)
  one <- function(...) unname(fix_stim_artifact(ep, ...)$data[1, 50:62, 1:2])

  lin <- cbind(c(-34.0590873126, -22.513246978, -6.5668911915, 4.1388007503, 14.8444926921, 25.5501846339, 36.2558765756, 46.9615685174, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942, -0.824619568),
               c(64.5084679766, 58.664130543, 46.9100578781, 33.2776427007, 19.6452275233, 6.0128123458, -7.6196028316, -21.252018009, -34.8844331864, -32.265626414, -23.2388112149, -9.0243712264, 8.455189235))
  expect_equal(one(), lin, tolerance = 1e-9)

  win <- cbind(c(-34.0590873126, -22.513246978, -6.5668911915, 2.905767263, 0.0, 0.0, 0.0, 14.9934490183, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942, -0.824619568),
               c(64.5084679766, 58.664130543, 46.9100578781, 7.7096155073, 0.0, 0.0, 0.0, -7.6855354434, -34.8844331864, -32.265626414, -23.2388112149, -9.0243712264, 8.455189235))
  expect_equal(one(mode = "window", tmax = 0.012), win, tolerance = 1e-9)

  con <- cbind(c(-34.0590873126, -22.513246978, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, -3.1900786844, 57.6672604592, 48.9186589309, 34.9134271167, 17.5484608942, -0.824619568),
               c(64.5084679766, 58.664130543, 16.4567862675, 16.4567862675, 16.4567862675, 16.4567862675, 16.4567862675, 16.4567862675, -34.8844331864, -32.265626414, -23.2388112149, -9.0243712264, 8.455189235))
  expect_equal(one(mode = "constant", baseline = c(-0.02, -0.005)), con, tolerance = 1e-9)
})

test_that("epochs: everything except the data stays as it was, and one line goes into the history", {
  ep  <- make_epochs(epoch_rec())
  out <- fix_stim_artifact(ep)
  for (f in setdiff(names(ep), c("data", "preprocessing_history"))) {
    expect_identical(out[[f]], ep[[f]], info = f)
  }
  expect_equal(length(out$preprocessing_history), length(ep$preprocessing_history) + 1L)
  expect_equal(out$preprocessing_history[[length(out$preprocessing_history)]],
               "Stimulation artifact fixed: 3 epoch(s), window 0 to 0.01 s, linear, 2 channel(s)")
})

test_that("epochs: 'events' must stay \"all\", the data must be loaded, and the window must fit", {
  ep <- make_epochs(epoch_rec())
  expect_error(fix_stim_artifact(ep, events = 1), "every epoch is fixed")
  no_data <- ep
  no_data$data <- NULL
  expect_error(fix_stim_artifact(no_data), "not loaded")
  no_times <- ep
  no_times$times <- no_times$times[-1]
  expect_error(fix_stim_artifact(no_times), "one time .* for each sample")
  expect_error(fix_stim_artifact(ep, tmax = 0.5), "does not fit inside the epochs")
  expect_error(fix_stim_artifact(ep, tmin = -0.2), "does not fit inside the epochs")
  expect_error(fix_stim_artifact(ep, mode = "constant", baseline = c(-0.5, -0.4)),
               "baseline .* does not fit inside the epochs")
})


# ============================================================================
#   7. WRONG INPUT
# ============================================================================

test_that("the arguments are checked", {
  x <- stim_eeg(1001L)
  expect_error(fix_stim_artifact(list(a = 1)), "object of class 'eeg'")
  expect_error(fix_stim_artifact(x, mode = "cubic"), "should be one of")
  expect_error(fix_stim_artifact(x, tmin = "0"), "'tmin' must be one number")
  expect_error(fix_stim_artifact(x, tmax = c(0.01, 0.02)), "'tmax' must be one number")
  expect_error(fix_stim_artifact(x, tmin = NA_real_), "'tmin' must be one number")
  expect_error(fix_stim_artifact(x, tmin = 0.01, tmax = 0.01), "must be less than tmax")
  expect_error(fix_stim_artifact(x, tmin = 0.02, tmax = 0.01), "must be less than tmax")
  expect_error(fix_stim_artifact(x, mode = "constant"), "needs 'baseline'")
  expect_error(fix_stim_artifact(x, mode = "constant", baseline = c(0, -0.01)), "start < end")
  expect_error(fix_stim_artifact(x, mode = "constant", baseline = -0.01), "start < end")
  expect_warning(fix_stim_artifact(x, events = 1, baseline = c(-0.02, -0.005)),
                 "only used with mode = \"constant\"")
})

test_that("data of the wrong kind are refused", {
  x <- stim_eeg(1001L)
  text <- x
  text$data <- matrix("a", 4, 2000)
  expect_error(fix_stim_artifact(text, events = 1), "numeric matrix")
  one_row <- x
  one_row$data <- x$data["Fz", ]
  expect_error(fix_stim_artifact(one_row, events = 1), "numeric matrix")

  ep <- make_epochs(epoch_rec())
  flat <- ep
  flat$data <- ep$data[, , 1]
  expect_error(fix_stim_artifact(flat), "numeric 3D array")
})

test_that("a window that is too short for the fill is refused, with the sample count", {
  x100 <- stim_eeg(101L, fs = 100)
  # the default 0.01 s at 100 Hz is 1 sample wide: a line would change nothing
  expect_error(fix_stim_artifact(x100, events = 1), "only 1 sample\\(s\\) wide at 100 Hz")
  expect_error(fix_stim_artifact(stim_eeg(1001L), events = 1, tmax = 0.005, mode = "window"),
               "only 3 sample\\(s\\) wide at 512 Hz, but mode = \"window\" needs at least 4")
  expect_error(fix_stim_artifact(stim_eeg(1001L), events = 1, tmin = 0.001, tmax = 0.0015,
                                 mode = "constant", baseline = c(-0.02, -0.005)),
               "only 0 sample\\(s\\) wide")
  # 200 Hz: 2 samples, one of them changes
  x200 <- stim_eeg(1001L, fs = 200)
  expect_equal(changed(fix_stim_artifact(x200, events = 1), x200), 1002L)
  # a baseline that holds no sample
  expect_error(fix_stim_artifact(stim_eeg(1001L), events = 1, mode = "constant",
                                 baseline = c(-0.0001, -0.00005)), "holds no sample")
})

test_that("a window that falls off the recording is refused and names the event", {
  expect_error(fix_stim_artifact(stim_eeg(1999L), events = 1),
               "event at sample 1999 does not fit inside the recording \\(samples 1 to 2000\\)")
  expect_error(fix_stim_artifact(stim_eeg(c(1001L, 1998L, 1999L)), events = 1),
               "event at sample 1998 .* 1 more event\\(s\\)")
  expect_error(fix_stim_artifact(stim_eeg(5L), events = 1, tmin = -0.02, tmax = 0.01),
               "event at sample 5 does not fit")
  expect_error(fix_stim_artifact(stim_eeg(5L), events = 1, mode = "constant",
                                 baseline = c(-0.02, -0.005)),
               "baseline .* event at sample 5 does not fit")
  # an event right at the last possible position is fine
  ok <- stim_eeg(1994L)
  expect_equal(changed(fix_stim_artifact(ok, events = 1), ok), 1995:1999)
})

test_that("samples that cannot be used as line ends or baseline are refused", {
  d <- stim_data()
  d["Fz", 1001] <- NA
  expect_error(fix_stim_artifact(stim_eeg(1001L, data = d), events = 1),
               "not a finite number")
  d <- stim_data()
  d["Fz", 1007] <- Inf
  expect_error(fix_stim_artifact(stim_eeg(1001L, data = d), events = 1),
               "not a finite number")
  # a missing value INSIDE the window is simply overwritten
  d <- stim_data()
  d["Fz", 1003] <- NaN
  out <- fix_stim_artifact(stim_eeg(1001L, data = d), events = 1)
  expect_false(anyNA(out$data["Fz", 1001:1007]))
  # ... but not in the baseline of the constant fill
  d <- stim_data()
  d["Fz", 995] <- NA
  expect_error(fix_stim_artifact(stim_eeg(1001L, data = d), events = 1, mode = "constant",
                                 baseline = c(-0.02, -0.005)), "baseline sample")
})


# ============================================================================
#   8. WHY IT HAS TO COME FIRST
# ============================================================================

# noise (about 1 microvolt) with a 5000 microvolt spike of 5 samples after each event
spiky_recording <- function(onsets, n, fs = 512) {
  set.seed(7)
  d <- rbind(Fz = rnorm(n), Cz = rnorm(n), VEOG = rnorm(n), Status = rep(0, n))
  for (o in onsets) d[c("Fz", "Cz"), (o + 1L):(o + 5L)] <- rep(c(5000, 4000, 3000, 2000, 1500), each = 2)
  stim_eeg(onsets, fs = fs, data = d)
}

test_that("the spike trips the amplitude limit of epoch_eeg() before the repair and not after", {
  rec <- spiky_recording(c(1001L, 2001L, 3001L), 4000L)
  epoch_it <- function(r) epoch_eeg(r, events = 1, tmin = -0.1, tmax = 0.3, baseline = NULL,
                                    reject_threshold = 200, verbose = FALSE)
  expect_error(epoch_it(rec), "rejected")
  fixed <- epoch_it(fix_stim_artifact(rec, events = 1))
  expect_equal(dim(fixed$data)[3], 3L)
})

test_that("filtering first smears the spike over a long time; repairing first does not", {
  n   <- 10240L
  rec <- spiky_recording(5121L, n)
  band <- function(r) eeg_bandpass(r, l_freq = 1, h_freq = 40, verbose = FALSE)
  # 0.1 to 0.3 s after the pulse: long after the spike itself (5 samples = 10 ms)
  far <- (5121L + 51L):(5121L + 154L)
  smeared <- max(abs(band(rec)$data["Fz", far]))
  clean   <- max(abs(band(fix_stim_artifact(rec, events = 1))$data["Fz", far]))
  expect_gt(smeared, 10 * clean)
})
