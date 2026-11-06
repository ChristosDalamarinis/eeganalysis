# ============================================================================
#                        STIMULATION-ARTIFACT REMOVAL
#            (replace the short spike that a stimulator leaves behind)
# ============================================================================
#
# When a stimulator fires (a TMS pulse, an electrical shock) the amplifier
# records a huge spike for a few milliseconds, far bigger than any brain
# signal. If it stays in the data, the next steps break: a filter smears it
# into a long wobble, ICA spends components on it, an average gets a bump at
# time zero. This file has ONE exported function, fix_stim_artifact(), that
# overwrites that short stretch with something harmless.
#
# HOW IT KNOWS WHERE THE STIMULATION WAS
#   It does not look for spikes in the data. It reads the events table
#   (eeg$events). The reader builds that table from the trigger (Status)
#   channel, so the experiment must send a trigger at the moment of the
#   stimulation; the user then says which trigger code(s) mean "stimulation".
#
# THE WINDOW AND THE THREE FILLS
#   For every stimulation event a window of a few samples (tmin to tmax, in
#   seconds, counted from the event) is overwritten. The fill is
#     "linear"    a straight line between the first and the last sample of
#                 the window (both are real data and are kept),
#     "window"    the window multiplied by 1, 0.25, 0, ..., 0, 0.25, 1,
#     "constant"  the average of a quiet stretch (the baseline) per channel.
#
# WHAT IT CHANGES
#   Only the data of the chosen channels, inside the windows. Events,
#   annotations, bads, montage, times and everything else stay as they were.
#   No "BAD_" annotation is added for the repaired stretches: epoch_eeg() throws
#   away every epoch that touches a "BAD_" annotation, so every trial would be
#   lost. One line in preprocessing_history records what was done.
#
# THE MATHS LIVES IN ONE PLACE
#   Continuous recordings and epochs both go through .stim_block(), which
#   works out the new samples for one window of one recording or one epoch;
#   the two paths of fix_stim_artifact() only cut the pieces out and put the
#   result back.
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================


# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

#' Turn the window and the baseline from seconds into samples (internal)
#'
#' The first and the last sample of the window are
#' \code{ceiling(rate * seconds)} away from the event. The samples that are
#' overwritten are the \code{n_win} ones from the first up to, but not
#' including, the last. Stops with a clear message when the window holds too
#' few samples for the chosen fill: \code{"linear"} needs 2 (it keeps the first
#' sample, so one sample would change nothing), \code{"window"} needs 4 and
#' \code{"constant"} needs 1.
#'
#' @param fs       Sampling rate in Hz.
#' @param tmin,tmax Window start and end in seconds, \code{tmin < tmax}.
#' @param mode     \code{"linear"}, \code{"window"} or \code{"constant"}.
#' @param baseline \code{c(start, end)} in seconds (only used for
#'   \code{"constant"}, otherwise ignored).
#' @return A list with \code{s_start}, \code{s_end} (sample offsets from the
#'   event), \code{n_win} (\code{s_end - s_start}), and \code{b_start},
#'   \code{b_end} (baseline offsets, \code{NA} unless \code{"constant"}; the
#'   baseline samples are \code{b_start} up to \code{b_end - 1}).
#' @keywords internal
.stim_window_samples <- function(fs, tmin, tmax, mode, baseline) {

  s_start <- ceiling(fs * tmin)
  s_end   <- ceiling(fs * tmax)
  n_win   <- s_end - s_start

  need <- c(linear = 2, window = 4, constant = 1)[[mode]]
  if (n_win < need) {
    stop("ERROR: the window (", format(tmin), " to ", format(tmax), " s) is ",
         "only ", n_win, " sample(s) wide at ", format(fs), " Hz, but mode = ",
         "\"", mode, "\" needs at least ", need, ". Make the window longer or ",
         "start it earlier (see 'tmin' and 'tmax').", call. = FALSE)
  }

  b_start <- b_end <- NA_real_
  if (mode == "constant") {
    b_start <- ceiling(fs * baseline[1])
    b_end   <- ceiling(fs * baseline[2])
    if (b_end - b_start < 1) {
      stop("ERROR: the baseline (", format(baseline[1]), " to ",
           format(baseline[2]), " s) holds no sample at ", format(fs),
           " Hz. Make it longer.", call. = FALSE)
    }
  }

  list(s_start = s_start, s_end = s_end, n_win = as.integer(n_win),
       b_start = b_start, b_end = b_end)
}


#' Which events are stimulations? (internal)
#'
#' Uses the same rule as \code{\link{epoch_eeg}}: only experimental triggers
#' count, that is events whose code, looked at through the lower 16 bits, is
#' not 0 (the upper bits of the trigger channel belong to the system).
#' \code{"all"} takes every one of them, a vector of codes takes only the
#' events with those codes.
#'
#' @param eeg    An \code{eeg} object.
#' @param events \code{"all"} or a vector of event codes.
#' @return Sorted, unique sample numbers (1-based) of the chosen events.
#' @keywords internal
.stim_event_samples <- function(eeg, events) {

  ev <- eeg$events
  if (!is.data.frame(ev) || nrow(ev) == 0L ||
      !all(c("onset", "type") %in% names(ev))) {
    stop("ERROR: the recording has no events. The times of the stimulation ",
         "come from the trigger events (see inspect_triggers()).",
         call. = FALSE)
  }

  code         <- suppressWarnings(bitwAnd(as.integer(ev$type), 65535L))
  experimental <- !is.na(code) & code != 0L

  if (identical(events, "all")) {
    keep <- experimental
    if (!any(keep)) {
      stop("ERROR: no experimental trigger found in the events (all of them ",
           "are system triggers).", call. = FALSE)
    }
  } else {
    if (!is.numeric(events) || length(events) < 1L || anyNA(events)) {
      stop("ERROR: 'events' must be \"all\" or a vector of event codes ",
           "(numbers).", call. = FALSE)
    }
    keep   <- experimental & code %in% events
    absent <- setdiff(events, code[experimental])
    if (!any(keep)) {
      stop("ERROR: no event has the code(s) ", paste(events, collapse = ", "),
           " (inspect_triggers() lists the codes in the recording).",
           call. = FALSE)
    }
    if (length(absent) > 0L) {
      warning("no event has the code(s) ", paste(absent, collapse = ", "),
              "; they are skipped.", call. = FALSE)
    }
  }

  sort(unique(as.integer(ev$onset[keep])))
}


#' The new samples for one window (internal)
#'
#' The one place where the fills are worked out. It gets the samples of one
#' window of one recording (or one epoch), all picked channels at once, and
#' returns what the first \code{n_win} of them become; the last sample of the
#' window is only read, never changed.
#'
#' @param seg     Numeric matrix, picked channels x (\code{n_win} + 1): the
#'   samples from the first up to and including the last sample of the window.
#' @param base    Numeric matrix, picked channels x baseline samples (only for
#'   \code{"constant"}, otherwise \code{NULL}).
#' @param mode    \code{"linear"}, \code{"window"} or \code{"constant"}.
#' @param n_win   Number of samples that are overwritten.
#' @param where   Text that says which event or epoch this is, for messages.
#' @return Numeric matrix, picked channels x \code{n_win}.
#' @keywords internal
.stim_block <- function(seg, base, mode, n_win, where) {

  if (mode == "linear") {
    ends <- seg[, c(1L, n_win + 1L), drop = FALSE]
    if (!all(is.finite(ends))) {
      stop("ERROR: the sample at the start or at the end of the window (",
           where, ") is not a finite number (NA, NaN or Inf), so no line can ",
           "be drawn.", call. = FALSE)
    }
    slope <- (ends[, 2L] - ends[, 1L]) / n_win
    return(outer(unname(slope), 0:(n_win - 1L)) + unname(ends[, 1L]))
  }

  if (mode == "window") {
    hann4  <- c(0, 0.75, 0.75, 0)
    weight <- 1 - c(hann4[1:2], rep(1, n_win - 4L), hann4[3:4])
    return(unname(seg[, seq_len(n_win), drop = FALSE]) *
             rep(weight, each = nrow(seg)))
  }

  if (!all(is.finite(base))) {
    stop("ERROR: a baseline sample (", where, ") is not a finite number ",
         "(NA, NaN or Inf), so no average can be taken.", call. = FALSE)
  }
  matrix(rowMeans(base), nrow = nrow(seg), ncol = n_win)
}


# ----------------------------------------------------------------------------
# fix_stim_artifact()
# ----------------------------------------------------------------------------
#' Replace the Short Spike a Stimulator Leaves in the Recording
#'
#' When a stimulator fires (a magnetic pulse, an electrical shock) the amplifier
#' records a huge spike for a few milliseconds, far larger than any brain
#' signal. This function overwrites that short stretch, called the
#' \emph{window}, after every stimulation event with a harmless fill, so that
#' later steps (filtering, ICA, averaging) are not damaged by the spike.
#'
#' \strong{How the stimulation is found.} The function does not look for
#' spikes in the data. It reads the events of the recording
#' (\code{eeg$events}), which the reader builds from the trigger channel. So
#' the experiment must send a trigger at the moment of the stimulation, and you
#' tell the function which code(s) mean "stimulation" with \code{events}. If
#' the recording has no trigger for the stimulation, the function cannot be
#' used. \code{events = "all"} uses every experimental trigger, which also
#' includes triggers that are not stimulations (a cue, a button press): the
#' same few milliseconds of real signal after each of them would be replaced by
#' a line. Give the stimulation code(s) instead; \code{\link{inspect_triggers}}
#' lists the codes.
#'
#' \strong{The window.} It runs from \code{tmin} to \code{tmax} seconds after
#' the event. In samples it starts \code{ceiling(rate * tmin)} and ends
#' \code{ceiling(rate * tmax)} samples after the event sample. The sample at
#' the start and the sample at the end are real data. The samples from the
#' start up to, but not including, the end are overwritten, in one of three
#' ways (\code{mode}):
#' \describe{
#'   \item{\code{"linear"} (default)}{A straight line from the sample at the
#'     start of the window to the sample at its end. Both of those are kept, so
#'     a window of \code{n} samples changes \code{n - 1} of them (at 512 Hz,
#'     0 to 0.01 s: 6 samples wide, 5 change). It needs at least 2 samples.}
#'   \item{\code{"window"}}{Each sample is multiplied by
#'     \code{1, 0.25, 0, ..., 0, 0.25, 1}, so the middle of the window becomes
#'     zero and the two edge samples stay as they are. It needs at least 4
#'     samples. Zero is only a harmless fill for data that are centred around
#'     zero (for example after filtering or re-referencing).}
#'   \item{\code{"constant"}}{Every sample of the window becomes the average of
#'     the channel over a quiet stretch, \code{baseline}, given in seconds
#'     relative to the event (the end of the baseline is not included).}
#' }
#' Two events that are closer together than the window are handled one after
#' the other, so the second line starts from the values the first one wrote.
#'
#' \strong{Where the spike starts matters.} The two ends of the line are real
#' samples. If the spike already sits on the first sample of the window (the
#' trigger sample), the line starts from a spike value and the repair fails.
#' Start the window a little earlier then, for example \code{tmin = -0.002}
#' (one sample at 512 Hz). If the trigger comes a fixed time before the
#' stimulation, move the window later; if that delay changes from trial to
#' trial, no fixed window fits. At low sampling rates the default window
#' (0.01 s) is only a few samples wide: at 100 Hz it is 1 sample and the
#' function stops with a message.
#'
#' \strong{When to use it.} On the continuous recording, \emph{before}
#' filtering, ICA, \code{\link{downsample}} and epoching. A filter turns a
#' spike of a few milliseconds into a wobble that lasts much longer and cannot
#' be cut out afterwards, and a spike of thousands of microvolts trips the
#' amplitude limit of \code{\link{epoch_eeg}}, which then drops the trial.
#'
#' \strong{Epochs.} For an \code{eeg_epochs} object the window is counted from
#' time zero of every epoch, and every epoch is fixed. That is only right when
#' the epochs are cut around the stimulation. If they are cut around another
#' event (a cue) and the stimulation comes later, fix the continuous recording
#' first (with \code{events = <stimulation code>}), then epoch.
#'
#' \strong{What it does not do.} Only the short spike is replaced. The slower
#' tail and the muscle twitch that a magnetic pulse also leaves are not
#' removed. No \code{"BAD_"} annotation is added for the repaired stretches
#' (\code{epoch_eeg} would reject every trial that touches one); the step is
#' recorded in \code{preprocessing_history}.
#'
#' @param eeg An object of class \code{eeg} (continuous) or \code{eeg_epochs}.
#' @param events Which events are stimulations. \code{"all"} (default) uses
#'   every experimental trigger (the same rule as \code{\link{epoch_eeg}}); a
#'   number or a vector of numbers uses only the events with those codes. For
#'   \code{eeg_epochs} every epoch is fixed and \code{events} must stay
#'   \code{"all"}.
#' @param tmin Start of the window in seconds, counted from the event (from
#'   time zero for epochs). Default 0; may be negative. Must be smaller than
#'   \code{tmax}.
#' @param tmax End of the window in seconds. Default 0.01.
#' @param mode \code{"linear"} (default), \code{"window"} or
#'   \code{"constant"}; see Details.
#' @param baseline \code{c(start, end)} in seconds, relative to the event.
#'   Needed for \code{mode = "constant"}; ignored (with a warning) otherwise.
#' @param picks Channels to fix: names or numbers. \code{NULL} (default) fixes
#'   every channel typed \code{"eeg"} (bad channels included); eye and trigger
#'   channels are left alone unless you name them.
#'
#' @return A new object of the same class with the windows overwritten and one
#'   line added to \code{preprocessing_history}. Only \code{data} changes. The
#'   input is left untouched.
#'
#' @examples
#' \dontrun{
#'   # trigger code 9 marks each pulse: repair the 10 ms after it
#'   eeg <- fix_stim_artifact(eeg, events = 9)
#'
#'   # the spike already starts on the trigger sample: start the window earlier
#'   eeg <- fix_stim_artifact(eeg, events = 9, tmin = -0.002, tmax = 0.012)
#'
#'   # the average of a quiet stretch before the pulse instead of a line
#'   eeg <- fix_stim_artifact(eeg, events = 9, mode = "constant",
#'                            baseline = c(-0.02, -0.005))
#' }
#'
#' @seealso \code{\link{inspect_triggers}}, \code{\link{epoch_eeg}},
#'   \code{\link{eeg_bandpass}}
#' @export
fix_stim_artifact <- function(eeg, events = "all", tmin = 0, tmax = 0.01,
                              mode = c("linear", "window", "constant"),
                              baseline = NULL, picks = NULL) {

  # ========== INPUT VALIDATION ==========
  kind <- .channels_kind(eeg, "fix_stim_artifact")
  mode <- match.arg(mode)

  one_number <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!one_number(tmin)) {
    stop("ERROR: 'tmin' must be one number (seconds).", call. = FALSE)
  }
  if (!one_number(tmax)) {
    stop("ERROR: 'tmax' must be one number (seconds).", call. = FALSE)
  }
  if (tmin >= tmax) {
    stop("ERROR: tmin (", format(tmin), ") must be less than tmax (",
         format(tmax), ").", call. = FALSE)
  }
  if (mode == "constant") {
    if (is.null(baseline)) {
      stop("ERROR: mode = \"constant\" needs 'baseline', c(start, end) in ",
           "seconds.", call. = FALSE)
    }
    if (!is.numeric(baseline) || length(baseline) != 2L ||
        !all(is.finite(baseline)) || baseline[1] >= baseline[2]) {
      stop("ERROR: 'baseline' must be c(start, end) in seconds with ",
           "start < end.", call. = FALSE)
    }
  } else if (!is.null(baseline)) {
    warning("'baseline' is only used with mode = \"constant\"; it is ignored.",
            call. = FALSE)
  }

  fs <- eeg$sampling_rate
  w  <- .stim_window_samples(fs, tmin, tmax, mode, baseline)

  # ---- which channels ------------------------------------------------------
  if (is.null(picks)) {
    rows <- which(.types_of(eeg) == "eeg")
    if (length(rows) == 0L) {
      stop("ERROR: no channel is typed \"eeg\"; name the channels to fix ",
           "with 'picks'.", call. = FALSE)
    }
  } else {
    rows <- match(.eog_resolve_picks(picks, eeg, "picks"), eeg$channels)
  }

  data <- eeg$data

  if (kind == "eeg") {

    # ========== CONTINUOUS RECORDING ==========
    if (!is.matrix(data) || !is.numeric(data)) {
      stop("ERROR: eeg$data must be a numeric matrix (channels x time).",
           call. = FALSE)
    }
    onsets <- .stim_event_samples(eeg, events)
    n_t    <- ncol(data)
    first  <- onsets + w$s_start
    last   <- onsets + w$s_end

    off <- which(first < 1 | last > n_t)
    if (length(off) > 0L) {
      stop("ERROR: the window (", format(tmin), " to ", format(tmax), " s) of ",
           "the event at sample ", onsets[off[1L]], " does not fit inside the ",
           "recording (samples 1 to ", n_t, ")",
           if (length(off) > 1L) {
             paste0("; ", length(off) - 1L, " more event(s) are too close to ",
                    "an edge as well")
           },
           ". Leave such events out ('events') or shorten the window.",
           call. = FALSE)
    }
    if (mode == "constant") {
      b_first <- onsets + w$b_start
      b_last  <- onsets + w$b_end - 1
      off <- which(b_first < 1 | b_last > n_t)
      if (length(off) > 0L) {
        stop("ERROR: the baseline (", format(baseline[1]), " to ",
             format(baseline[2]), " s) of the event at sample ",
             onsets[off[1L]], " does not fit inside the recording (samples 1 ",
             "to ", n_t, ").", call. = FALSE)
      }
    }

    for (k in seq_along(onsets)) {
      seg  <- data[rows, first[k]:last[k], drop = FALSE]
      base <- if (mode == "constant") {
        data[rows, b_first[k]:b_last[k], drop = FALSE]
      }
      data[rows, first[k]:(last[k] - 1)] <-
        .stim_block(seg, base, mode, w$n_win,
                    paste0("the event at sample ", onsets[k]))
    }
    what <- paste0(length(onsets), " event(s)",
                   if (!identical(events, "all")) {
                     paste0(" (code ", paste(events, collapse = ", "), ")")
                   })

  } else {

    # ========== EPOCHS ==========
    if (!identical(events, "all")) {
      stop("ERROR: for epochs every epoch is fixed, so 'events' must stay ",
           "\"all\".", call. = FALSE)
    }
    if (is.null(data)) {
      stop("ERROR: epochs$data is not loaded - re-run epoch_eeg() with ",
           "preload = TRUE.", call. = FALSE)
    }
    if (!is.numeric(data) || length(dim(data)) != 3L) {
      stop("ERROR: epochs$data must be a numeric 3D array (channels x time x ",
           "epochs).", call. = FALSE)
    }
    n_t <- dim(data)[2L]
    n_e <- dim(data)[3L]
    if (!is.numeric(eeg$times) || length(eeg$times) != n_t) {
      stop("ERROR: epochs$times must hold one time (in seconds) for each ",
           "sample of the epochs.", call. = FALSE)
    }

    # the column of time zero: epoch_eeg() starts every epoch round(tmin * fs)
    # samples before the event
    zero  <- 1 - round(eeg$times[1L] * fs)
    first <- zero + w$s_start
    last  <- zero + w$s_end
    if (first < 1 || last > n_t) {
      stop("ERROR: the window (", format(tmin), " to ", format(tmax), " s, ",
           "counted from time zero) does not fit inside the epochs (",
           format(round(eeg$times[1L], 4)), " to ",
           format(round(eeg$times[n_t], 4)), " s).", call. = FALSE)
    }
    if (mode == "constant") {
      b_first <- zero + w$b_start
      b_last  <- zero + w$b_end - 1
      if (b_first < 1 || b_last > n_t) {
        stop("ERROR: the baseline (", format(baseline[1]), " to ",
             format(baseline[2]), " s, counted from time zero) does not fit ",
             "inside the epochs (", format(round(eeg$times[1L], 4)), " to ",
             format(round(eeg$times[n_t], 4)), " s).", call. = FALSE)
      }
    }

    for (e in seq_len(n_e)) {
      seg  <- matrix(data[rows, first:last, e], nrow = length(rows))
      base <- if (mode == "constant") {
        matrix(data[rows, b_first:b_last, e], nrow = length(rows))
      }
      data[rows, first:(last - 1), e] <-
        .stim_block(seg, base, mode, w$n_win, paste0("epoch ", e))
    }
    what <- paste0(n_e, " epoch(s)")
  }

  # ========== RETURN ==========
  out <- eeg
  out$data <- data
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Stimulation artifact fixed: ", what, ", window ",
                format(tmin), " to ", format(tmax), " s, ", mode,
                if (mode == "constant") {
                  paste0(" (baseline ", format(baseline[1]), " to ",
                         format(baseline[2]), " s)")
                }, ", ", length(rows), " channel(s)")))
  out
}
