# ============================================================================
#                         ADD EVENTS BY HAND
#          (write events into the events table of a recording)
# ============================================================================
#
# The events table (eeg$events) says "something happened at this sample". The
# readers fill it from the trigger channel. If the experiment never sent a
# trigger for something (the stimulation, a cue, a button press), or the
# trigger was lost, there is no row for it, and every step that needs events
# (epoch_eeg(), fix_stim_artifact(), ...) cannot see it. This file has ONE
# exported function, add_events(), that writes such rows from times the user
# gives.
#
# WHAT IT DOES AND WHAT IT DOES NOT
#   It does not look at the data. It cannot see when something happened: it
#   writes the time the user gives it, turned into the nearest sample. The
#   user must know the times from somewhere else (a log file of the
#   stimulator, a regular schedule), and they must be on the clock of the
#   recording: seconds from the start of the recording.
#
# THE ROWS
#   One row per event, with the same four columns the reader writes:
#     onset        the sample number (1-based, integer),
#     onset_time   (onset - 1) / sampling rate, in seconds,
#     type         the trigger code (integer, 1 to 65535),
#     description  text, by default "Trigger: <code>" like the reader.
#   The old events stay (with any extra column, for example "label", filled
#   with NA for the new rows) unless replace = TRUE. The table is sorted by
#   onset. Two events on one sample are refused, because epoch_eeg() keeps
#   only the first of them.
#
# THE ONE TRICKY PART
#   The empty events table that new_eeg() makes when it is given no events
#   stores "type" as text. Adding a number to it the plain way would turn the
#   whole column into text and break matching such as events = 9. The helper
#   .events_existing() turns "type" into integers first.
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================


# ----------------------------------------------------------------------------
# Internal helper
# ----------------------------------------------------------------------------

#' The old events table, ready to take new rows (internal)
#'
#' A missing table (\code{NULL}) becomes an empty table with the four columns
#' the reader writes. Otherwise the table must have the columns \code{onset},
#' \code{onset_time}, \code{type} and \code{description}; their classes are
#' made integer, numeric, integer and character (so the empty table that
#' \code{\link{new_eeg}} makes, whose \code{type} is text, takes numbers
#' without turning the column into text). Any other column, for example
#' \code{label}, is kept as it is.
#'
#' @param eeg A continuous \code{eeg} object.
#' @return A data frame. Stops with a clear message when the columns are
#'   missing or when \code{type} holds text that is not a number.
#' @keywords internal
.events_existing <- function(eeg) {

  ev <- eeg$events
  if (is.null(ev)) {
    return(data.frame(onset = integer(0), onset_time = numeric(0),
                      type = integer(0), description = character(0),
                      stringsAsFactors = FALSE))
  }

  needed <- c("onset", "onset_time", "type", "description")
  if (!is.data.frame(ev) || !all(needed %in% names(ev))) {
    stop("ERROR: eeg$events must be a data frame with the columns onset, ",
         "onset_time, type and description.", call. = FALSE)
  }

  type <- suppressWarnings(as.integer(ev$type))
  bad  <- is.na(type) & !is.na(ev$type)
  if (any(bad)) {
    shown <- unique(as.character(ev$type[bad]))
    stop("ERROR: eeg$events$type holds values that are not numbers (",
         paste(shown[seq_len(min(3L, length(shown)))], collapse = ", "),
         "), so new events cannot be added to it.", call. = FALSE)
  }

  ev$onset       <- as.integer(ev$onset)
  ev$onset_time  <- as.numeric(ev$onset_time)
  ev$type        <- type
  ev$description <- as.character(ev$description)
  ev
}


# ----------------------------------------------------------------------------
# add_events()
# ----------------------------------------------------------------------------
#' Add Events to a Recording by Hand
#'
#' Writes events, a time and a code each, into the events table of a
#' continuous recording. Use it when the experiment never sent a trigger for
#' something, or the trigger was lost: the stimulation, a cue, a button press,
#' or a regular schedule such as one event every 2 seconds for resting-state
#' data. The new events are used like any trigger the reader made, for example
#' by \code{\link{epoch_eeg}} and \code{\link{fix_stim_artifact}}.
#'
#' \strong{The function does not look at the data.} It cannot see when
#' something happened. It writes the time you give it, so you must know the
#' times from somewhere else (for example a log file of the stimulator) and
#' they must be on the clock of the recording, in seconds from its start. If
#' your log uses the time of day, subtract the time of day at which the
#' recording started. If the log and the recording use different clocks, match
#' them first; the function cannot detect a mismatch. Looking at the data
#' around the new events is a good check: if every spike is off by the same
#' amount, add that amount to all times and run it again.
#'
#' \strong{From seconds to a sample.} A time becomes a sample number with
#' \code{onset = round(time * rate) + 1} (sample 1 is time 0), and the stored
#' \code{onset_time} is \code{(onset - 1) / rate}, the time of that sample.
#' So a time is moved to the nearest sample, by at most half a sample (about 1
#' ms at 512 Hz). The sample must lie inside the recording.
#'
#' \strong{The rows.} Each new event gets the four columns the reader writes:
#' \code{onset}, \code{onset_time}, \code{type} (the code) and
#' \code{description}. The old events stay, with all their columns (a column
#' such as \code{label} is filled with \code{NA} for the new rows), unless you
#' use \code{replace = TRUE}. The table is sorted by \code{onset}. Two events
#' on the same sample are refused, because \code{\link{epoch_eeg}} would keep
#' only the first of them.
#'
#' @param eeg An object of class \code{eeg} (continuous data). An
#'   \code{eeg_epochs} object is refused: its events are fixed.
#' @param times Numeric vector: the times of the new events in seconds from
#'   the start of the recording. Must not be negative.
#' @param code The trigger code: whole numbers from 1 to 65535, one number for
#'   all events or one per time. (Codes whose lower 16 bits are 0 would be
#'   ignored by \code{\link{epoch_eeg}} and \code{\link{fix_stim_artifact}}.)
#' @param description Text for the \code{description} column: one for all
#'   events or one per time. \code{NULL} (default) writes
#'   \code{"Trigger: <code>"}, the same text the reader writes.
#' @param replace Logical. \code{FALSE} (default) keeps the old events and adds
#'   the new ones. \code{TRUE} throws the old events away and keeps only the
#'   new ones.
#'
#' @return A new \code{eeg} object. Only \code{events} changes, and one line
#'   is added to \code{preprocessing_history}. The input is left untouched.
#'
#' @examples
#' \dontrun{
#'   # pulse times from a stimulator log (seconds after the recording started)
#'   eeg <- add_events(eeg, times = c(12.5, 17.5, 22.5), code = 9)
#'   eeg <- fix_stim_artifact(eeg, events = 9)
#'
#'   # resting state: one event every 2 seconds, then equal windows
#'   eeg <- add_events(eeg, times = seq(0, 298, by = 2), code = 1)
#'   epochs <- epoch_eeg(eeg, events = 1, tmin = 0, tmax = 2, baseline = NULL)
#'
#'   # throw the old events away and keep only the new ones
#'   eeg <- add_events(eeg, times = c(10, 20), code = 1, replace = TRUE)
#' }
#'
#' @seealso \code{\link{inspect_triggers}}, \code{\link{epoch_eeg}},
#'   \code{\link{fix_stim_artifact}}, \code{\link{new_eeg}}
#' @export
add_events <- function(eeg, times, code, description = NULL, replace = FALSE) {

  # ========== INPUT VALIDATION ==========
  .channels_kind(eeg, "add_events", epochs_ok = FALSE)

  if (!is.numeric(times) || length(times) < 1L) {
    stop("ERROR: 'times' must be a non-empty numeric vector (seconds from ",
         "the start of the recording).", call. = FALSE)
  }
  times <- as.numeric(unname(times))
  if (!all(is.finite(times))) {
    stop("ERROR: 'times' must not hold NA, NaN or Inf.", call. = FALSE)
  }
  if (any(times < 0)) {
    stop("ERROR: 'times' must not be negative (first one: ",
         format(times[times < 0][1L]), " s).", call. = FALSE)
  }
  n_new <- length(times)

  if (!is.numeric(code) || length(code) < 1L || !all(is.finite(code)) ||
      any(code != round(code)) || any(code < 1) || any(code > 65535)) {
    stop("ERROR: 'code' must be whole numbers from 1 to 65535.",
         call. = FALSE)
  }
  if (!(length(code) %in% c(1L, n_new))) {
    stop("ERROR: 'code' must be one number or one number per time (",
         n_new, " time(s), ", length(code), " code(s)).", call. = FALSE)
  }
  code <- as.integer(code)

  if (is.null(description)) {
    description <- paste0("Trigger: ", code)
  } else if (!is.character(description) || anyNA(description) ||
             !(length(description) %in% c(1L, n_new))) {
    stop("ERROR: 'description' must be text: one for all events or one per ",
         "time.", call. = FALSE)
  }

  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("ERROR: 'replace' must be TRUE or FALSE.", call. = FALSE)
  }

  # ========== TIME TO SAMPLE ==========
  fs   <- eeg$sampling_rate
  n    <- ncol(eeg$data)
  samp <- round(times * fs)                  # samples after the first one
  outside <- which(samp > n - 1)
  if (length(outside) > 0L) {
    stop("ERROR: the time ", format(times[outside[1L]]), " s is outside the ",
         "recording (0 to ", format(round((n - 1) / fs, 4)), " s)",
         if (length(outside) > 1L) {
           paste0("; ", length(outside) - 1L, " more time(s) are outside as ",
                  "well")
         }, ".", call. = FALSE)
  }
  onset <- as.integer(samp) + 1L

  dup <- which(duplicated(onset))
  if (length(dup) > 0L) {
    stop("ERROR: two of the new events fall on the same sample (sample ",
         onset[dup[1L]], "). Give them different times.", call. = FALSE)
  }

  # ========== THE OLD EVENTS ==========
  n_old <- if (is.data.frame(eeg$events)) nrow(eeg$events) else 0L
  old   <- if (replace) NULL else .events_existing(eeg)
  if (!is.null(old)) {
    clash <- which(onset %in% old$onset)
    if (length(clash) > 0L) {
      stop("ERROR: an event already sits at sample ", onset[clash[1L]],
           " (time ", format((onset[clash[1L]] - 1) / fs), " s). Remove it ",
           "first or use replace = TRUE.", call. = FALSE)
    }
  }

  # ========== THE NEW TABLE ==========
  new <- data.frame(onset       = onset,
                    onset_time  = (onset - 1) / fs,
                    type        = rep_len(code, n_new),
                    description = rep_len(description, n_new),
                    stringsAsFactors = FALSE)
  if (is.null(old)) {
    events <- new
  } else {
    # a column the new rows do not have (for example "label") gets NA of the
    # same kind
    for (col in setdiff(names(old), names(new))) {
      new[[col]] <- old[[col]][rep(NA_integer_, n_new)]
    }
    events <- rbind(old, new[names(old)])
  }
  events <- events[order(events$onset), , drop = FALSE]
  rownames(events) <- NULL

  # ========== RETURN ==========
  codes <- paste(sort(unique(code)), collapse = ", ")
  note <- if (replace) {
    paste0("Events replaced: ", n_old, " old event(s) removed, ", n_new,
           " new event(s), code(s) ", codes)
  } else {
    paste0("Events added: ", n_new, " event(s), code(s) ", codes)
  }

  out <- eeg
  out$events <- events
  out$preprocessing_history <- c(out$preprocessing_history, list(note))
  out
}
