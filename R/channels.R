# ============================================================================
#              CHANNEL AND DATA SURGERY (choose, drop, rename, add,
#                      combine channels; trim and join recordings)
# ============================================================================
#
# Small tools for cutting a recording up and putting it back together.
#
#   drop_channels()     remove the channels you name
#   pick_channels()     keep the channels you name (in your order)
#   pick_types()        keep the channels of a kind (eeg, eog, ...)
#   rename_channels()   change channel names
#   add_channels()      put extra channels next to the existing ones
#   combine_channels()  merge groups of channels into one channel each
#   eeg_crop()          cut time off the start and the end
#   eeg_concatenate()   join recordings one after the other
#
# (set_channel_types(), which says what KIND a channel is, lives in
# channel_types.R.)
#
# WHY ONE SHARED HELPER
#   An eeg object is not just a data matrix. After any change to the channels
#   all of these must still agree: data, channels, channel_types, bads, the
#   montage, the channel column of the annotations (and, for epochs, the
#   dimnames of the data array). Epoching, for example, lets only the channels
#   typed "eeg" that are not in bads take part in rejecting a trial, so a
#   channel that is gone from the data but still listed in channel_types
#   would make the wrong channels vote. drop_channels(), pick_channels(),
#   pick_types() and rename_channels() therefore all go through ONE function,
#   .select_channels(), which updates every field in one place, and that one
#   place is what the tests look at most closely.
#
# WHAT IS NOT TOUCHED
#   A fitted ICA or EOG regression model remembers the channel names it was
#   fitted on. These functions do not change such objects. If a channel they
#   need is no longer there, applying them stops with an error that lists the
#   missing channels (see apply_eog_regression() and the ICA functions), so
#   the order is: use the eye channels first (find_bads_eog(),
#   fit_eog_regression()), drop them afterwards.
#
# Author: Christos Dalamarinis
# Date: Nov - 2026
# ============================================================================


# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

#' Which kind of recording object is this? (internal)
#'
#' Stops with a clear message unless \code{x} is a continuous \code{eeg}
#' object or, when \code{epochs_ok} is TRUE, an \code{eeg_epochs} object.
#'
#' @param x         The object given by the user.
#' @param fn        Name of the calling function, for the message.
#' @param epochs_ok Logical. May \code{x} also be an \code{eeg_epochs} object?
#' @return \code{"eeg"} or \code{"eeg_epochs"}.
#' @keywords internal
.channels_kind <- function(x, fn, epochs_ok = TRUE) {
  if (inherits(x, "eeg")) return("eeg")
  if (inherits(x, "eeg_epochs")) {
    if (epochs_ok) return("eeg_epochs")
    stop("ERROR: ", fn, "() works on continuous 'eeg' objects only, not on ",
         "epochs.", call. = FALSE)
  }
  stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg())",
       if (epochs_ok) " or 'eeg_epochs'", ".", call. = FALSE)
}

#' Channel types of an object, with the name-based default (internal)
#'
#' Hand-built objects may have no \code{channel_types}; then the default of
#' \code{\link{new_eeg}} (see \code{\link{classify_channels}}) is used.
#'
#' @param x An \code{eeg} or \code{eeg_epochs} object.
#' @return Character vector, one type per channel.
#' @keywords internal
.types_of <- function(x) {
  if (is.null(x$channel_types)) {
    classify_channels(as.character(x$channels))
  } else {
    x$channel_types
  }
}

#' Check a vector of channel names given by the user (internal)
#'
#' @param v   The value given by the user.
#' @param arg Argument name, for the message.
#' @return \code{v}, unchanged, if it is a non-empty character vector of
#'   non-missing, non-empty names without repeats.
#' @keywords internal
.check_name_vector <- function(v, arg) {
  if (!is.character(v) || length(v) < 1L || anyNA(v) || any(!nzchar(v))) {
    stop("ERROR: '", arg, "' must be a non-empty character vector of ",
         "non-missing, non-empty channel names.", call. = FALSE)
  }
  repeated <- unique(v[duplicated(v)])
  if (length(repeated) > 0L) {
    stop("ERROR: '", arg, "' repeats: ", paste(repeated, collapse = ", "),
         ".", call. = FALSE)
  }
  v
}

#' A short list of names for a history line (internal)
#'
#' @param v   Character vector.
#' @param max How many names to show before the rest is summarised.
#' @return One string, e.g. \code{"Pz, O1"} or \code{"A, B, C ... (7 more)"}.
#' @keywords internal
.name_list <- function(v, max = 10L) {
  if (length(v) <= max) return(paste(v, collapse = ", "))
  paste0(paste(v[seq_len(max)], collapse = ", "), " ... (",
         length(v) - max, " more)")
}

#' Keep only some channels of a montage (internal)
#'
#' @param montage An object of class \code{montage}, or NULL.
#' @param keep_names Character vector of channel names to keep.
#' @return The montage with only those channels, or NULL if none are left
#'   (or if there was no montage).
#' @keywords internal
.montage_keep <- function(montage, keep_names) {
  if (is.null(montage)) return(NULL)
  idx <- which(montage$channels %in% keep_names)
  if (length(idx) == 0L) return(NULL)
  out <- montage
  out$channels  <- montage$channels[idx]
  out$positions <- montage$positions[idx, , drop = FALSE]
  out
}

#' Rename channels inside a montage (internal)
#'
#' @param montage An object of class \code{montage}, or NULL.
#' @param old,new Character vectors of the same length: old names and the new
#'   name of each.
#' @return The montage with the channels renamed (positions untouched).
#' @keywords internal
.montage_rename <- function(montage, old, new) {
  if (is.null(montage)) return(NULL)
  m <- match(montage$channels, old)
  names_new <- montage$channels
  names_new[!is.na(m)] <- new[m[!is.na(m)]]
  out <- montage
  out$channels <- names_new
  if (is.data.frame(out$positions) && "channel" %in% names(out$positions)) {
    out$positions$channel <- names_new
  }
  out
}

#' Remove annotation rows aimed at removed channels (internal)
#'
#' An annotation row can be general (\code{channel} is NA) or aimed at one
#' channel. A row aimed at a channel that no longer exists means nothing any
#' more, so it is deleted; general rows stay.
#'
#' @param anns    The annotations data frame, or NULL.
#' @param removed Character vector of channel names that are gone.
#' @return The annotations without those rows.
#' @keywords internal
.annotations_without_channels <- function(anns, removed) {
  if (is.null(anns) || nrow(anns) == 0L || length(removed) == 0L ||
      !("channel" %in% names(anns))) {
    return(anns)
  }
  gone <- !is.na(anns$channel) & anns$channel %in% removed
  out <- anns[!gone, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Rebuild an eeg or eeg_epochs object with some of its channels (internal)
#'
#' The one place where a channel change is carried out on every field that
#' mentions channels: the data (rows, and for epochs the dimnames), the names,
#' \code{channel_types}, \code{bads}, the montage and the \code{channel}
#' column of the annotations. Used by \code{\link{drop_channels}},
#' \code{\link{pick_channels}}, \code{\link{pick_types}} and
#' \code{\link{rename_channels}}.
#'
#' @param x         An \code{eeg} or \code{eeg_epochs} object.
#' @param keep      Integer positions of the channels to keep, in the order
#'   they should have in the result.
#' @param new_names Character vector (same length as \code{keep}) with the
#'   new name of each kept channel, or NULL to keep the names.
#' @param note      One line for \code{preprocessing_history}.
#' @return The new object (the input is left untouched).
#' @keywords internal
.select_channels <- function(x, keep, new_names = NULL, note) {

  old_names <- as.character(x$channels)
  kept_old  <- old_names[keep]
  kept_new  <- if (is.null(new_names)) kept_old else new_names
  out <- x

  # ---- the data -----------------------------------------------------------
  if (!is.null(x$data)) {
    if (inherits(x, "eeg_epochs")) {              # channels x samples x epochs
      out$data <- x$data[keep, , , drop = FALSE]
      dn <- dimnames(out$data)
      if (!is.null(dn) && !is.null(dn[[1L]])) {
        dn[[1L]] <- kept_new
        dimnames(out$data) <- dn
      }
    } else {                                      # channels x samples
      out$data <- x$data[keep, , drop = FALSE]
      if (!is.null(rownames(out$data))) rownames(out$data) <- kept_new
    }
  }

  # ---- names, types, bad channels -----------------------------------------
  out$channels <- kept_new
  if (!is.null(x$channel_types)) out$channel_types <- x$channel_types[keep]

  bads_old <- as.character(x$bads)
  bads_old <- bads_old[bads_old %in% kept_old]
  out$bads <- kept_new[match(bads_old, kept_old)]

  # ---- the montage ----------------------------------------------------------
  if (!is.null(x$montage)) {
    m <- .montage_keep(x$montage, kept_old)
    if (!is.null(new_names)) m <- .montage_rename(m, kept_old, kept_new)
    out["montage"] <- list(m)                     # keeps the field even if NULL
  }

  # ---- annotations aimed at a channel (continuous data only) -------------
  if (!is.null(x$annotations)) {
    anns <- .annotations_without_channels(x$annotations,
                                          setdiff(old_names, kept_old))
    if (!is.null(new_names) && nrow(anns) > 0L && "channel" %in% names(anns)) {
      has_ch <- !is.na(anns$channel)
      anns$channel[has_ch] <- kept_new[match(anns$channel[has_ch], kept_old)]
    }
    out$annotations <- anns
  }

  out$preprocessing_history <- c(out$preprocessing_history, list(note))
  out
}


# ----------------------------------------------------------------------------
# drop_channels()
# ----------------------------------------------------------------------------
#' Drop Channels
#'
#' Removes the channels you name and keeps everything else. The channels are
#' also removed from \code{bads}, from the montage, and from the annotations
#' that were aimed at them (an annotation row aimed only at a dropped channel
#' is deleted; general rows stay).
#'
#' Use it, for example, to get rid of the eye, heart or status channels once
#' the eye-based steps (\code{\link{find_bads_eog}},
#' \code{\link{fit_eog_regression}}) are done. A fitted ICA or EOG regression
#' model is not changed by this: if you drop a channel such a model was
#' fitted on, applying it later stops with an error that names the missing
#' channel.
#'
#' @param eeg An object of class \code{eeg} or \code{eeg_epochs}.
#' @param channels Character vector of channel names to remove.
#' @param on_missing \code{"raise"} (default): a name that is not a channel of
#'   \code{eeg} is an error and nothing is dropped. \code{"ignore"}: such
#'   names are skipped.
#'
#' @return A new object of the same class without those channels, with a note
#'   appended to \code{preprocessing_history}. The input is left untouched, so
#'   reassign the result (\code{eeg <- drop_channels(eeg, "EXG7")}). Dropping
#'   every channel is an error.
#'
#' @examples
#' \dontrun{
#'   eeg <- drop_channels(eeg, c("VEOG", "HEOG", "EXG7"))
#' }
#'
#' @seealso \code{\link{pick_channels}}, \code{\link{pick_types}},
#'   \code{\link{set_channel_types}}
#' @export
drop_channels <- function(eeg, channels, on_missing = c("raise", "ignore")) {

  .channels_kind(eeg, "drop_channels")
  on_missing <- match.arg(on_missing)
  channels   <- .check_name_vector(channels, "channels")

  missing_ch <- setdiff(channels, eeg$channels)
  if (length(missing_ch) > 0L && on_missing == "raise") {
    stop("ERROR: Channel(s) ", paste(missing_ch, collapse = ", "),
         " not found, nothing dropped.", call. = FALSE)
  }
  to_drop <- intersect(channels, eeg$channels)
  if (length(to_drop) == 0L) return(eeg)          # nothing to do ("ignore")
  if (length(to_drop) == length(eeg$channels)) {
    stop("ERROR: All channels would be dropped.", call. = FALSE)
  }

  keep <- which(!(eeg$channels %in% to_drop))
  .select_channels(eeg, keep,
                   note = paste0("Channels dropped: ", .name_list(to_drop),
                                 " (", length(to_drop), " of ",
                                 length(eeg$channels), ")"))
}


# ----------------------------------------------------------------------------
# pick_channels()
# ----------------------------------------------------------------------------
#' Pick Channels by Name
#'
#' Keeps only the channels you name and removes all the others.
#'
#' @param eeg An object of class \code{eeg} or \code{eeg_epochs}.
#' @param channels Character vector of channel names to keep.
#' @param ordered Logical. \code{TRUE} (default): the result has the channels
#'   in the order you gave them (\code{c("Cz", "Fz")} gives Cz, then Fz).
#'   \code{FALSE}: the channels keep their original order.
#' @param on_missing \code{"raise"} (default): a name that is not a channel of
#'   \code{eeg} is an error. \code{"ignore"}: such names are skipped. (A
#'   misspelled name is an error here even with \code{ordered = FALSE}; it is
#'   never skipped silently.)
#'
#' @return A new object of the same class with only those channels.
#'   \code{bads} keeps the bad channels among them. A note is appended to
#'   \code{preprocessing_history}. Picking no channel at all is an error.
#'
#' @examples
#' \dontrun{
#'   eeg <- pick_channels(eeg, c("Fz", "Cz", "Pz"))
#' }
#'
#' @seealso \code{\link{drop_channels}}, \code{\link{pick_types}}
#' @export
pick_channels <- function(eeg, channels, ordered = TRUE,
                          on_missing = c("raise", "ignore")) {

  .channels_kind(eeg, "pick_channels")
  on_missing <- match.arg(on_missing)
  channels   <- .check_name_vector(channels, "channels")
  if (!is.logical(ordered) || length(ordered) != 1L || is.na(ordered)) {
    stop("ERROR: 'ordered' must be TRUE or FALSE.", call. = FALSE)
  }

  missing_ch <- setdiff(channels, eeg$channels)
  if (length(missing_ch) > 0L && on_missing == "raise") {
    stop("ERROR: Missing channel(s) required by 'channels': ",
         paste(missing_ch, collapse = ", "), ".", call. = FALSE)
  }
  present <- intersect(channels, eeg$channels)
  keep <- if (ordered) match(present, eeg$channels) else
    which(eeg$channels %in% present)
  if (length(keep) == 0L) {
    stop("ERROR: No channels match the selection.", call. = FALSE)
  }

  .select_channels(eeg, keep,
                   note = paste0("Channels picked: ",
                                 .name_list(eeg$channels[keep]), " (",
                                 length(keep), " of ", length(eeg$channels),
                                 ")"))
}


# ----------------------------------------------------------------------------
# pick_types()
# ----------------------------------------------------------------------------
#' Pick Channels by Type
#'
#' Keeps only the channels of the kinds you name (for example all scalp EEG
#' channels) and removes the rest. The kind of every channel is
#' \code{eeg$channel_types}; set it with \code{\link{set_channel_types}}.
#'
#' By default the bad channels are left out as well (\code{exclude =
#' "bads"}), and they disappear from \code{bads} because they are no longer in
#' the data. Pass \code{exclude = character(0)} to keep them.
#'
#' @param eeg An object of class \code{eeg} or \code{eeg_epochs}.
#' @param types Character vector of the kinds to keep, from \code{"eeg"},
#'   \code{"eog"}, \code{"ecg"}, \code{"emg"}, \code{"resp"}, \code{"gsr"},
#'   \code{"temp"}, \code{"bio"}, \code{"misc"} and \code{"status"}. (One
#'   argument instead of a TRUE/FALSE flag per kind, because the first
#'   argument is already called \code{eeg}.)
#' @param exclude Character vector of channel names to leave out even though
#'   their type is wanted. The special value \code{"bads"} (default) stands for
#'   the channels in \code{eeg$bads}. Use \code{character(0)} (or NULL) to
#'   leave out nothing.
#'
#' @return A new object of the same class with only those channels, in their
#'   original order, with a note appended to \code{preprocessing_history}. If
#'   no channel matches, an error.
#'
#' @examples
#' \dontrun{
#'   # keep every scalp channel, bad ones included, drop eye/heart/status
#'   eeg <- pick_types(eeg, types = "eeg", exclude = character(0))
#'
#'   # keep the good scalp channels only (the bad ones are removed)
#'   eeg <- pick_types(eeg, types = "eeg")
#' }
#'
#' @seealso \code{\link{pick_channels}}, \code{\link{drop_channels}},
#'   \code{\link{set_channel_types}}
#' @export
pick_types <- function(eeg, types, exclude = "bads") {

  .channels_kind(eeg, "pick_types")

  if (missing(types) || !is.character(types) || length(types) < 1L ||
      anyNA(types)) {
    stop("ERROR: 'types' must be a non-empty character vector, e.g. ",
         "types = \"eeg\".", call. = FALSE)
  }
  unknown <- setdiff(unique(types), .valid_channel_types())
  if (length(unknown) > 0L) {
    stop("ERROR: 'types' must be from: ",
         paste(.valid_channel_types(), collapse = ", "), "; got: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (!is.null(exclude) && (!is.character(exclude) || anyNA(exclude))) {
    stop("ERROR: 'exclude' must be \"bads\", a character vector of channel ",
         "names, character(0) or NULL.", call. = FALSE)
  }

  excluded <- if (is.null(exclude)) character(0) else exclude
  if ("bads" %in% excluded) {
    excluded <- c(setdiff(excluded, "bads"), as.character(eeg$bads))
  }

  wanted_type <- .types_of(eeg) %in% types
  keep <- which(wanted_type & !(eeg$channels %in% excluded))
  if (length(keep) == 0L) {
    stop("ERROR: No channels match the selection.", call. = FALSE)
  }

  # channels of a wanted type that were left out only because of 'exclude'
  left_out <- as.character(eeg$channels[wanted_type & eeg$channels %in% excluded])

  .select_channels(eeg, keep,
                   note = paste0("Channels picked by type (",
                                 paste(unique(types), collapse = ", "), "): ",
                                 length(keep), " of ", length(eeg$channels),
                                 " kept",
                                 if (length(left_out) > 0L)
                                   paste0("; left out (excluded): ",
                                          .name_list(left_out))
                                 else ""))
}


# ----------------------------------------------------------------------------
# rename_channels()
# ----------------------------------------------------------------------------
#' Rename Channels
#'
#' Changes channel names. The new names are carried into \code{bads}, the
#' montage, the annotations aimed at a channel, and (for epochs) the dimnames
#' of the data. The data themselves do not move: a channel keeps its signal
#' and gets another name.
#'
#' A fitted ICA or EOG regression model is not renamed. If you rename channels
#' after fitting one, applying it stops with an error that names the channels
#' it can no longer find.
#'
#' @param eeg An object of class \code{eeg} or \code{eeg_epochs}.
#' @param mapping Either a named character vector, with the old names as names
#'   and the new names as values (\code{c(Fz = "F_z", EXG1 = "VEOG_up")}), or a
#'   function that is applied to the vector of all channel names and returns
#'   the new names (for example \code{toupper}).
#'
#' @return A new object of the same class with the channels renamed, with a
#'   note appended to \code{preprocessing_history}. It is an error to name a
#'   channel that does not exist, or to end up with two channels of the same
#'   name. Swapping two names (\code{c(Fz = "Cz", Cz = "Fz")}) is allowed:
#'   each channel keeps its own signal and just gets the other name.
#'
#' @examples
#' \dontrun{
#'   eeg <- rename_channels(eeg, c(EXG1 = "VEOG_up", EXG2 = "VEOG_down"))
#'   eeg <- rename_channels(eeg, toupper)
#' }
#'
#' @seealso \code{\link{set_channel_types}}, \code{\link{drop_channels}}
#' @export
rename_channels <- function(eeg, mapping) {

  .channels_kind(eeg, "rename_channels")
  old_all <- as.character(eeg$channels)

  if (is.function(mapping)) {
    new_all <- mapping(old_all)
    if (!is.character(new_all) || length(new_all) != length(old_all) ||
        anyNA(new_all) || any(!nzchar(new_all))) {
      stop("ERROR: 'mapping', a function, must return one non-empty ",
           "character name for every channel.", call. = FALSE)
    }
    how <- "with a function"
  } else {
    if (!is.character(mapping) || length(mapping) < 1L ||
        is.null(names(mapping)) || anyNA(names(mapping)) ||
        any(!nzchar(names(mapping))) || anyNA(mapping) ||
        any(!nzchar(mapping))) {
      stop("ERROR: 'mapping' must be a named character vector, e.g. ",
           "c(Fz = \"F_z\"), or a function.", call. = FALSE)
    }
    if (anyDuplicated(names(mapping)) > 0L) {
      stop("ERROR: 'mapping' names the same channel more than once.",
           call. = FALSE)
    }
    unknown <- setdiff(names(mapping), old_all)
    if (length(unknown) > 0L) {
      stop("ERROR: Invalid channel name(s) ", paste(unknown, collapse = ", "),
           " - not present in the channels.", call. = FALSE)
    }
    new_all <- old_all
    new_all[match(names(mapping), old_all)] <- unname(mapping)
    how <- NULL
  }

  if (anyDuplicated(new_all) > 0L) {
    stop("ERROR: New channel names are not unique, renaming failed.",
         call. = FALSE)
  }

  changed <- which(new_all != old_all)
  if (length(changed) == 0L) return(eeg)           # nothing to rename

  pairs <- paste0(old_all[changed], " -> ", new_all[changed])
  note  <- if (is.null(how)) {
    paste0("Channels renamed: ", .name_list(pairs))
  } else {
    paste0("Channels renamed ", how, ": ", .name_list(pairs))
  }
  .select_channels(eeg, seq_along(old_all), new_names = new_all, note = note)
}


# ----------------------------------------------------------------------------
# add_channels()
# ----------------------------------------------------------------------------
#' Add Channels From Another Recording Object
#'
#' Puts the channels of \code{new} after the channels of \code{eeg}. Typical
#' use: a signal that came from another device or that you computed
#' yourself (heart rate, a derived channel), built into an \code{eeg} object
#' with \code{\link{new_eeg}} and then joined to the main recording.
#'
#' Both objects must have the same sampling rate and the same number of
#' samples (sample for sample the same stretch of time), and no channel name
#' may exist in both.
#'
#' The new channels keep the types they have in \code{new}
#' (\code{new$channel_types}; set them with \code{\link{set_channel_types}}),
#' and their bad channels, if any, are added to \code{bads}. The montage is
#' left as it is: the new channels have no scalp position, the same as the
#' channel added by \code{\link{add_reference_channels}}; use
#' \code{\link{set_montage}} afterwards if you need one. Events and
#' annotations of \code{new} are not copied.
#'
#' @param eeg An object of class \code{eeg} (continuous data).
#' @param new An object of class \code{eeg} holding the channels to add.
#'
#' @return A new \code{eeg} object with the channels of \code{new} appended at
#'   the end, with a note appended to \code{preprocessing_history}.
#'
#' @examples
#' \dontrun{
#'   heart <- new_eeg(data = ecg_matrix, channels = "ECG",
#'                    sampling_rate = eeg$sampling_rate,
#'                    channel_types = "ecg")
#'   eeg <- add_channels(eeg, heart)
#' }
#'
#' @seealso \code{\link{add_reference_channels}} (adds a flat channel for the
#'   recording reference), \code{\link{set_bipolar_reference}} (adds a channel
#'   derived from two others), \code{\link{drop_channels}}
#' @export
add_channels <- function(eeg, new) {

  .channels_kind(eeg, "add_channels", epochs_ok = FALSE)
  if (!inherits(new, "eeg")) {
    stop("ERROR: 'new' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  if (!isTRUE(abs(eeg$sampling_rate - new$sampling_rate) <
              1e-8 * max(eeg$sampling_rate, 1))) {
    stop("ERROR: The sampling rates differ (", eeg$sampling_rate, " Hz in ",
         "'eeg', ", new$sampling_rate, " Hz in 'new').", call. = FALSE)
  }
  if (ncol(eeg$data) != ncol(new$data)) {
    stop("ERROR: All data dimensions except channels must match, got [",
         ncol(eeg$data), "] != [", ncol(new$data), "].", call. = FALSE)
  }
  repeated <- unique(new$channels[duplicated(new$channels)])
  if (length(repeated) > 0L) {
    stop("ERROR: 'new' repeats channel name(s): ",
         paste(repeated, collapse = ", "), ".", call. = FALSE)
  }
  clash <- intersect(new$channels, eeg$channels)
  if (length(clash) > 0L) {
    stop("ERROR: channel(s) already in eeg$channels: ",
         paste(clash, collapse = ", "), ".", call. = FALSE)
  }

  out <- eeg
  data <- rbind(unname(eeg$data), unname(new$data))
  if (!is.null(rownames(eeg$data))) {
    rownames(data) <- c(eeg$channels, new$channels)
  }
  out$data          <- data
  out$channels      <- c(as.character(eeg$channels), as.character(new$channels))
  out$channel_types <- c(.types_of(eeg), .types_of(new))
  out$bads          <- union(as.character(eeg$bads), as.character(new$bads))

  ignored <- character(0)
  if (!is.null(new$events) && nrow(new$events) > 0L) ignored <- c(ignored, "events")
  if (!is.null(new$annotations) && nrow(new$annotations) > 0L) {
    ignored <- c(ignored, "annotations")
  }
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Channels added: ", .name_list(new$channels), " (",
                length(new$channels), ")",
                if (length(ignored) > 0L)
                  paste0("; the ", paste(ignored, collapse = " and "),
                         " of the added object were not copied") else "")))
  out
}


# ----------------------------------------------------------------------------
# combine_channels()
# ----------------------------------------------------------------------------
#' Combine Groups of Channels Into One Channel Each
#'
#' Merges every group of channels into a single new channel, sample by sample,
#' with the mean, the median, the standard deviation, or your own function.
#' A typical use is a region of interest: the mean of the frontal electrodes.
#'
#' The result holds only the combined channels (plus the status channels if
#' you ask for them). The channels that went into a group must all be of the
#' same type, and the combined channel gets that type. The combined channels
#' have no scalp position, so the montage is removed; annotation rows aimed at
#' a channel that is no longer there are deleted; \code{bads} is empty.
#'
#' @param eeg An object of class \code{eeg} (continuous data).
#' @param groups A named list. Each name becomes the name of a new channel, and
#'   the element lists the channels to combine, as channel names or as channel
#'   numbers: \code{list(Frontal = c("Fz", "F3", "F4"), Parietal = c("Pz",
#'   "P3"))}.
#' @param method \code{"mean"} (default), \code{"median"}, \code{"std"} (the
#'   population standard deviation: the sum of squares is divided by the
#'   number of channels, not by that number minus one), or a function that
#'   receives the group as a matrix (channels in rows, samples in columns) and
#'   returns one value per sample.
#' @param keep_status Logical. Keep the status channel(s) in the result? They
#'   come first, then the combined channels. Default: \code{FALSE}.
#' @param drop_bad Logical. Leave the channels in \code{eeg$bads} out of the
#'   groups? Default: \code{FALSE}, which includes them.
#'
#' @return A new \code{eeg} object with the combined channels, with a note
#'   appended to \code{preprocessing_history}.
#'
#' @examples
#' \dontrun{
#'   roi <- combine_channels(eeg, list(Frontal = c("Fz", "F3", "F4")),
#'                           drop_bad = TRUE)
#' }
#'
#' @seealso \code{\link{pick_channels}}, \code{\link{set_bipolar_reference}}
#' @export
combine_channels <- function(eeg, groups, method = "mean",
                             keep_status = FALSE, drop_bad = FALSE) {

  .channels_kind(eeg, "combine_channels", epochs_ok = FALSE)

  if (!is.list(groups) || length(groups) < 1L || is.null(names(groups)) ||
      anyNA(names(groups)) || any(!nzchar(names(groups)))) {
    stop("ERROR: 'groups' must be a named list, e.g. ",
         "list(Frontal = c(\"Fz\", \"F3\")).", call. = FALSE)
  }
  if (anyDuplicated(names(groups)) > 0L) {
    stop("ERROR: 'groups' has a name more than once.", call. = FALSE)
  }
  if (!is.function(method) &&
      !(is.character(method) && length(method) == 1L &&
        method %in% c("mean", "median", "std"))) {
    stop("ERROR: 'method' must be a function, or one of \"mean\", ",
         "\"median\", \"std\"; got ",
         if (is.character(method)) paste0("\"", method[1L], "\"") else
           class(method)[1L], ".", call. = FALSE)
  }
  for (flag in c("keep_status", "drop_bad")) {
    v <- get(flag)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("ERROR: '", flag, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }

  n_ch  <- nrow(eeg$data)
  types <- .types_of(eeg)

  # ---- which channels are in each group ------------------------------------
  members <- lapply(names(groups), function(g) {
    el <- groups[[g]]
    if (is.numeric(el) && length(el) >= 1L && !anyNA(el)) {
      idx <- as.integer(el)
      if (any(idx < 1L | idx > n_ch)) {
        stop("ERROR: Some channel indices are out of bounds (group '", g,
             "').", call. = FALSE)
      }
    } else if (is.character(el) && length(el) >= 1L && !anyNA(el)) {
      idx <- match(el, eeg$channels)
      if (anyNA(idx)) {
        stop("ERROR: Channel(s) ", paste(el[is.na(idx)], collapse = ", "),
             " of group '", g, "' not found.", call. = FALSE)
      }
    } else {
      stop("ERROR: Group '", g, "' must be a non-empty vector of channel ",
           "names or channel numbers.", call. = FALSE)
    }
    if (anyDuplicated(idx) > 0L) idx <- unique(idx)
    if (isTRUE(drop_bad)) idx <- idx[!(eeg$channels[idx] %in% eeg$bads)]
    if (length(idx) == 0L) {
      stop("ERROR: Group '", g, "' has no channels left after leaving out ",
           "the bad channels.", call. = FALSE)
    }
    gt <- unique(types[idx])
    if (length(gt) > 1L) {
      stop("ERROR: Cannot combine sensors of different types; \"", g,
           "\" contains types ", paste(gt, collapse = ", "), ".",
           call. = FALSE)
    }
    idx
  })
  names(members) <- names(groups)
  group_types <- vapply(members, function(i) types[i[1L]], "")

  # ---- combine, sample by sample -------------------------------------------
  combined <- matrix(0, nrow = length(members), ncol = ncol(eeg$data))
  for (g in seq_along(members)) {
    block <- eeg$data[members[[g]], , drop = FALSE]
    k <- nrow(block)
    value <- if (is.function(method)) {
      method(block)
    } else if (method == "mean") {
      colMeans(block)
    } else if (method == "median") {
      apply(block, 2L, stats::median)
    } else {                                      # "std": divide by k, not k - 1
      sqrt(colMeans((block - rep(colMeans(block), each = k))^2))
    }
    if (!is.numeric(value) || length(value) != ncol(eeg$data)) {
      stop("ERROR: 'method' must return one number per sample (",
           ncol(eeg$data), "); group '", names(members)[g], "' got ",
           length(value), ".", call. = FALSE)
    }
    combined[g, ] <- value
  }

  status_idx <- if (isTRUE(keep_status)) which(types == "status") else integer(0)
  new_names  <- c(as.character(eeg$channels[status_idx]), names(members))
  clash <- new_names[duplicated(new_names)]
  if (length(clash) > 0L) {
    stop("ERROR: The new channel name(s) ", paste(unique(clash), collapse = ", "),
         " would exist twice (a group name equals a status channel).",
         call. = FALSE)
  }

  out <- eeg
  data <- rbind(eeg$data[status_idx, , drop = FALSE], combined)
  if (!is.null(rownames(eeg$data))) rownames(data) <- new_names else
    rownames(data) <- NULL
  out$data          <- data
  out$channels      <- new_names
  out$channel_types <- c(types[status_idx], unname(group_types))
  out$bads          <- character(0)
  out["montage"]    <- list(NULL)
  if (!is.null(eeg$annotations)) {
    out$annotations <- .annotations_without_channels(
      eeg$annotations, setdiff(as.character(eeg$channels),
                               as.character(eeg$channels[status_idx])))
  }

  method_txt <- if (is.function(method)) "custom function" else method
  parts <- vapply(names(members), function(g) {
    paste0(g, " = ", paste(eeg$channels[members[[g]]], collapse = " + "))
  }, "")
  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Channels combined (", method_txt, "): ",
                .name_list(parts, max = 5L),
                "; the original channels were not kept",
                if (length(status_idx) > 0L) " (status kept)" else "")))
  out
}


# ----------------------------------------------------------------------------
# eeg_crop()
# ----------------------------------------------------------------------------
#' Cut Time Off the Start and the End of a Recording
#'
#' Keeps the part of a continuous recording between \code{tmin} and
#' \code{tmax}. The first kept sample is the one nearest to \code{tmin}, the
#' last one the one nearest to \code{tmax} (both are kept; see
#' \code{include_tmax}). After cropping, \code{times} start again at 0.
#'
#' Everything that is tied to time moves with the data:
#' \describe{
#'   \item{events}{Events outside the kept part are removed; the others are
#'     shifted, in both \code{onset} (sample number) and \code{onset_time}
#'     (seconds).}
#'   \item{annotations}{Shifted by the start time and cut to the kept part. An
#'     annotation that only partly lies in the kept part is shortened; one
#'     completely outside is removed. (An annotation that ends exactly at the
#'     new start stays, with length 0, and one that starts at the last kept
#'     sample stays with a length of one sample.)}
#' }
#'
#' @param eeg An object of class \code{eeg} (continuous data).
#' @param tmin Start of the part to keep, in seconds from the start of the
#'   recording. Default: 0. Must not be negative.
#' @param tmax End of the part to keep, in seconds. \code{NULL} (default) means
#'   the end of the recording. Must not be larger than the time of the last
#'   sample, and not smaller than \code{tmin}.
#' @param include_tmax Logical. Keep the sample at \code{tmax}? Default
#'   \code{TRUE}. With \code{FALSE} the last sample is left out (also when
#'   \code{tmax} is \code{NULL}).
#'
#' @return A new \code{eeg} object with the kept part, with a note appended to
#'   \code{preprocessing_history} that says what was cut. The input is left
#'   untouched.
#'
#' @examples
#' \dontrun{
#'   # drop the first 30 s and everything after minute 10
#'   eeg <- eeg_crop(eeg, tmin = 30, tmax = 600)
#' }
#'
#' @seealso \code{\link{eeg_concatenate}}
#' @export
eeg_crop <- function(eeg, tmin = 0, tmax = NULL, include_tmax = TRUE) {

  .channels_kind(eeg, "eeg_crop", epochs_ok = FALSE)

  sf        <- eeg$sampling_rate
  n         <- ncol(eeg$data)
  last_time <- (n - 1) / sf

  one_number <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!one_number(tmin)) {
    stop("ERROR: 'tmin' must be one number (seconds).", call. = FALSE)
  }
  if (is.null(tmax)) tmax <- last_time
  if (!one_number(tmax)) {
    stop("ERROR: 'tmax' must be one number (seconds) or NULL.", call. = FALSE)
  }
  if (!is.logical(include_tmax) || length(include_tmax) != 1L ||
      is.na(include_tmax)) {
    stop("ERROR: 'include_tmax' must be TRUE or FALSE.", call. = FALSE)
  }

  if (tmin < 0) {
    stop("ERROR: tmin (", format(tmin), ") must be >= 0.", call. = FALSE)
  }
  if (tmax > last_time + 1e-9) {
    stop("ERROR: tmax (", format(tmax), ") must be less than or equal to the ",
         "max time (", format(round(last_time, 4)), " s).", call. = FALSE)
  }
  if (tmin > tmax) {
    stop("ERROR: tmin (", format(tmin), ") must not be greater than tmax (",
         format(tmax), ").", call. = FALSE)
  }

  start0 <- as.integer(round(tmin * sf))          # 0-based sample numbers
  stop0  <- as.integer(round(tmax * sf))
  if (!include_tmax) stop0 <- stop0 - 1L
  if (stop0 < start0) {
    stop("ERROR: No samples would be left after cropping.", call. = FALSE)
  }
  n_new <- stop0 - start0 + 1L

  out <- eeg
  out$data  <- eeg$data[, (start0 + 1L):(stop0 + 1L), drop = FALSE]
  out$times <- (0:(n_new - 1L)) / sf

  # ---- events: keep the ones inside, shift them ---------------------------
  ev <- eeg$events
  if (is.data.frame(ev) && nrow(ev) > 0L && "onset" %in% names(ev)) {
    inside <- ev$onset >= (start0 + 1L) & ev$onset <= (stop0 + 1L)
    ev <- ev[inside, , drop = FALSE]
    ev$onset <- ev$onset - start0
    if ("onset_time" %in% names(ev)) ev$onset_time <- (ev$onset - 1) / sf
    rownames(ev) <- NULL
    out$events <- ev
  }

  # ---- annotations: shift, cut to the kept part ------------------------------
  an <- eeg$annotations
  if (is.data.frame(an) && nrow(an) > 0L) {
    tol      <- 1e-9
    t0       <- start0 / sf
    span_end <- n_new / sf                         # one sample past the last one
    last_t   <- (n_new - 1) / sf
    new_on   <- an$onset - t0
    new_end  <- new_on + an$duration
    inside   <- new_end >= -tol & new_on <= last_t + tol
    an <- an[inside, , drop = FALSE]
    new_on  <- new_on[inside]
    new_end <- new_end[inside]
    an$onset    <- pmax(new_on, 0)
    an$duration <- pmax(pmin(new_end, span_end) - an$onset, 0)
    rownames(an) <- NULL
    out$annotations <- an
  }

  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Cropped: samples ", start0 + 1L, " to ", stop0 + 1L, " of ",
                n, " kept (", format(round(start0 / sf, 4)), " s to ",
                format(round(stop0 / sf, 4)), " s of the original; ", n_new,
                " samples)")))
  out
}


# ----------------------------------------------------------------------------
# eeg_concatenate()
# ----------------------------------------------------------------------------
#' Join Recordings One After the Other
#'
#' Joins two or more continuous recordings in the order given, for example the
#' runs of one session. The recordings must belong together: the same channels
#' in the same order, the same types, the same bad channels, the same sampling
#' rate and the same reference. Anything else is an error that says what
#' differs. (If the bad channels differ, mark the same channels as bad in every
#' recording first, for example \code{eeg$bads <- union(a$bads, b$bads)}.)
#'
#' The result has one continuous \code{times} axis that starts at 0. Events
#' and annotations of later recordings are moved by the length of the
#' recordings before them. At every join one annotation row is added,
#' \code{"BAD_boundary"} with length 0 and no channel, so that epochs that
#' stretch across a join are rejected like any epoch that overlaps a bad
#' stretch. \code{\link{epoch_eeg}} does this with its default
#' (\code{reject_by_annotation = TRUE}); if you give that argument a list of
#' descriptions, put \code{"BAD_boundary"} in the list. Metadata, montage and
#' the reference come from the first recording.
#'
#' Filter each recording before joining them: filtering after the join creates
#' edge effects at every join.
#'
#' @param ... Two or more objects of class \code{eeg}, or one list of them.
#'
#' @return A new \code{eeg} object. A note in \code{preprocessing_history} gives
#'   the number of recordings and their lengths.
#'
#' @examples
#' \dontrun{
#'   eeg <- eeg_concatenate(run1, run2, run3)
#'   eeg <- eeg_concatenate(list(run1, run2))
#' }
#'
#' @seealso \code{\link{eeg_crop}}
#' @export
eeg_concatenate <- function(...) {

  recs <- list(...)
  if (length(recs) == 1L && is.list(recs[[1L]]) &&
      !inherits(recs[[1L]], "eeg")) {
    recs <- recs[[1L]]
  }
  recs <- unname(recs)
  if (length(recs) < 2L) {
    stop("ERROR: Give at least two recordings to join.", call. = FALSE)
  }
  for (i in seq_along(recs)) {
    if (!inherits(recs[[i]], "eeg")) {
      stop("ERROR: Recording ", i, " is not an object of class 'eeg' (see ",
           "new_eeg()); epochs cannot be joined with eeg_concatenate().",
           call. = FALSE)
    }
  }

  first <- recs[[1L]]
  for (i in 2:length(recs)) {
    r <- recs[[i]]
    where <- paste0("Recording ", i, " does not match recording 1: ")
    if (!identical(as.character(r$channels), as.character(first$channels))) {
      if (length(r$channels) != length(first$channels)) {
        stop("ERROR: ", where, "it has ", length(r$channels), " channels, ",
             "recording 1 has ", length(first$channels), ".", call. = FALSE)
      }
      pos <- which(as.character(r$channels) != as.character(first$channels))[1L]
      stop("ERROR: ", where, "the channel names or their order differ (first ",
           "difference at channel ", pos, ": '", r$channels[pos], "' versus '",
           first$channels[pos], "'). Rename or reorder the channels first.",
           call. = FALSE)
    }
    if (!isTRUE(abs(r$sampling_rate - first$sampling_rate) <
                1e-8 * max(first$sampling_rate, 1))) {
      stop("ERROR: ", where, "the sampling rates differ (", r$sampling_rate,
           " Hz versus ", first$sampling_rate, " Hz).", call. = FALSE)
    }
    if (!identical(as.character(.types_of(r)), as.character(.types_of(first)))) {
      stop("ERROR: ", where, "the channel types differ.", call. = FALSE)
    }
    if (!setequal(as.character(r$bads), as.character(first$bads))) {
      stop("ERROR: ", where, "the bad channels differ (recording 1: ",
           if (length(first$bads)) paste(first$bads, collapse = ", ") else
             "none", "; recording ", i, ": ",
           if (length(r$bads)) paste(r$bads, collapse = ", ") else "none",
           "). Mark the same channels as bad in every recording first.",
           call. = FALSE)
    }
    if (!identical(as.character(r$reference), as.character(first$reference))) {
      stop("ERROR: ", where, "the reference differs ('", r$reference,
           "' versus '", first$reference, "').", call. = FALSE)
    }
  }

  sf      <- first$sampling_rate
  n_each  <- vapply(recs, function(r) ncol(r$data), 1L)
  offsets <- c(0L, cumsum(n_each))[seq_along(recs)]        # samples before each
  n_total <- sum(n_each)

  out <- first
  data <- do.call(cbind, lapply(recs, function(r) unname(r$data)))
  if (!is.null(rownames(first$data))) rownames(data) <- rownames(first$data)
  out$data  <- data
  out$times <- (0:(n_total - 1L)) / sf

  # ---- events: shift by the samples before each recording ------------------
  ev_parts <- lapply(seq_along(recs), function(i) {
    ev <- recs[[i]]$events
    if (!is.data.frame(ev) || nrow(ev) == 0L || !("onset" %in% names(ev))) {
      return(NULL)
    }
    ev$onset <- ev$onset + offsets[i]
    if ("onset_time" %in% names(ev)) ev$onset_time <- (ev$onset - 1) / sf
    ev
  })
  ev_parts <- Filter(Negate(is.null), ev_parts)
  if (length(ev_parts) > 0L) {
    ev_names <- lapply(ev_parts, names)
    if (!all(vapply(ev_names, function(nm) setequal(nm, ev_names[[1L]]), TRUE))) {
      stop("ERROR: The events tables of the recordings have different ",
           "columns.", call. = FALSE)
    }
    ev_all <- do.call(rbind, lapply(ev_parts, function(e) e[, ev_names[[1L]],
                                                            drop = FALSE]))
    rownames(ev_all) <- NULL
    out$events <- ev_all
  }

  # ---- annotations: shift, add one boundary row at every join -----------------
  an_parts <- lapply(seq_along(recs), function(i) {
    an <- recs[[i]]$annotations
    if (is.null(an)) an <- .empty_annotations()
    if (nrow(an) > 0L) an$onset <- an$onset + offsets[i] / sf
    an
  })
  boundaries <- lapply(2:length(recs), function(i) {
    data.frame(onset = offsets[i] / sf, duration = 0,
               description = "BAD_boundary", channel = NA_character_,
               stringsAsFactors = FALSE)
  })
  an_all <- do.call(rbind, c(an_parts, boundaries))
  an_all <- an_all[order(an_all$onset), , drop = FALSE]
  rownames(an_all) <- NULL
  out$annotations <- an_all

  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Concatenated ", length(recs), " recordings (",
                paste(n_each, collapse = " + "), " samples)")))
  out
}
