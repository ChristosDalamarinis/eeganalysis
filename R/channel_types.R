#' ============================================================================
#'              Channel Types: Tell the Package What a Channel Is
#' ============================================================================
#'
#' Every eeg object carries eeg$channel_types - "eeg", "status", or one of a
#' handful of specific physiological roles ("eog", "ecg", "emg", "resp",
#' "gsr", "temp", "bio", "misc") - that later steps read to tell scalp EEG
#' channels (bad-channel detection, ICA, interpolation, ...) apart from eye/
#' heart/muscle channels (EOG regression, find_bads_eog(), ...).
#'
#' new_eeg() only ever works out "eeg" vs "status" from the channel NAME (see
#' classify_channels(), R/eeg_class.R) - nothing more specific is ever
#' guessed. set_channel_types() is where a channel's real role is stated
#' after the fact: at load time, through a reader's own eog=/misc=
#' arguments (see read_bdf_native()), or any time afterward with this
#' function.
#'
#' Author: Christos Dalamarinis
#' Date: Oct - 2026
#' ============================================================================

# ----------------------------------------------------------------------------
# set_channel_types() - change the type of the channels you name
# ----------------------------------------------------------------------------
#' Set the Type of Individual Channels
#'
#' Changes the type of the channels you name and leaves every other channel
#' exactly as it is. Use it when a channel's real role isn't - and can't be -
#' worked out from its name: an eye, heart, or muscle channel that
#' \code{\link{new_eeg}} would otherwise leave typed \code{"eeg"} by default.
#' A port of the idea behind MNE-Python's \code{raw.set_channel_types()}.
#'
#' @param eeg An object of class \code{eeg} (see \code{\link{new_eeg}}).
#' @param types Named character vector: the names are channel names from
#'   \code{eeg$channels}, the values are the new type for each one - one of
#'   \code{"eeg"}, \code{"eog"}, \code{"ecg"}, \code{"emg"}, \code{"resp"},
#'   \code{"gsr"}, \code{"temp"}, \code{"bio"}, \code{"misc"}, or
#'   \code{"status"} (see \code{\link{new_eeg}}'s \code{channel_types}
#'   argument for what each one means). E.g.
#'   \code{c(VEOG = "eog", EXG5 = "gsr")}.
#'
#' @return A new object of class \code{eeg} with the named channels' types
#'   changed and a note appended to \code{preprocessing_history}. Since R
#'   does not change arguments in place, reassign the result
#'   (\code{eeg <- set_channel_types(eeg, ...)}); the input object itself is
#'   left untouched.
#'
#' @examples
#' \dontrun{
#'   eeg <- read_bdf_native("sub_0_ses_1.bdf")
#'   eeg <- set_channel_types(eeg, c(VEOG = "eog", HEOG = "eog", EXG5 = "gsr"))
#' }
#'
#' @seealso \code{\link{new_eeg}}, \code{\link{read_bdf_native}}
#'
#' @export
set_channel_types <- function(eeg, types) {

  # ========== VALIDATE eeg ==========

  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  # ========== VALIDATE types ==========

  if (!is.character(types) || length(types) == 0 ||
      is.null(names(types)) || any(is.na(names(types)) | names(types) == "")) {
    stop("ERROR: 'types' must be a named character vector, e.g. ",
         "c(VEOG = \"eog\").", call. = FALSE)
  }

  unknown_channels <- setdiff(names(types), eeg$channels)
  if (length(unknown_channels) > 0) {
    stop("ERROR: channel(s) not found in eeg$channels: ",
         paste(unknown_channels, collapse = ", "), call. = FALSE)
  }

  unknown_types <- setdiff(unique(unname(types)), .valid_channel_types())
  if (length(unknown_types) > 0) {
    stop("ERROR: types must be one of: ",
         paste(.valid_channel_types(), collapse = ", "), "; got: ",
         paste(unknown_types, collapse = ", "), call. = FALSE)
  }

  # ========== SET THE TYPES ==========

  out <- eeg

  # An old or hand-built object may have no types yet: start from the
  # name-based default so the channels that were not named still get one.
  if (is.null(out$channel_types)) {
    out$channel_types <- classify_channels(as.character(out$channels))
  }

  out$channel_types[match(names(types), out$channels)] <- unname(types)

  # ========== PREPROCESSING HISTORY ==========

  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Channel types set: ",
                paste0(names(types), " -> ", unname(types), collapse = ", ")))
  )

  out
}
