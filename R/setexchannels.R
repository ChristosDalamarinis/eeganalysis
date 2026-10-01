#' ============================================================================
#'                    External Channel Detection (Database-Driven)
#' ============================================================================
#'
#' A single question: "does this channel name match a known BioSemi
#' auxiliary port?" - looked up directly in the electrode database, never
#' guessed or inferred. This is a suggestion tool only: it does not read or
#' write channel_types, and does not rename anything.
#'
#' This used to be step one of a three-step workflow (detect, then
#' interactively label, then rename the channel to bake the label into its
#' name) built before channel_types existed on the eeg object. The other two
#' steps - identify_external_channels() and apply_external_labels() - are
#' retired: the fact they used to encode by renaming a channel (e.g. "EXG1"
#' -> "EOG_L (EXG1)") now lives directly in eeg$channel_types instead, set
#' with \code{\link{set_channel_types}} (or a reader's own eog=/misc=
#' arguments, see \code{\link{read_bdf_native}}). Nothing needs to be renamed
#' for the package to know a channel's role any more.
#'
#' detect_external_channels() still has a job: suggesting candidates before
#' you decide what to tell set_channel_types(), e.g.
#' \code{detect_external_channels(eeg$channels)} to see which of your
#' channels look like known auxiliary ports.
#'
#' Author: Christos Dalamarinis
#' Date: Jan - 2026
#' ============================================================================
#'
#' Detect External Channels (Non-Interactive, Database-Driven)
#'
#' @description
#' Detects external channels using the electrode database without user interaction.
#' Useful for scripting and automation.
#'
#' @param data A data frame, matrix, list, or character vector containing channel names
#' @param channel_col Character string specifying column with channel names (for data frames)
#'
#' @return Character vector of detected external channel names
#'
#' @export
detect_external_channels <- function(data, channel_col = NULL) {

  # Extract channel names
  if (is.character(data)) {
    channel_names <- trimws(data)
  } else if (is.data.frame(data) || is.matrix(data)) {
    if (is.null(channel_col)) {
      channel_names <- colnames(data)
    } else {
      if (!channel_col %in% colnames(data)) {
        stop("Specified channel_col '", channel_col, "' not found in data.")
      }
      channel_names <- unique(data[[channel_col]])
    }
  } else if (is.list(data)) {
    if (!is.null(data$channels)) {
      channel_names <- trimws(data$channels)
    } else if (!is.null(data$channel_names)) {
      channel_names <- trimws(data$channel_names)
    } else {
      stop("Cannot find channel names in list.")
    }
  } else {
    stop("Invalid input type.")
  }

  # Get electrode database
  electrode_db <- get_electrode_database()

  # Identify external channels
  external_channels <- character()

  for (ch_name in channel_names) {
    ch_lower <- tolower(ch_name)

    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]

      if (electrode_info$position_type %in% c("External", "GSR", "Ergo/AUX",
                                              "Respiration", "Plethysmograph",
                                              "Temperature")) {
        external_channels <- c(external_channels, ch_name)
      }
    }
  }

  return(external_channels)
}
