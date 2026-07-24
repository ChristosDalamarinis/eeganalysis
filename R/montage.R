#' ============================================================================
#'                     EEG Montage: Channel Scalp Positions
#' ============================================================================
#'
#' A montage stores the spatial layout of a set of channels on the scalp
#' (Cartesian x/y/z and spherical theta/phi/radius coordinates), built from
#' the built-in BioSemi 64-channel electrode database (see
#' get_electrode_database(), R/channel_info2.R) and attached to an eeg object
#' via set_montage() - mirroring MNE-Python's raw.set_montage() workflow.
#'
#' Author: Christos Dalamarinis
#' Date: July 2026
#' ============================================================================

#' Create a New Montage Object (low-level constructor)
#'
#' Low-level constructor for the 'montage' S3 class. Most users should call
#' \code{\link{create_montage}} instead, which builds a montage from the
#' package's built-in BioSemi electrode template. Use \code{new_montage()}
#' directly only when supplying custom channel positions of your own.
#'
#' @param channels Character vector of channel names.
#' @param positions Data frame with one row per channel (in the same order as
#'   \code{channels}) and columns \code{channel, x, y, z, theta, phi, radius}.
#' @param coord_frame Character string identifying the coordinate template or
#'   source this montage was built from (e.g. \code{"biosemi64"}).
#'
#' @return An object of class \code{'montage'}: a list with \code{channels},
#'   \code{positions}, and \code{coord_frame}.
#'
#' @export
new_montage <- function(channels, positions, coord_frame = "biosemi64") {
  structure(
    list(
      channels = as.character(channels),
      positions = positions,
      coord_frame = as.character(coord_frame)
    ),
    class = "montage"
  )
}

#' Create a Montage From the Built-in BioSemi Electrode Template
#'
#' Builds a \code{montage} object giving scalp positions (Cartesian
#' x/y/z and spherical theta/phi/radius) for a set of channel names, drawn
#' from the package's built-in 64-channel BioSemi electrode database
#' (\code{get_electrode_database}). Channel names may be given in
#' either standard 10-20/10-10 notation (e.g. \code{"Cz"}) or BioSemi internal
#' notation (e.g. \code{"B16"}); matching is case-insensitive.
#'
#' @param channels Character vector of channel names to include. If
#'   \code{NULL} (default), all 64 standard 10-20/10-10 electrode names are
#'   used.
#' @param template Character string naming the coordinate template to use.
#'   Currently only \code{"biosemi64"} (the default) is supported.
#'
#' @return An object of class \code{'montage'}. See \code{new_montage}.
#'
#' @details
#' Requested channels that have no scalp position in the template - either
#' because they are not in the database at all, or because they are external
#' channels (EOG/ECG/EMG/GSR/etc., which have no coordinates) - are dropped
#' with a warning naming exactly which channels were excluded. Use
#' \code{detect_electrode_naming_system} first if you are unsure which
#' naming convention a dataset uses.
#'
#' @examples
#' \dontrun{
#'   # All 64 standard electrodes
#'   m <- create_montage()
#'   print(m)
#'
#'   # Just a subset, mixing standard and BioSemi names
#'   m <- create_montage(c("Cz", "Fz", "Pz", "B1"))
#' }
#'
#' @seealso \code{\link{set_montage}}, \code{\link{get_electrode_database}}
#'
#' @export
create_montage <- function(channels = NULL, template = "biosemi64") {

  if (!identical(template, "biosemi64")) {
    stop("ERROR: Unknown montage template '", template,
         "'. Only 'biosemi64' is currently supported.", call. = FALSE)
  }

  electrode_db <- get_electrode_database()

  if (is.null(channels)) {
    is_scalp <- vapply(electrode_db, function(e) {
      isTRUE(e$standard_system %in% c("10-20", "10-10"))
    }, logical(1))
    channels <- unique(vapply(electrode_db[is_scalp],
                               function(e) e$standard_name, character(1)))
  }

  channels <- as.character(channels)
  lookup_keys <- tolower(channels)

  in_db <- lookup_keys %in% names(electrode_db)

  has_coords <- rep(FALSE, length(channels))
  if (any(in_db)) {
    has_coords[in_db] <- !is.na(vapply(
      electrode_db[lookup_keys[in_db]],
      function(e) e$cartesian_coords$x, numeric(1)
    ))
  }

  ok <- in_db & has_coords
  if (any(!ok)) {
    warning("create_montage(): ", sum(!ok),
            " channel(s) have no scalp position in the '", template,
            "' template and will be dropped: ",
            paste(channels[!ok], collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }
  if (!any(ok)) {
    stop("ERROR: None of the requested channels have a scalp position in ",
         "the '", template, "' template.", call. = FALSE)
  }

  channels_ok <- channels[ok]
  entries <- electrode_db[lookup_keys[ok]]

  positions <- data.frame(
    channel = channels_ok,
    x      = vapply(entries, function(e) e$cartesian_coords$x, numeric(1)),
    y      = vapply(entries, function(e) e$cartesian_coords$y, numeric(1)),
    z      = vapply(entries, function(e) e$cartesian_coords$z, numeric(1)),
    theta  = vapply(entries, function(e) e$spherical_coords$theta, numeric(1)),
    phi    = vapply(entries, function(e) e$spherical_coords$phi, numeric(1)),
    radius = vapply(entries, function(e) e$spherical_coords$radius, numeric(1)),
    stringsAsFactors = FALSE
  )

  new_montage(channels = channels_ok, positions = positions, coord_frame = template)
}

#' Attach a Montage to an EEG Object
#'
#' Attaches channel scalp positions to an \code{eeg} object, mirroring
#' MNE-Python's \code{raw.set_montage()}. Only channels classified as
#' \code{"eeg"} in \code{eeg_obj$channel_types} (see the internal
#' \code{classify_channels()}, R/eeg_class.R) receive positions - external
#' (EOG/ECG/EMG/GSR/etc.) and status channels are never spatial and are left
#' out of the attached montage.
#'
#' @param eeg_obj An object of class \code{'eeg'}.
#' @param montage An object of class \code{'montage'}, typically built with
#'   \code{\link{create_montage}}.
#'
#' @return The input \code{eeg} object with \code{$montage} set to a montage
#'   restricted to the channels actually present in \code{eeg_obj}, and an
#'   entry appended to \code{$preprocessing_history}.
#'
#' @details
#' Channels of type \code{"eeg"} in \code{eeg_obj} that have no matching
#' entry in \code{montage}, and channels in \code{montage} that are not
#' present in \code{eeg_obj}, both produce a warning (not an error) - such
#' channels simply have no scalp position to plot.
#'
#' @examples
#' \dontrun{
#'   eeg <- read_bdf_native("recording.bdf")
#'   eeg <- set_montage(eeg, create_montage())
#' }
#'
#' @seealso \code{\link{create_montage}}, \code{\link{plot_topography}}
#'
#' @export
set_montage <- function(eeg_obj, montage) {

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!inherits(montage, "montage")) {
    stop("ERROR: 'montage' must be an object of class 'montage' ",
         "(see create_montage()).", call. = FALSE)
  }

  eeg_channels <- eeg_obj$channels[eeg_obj$channel_types == "eeg"]

  missing_in_montage <- setdiff(eeg_channels, montage$channels)
  if (length(missing_in_montage) > 0) {
    warning("set_montage(): ", length(missing_in_montage),
            " EEG channel(s) in the data have no position in this montage: ",
            paste(missing_in_montage, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  missing_in_data <- setdiff(montage$channels, eeg_obj$channels)
  if (length(missing_in_data) > 0) {
    warning("set_montage(): ", length(missing_in_data),
            " montage channel(s) are not present in the data and will be ",
            "ignored: ", paste(missing_in_data, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  kept <- intersect(montage$channels, eeg_channels)
  if (length(kept) == 0) {
    stop("ERROR: No overlap between montage channels and EEG channels in ",
         "this object.", call. = FALSE)
  }

  eeg_obj$montage <- new_montage(
    channels = kept,
    positions = montage$positions[montage$positions$channel %in% kept, , drop = FALSE],
    coord_frame = montage$coord_frame
  )

  eeg_obj$preprocessing_history <- c(
    eeg_obj$preprocessing_history,
    list(paste0("Montage attached: ", montage$coord_frame, " (",
                length(kept), " channel(s))"))
  )

  eeg_obj
}

#' Print Method for Montage Objects
#'
#' Displays a formatted summary of a montage: coordinate template, channel
#' count, and channel names. Mirrors the style of \code{\link{print.eeg}}.
#'
#' @param x An object of class \code{'montage'}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.montage <- function(x, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Montage Object Summary\n")
  cat(strrep("=", 70), "\n\n")

  cat("  Coordinate template: ", x$coord_frame, "\n")
  cat("  Channels:            ", length(x$channels), "\n")

  channel_display <- paste(x$channels[seq_len(min(5, length(x$channels)))],
                            collapse = ", ")
  if (length(x$channels) > 5) {
    channel_display <- paste0(channel_display, " ... (+",
                               length(x$channels) - 5, " more)")
  }
  cat("    List:              ", channel_display, "\n")

  cat("\n", strrep("=", 70), "\n\n", sep = "")

  invisible(x)
}
