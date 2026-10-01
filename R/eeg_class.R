#' ============================================================================
#'                        EEG S3 Class Definition
#' ============================================================================
#' 
#' This file defines the core EEG data structure as an S3 class in R.
#' The eeg class is the fundamental data container for all EEG analysis.
#' 
#' Author: Christos Dalamarinis
#' Date: Dec - 2025
#' ============================================================================
#'
#' Create a New EEG Object
#'
#' This function creates an S3 object of class 'eeg' that stores EEG data,
#' metadata, and processing history. This is the core data structure for
#' all eeganalysis functions.
#'
#' @param data Numeric matrix of EEG signal values
#'            Dimensions: rows = channels, columns = time points
#'            Units: typically microvolts (microV)
#'
#' @param channels Character vector of channel names (e.g., "Cz", "Pz", "Oz")
#'                 Length must match nrow(data)
#'
#' @param bads Character vector of channel names to mark as bad (optional).
#'            Default: NULL (no channels marked bad). Every name must also
#'            appear in \code{channels}. This is a single shared list, the
#'            same idea as MNE's \code{raw.info['bads']}: functions that
#'            touch channels (re-referencing, ICA, interpolation) should
#'            read/write \code{eeg$bads} instead of each one taking its own
#'            exclude list, so marking a channel bad once is enough for it
#'            to stay excluded everywhere downstream.
#'
#' @param annotations Data frame of bad time ranges (optional). Default:
#'            NULL, which creates a zero-row frame with columns
#'            \code{onset} (numeric, seconds), \code{duration} (numeric,
#'            seconds), \code{description} (character, e.g.
#'            \code{"BAD_muscle"}), and \code{channel} (character, NA
#'            unless the row is channel-specific). The time-domain sibling
#'            of \code{bads}: where \code{bads} disqualifies a whole
#'            channel, a row here disqualifies only a stretch of time, for
#'            downstream steps (epoching, ICA fitting) to check against.
#'            Written by \code{\link{annotate_amplitude}},
#'            \code{\link{annotate_muscle}}, \code{\link{annotate_nan}},
#'            and \code{\link{annotate_break}} (see R/annotations.R).
#'
#' @param sampling_rate Numeric value - sampling rate in Hz
#'                      Common values: 256, 512, 1024, 2048 Hz
#'
#' @param times Numeric vector of time points in seconds (optional)
#'             If NULL, automatically created from sampling_rate
#'             Length must equal ncol(data)
#'
#' @param events Data frame with event/trigger information (optional)
#'              Columns should include: onset (sample index), 
#'              type (trigger code), description
#'
#' @param metadata List containing experiment metadata (optional)
#'                Examples: subject_id, session, date, device, etc.
#'
#' @param reference Character string indicating reference scheme
#'                 Default: "original" (no change)
#'                 Options: "average", "linked_mastoids", "CMS/DRL", etc.
#'
#' @param preprocessing_history List tracking all preprocessing steps applied
#'                             Useful for reproducibility and auditing
#'
#' @param montage An object of class 'montage' giving channel scalp positions
#'                (optional). \code{NULL} until attached via
#'                \code{\link{set_montage}}. See \code{\link{create_montage}}.
#'
#' @param channel_types Character vector with the type of each channel
#'            (optional), same length and order as \code{channels}. One of
#'            \code{"eeg"}, \code{"eog"}, \code{"ecg"}, \code{"emg"},
#'            \code{"resp"}, \code{"gsr"}, \code{"temp"}, \code{"bio"},
#'            \code{"misc"}, or \code{"status"} per channel (see
#'            \code{\link{classify_channels}} for what each one means).
#'            Default: \code{NULL}, which assumes every channel is
#'            \code{"eeg"} except one literally named \code{"status"} or
#'            \code{"trigger"} - nothing else is guessed. Pass it when a
#'            channel's real role is known: at load time (e.g. a reader's own
#'            \code{eog =}/\code{misc =} arguments), or when a function
#'            rebuilds an eeg object and must carry a type forward that the
#'            channel's name alone would not reveal (e.g. a channel built by
#'            \code{\link{set_bipolar_reference}}, typed \code{"eog"} but
#'            named just \code{"VEOG"}).
#'
#' @return An object of class 'eeg' containing:
#'  \describe{
#'    \item{data}{Numeric matrix of EEG values (channels x time points)}
#'    \item{channels}{Character vector of channel names}
#'    \item{channel_types}{Character vector, same length and order as
#'      \code{channels} - see the \code{channel_types} argument above for
#'      the allowed values and how this field is filled in. Downstream code
#'      should read this field rather than re-deriving it from
#'      \code{channels}.}
#'    \item{bads}{Character vector of channel names marked as bad (empty
#'      character vector if none). A shared, persistent list - other
#'      functions should read/write this instead of taking their own
#'      per-call exclude list.}
#'    \item{annotations}{Data frame of bad time ranges, columns
#'      \code{onset, duration, description, channel} (zero rows if none).
#'      The time-domain sibling of \code{bads} - see
#'      \code{\link{annotate_amplitude}} and friends in R/annotations.R.}
#'    \item{sampling_rate}{Numeric sampling rate}
#'    \item{times}{Numeric time vector}
#'    \item{events}{Data frame with event information}
#'    \item{metadata}{List with experiment metadata}
#'    \item{reference}{Reference scheme used}
#'    \item{preprocessing_history}{Processing log}
#'    \item{montage}{Object of class 'montage' with channel scalp positions,
#'      or \code{NULL} if none has been attached yet}
#'  }
#'
#' @examples
#' \dontrun{
#'   # Create EEG object from raw data
#'   eeg <- new_eeg(
#'     data = eeg_matrix,
#'     channels = c("Cz", "Pz", "Oz"),
#'     sampling_rate = 2048,
#'     metadata = list(subject = "S01", date = "2025-01-15")
#'   )
#' }
#'
#' @export
new_eeg <- function(data,
                    channels,
                    sampling_rate,
                    times = NULL,
                    events = NULL,
                    metadata = NULL,
                    reference = "original",
                    preprocessing_history = NULL,
                    montage = NULL,
                    bads = NULL,
                    annotations = NULL,
                    channel_types = NULL) {
  
  # ========== INPUT VALIDATION ==========
  
  # Convert data to matrix if needed
  if (!is.matrix(data) && !is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Ensure data is numeric (double) type 16/02/2026
  storage.mode(data) <- "double"
  
  # Validate channel count matches data dimensions
  if (nrow(data) != length(channels)) {
    stop("ERROR: Number of channels (", length(channels), 
         ") does not match number of columns in data (", 
         nrow(data), ")")
  }
  
  # ========== CREATE TIME VECTOR ==========
  
  if (is.null(times)) {
    # Auto-generate time vector from sampling rate
    times <- (0:(ncol(data) - 1)) / sampling_rate
  } else {
    # Validate provided time vector
    if (length(times) != ncol(data)) {
      stop("ERROR: Length of times (", length(times), 
           ") does not match number of rows in data (", 
           ncol(data), ")")
    }
  }
  
  # ========== CLASSIFY CHANNELS ==========

  # NULL (the default) works the type out for every channel from its name -
  # see classify_channels(): only a channel literally named "status" or
  # "trigger" is special-cased, everything else defaults to "eeg". A caller
  # that already knows a channel's real role (a reader's own eog=/misc=
  # arguments, or a function rebuilding an eeg object that must carry a type
  # forward) hands it in here instead, and it is kept exactly as given.
  if (is.null(channel_types)) {
    channel_types <- classify_channels(as.character(channels))
  } else {
    channel_types <- as.character(channel_types)
    if (length(channel_types) != length(channels)) {
      stop("ERROR: Length of channel_types (", length(channel_types),
           ") does not match number of channels (", length(channels), ")")
    }
    unknown_types <- setdiff(unique(channel_types), .valid_channel_types())
    if (length(unknown_types) > 0) {
      stop("ERROR: 'channel_types' must be one of: ",
           paste(.valid_channel_types(), collapse = ", "), "; got: ",
           paste(unknown_types, collapse = ", "))
    }
  }

  # ========== VALIDATE BAD CHANNELS ==========

  if (is.null(bads)) {
    bads <- character(0)
  } else {
    bads <- as.character(bads)
    unknown_bads <- setdiff(bads, channels)
    if (length(unknown_bads) > 0) {
      stop("ERROR: 'bads' contains channel name(s) not found in 'channels': ",
           paste(unknown_bads, collapse = ", "))
    }
  }

  # ========== DEFAULT ANNOTATIONS ==========

  # Time-domain sibling of `bads` (see R/annotations.R): a zero-row frame
  # by default so annotate_*() functions never have to fall back on their
  # own .ensure_annotations() check for a NULL field.
  if (is.null(annotations)) {
    annotations <- .empty_annotations()
  }

  # ========== CREATE EVENTS DATAFRAME ==========
  
  if (is.null(events)) {
    # Create empty events data frame with proper structure
    events <- data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0)
    )
  }
  
  # ========== CREATE METADATA LIST ==========
  
  if (is.null(metadata)) {
    metadata <- list()
  }
  
  # ========== CREATE PREPROCESSING HISTORY ==========
  
  if (is.null(preprocessing_history)) {
    preprocessing_history <- list()
  }
  
  # ========== CONSTRUCT EEG OBJECT ==========
  
  # Create the S3 object with all components
  eeg_object <- structure(
    list(
      data = as.matrix(data),
      channels = as.character(channels),
      channel_types = channel_types,
      bads = bads,
      annotations = annotations,
      sampling_rate = as.numeric(sampling_rate),
      times = as.numeric(times),
      events = events,
      metadata = metadata,
      reference = as.character(reference),
      preprocessing_history = preprocessing_history,
      montage = montage
    ),
    class = "eeg"
  )
  
  return(eeg_object)
}

#' Valid Channel Types (internal)
#'
#' The complete set of values \code{new_eeg()}'s \code{channel_types}
#' argument accepts, and that \code{\link{set_channel_types}} accepts for
#' relabeling a channel after loading. Defined once here so every function
#' that validates or documents channel types stays in sync.
#'
#' \code{"eeg"} is a scalp electrode. \code{"status"} is the BioSemi
#' status/trigger channel. The rest name a channel's physiological role and
#' are never guessed - they must be stated explicitly, because that fact
#' comes from how the recording was set up, not from the channel's name:
#' \code{"eog"} (eye movement), \code{"ecg"} (heart), \code{"emg"} (muscle),
#' \code{"resp"} (respiration belt), \code{"gsr"} (skin conductance),
#' \code{"temp"} (temperature), \code{"misc"} (anything else not meant to be
#' analysed as a signal), or \code{"bio"} (any other physiological signal
#' with no more specific type here).
#'
#' @return Character vector of the allowed \code{channel_types} values.
#' @keywords internal
.valid_channel_types <- function() {
  c("eeg", "eog", "ecg", "emg", "resp", "gsr", "temp", "bio", "misc", "status")
}

#' Default Channel Types From Channel Names (internal)
#'
#' Works out a starting \code{channel_types} value for a vector of channel
#' names, used by \code{new_eeg()} only when \code{channel_types} is not
#' supplied. Every channel defaults to \code{"eeg"}, except one literally
#' named \code{"status"} or \code{"trigger"} (case-insensitive), which is
#' \code{"status"} - the same name-based rule MNE-Python's readers use for
#' their STIM channel (\code{stim_channel = "auto"}).
#'
#' Nothing else is guessed. A channel's real physiological role - eog, ecg,
#' emg, resp, gsr, temp, bio, or misc - is a fact about how the recording
#' was set up, not something derivable from its name, so it must be stated
#' explicitly: via a reader's own arguments at load time, or via
#' \code{\link{set_channel_types}} afterward.
#'
#' @param channels Character vector of channel names.
#' @return Character vector the same length as \code{channels}, with values
#'   \code{"eeg"} or \code{"status"}, in the same order as \code{channels}.
#' @keywords internal
classify_channels <- function(channels) {
  ch    <- channels
  types <- rep("eeg", length(ch))

  status_idx <- which(tolower(ch) %in% c("status", "trigger"))
  types[status_idx] <- "status"

  types
}

#' Print Method for EEG Objects
#'
#' Custom print method that displays a formatted summary of an EEG object.
#' Shows key information about the recording without overwhelming the console.
#'
#' @param x An object of class 'eeg'
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns x (following R print method convention)
#'
#' @examples
#' \dontrun{
#'   eeg <- read_biosemi("data.bdf")
#'   print(eeg)  # Calls this method automatically
#' }
#'
#' @export
print.eeg <- function(x, ...) {
  
  # ========== HEADER ==========
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("EEG Object Summary\n")
  cat(strrep("=", 70), "\n\n")
  
  # ========== BASIC INFORMATION ==========
  cat("RECORDING INFORMATION:\n")
  cat("  Channels:        ", length(x$channels), "\n")
  
  # Show first 5 channel names, then "... (+N more)" if necessary
  channel_display <- paste(x$channels[1:min(5, length(x$channels))], collapse = ", ")
  if (length(x$channels) > 5) {
    channel_display <- paste0(channel_display, " ... (+", length(x$channels) - 5, " more)")
  }
  cat("    List:          ", channel_display, "\n")
  
  cat("  Time points:     ", ncol(x$data), "\n")
  cat("  Duration:        ", sprintf("%.2f", ncol(x$data) / x$sampling_rate), " seconds\n")
  cat("  Sampling rate:   ", x$sampling_rate, " Hz\n")
  
  # ========== DATA STATISTICS ==========
  # Stats are scoped to EEG channels only, excluding status, external
  # (EXG/EOG/ECG/EMG/GSR/etc.), and bad channels. Classification is read
  # from x$channel_types, computed once by classify_channels() at
  # construction in new_eeg() - not re-derived here - and bad channels are
  # read from x$bads, the shared exclude list other functions should also
  # use instead of taking their own per-call exclude argument.

  .eeg_idx <- which(x$channel_types == "eeg" & !(x$channels %in% x$bads))

  .eeg_data  <- x$data[.eeg_idx, , drop = FALSE]
  
  cat("\nDATA STATISTICS (EEG channels only):\n")
  cat("  Amplitude range: ", round(min(.eeg_data), 2), " to ",
      round(max(.eeg_data), 2), " microV\n", sep = "")
  cat("  Mean amplitude:  ", round(mean(.eeg_data), 2), " microV\n")
  cat("  Std deviation:   ", round(sd(.eeg_data), 2), " microV\n")
  
  # ========== REFERENCE INFORMATION ==========
  cat("\nREFERENCE:\n")
  cat("  Scheme:          ", x$reference, "\n")

  # ========== BAD CHANNELS ==========
  cat("\nBAD CHANNELS:\n")
  if (length(x$bads) > 0) {
    cat("  Marked bad:      ", paste(x$bads, collapse = ", "),
        " (", length(x$bads), ")\n", sep = "")
  } else {
    cat("  Marked bad:       None\n")
  }

  # ========== ANNOTATIONS ==========
  cat("\nANNOTATIONS:\n")
  cat("  Total annotations:", nrow(x$annotations), "\n")
  if (nrow(x$annotations) > 0) {
    annotation_types <- unique(x$annotations$description)
    cat("  Types:            ", paste(annotation_types, collapse = ", "), "\n")
  }

  # ========== EVENT INFORMATION ==========
  cat("\nEVENTS:\n")
  cat("  Total events:    ", nrow(x$events), "\n")
  
  if (nrow(x$events) > 0) {
    event_types <- unique(x$events$type)
    cat("  Event types:     ", paste(event_types, collapse = ", "), "\n")
  }
  
  # ========== METADATA DISPLAY ==========
  if (length(x$metadata) > 0) {
    cat("\nMETADATA:\n")
    for (key in names(x$metadata)) {
      val <- x$metadata[[key]]
      
      # Ensure val is always coercible to character
      if (!is.atomic(val)) {
        val <- paste0("[", class(val)[1], " object]")
      } else if (is.character(val) && length(val) > 1) {
        val <- paste(val[1:min(3, length(val))], collapse = ", ")
      } else {
        val <- paste(as.character(val), collapse = ", ")
      }
      
      if (nchar(val) > 40) {
        val <- paste0(substr(val, 1, 37), "...")
      }
      
      cat(" ", key, ": ", val, "\n", sep = "")
    }
  }
  
  # ========== PREPROCESSING HISTORY ==========
  if (length(x$preprocessing_history) > 0) {
    cat("\nPREPROCESSING HISTORY:\n")
    for (i in seq_along(x$preprocessing_history)) {
      cat("  ", i, ". ", x$preprocessing_history[[i]], "\n", sep = "")
    }
  }
  
  # ========== FOOTER ==========
  cat(strrep("=", 70), "\n\n")
  
  # Return invisibly (standard R convention)
  invisible(x)
}
