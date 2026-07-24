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
#' @return An object of class 'eeg' containing:
#'  \describe{
#'    \item{data}{Numeric matrix of EEG values (channels x time points)}
#'    \item{channels}{Character vector of channel names}
#'    \item{channel_types}{Character vector, same length and order as
#'      \code{channels}, classifying each channel as \code{"eeg"},
#'      \code{"external"} (EOG/ECG/EMG/GSR/etc.), or \code{"status"}.
#'      Computed once at construction via \code{classify_channels()} - see
#'      that function for the classification rules - so downstream code
#'      should read this field rather than re-deriving it from
#'      \code{channels}.}
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
                    montage = NULL) {
  
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

  # Computed once here so downstream functions (print.eeg(), etc.) read
  # channel_types instead of re-deriving the eeg/external/status split.
  channel_types <- classify_channels(as.character(channels))

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

#' Classify Channels as EEG, External, or Status (internal)
#'
#' Computes the eeg/external/status classification for a vector of channel
#' names, once, so that \code{new_eeg()} can store the result on the
#' \code{eeg} object as \code{channel_types} instead of every downstream
#' function re-deriving it. Combines the BioSemi status-channel check with
#' the same two-pass external-channel detection used elsewhere in the
#' package: pass 1 is the electrode-database lookup
#' (\code{detect_external_channels()}, see R/setexchannels.R), pass 2 is a
#' regex fallback for channels renamed with the original name kept in
#' parentheses (e.g. \code{"MASTOID LEFT (EXG5)"}).
#'
#' This only distinguishes broad channel categories (eeg/external/status),
#' not the physiological sub-type of an external channel (EOG vs ECG vs
#' EMG vs GSR). BioSemi's EXG1-EXG8 ports are generic in the electrode
#' database - which physiological signal is wired to a given port is a
#' fact about how that specific recording session was set up, not
#' something derivable from the channel name. That labeling is handled
#' separately by \code{identify_external_channels()} /
#' \code{apply_external_labels()}.
#'
#' @param channels Character vector of channel names.
#' @return Character vector the same length as \code{channels}, with values
#'   \code{"eeg"}, \code{"external"}, or \code{"status"}, in the same order
#'   as \code{channels}.
#' @keywords internal
classify_channels <- function(channels) {
  ch    <- channels
  types <- rep("eeg", length(ch))

  status_idx <- which(grepl("^status$", tolower(ch)))
  types[status_idx] <- "status"

  exg_pattern <- "\\((EXG[1-8]|GSR[12]|Plet|Temp|Resp|Erg[12])\\)"
  exg_pass1 <- detect_external_channels(ch)
  exg_pass2 <- ch[grepl(exg_pattern, ch, ignore.case = TRUE)]
  exg_idx   <- setdiff(which(ch %in% unique(c(exg_pass1, exg_pass2))), status_idx)
  types[exg_idx] <- "external"

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
  # Stats are scoped to EEG channels only, excluding status and external
  # (EXG/EOG/ECG/EMG/GSR/etc.) channels. Classification is read from
  # x$channel_types, computed once by classify_channels() at construction
  # in new_eeg() - not re-derived here.

  .eeg_idx <- which(x$channel_types == "eeg")

  .eeg_data  <- x$data[.eeg_idx, , drop = FALSE]
  
  cat("\nDATA STATISTICS (EEG channels only):\n")
  cat("  Amplitude range: ", round(min(.eeg_data), 2), " to ",
      round(max(.eeg_data), 2), " microV\n", sep = "")
  cat("  Mean amplitude:  ", round(mean(.eeg_data), 2), " microV\n")
  cat("  Std deviation:   ", round(sd(.eeg_data), 2), " microV\n")
  
  # ========== REFERENCE INFORMATION ==========
  cat("\nREFERENCE:\n")
  cat("  Scheme:          ", x$reference, "\n")
  
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
