#' ============================================================================
#' BioSemi BDF File Reader Functions
#' ============================================================================
#'
#' This module provides functions to read and import EEG data recorded from
#' BioSemi systems (BDF format). BDF (BioSemi Data Format) is a binary format
#' for storing EEG data with metadata and event markers.
#'
#' BioSemi systems typically use:
#'  - CMS/DRL reference scheme (Common Mode Sense / Driven Right Leg)
#'  - High sampling rates (256-2048 Hz typical)
#'  - 64, 128, or 256 electrode channels
#'  - Status channel for event markers/triggers
#'
#' Author: Christos Dalamarinis
#' Date: 2025-03-15
#' ============================================================================

#' Read BioSemi BDF File
#'
#' This is the main function to import EEG data from a BioSemi .bdf file.
#' It reads the binary file, extracts all EEG channels and metadata,
#' parses event triggers, and returns a properly structured eeg object
#' ready for analysis.
#'
#' @param file_path Character string with full path to the .bdf file
#'                  Example: "/Users/user/data/subject_01.bdf"
#'
#' @details
#' The function performs the following steps:
#' 1. Validates the file exists and has .bdf extension
#' 2. Reads the BDF header using edfReader::readEdfHeader()
#' 3. Reads all signal channels using edfReader::readEdfSignals()
#' 4. Extracts channel names, sampling rate, and recording metadata
#' 5. Separates EEG channels from status/trigger channel
#' 6. Extracts and parses event markers from status channel
#' 7. Returns an eeg object with all data and metadata
#'
#' The status channel is typically:
#'  - Named "Status" in BioSemi recordings
#'  - Contains trigger codes corresponding to experimental events
#'  - Bit values encode different triggers (0, 1, 2, 4, 8, etc.)
#'  - Used to epoch data around events of interest
#'
#' @return An object of class 'eeg' containing:
#'  \describe{
#'    \item{data}{Matrix of EEG signal values (channels × time points)}
#'    \item{channels}{Character vector of EEG channel names}
#'    \item{sampling_rate}{Numeric - sampling rate in Hz}
#'    \item{times}{Numeric vector - time points in seconds}
#'    \item{events}{Data frame - event markers (onset, type, description)}
#'    \item{metadata}{List - recording information and file metadata}
#'    \item{reference}{Character - reference scheme (Biosemi CMS/DRL)}
#'    \item{preprocessing_history}{List - marks file as imported}
#'  }
#'
#' @examples
#' \dontrun{
#'   # Load a BioSemi file
#'   eeg <- read_biosemi("~/data/experiment/subject_01.bdf")
#'   
#'   # View summary information
#'   print(eeg)
#'   
#'   # Access components
#'   head(eeg$data)           # First rows of EEG data
#'   eeg$channels             # Channel names
#'   eeg$sampling_rate        # Sampling rate
#'   eeg$events               # Event table
#'   eeg$metadata             # Recording metadata
#' }
#'
#' @importFrom edfReader readEdfHeader readEdfSignals
#' @export
read_biosemi <- function(file_path) {
  
  # ========== INPUT VALIDATION ==========
  
  # Check file exists
  if (!file.exists(file_path)) {
    stop("ERROR: File not found at path: ", file_path)
  }
  
  # Check file has .bdf extension (case insensitive)
  if (!grepl("\\.bdf$", file_path, ignore.case = TRUE)) {
    warning("WARNING: File does not have .bdf extension. ",
            "Expected BioSemi BDF format. Attempting to read anyway...")
  }
  
  # ========== READ BDF FILE HEADER ==========
  
  cat("[1/4] Reading BDF header...")  # Progress indicator
  
  # Use edfReader to read the BDF file header
  # The header contains metadata about channels, sampling rate, duration, etc.
  header <- tryCatch({
    edfReader::readEdfHeader(file_path)
  }, error = function(e) {
    stop("ERROR reading BDF file header: ", e$message, 
         "\nCheck that the file is a valid BioSemi BDF file.")
  })
  
  cat(" Done!\n")
  
  # ========== READ BDF SIGNAL DATA ==========
  
  cat("[2/4] Reading EEG signals...")  # Progress indicator
  
  # Read all signal channels from the BDF file
  # Each signal contains the time series data for one channel
  signals <- tryCatch({
    edfReader::readEdfSignals(header)
  }, error = function(e) {
    stop("ERROR reading BDF signals: ", e$message)
  })
  
  cat(" Done!\n")
  
  # ========== EXTRACT CHANNEL INFORMATION ==========
  
  # Get channel names from header
  channels_raw <- header$sHeaders$label
  channels_clean <- trimws(channels_raw)  # Remove leading/trailing whitespace
  
  # Get sampling rate from header (same for all channels in BioSemi)
  sampling_rate <- header$sHeaders$sRate[1]
  
  cat("  Found ", length(channels_clean), " channels\n", sep = "")
  cat("  Sampling rate: ", sampling_rate, " Hz\n", sep = "")
  
  # ========== CONVERT SIGNALS TO MATRIX ==========
  
  # Count number of channels and samples
  n_channels <- length(signals)
  
  #Check that all channels have the same number of samples
  sample_lengths <- c()
  for (i in 1:n_channels) {
    sample_lengths[i] <- length(signals[[i]]$signal)
  }
  
  if (length(unique(sample_lengths)) > 1) {
    stop("ERROR: Channels have different numbers of samples!\n",
         "  Sample lengths: ", paste(unique(sample_lengths), collapse = ", "), "\n",
         "  This file may be corrupted or have variable sampling rates.")
  }
  
  n_samples <- sample_lengths[1]
  
  # Check that all channels have the same number of samples
  sample_lengths <- c()
  for (i in 1:n_channels) {
    sample_lengths[i] <- length(signals[[i]]$signal)
  }
  
  # DIAGNOSTIC: Show sample lengths
  cat("  Sample lengths per channel:\n")
  cat("    Min: ", min(sample_lengths), "\n", sep = "")
  cat("    Max: ", max(sample_lengths), "\n", sep = "")
  cat("    Unique values: ", paste(unique(sample_lengths), collapse = ", "), "\n", sep = "")
  
  if (length(unique(sample_lengths)) > 1) {
    stop("ERROR: Channels have different numbers of samples!\n",
         "  Sample lengths: ", paste(unique(sample_lengths), collapse = ", "), "\n",
         "  This file may be corrupted or have variable sampling rates.")
  }
  
  # Create matrix to hold all channel data
  # Rows = channels, Columns = time points
  eeg_matrix <- matrix(0, nrow = n_channels, ncol = n_samples)
  
  # Fill matrix with signal data from each channel
  # edfReader returns signals as a list, each element is a list with $signal
  # Fill matrix with signal data from each channel
  for (i in 1:n_channels) {
    signal_data <- signals[[i]]$signal
    eeg_matrix[i, ] <- signal_data
  }
  
  # ========== IDENTIFY AND SEPARATE STATUS CHANNEL ==========
  
  # BioSemi recordings include a status channel for event markers/triggers
  # This is typically the last channel and named "Status"
  
  # Look for channel named "Status"
  status_channel_idx <- which(tolower(channels_clean) == "status")
  
  # If not found, assume last channel is status
  if (length(status_channel_idx) == 0) {
    status_channel_idx <- n_channels
    cat("  Status channel: Using last channel (Ch", status_channel_idx, ")\n", sep = "")
  } else {
    cat("  Status channel: ", channels_clean[status_channel_idx[1]], "\n")
  }
  
  # Extract status channel for event parsing
  status_signal <- as.numeric(eeg_matrix[status_channel_idx[1], ])
  
  # Separate EEG channels from status channel
  eeg_channels_idx <- setdiff(1:n_channels, status_channel_idx)
  eeg_data <- eeg_matrix[eeg_channels_idx, , drop = FALSE]
  channels_final <- channels_clean[eeg_channels_idx]
  
  cat("  EEG channels: ", length(eeg_channels_idx), "\n", sep = "")
  
  # ========== CREATE TIME VECTOR ==========
  
  cat("[3/4] Creating time vector...")  # Progress indicator
  
  # Time vector in seconds: 0, 1/srate, 2/srate, ..., (n_samples-1)/srate
  times <- (0:(n_samples - 1)) / sampling_rate
  
  cat(" Done!\n")
  
  # ========== EXTRACT EVENTS FROM STATUS CHANNEL ==========
  
  cat("[4/4] Extracting events from status channel...")  # Progress indicator
  
  # Parse event markers from the status signal
  # The status channel contains trigger codes when events occur
  events_df <- extract_biosemi_events(status_signal, sampling_rate)
  
  cat(" Done!\n")
  cat("  Events found: ", nrow(events_df), "\n", sep = "")
  
  # ========== PREPARE METADATA ==========
  
  # Collect all metadata from the BDF file header
  metadata <- list(
    # File information
    file_path = file_path,
    file_name = basename(file_path),
    
    # Recording date/time
    date_recorded = tryCatch(header$startDate, error = function(e) NA),
    time_recorded = tryCatch(header$startTime, error = function(e) NA),
    
    # Device and recording information
    device = "BioSemi",
    n_channels_recorded = n_channels,
    n_eeg_channels = length(eeg_channels_idx),
    status_channel_name = ifelse(length(status_channel_idx) > 0, 
                                 channels_raw[status_channel_idx[1]], 
                                 NA),
    
    # Reference information
    reference_scheme = "Biosemi CMS/DRL (Common Mode Sense / Driven Right Leg)",
    
    # Recording duration
    duration_seconds = max(times),
    n_samples = n_samples,
    n_events = nrow(events_df)
  )
  
  # ========== CREATE EEG OBJECT ==========
  
  # Construct the eeg object using new_eeg() function
  eeg_obj <- new_eeg(
    data = eeg_data,
    channels = channels_final,
    sampling_rate = sampling_rate,
    times = times,
    events = events_df,
    metadata = metadata,
    reference = "Biosemi CMS/DRL",
    preprocessing_history = list(
      paste0("Imported from BDF file: ", basename(file_path))
    )
  )
  
  # ========== RETURN RESULT ==========
  
  cat("\n✓ Successfully imported BioSemi BDF file!\n\n")
  
  return(eeg_obj)
}


#' Extract Events from BioSemi Status Channel
#'
#' This internal function parses event markers from the BioSemi status channel.
#' BioSemi systems encode experiment triggers as numeric codes in the status
#' channel. This function identifies when these values change (indicating events)
#' and creates an event table.
#'
#' @param status_signal Numeric vector containing the status channel values
#'                      at each time point
#'
#' @param sampling_rate Numeric - sampling rate in Hz
#'                      Used to convert sample indices to time (seconds)
#'
#' @return Data frame with event information:
#'  \describe{
#'    \item{onset}{Integer - event onset in samples (0-indexed)}
#'    \item{onset_time}{Numeric - event onset in seconds}
#'    \item{type}{Character - trigger code as string (e.g., "255")}
#'    \item{description}{Character - human-readable description}
#'  }
#'
#' @details
#' Events are detected by identifying non-zero values in the status channel.
#' BioSemi systems can encode multiple triggers using bit flags:
#'  - Trigger 1 = value 1 (2^0)
#'  - Trigger 2 = value 2 (2^1)
#'  - Trigger 3 = value 4 (2^2)
#'  - Trigger 4 = value 8 (2^3)
#'  - etc.
#'
#' Combined triggers are the sum (e.g., triggers 1+2 = value 3)
#'
#' @keywords internal
extract_biosemi_events <- function(status_signal, sampling_rate) {
  
  # ========== HANDLE EMPTY STATUS ==========
  
  # If status signal is NULL or empty, return empty events data frame
  if (is.null(status_signal) || length(status_signal) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0)
    ))
  }
  
  # ========== DETECT EVENT TRANSITIONS ==========
  
  # Convert to numeric (ensure numeric type)
  status_numeric <- as.numeric(status_signal)
  
  # Find where status has non-zero values (events)
  event_idx <- which(status_numeric != 0)
  
  # ========== HANDLE CASE WITH NO EVENTS ==========
  
  if (length(event_idx) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0)
    ))
  }
  
  # ========== CREATE EVENTS DATA FRAME ==========
  
  # Get unique trigger values by finding where transitions occur
  # Only keep first occurrence of each trigger value
  status_changes <- c(TRUE, diff(status_numeric) != 0)
  event_idx_unique <- which(status_changes & status_numeric != 0)
  
  # Create the events data frame
  events_df <- data.frame(
    onset = event_idx_unique,  # Sample index (1-indexed in R)
    onset_time = (event_idx_unique - 1) / sampling_rate,  # Convert to seconds (0-indexed)
    type = as.character(status_numeric[event_idx_unique]),  # Trigger code as string
    description = paste0("Trigger code: ", 
                         as.character(status_numeric[event_idx_unique])),
    stringsAsFactors = FALSE
  )
  
  return(events_df)
}


#' Generate Quality Summary of BioSemi Import
#'
#' After importing a BioSemi file, this function provides a comprehensive
#' quality report including data statistics, channel information, event counts,
#' and any issues detected during import.
#'
#' @param eeg An eeg object imported via read_biosemi()
#'
#' @return Invisibly returns a list with quality metrics:
#'  \describe{
#'    \item{n_channels}{Number of EEG channels}
#'    \item{n_samples}{Total time points per channel}
#'    \item{duration_sec}{Recording duration in seconds}
#'    \item{sampling_rate}{Sampling rate in Hz}
#'    \item{channel_vars}{Variance per channel}
#'    \item{n_events}{Number of events detected}
#'    \item{event_types}{Unique trigger codes}
#'  }
#'
#' @details
#' This function prints a formatted report and silently returns metrics.
#' The report includes:
#'  - Channel count and data dimensions
#'  - Data quality statistics (variance, amplitude range)
#'  - Event/trigger information
#'  - Potential issues (no events, flat channels, etc.)
#'
#' @examples
#' \dontrun{
#'   eeg <- read_biosemi("data.bdf")
#'   quality_metrics <- summarize_biosemi_import(eeg)
#' }
#'
#' @export
summarize_biosemi_import <- function(eeg) {
  
  # ========== HEADER ==========
  cat("\n")
  cat(strrep("=", 75), "\n")
  cat("BioSemi Import Quality Summary\n")
  cat(strrep("=", 75), "\n\n")
  
  # ========== DATA STATISTICS ==========
  cat("DATA STATISTICS:\n")
  cat("  Channels:        ", length(eeg$channels), "\n")
  cat("  Time points:     ", ncol(eeg$data), "\n")
  cat("  Duration:        ", round(max(eeg$times), 3), " seconds\n")
  cat("  Sampling rate:   ", eeg$sampling_rate, " Hz\n")
  cat("  Data type:       numeric matrix\n")
  
  # ========== CHANNEL INFORMATION ==========
  cat("\nCHANNEL INFORMATION:\n")
  
  # Display first 8 channels, then "... (+N more)" if needed
  channels_display <- paste(eeg$channels[1:min(8, length(eeg$channels))], collapse = ", ")
  if (length(eeg$channels) > 8) {
    channels_display <- paste0(channels_display, " ... (+", 
                               length(eeg$channels) - 8, " more)")
  }
  cat("  Channels: ", channels_display, "\n")
  
  # ========== DATA QUALITY ANALYSIS ==========
  cat("\nDATA QUALITY METRICS:\n")
  
  # Calculate variance for each channel
  channel_vars <- apply(eeg$data, 1, var, na.rm = TRUE)
  channel_means <- apply(eeg$data, 1, mean, na.rm = TRUE)
  channel_sds <- apply(eeg$data, 1, sd, na.rm = TRUE)
  
  # Report amplitude statistics
  cat("  Amplitude range:       [", round(min(eeg$data), 2), ", ",
      round(max(eeg$data), 2), "] µV\n", sep = "")
  cat("  Mean across channels:  ", round(mean(channel_means), 2), " µV\n")
  cat("  Variance range:        [", round(min(channel_vars), 2), ", ",
      round(max(channel_vars), 2), "]\n", sep = "")
  
  # Identify potentially problematic channels
  flat_channels <- which(channel_vars < 0.01)
  if (length(flat_channels) > 0) {
    cat("  WARNING: ", length(flat_channels), " potentially flat channels: ",
        paste(eeg$channels[flat_channels], collapse = ", "), "\n", sep = "")
  }
  
  high_noise_channels <- which(channel_vars > quantile(channel_vars, 0.95))
  if (length(high_noise_channels) > 0) {
    cat("  NOTE: ", length(high_noise_channels), " channels with high variance\n", sep = "")
  }
  
  # ========== EVENT INFORMATION ==========
  cat("\nEVENT INFORMATION:\n")
  cat("  Total events:    ", nrow(eeg$events), "\n")
  
  if (nrow(eeg$events) > 0) {
    event_types <- unique(eeg$events$type)
    cat("  Event types:     ", paste(event_types, collapse = ", "), "\n")
    cat("  Event counts:\n")
    
    # Count occurrences of each trigger type
    for (et in event_types) {
      count <- sum(eeg$events$type == et)
      cat("    Trigger ", et, ": ", count, " occurrences\n", sep = "")
    }
  } else {
    cat("  WARNING: No events detected in status channel!\n")
    cat("    Check that triggers were properly recorded.\n")
  }
  
  # ========== REFERENCE INFORMATION ==========
  cat("\nREFERENCE INFORMATION:\n")
  cat("  Reference:      ", eeg$reference, "\n")
  cat("  Description:    CMS (Common Mode Sense) at midline,\n")
  cat("                  DRL (Driven Right Leg) at right leg\n")
  
  # ========== FOOTER ==========
  cat(strrep("=", 75), "\n\n")
  
  # ========== RETURN METRICS INVISIBLY ==========
  
  # Create a list of quality metrics
  quality_metrics <- list(
    n_channels = length(eeg$channels),
    n_samples = nrow(eeg$data),
    duration_sec = max(eeg$times),
    sampling_rate = eeg$sampling_rate,
    channel_vars = channel_vars,
    n_events = nrow(eeg$events),
    event_types = if (nrow(eeg$events) > 0) unique(eeg$events$type) else character(0),
    flat_channels = if (length(flat_channels) > 0) eeg$channels[flat_channels] else character(0),
    high_noise_channels = if (length(high_noise_channels) > 0) 
      eeg$channels[high_noise_channels] else character(0)
  )
  
  # Return invisibly (won't print when called without assignment)
  invisible(quality_metrics)
}

# End of read_biosemi.r