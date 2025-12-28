#' ============================================================================
#' BioSemi Channel Information and Inspection Utilities
#' ============================================================================
#'
#' This module provides functions to inspect and analyze BioSemi BDF files
#' for channel information without needing to import the full dataset.
#' Useful for quick file preview, metadata extraction, and electrode layout verification.
#'
#' Author: Christos Dalamarinis
#' Date: 2025
#' ============================================================================

#' Get Standard Electrode Position Database
#'
#' Returns a comprehensive lookup table of standard electrode positions
#' used in EEG research (10-20, 10-10 systems) and common external channels
#' (EOG, ECG, EMG, GSR).
#'
#' @return Named list where names are electrode names and values are
#' position information lists with: position_name, position_type, region
#'
#' @details
#' This database includes:
#' - Standard 10-20 system: Fp, F, C, P, O, T positions
#' - Extended 10-10 system: Additional positions for higher density
#' - External channels: EOG (eyes), ECG (heart), EMG (muscle), GSR (skin)
#' - Reference electrodes: Cz, CPz, M1/M2 (mastoids)
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' electrode_db <- get_electrode_database()
#' electrode_db$Cz  # Get info about Cz electrode
#' }
#'
#' @export
get_electrode_database <- function() {
  
  # ========== STANDARD 10-20 ELECTRODE SYSTEM ==========
  # Primary electrode positions
  
  electrode_db <- list(
    # ===== PREFRONTAL (Fp) =====
    "Fp1" = list(
      position_name = "Prefrontal left",
      position_type = "Prefrontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "Fpz" = list(
      position_name = "Prefrontal midline",
      position_type = "Prefrontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "Fp2" = list(
      position_name = "Prefrontal right",
      position_type = "Prefrontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    
    # ===== ANTERIOR FRONTAL (AF) - 10-10 System =====
    "AF7" = list(
      position_name = "Anterior frontal left",
      position_type = "Anterior Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "AF3" = list(
      position_name = "Anterior frontal left-center",
      position_type = "Anterior Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "AFz" = list(
      position_name = "Anterior frontal midline",
      position_type = "Anterior Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "AF4" = list(
      position_name = "Anterior frontal right-center",
      position_type = "Anterior Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "AF8" = list(
      position_name = "Anterior frontal right",
      position_type = "Anterior Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    
    # ===== FRONTAL (F) =====
    "F7" = list(
      position_name = "Frontal left temporal",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "F5" = list(
      position_name = "Frontal left-center",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "F3" = list(
      position_name = "Frontal left",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "F1" = list(
      position_name = "Frontal left-midline",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "Fz" = list(
      position_name = "Frontal midline",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "F2" = list(
      position_name = "Frontal right-midline",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "F4" = list(
      position_name = "Frontal right",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    "F6" = list(
      position_name = "Frontal right-center",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-10"
    ),
    "F8" = list(
      position_name = "Frontal right temporal",
      position_type = "Frontal",
      region = "Frontal",
      standard_system = "10-20"
    ),
    
    # ===== FRONTO-TEMPORAL (FT) - 10-10 System =====
    "FT7" = list(
      position_name = "Fronto-temporal left",
      position_type = "Fronto-Temporal",
      region = "Temporal",
      standard_system = "10-10"
    ),
    "FT8" = list(
      position_name = "Fronto-temporal right",
      position_type = "Fronto-Temporal",
      region = "Temporal",
      standard_system = "10-10"
    ),
    
    # ===== FRONTO-CENTRAL (FC) - 10-10 System =====
    "FC5" = list(
      position_name = "Fronto-central left",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FC3" = list(
      position_name = "Fronto-central left-center",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FC1" = list(
      position_name = "Fronto-central left-midline",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FCz" = list(
      position_name = "Fronto-central midline",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FC2" = list(
      position_name = "Fronto-central right-midline",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FC4" = list(
      position_name = "Fronto-central right-center",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "FC6" = list(
      position_name = "Fronto-central right",
      position_type = "Fronto-Central",
      region = "Central",
      standard_system = "10-10"
    ),
    
    # ===== TEMPORAL (T) =====
    "T7" = list(
      position_name = "Temporal left",
      position_type = "Temporal",
      region = "Temporal",
      standard_system = "10-20"
    ),
    "T8" = list(
      position_name = "Temporal right",
      position_type = "Temporal",
      region = "Temporal",
      standard_system = "10-20"
    ),
    
    # ===== CENTRAL (C) =====
    "C5" = list(
      position_name = "Central left",
      position_type = "Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "C3" = list(
      position_name = "Central left",
      position_type = "Central",
      region = "Central",
      standard_system = "10-20"
    ),
    "C1" = list(
      position_name = "Central left-midline",
      position_type = "Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "Cz" = list(
      position_name = "Central midline",
      position_type = "Central",
      region = "Central",
      standard_system = "10-20"
    ),
    "C2" = list(
      position_name = "Central right-midline",
      position_type = "Central",
      region = "Central",
      standard_system = "10-10"
    ),
    "C4" = list(
      position_name = "Central right",
      position_type = "Central",
      region = "Central",
      standard_system = "10-20"
    ),
    "C6" = list(
      position_name = "Central right",
      position_type = "Central",
      region = "Central",
      standard_system = "10-10"
    ),
    
    # ===== CENTRO-TEMPORAL (CT/TP) - 10-10 System =====
    "TP7" = list(
      position_name = "Centro-temporal left",
      position_type = "Centro-Temporal",
      region = "Temporal",
      standard_system = "10-10"
    ),
    "TP8" = list(
      position_name = "Centro-temporal right",
      position_type = "Centro-Temporal",
      region = "Temporal",
      standard_system = "10-10"
    ),
    
    # ===== CENTRO-PARIETAL (CP) - 10-10 System =====
    "CP5" = list(
      position_name = "Centro-parietal left",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CP3" = list(
      position_name = "Centro-parietal left-center",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CP1" = list(
      position_name = "Centro-parietal left-midline",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CPz" = list(
      position_name = "Centro-parietal midline",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CP2" = list(
      position_name = "Centro-parietal right-midline",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CP4" = list(
      position_name = "Centro-parietal right-center",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "CP6" = list(
      position_name = "Centro-parietal right",
      position_type = "Centro-Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    
    # ===== PARIETAL (P) =====
    "P7" = list(
      position_name = "Parietal left temporal",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-20"
    ),
    "P5" = list(
      position_name = "Parietal left-center",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "P3" = list(
      position_name = "Parietal left",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-20"
    ),
    "P1" = list(
      position_name = "Parietal left-midline",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "Pz" = list(
      position_name = "Parietal midline",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-20"
    ),
    "P2" = list(
      position_name = "Parietal right-midline",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "P4" = list(
      position_name = "Parietal right",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-20"
    ),
    "P6" = list(
      position_name = "Parietal right-center",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-10"
    ),
    "P8" = list(
      position_name = "Parietal right temporal",
      position_type = "Parietal",
      region = "Parietal",
      standard_system = "10-20"
    ),
    
    # ===== PARIETO-OCCIPITAL (PO) - 10-10 System =====
    "PO7" = list(
      position_name = "Parieto-occipital left",
      position_type = "Parieto-Occipital",
      region = "Occipital",
      standard_system = "10-10"
    ),
    "PO3" = list(
      position_name = "Parieto-occipital left-center",
      position_type = "Parieto-Occipital",
      region = "Occipital",
      standard_system = "10-10"
    ),
    "POz" = list(
      position_name = "Parieto-occipital midline",
      position_type = "Parieto-Occipital",
      region = "Occipital",
      standard_system = "10-10"
    ),
    "PO4" = list(
      position_name = "Parieto-occipital right-center",
      position_type = "Parieto-Occipital",
      region = "Occipital",
      standard_system = "10-10"
    ),
    "PO8" = list(
      position_name = "Parieto-occipital right",
      position_type = "Parieto-Occipital",
      region = "Occipital",
      standard_system = "10-10"
    ),
    
    # ===== OCCIPITAL (O) =====
    "O1" = list(
      position_name = "Occipital left",
      position_type = "Occipital",
      region = "Occipital",
      standard_system = "10-20"
    ),
    "Oz" = list(
      position_name = "Occipital midline",
      position_type = "Occipital",
      region = "Occipital",
      standard_system = "10-20"
    ),
    "O2" = list(
      position_name = "Occipital right",
      position_type = "Occipital",
      region = "Occipital",
      standard_system = "10-20"
    ),
    
    # ===== MASTOID REFERENCE (M) =====
    "M1" = list(
      position_name = "Mastoid left",
      position_type = "Mastoid Reference",
      region = "Reference",
      standard_system = "Reference"
    ),
    "M2" = list(
      position_name = "Mastoid right",
      position_type = "Mastoid Reference",
      region = "Reference",
      standard_system = "Reference"
    ),
    
    # ===== EXTERNAL CHANNELS - EOG (Electrooculogram) =====
    "EOG_L" = list(
      position_name = "Electrooculogram left",
      position_type = "EOG",
      region = "External",
      standard_system = "EOG"
    ),
    "EOG_R" = list(
      position_name = "Electrooculogram right",
      position_type = "EOG",
      region = "External",
      standard_system = "EOG"
    ),
    "HEOG_L" = list(
      position_name = "Horizontal EOG left",
      position_type = "EOG",
      region = "External",
      standard_system = "EOG"
    ),
    "HEOG_R" = list(
      position_name = "Horizontal EOG right",
      position_type = "EOG",
      region = "External",
      standard_system = "EOG"
    ),
    "VEOG" = list(
      position_name = "Vertical EOG",
      position_type = "EOG",
      region = "External",
      standard_system = "EOG"
    ),
    
    # ===== EXTERNAL CHANNELS - ECG (Electrocardiogram) =====
    "ECG" = list(
      position_name = "Electrocardiogram (heart rate)",
      position_type = "ECG",
      region = "External",
      standard_system = "ECG"
    ),
    "ECG1" = list(
      position_name = "Electrocardiogram channel 1",
      position_type = "ECG",
      region = "External",
      standard_system = "ECG"
    ),
    "ECG2" = list(
      position_name = "Electrocardiogram channel 2",
      position_type = "ECG",
      region = "External",
      standard_system = "ECG"
    ),
    
    # ===== EXTERNAL CHANNELS - EMG (Electromyogram) =====
    "EMG" = list(
      position_name = "Electromyogram (muscle activity)",
      position_type = "EMG",
      region = "External",
      standard_system = "EMG"
    ),
    "EMG_Jaw" = list(
      position_name = "Electromyogram jaw",
      position_type = "EMG",
      region = "External",
      standard_system = "EMG"
    ),
    "EMG_Chin" = list(
      position_name = "Electromyogram chin",
      position_type = "EMG",
      region = "External",
      standard_system = "EMG"
    ),
    
    # ===== EXTERNAL CHANNELS - GSR (Galvanic Skin Response) =====
    "GSR" = list(
      position_name = "Galvanic skin response (skin conductance)",
      position_type = "GSR",
      region = "External",
      standard_system = "GSR"
    ),
    "GSR1" = list(
      position_name = "Galvanic skin response channel 1",
      position_type = "GSR",
      region = "External",
      standard_system = "GSR"
    ),
    "GSR2" = list(
      position_name = "Galvanic skin response channel 2",
      position_type = "GSR",
      region = "External",
      standard_system = "GSR"
    ),
    
    # ===== STATUS CHANNEL =====
    "Status" = list(
      position_name = "Status channel (event markers/triggers)",
      position_type = "Status",
      region = "Metadata",
      standard_system = "BioSemi"
    )
  )
  
  return(electrode_db)
}

#' Get Electrode Position Information
#'
#' Looks up standard position information for a given electrode name.
#' Returns details about region, position type, and standard system.
#'
#' @param channel_name Character - name of the electrode (e.g., "Cz", "Fp1", "EOG_L")
#'
#' @return List with position information:
#' \describe{
#'   \item{position_name}{Full descriptive name of the electrode}
#'   \item{position_type}{Category (Frontal, Central, Parietal, EOG, ECG, etc.)}
#'   \item{region}{Brain region or external location}
#'   \item{standard_system}{10-20, 10-10, or external channel type}
#' }
#' Returns NULL with "Unknown" status if electrode not found.
#'
#' @details
#' If the electrode name is not in the standard database, returns a list
#' with "Unknown" values. Case-sensitive matching.
#'
#' @examples
#' \dontrun{
#' get_electrode_position("Cz")    # Central midline
#' get_electrode_position("EOG_L") # Left EOG
#' get_electrode_position("Ch99")  # Unknown channel
#' }
#'
#' @export
get_electrode_position <- function(channel_name) {
  
  # Get the electrode database
  electrode_db <- get_electrode_database()
  
  # Look up the channel name
  if (channel_name %in% names(electrode_db)) {
    return(electrode_db[[channel_name]])
  } else {
    # Return unknown electrode info
    return(list(
      position_name = paste("Unknown electrode:", channel_name),
      position_type = "Unknown",
      region = "Unknown",
      standard_system = "Unknown"
    ))
  }
}

#' Scan BioSemi BDF File for Channel Information
#'
#' Reads a BioSemi BDF file header and extracts channel information without
#' loading the full dataset. Provides quick inspection of file contents.
#'
#' @param file_path Character - full path to the .bdf file
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{channel_number}{Sequential channel number (1-indexed)}
#'   \item{channel_name}{Channel label from BDF header}
#'   \item{position_type}{Standard position category}
#'   \item{region}{Brain region or external location}
#'   \item{standard_system}{10-20, 10-10, EOG, ECG, etc.}
#' }
#'
#' @details
#' This function performs minimal file I/O - only reads the BDF header
#' to extract channel metadata. Useful for:
#' - Quick file preview before full import
#' - Verifying expected electrode layout
#' - Extracting channel names programmatically
#' - Quality control checks
#'
#' Uses the edfReader package (same as read_biosemi).
#'
#' @examples
#' \dontrun{
#' # Quick preview of a BDF file
#' channels <- scan_biosemi_channels("data/subject_01.bdf")
#' print(channels)
#'
#' # Check specific channels
#' eeg_only <- channels[channels$standard_system != "Status", ]
#' external <- channels[channels$region == "External", ]
#' }
#'
#' @importFrom edfReader readEdfHeader
#' @export
scan_biosemi_channels <- function(file_path) {
  
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
  
  # ========== READ BDF HEADER ==========
  
  # Use edfReader to read only the header (no signal data)
  header <- tryCatch({
    edfReader::readEdfHeader(file_path)
  }, error = function(e) {
    stop("ERROR reading BDF file header: ", e$message,
         "\nCheck that the file is a valid BioSemi BDF file.")
  })
  
  # ========== EXTRACT CHANNEL INFORMATION ==========
  
  # Get channel names and clean whitespace
  channels_raw <- header$sHeaders$label
  channels_clean <- trimws(channels_raw)
  
  # Total number of channels
  n_channels <- length(channels_clean)
  
  # ========== CREATE CHANNEL INFO DATA FRAME ==========
  
  # Initialize data frame
  channel_info_df <- data.frame(
    channel_number = 1:n_channels,
    channel_name = channels_clean,
    position_type = character(n_channels),
    region = character(n_channels),
    standard_system = character(n_channels),
    stringsAsFactors = FALSE
  )
  
  # ========== POPULATE POSITION INFORMATION ==========
  
  # For each channel, look up position information
  for (i in 1:n_channels) {
    channel_name <- channels_clean[i]
    position_info <- get_electrode_position(channel_name)
    
    channel_info_df$position_type[i] <- position_info$position_type
    channel_info_df$region[i] <- position_info$region
    channel_info_df$standard_system[i] <- position_info$standard_system
  }
  
  return(channel_info_df)
}

#' Inspect BioSemi BDF File
#'
#' Provides a formatted summary of channel information from a BioSemi BDF file.
#' Prints a comprehensive report without loading the full dataset.
#'
#' @param file_path Character - full path to the .bdf file
#'
#' @return Invisibly returns the channel info data frame from scan_biosemi_channels()
#'
#' @details
#' Prints a formatted report including:
#' - File information
#' - Total channel count
#' - Channel breakdown (EEG, external, status)
#' - Channel listing with positions and regions
#' - Summary statistics
#'
#' @examples
#' \dontrun{
#' # Inspect a BDF file before import
#' inspect_biosemi_file("data/subject_01.bdf")
#' }
#'
#' @export
inspect_biosemi_file <- function(file_path) {
  
  # ========== FILE VALIDATION ==========
  
  if (!file.exists(file_path)) {
    stop("ERROR: File not found at path: ", file_path)
  }
  
  # ========== GET CHANNEL INFORMATION ==========
  
  channel_info <- scan_biosemi_channels(file_path)
  
  # ========== CALCULATE SUMMARY STATISTICS ==========
  
  total_channels <- nrow(channel_info)
  status_channels <- sum(channel_info$standard_system == "Status")
  eeg_channels <- sum(channel_info$region == "Parietal" | 
                        channel_info$region == "Frontal" |
                        channel_info$region == "Central" |
                        channel_info$region == "Temporal" |
                        channel_info$region == "Occipital")
  external_channels <- sum(channel_info$region == "External")
  reference_channels <- sum(channel_info$region == "Reference")
  
  # Count by position type
  position_types <- unique(channel_info$position_type)
  
  # ========== PRINT HEADER ==========
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("BioSemi BDF File Channel Inspection\n")
  cat(strrep("=", 80), "\n\n")
  
  # ========== PRINT FILE INFORMATION ==========
  
  cat("FILE INFORMATION:\n")
  cat(" Path: ", file_path, "\n", sep = "")
  cat(" File name: ", basename(file_path), "\n", sep = "")
  cat(" File size: ", format(file.size(file_path), units = "auto"), "\n", sep = "")
  
  # ========== PRINT CHANNEL SUMMARY ==========
  
  cat("\nCHANNEL SUMMARY:\n")
  cat(" Total channels: ", total_channels, "\n", sep = "")
  cat(" EEG channels: ", eeg_channels, "\n", sep = "")
  cat(" External channels: ", external_channels, " (EOG, ECG, EMG, GSR, etc.)\n", sep = "")
  cat(" Reference channels: ", reference_channels, "\n", sep = "")
  cat(" Status/Trigger channel: ", status_channels, "\n", sep = "")
  
  # ========== PRINT CHANNEL BREAKDOWN BY TYPE ==========
  
  cat("\nCHANNEL BREAKDOWN BY TYPE:\n")
  for (pos_type in sort(unique(channel_info$position_type))) {
    count <- sum(channel_info$position_type == pos_type)
    cat(" ", pos_type, ": ", count, "\n", sep = "")
  }
  
  # ========== PRINT DETAILED CHANNEL LISTING ==========
  
  cat("\nDETAILED CHANNEL LISTING:\n")
  cat(strrep("-", 80), "\n")
  
  # Print header
  cat(sprintf("%-5s %-15s %-20s %-15s %-15s\n",
              "Ch#", "Channel Name", "Position Type", "Region", "System"))
  cat(strrep("-", 80), "\n")
  
  # Print each channel
  for (i in 1:nrow(channel_info)) {
    cat(sprintf("%-5d %-15s %-20s %-15s %-15s\n",
                channel_info$channel_number[i],
                channel_info$channel_name[i],
                channel_info$position_type[i],
                channel_info$region[i],
                channel_info$standard_system[i]))
  }
  
  cat(strrep("-", 80), "\n")
  
  # ========== PRINT CHANNEL LISTING BY REGION ==========
  
  cat("\nCHANNELS BY REGION:\n")
  regions <- sort(unique(channel_info$region))
  for (region in regions) {
    channels_in_region <- channel_info$channel_name[channel_info$region == region]
    channels_str <- paste(channels_in_region, collapse = ", ")
    
    # Truncate long lists
    if (nchar(channels_str) > 70) {
      channels_str <- paste0(substr(channels_str, 1, 67), "...")
    }
    
    cat(" ", region, ": ", channels_str, "\n", sep = "")
  }
  
  # ========== PRINT FOOTER ==========
  
  cat(strrep("=", 80), "\n\n")
  
  # ========== RETURN INVISIBLY ==========
  
  invisible(channel_info)
}


#' ============================================================================
#' Detect Electrode Naming System in EEG Dataset
#'
#' Automatically identifies whether electrodes in a BioSemi BDF file (or loaded
#' dataset) are coded using standard international notation (10-20/10-10 system)
#' or BioSemi's internal classification system (A1-A32, B1-B32).
#'
#' @param file_path Character - full path to the BioSemi .bdf file, OR NULL if
#'   using channel_names parameter
#' @param channel_names Character vector - channel names from already-loaded dataset.
#'   Use this if your data is already in R environment. If provided, file_path
#'   is ignored.
#'
#' @return List with the following components:
#' \describe{
#'   \item{naming_system}{Character. One of:
#'     - "10-20/10-10" = Standard international electrode naming
#'     - "BioSemi_AB" = BioSemi internal A1-A32 + B1-B32 naming (64 electrodes)
#'     - "Mixed" = Both standard and BioSemi naming detected
#'     - "Unknown" = Unrecognized naming convention
#'   }
#'   \item{confidence}{Numeric (0-1). Confidence score based on percentage
#'     of recognized electrodes. 1.0 = 100% match to identified system.
#'   }
#'   \item{channels_detected}{Data frame containing all channel information
#'     with columns:
#'     - channel_number: Channel index
#'     - channel_name: Original channel name from BDF
#'     - naming_type: Classification (Standard 10-20/10-10, BioSemi A/B, External, etc.)
#'     - system_category: Primary system category
#'     - standard_equivalent: 10-20/10-10 name if BioSemi detected (or NA)
#'   }
#'   \item{summary}{Character vector with human-readable summary of findings:
#'     - Number of channels in each system
#'     - Unrecognized channels (if any)
#'     - Recommendation for next steps
#'   }
#'   \item{detailed_breakdown}{List with counts by system type:
#'     \describe{
#'       \item{standard_10_20}{Integer. Count of recognized 10-20 electrodes}
#'       \item{standard_10_10}{Integer. Count of recognized 10-10 electrodes}
#'       \item{biosemi_a}{Integer. Count of BioSemi A1-A32 electrodes}
#'       \item{biosemi_b}{Integer. Count of BioSemi B1-B32 electrodes}
#'       \item{external}{Integer. Count of external channels (EOG, ECG, etc.)}
#'       \item{unrecognized}{Integer. Count of unrecognized channels}
#'     }
#'   }
#' }
#'
#' @details
#' This function performs automatic detection of electrode naming conventions
#' by analyzing channel names from either:
#' 1. A BDF file header (if file_path is provided)
#' 2. A character vector of channel names (if channel_names is provided)
#'
#' **Two Usage Modes:**
#'
#' **Mode 1: File Path (reads from disk)**
#' ```
#' result <- detect_electrode_naming_system(file_path = "data/subject_01.bdf")
#' ```
#'
#' **Mode 2: Loaded Dataset (from R environment)**
#' ```
#' my_data <- read_biosemi("data/subject_01.bdf")
#' result <- detect_electrode_naming_system(channel_names = colnames(my_data))
#' # OR
#' result <- detect_electrode_naming_system(channel_names = my_data$channel_names)
#' ```
#'
#' The function then:
#'
#' 1. **Classifies each channel** against known naming patterns:
#'    - Standard 10-20: Fp1, Fpz, Cz, Pz, O1, O2, etc.
#'    - Standard 10-10: AF3, FC5, CP3, PO3, etc.
#'    - BioSemi A-series: A1-A32 (64-channel cap)
#'    - BioSemi B-series: B1-B32 (64-channel cap)
#'    - External: EXG1-EXG8 (reference electrodes shown as green in BioSemi layout)
#'    - Status: Status/trigger channel
#'
#' 2. **Provides BioSemi-to-Standard mapping** using official BioSemi 64-electrode
#'    cap layout (A1-A32 + B1-B32 mapped to 10-20/10-10 positions)
#'
#' 3. **Calculates confidence** based on percentage of channels correctly
#'    classified within detected system
#'
#' 4. **Provides detailed breakdown** of findings with actionable recommendations
#'
#' **BioSemi 64-Channel Layout (from official documentation):**
#'
#' The function uses the official BioSemi mapping:
#' - A1-A32: First 32 electrodes (mixed 10-20 and 10-10 positions)
#' - B1-B32: Second 32 electrodes (mixed 10-20 and 10-10 positions)
#' - External references (green circles in layout): EXG1-EXG8
#'
#' **Why This Matters:**
#'
#' BioSemi hardware records with internal A1-A32/B1-B32 labeling, but many
#' analysis pipelines expect standard 10-20/10-10 names. This function helps you:
#'
#' - Identify if channel mapping is needed
#' - Get the exact 10-20/10-10 equivalent for each BioSemi channel
#' - Choose appropriate preprocessing pipeline
#' - Verify expected electrode layout
#' - Document your data naming convention
#'
#' **Important Notes:**
#'
#' - Matching is **case-sensitive** (e.g., "Cz" recognized, "cz" is not)
#' - When using file_path, reads **header only** (very fast, ~50-100 ms)
#' - Works with any BioSemi BDF file format
#' - External channels (EXG1-EXG8) are recognized as reference electrodes
#'
#' @examples
#' \dontrun{
#'
#'   # Example 1: Check a BioSemi file from disk
#'   result <- detect_electrode_naming_system(file_path = "data/subject_01.bdf")
#'   print(result$naming_system)      # "10-20/10-10"
#'   print(result$confidence)         # 1.0
#'
#'   # Example 2: Check already-loaded dataset
#'   my_data <- read_biosemi("data/subject_01.bdf")
#'   result <- detect_electrode_naming_system(channel_names = colnames(my_data))
#'   print(result$naming_system)      # "BioSemi_AB"
#'   
#'   # Example 3: With eegUtils or other packages
#'   library(eegUtils)
#'   eeg_data <- import_raw("data/subject_01.bdf")
#'   result <- detect_electrode_naming_system(channel_names = channel_names(eeg_data))
#'   
#'   # Example 4: See the mapping
#'   View(result$channels_detected[, c("channel_name", "standard_equivalent")])
#'   
#'   # Example 5: Export mapping for reference
#'   write.csv(result$channels_detected, "biosemi_mapping.csv", row.names = FALSE)
#'
#'   # Example 6: Batch processing from files
#'   files <- list.files("data/", pattern = ".bdf$", full.names = TRUE)
#'   for (file in files) {
#'     result <- detect_electrode_naming_system(file_path = file)
#'     cat("File:", basename(file), "System:", result$naming_system, "\n")
#'   }
#'
#'   # Example 7: Batch processing from loaded datasets
#'   dataset_list <- list(data1, data2, data3)
#'   for (i in seq_along(dataset_list)) {
#'     result <- detect_electrode_naming_system(channel_names = colnames(dataset_list[[i]]))
#'     cat("Dataset", i, "System:", result$naming_system, "\n")
#'   }
#' }
#'
#' @keywords internal
#'
#' @importFrom edfReader readEdfHeader
#' @export
#'
detect_electrode_naming_system <- function(file_path = NULL, channel_names = NULL) {
  
  # ========================================================================
  # INPUT VALIDATION
  # ========================================================================
  
  # Must provide either file_path OR channel_names (but not both or neither)
  if (is.null(file_path) && is.null(channel_names)) {
    stop("ERROR: Must provide either 'file_path' OR 'channel_names' parameter.\n",
         "  Usage:\n",
         "  - For file: detect_electrode_naming_system(file_path = 'data.bdf')\n",
         "  - For loaded data: detect_electrode_naming_system(channel_names = colnames(my_data))")
  }
  
  # If both provided, use channel_names and warn
  if (!is.null(file_path) && !is.null(channel_names)) {
    warning("Both 'file_path' and 'channel_names' provided. Using 'channel_names' and ignoring 'file_path'.")
    file_path <- NULL
  }
  
  # ========================================================================
  # GET CHANNEL NAMES (from file OR from parameter)
  # ========================================================================
  
  if (!is.null(channel_names)) {
    # MODE 1: Channel names provided directly (loaded dataset)
    
    # Validate channel_names is a character vector
    if (!is.character(channel_names)) {
      stop("ERROR: 'channel_names' must be a character vector.\n",
           "  Example: channel_names = c('A1', 'A2', 'A3', ...)")
    }
    
    # Clean channel names
    channels_clean <- trimws(channel_names)
    source_info <- "from provided channel names"
    
  } else {
    # MODE 2: File path provided (read from disk)
    
    # Check file exists
    if (!file.exists(file_path)) {
      stop("ERROR: File not found at path: ", file_path)
    }
    
    # Check file has .bdf extension
    if (!grepl("\\.bdf$", file_path, ignore.case = TRUE)) {
      warning("WARNING: File does not have .bdf extension. ",
              "Expected BioSemi BDF format. Attempting to read anyway...")
    }
    
    # Read BDF header
    header <- tryCatch({
      edfReader::readEdfHeader(file_path)
    }, error = function(e) {
      stop("ERROR reading BDF file header: ", e$message,
           "\nCheck that the file is a valid BioSemi BDF file.")
    })
    
    # Extract and clean channel names
    channels_raw <- header$sHeaders$label
    channels_clean <- trimws(channels_raw)
    source_info <- basename(file_path)
  }
  
  n_channels <- length(channels_clean)
  
  # ========================================================================
  # DEFINE BIOSEMI 64-CHANNEL MAPPING (FROM OFFICIAL LAYOUT)
  # ========================================================================
  # Based on the BioSemi 64+2 electrode cap layout diagram
  # Maps BioSemi internal naming (A1-A32, B1-B32) to 10-20/10-10 standard
  
  biosemi_to_standard <- list(
    # A-series electrodes (A1-A32)
    "A1"  = "Fp1",
    "A2"  = "AF7",
    "A3"  = "AF3",
    "A4"  = "F1",
    "A5"  = "F3",
    "A6"  = "F5",
    "A7"  = "F7",
    "A8"  = "FT7",
    "A9"  = "FC5",
    "A10" = "FC3",
    "A11" = "FC1",
    "A12" = "C1",
    "A13" = "C3",
    "A14" = "C5",
    "A15" = "T7",
    "A16" = "TP7",
    "A17" = "CP5",
    "A18" = "CP3",
    "A19" = "CP1",
    "A20" = "P1",
    "A21" = "P3",
    "A22" = "P5",
    "A23" = "P7",
    "A24" = "P9",
    "A25" = "PO7",
    "A26" = "PO3",
    "A27" = "O1",
    "A28" = "Iz",
    "A29" = "Oz",
    "A30" = "POz",
    "A31" = "Pz",
    "A32" = "CPz",
    
    # B-series electrodes (B1-B32)
    "B1"  = "Fpz",
    "B2"  = "Fp2",
    "B3"  = "AF8",
    "B4"  = "AF4",
    "B5"  = "AFz",
    "B6"  = "Fz",
    "B7"  = "F2",
    "B8"  = "F4",
    "B9"  = "F6",
    "B10" = "F8",
    "B11" = "FT8",
    "B12" = "FC6",
    "B13" = "FC4",
    "B14" = "FC2",
    "B15" = "FCz",
    "B16" = "Cz",
    "B17" = "C2",
    "B18" = "C4",
    "B19" = "C6",
    "B20" = "T8",
    "B21" = "TP8",
    "B22" = "CP6",
    "B23" = "CP4",
    "B24" = "CP2",
    "B25" = "P2",
    "B26" = "P4",
    "B27" = "P6",
    "B28" = "P8",
    "B29" = "P10",
    "B30" = "PO8",
    "B31" = "PO4",
    "B32" = "O2"
  )
  
  # ========================================================================
  # DEFINE NAMING PATTERNS
  # ========================================================================
  
  # Standard 10-20 system electrodes
  standard_10_20 <- c(
    # Prefrontal
    "Fp1", "Fpz", "Fp2",
    # Frontal
    "F7", "F3", "Fz", "F4", "F8",
    # Temporal
    "T7", "T8",
    # Central
    "C3", "Cz", "C4",
    # Parietal
    "P7", "P3", "Pz", "P4", "P8",
    # Occipital
    "O1", "Oz", "O2",
    # Additional temporal
    "P9", "P10",
    # Inion
    "Iz"
  )
  
  # Standard 10-10 system electrodes (extended)
  standard_10_10 <- c(
    # Anterior Frontal
    "AF7", "AF3", "AFz", "AF4", "AF8",
    # Fronto-Temporal
    "FT7", "FT8",
    # Fronto-Central
    "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6",
    # Central extended
    "C5", "C1", "C2", "C6",
    # Centro-Temporal
    "TP7", "TP8",
    # Centro-Parietal
    "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6",
    # Parietal extended
    "P5", "P1", "P2", "P6",
    # Parieto-Occipital
    "PO7", "PO3", "POz", "PO4", "PO8"
  )
  
  # Combine all standard electrodes
  all_standard <- c(standard_10_20, standard_10_10)
  
  # BioSemi internal naming patterns
  biosemi_a_pattern <- "^A([0-9]{1,2})$"  # A1 to A32
  biosemi_b_pattern <- "^B([0-9]{1,2})$"  # B1 to B32
  
  # External channel patterns (reference electrodes - shown as green in layout)
  # BioSemi uses EXG1-EXG8 for external/reference electrodes
  external_patterns <- c(
    "^EXG[0-9]",                       # BioSemi external (EXG1-EXG8)
    "^EOG", "^ECG", "^EMG", "^GSR",    # Standard external
    "^HEOG", "^VEOG",                  # EOG variants
    "Status", "Trigger", "Marker"      # Event channels
  )
  
  # ========================================================================
  # CLASSIFY EACH CHANNEL
  # ========================================================================
  
  channel_classification <- data.frame(
    channel_number = 1:n_channels,
    channel_name = channels_clean,
    naming_type = character(n_channels),
    system_category = character(n_channels),
    standard_equivalent = character(n_channels),  # NEW: for BioSemi mapping
    stringsAsFactors = FALSE
  )
  
  # Initialize counters
  count_standard_10_20 <- 0
  count_standard_10_10 <- 0
  count_biosemi_a <- 0
  count_biosemi_b <- 0
  count_external <- 0
  count_unrecognized <- 0
  
  # Classify each channel
  for (i in 1:n_channels) {
    ch_name <- channels_clean[i]
    
    # Check if external channel
    is_external <- any(sapply(external_patterns, function(pattern) {
      grepl(pattern, ch_name, ignore.case = TRUE)
    }))
    
    if (is_external) {
      channel_classification$naming_type[i] <- "External"
      channel_classification$system_category[i] <- "External"
      channel_classification$standard_equivalent[i] <- NA
      count_external <- count_external + 1
    }
    # Check if standard 10-20
    else if (ch_name %in% standard_10_20) {
      channel_classification$naming_type[i] <- "Standard 10-20"
      channel_classification$system_category[i] <- "10-20/10-10"
      channel_classification$standard_equivalent[i] <- ch_name  # Already standard
      count_standard_10_20 <- count_standard_10_20 + 1
    }
    # Check if standard 10-10
    else if (ch_name %in% standard_10_10) {
      channel_classification$naming_type[i] <- "Standard 10-10"
      channel_classification$system_category[i] <- "10-20/10-10"
      channel_classification$standard_equivalent[i] <- ch_name  # Already standard
      count_standard_10_10 <- count_standard_10_10 + 1
    }
    # Check if BioSemi A-series
    else if (grepl(biosemi_a_pattern, ch_name)) {
      channel_classification$naming_type[i] <- "BioSemi A-series"
      channel_classification$system_category[i] <- "BioSemi_AB"
      # Map to standard equivalent
      if (ch_name %in% names(biosemi_to_standard)) {
        channel_classification$standard_equivalent[i] <- biosemi_to_standard[[ch_name]]
      } else {
        channel_classification$standard_equivalent[i] <- NA
      }
      count_biosemi_a <- count_biosemi_a + 1
    }
    # Check if BioSemi B-series
    else if (grepl(biosemi_b_pattern, ch_name)) {
      channel_classification$naming_type[i] <- "BioSemi B-series"
      channel_classification$system_category[i] <- "BioSemi_AB"
      # Map to standard equivalent
      if (ch_name %in% names(biosemi_to_standard)) {
        channel_classification$standard_equivalent[i] <- biosemi_to_standard[[ch_name]]
      } else {
        channel_classification$standard_equivalent[i] <- NA
      }
      count_biosemi_b <- count_biosemi_b + 1
    }
    # Unrecognized
    else {
      channel_classification$naming_type[i] <- "Unrecognized"
      channel_classification$system_category[i] <- "Unknown"
      channel_classification$standard_equivalent[i] <- NA
      count_unrecognized <- count_unrecognized + 1
    }
  }
  
  # ========================================================================
  # DETERMINE OVERALL NAMING SYSTEM
  # ========================================================================
  
  # Count by primary system
  standard_total <- count_standard_10_20 + count_standard_10_10
  biosemi_total <- count_biosemi_a + count_biosemi_b
  
  # Exclude external and status channels from system determination
  non_external <- n_channels - count_external
  non_external <- ifelse(non_external < 1, 1, non_external)  # Avoid division by zero
  
  # Determine primary naming system
  naming_system <- "Unknown"
  confidence <- 0
  
  if (standard_total > 0 && biosemi_total == 0 && count_unrecognized == 0) {
    naming_system <- "10-20/10-10"
    confidence <- standard_total / non_external
  } else if (biosemi_total > 0 && standard_total == 0 && count_unrecognized == 0) {
    naming_system <- "BioSemi_AB"
    confidence <- biosemi_total / non_external
  } else if (standard_total > 0 && biosemi_total > 0) {
    naming_system <- "Mixed"
    confidence <- max(standard_total, biosemi_total) / non_external
  } else if (count_unrecognized > 0) {
    naming_system <- "Unknown"
    confidence <- (non_external - count_unrecognized) / non_external
  }
  
  # ========================================================================
  # GENERATE SUMMARY
  # ========================================================================
  
  summary_lines <- c()
  
  # Header
  summary_lines <- c(summary_lines,
                     paste("Electrode Naming System Detection for:", source_info))
  summary_lines <- c(summary_lines,
                     strrep("=", 70))
  summary_lines <- c(summary_lines, "")
  
  # Main findings
  summary_lines <- c(summary_lines,
                     paste("PRIMARY NAMING SYSTEM DETECTED:", naming_system))
  summary_lines <- c(summary_lines,
                     paste("Confidence Score:", round(confidence * 100, 1), "%"))
  summary_lines <- c(summary_lines, "")
  
  # Detailed breakdown
  summary_lines <- c(summary_lines, "CHANNEL BREAKDOWN:")
  summary_lines <- c(summary_lines,
                     paste("  Standard 10-20 electrodes:", count_standard_10_20))
  summary_lines <- c(summary_lines,
                     paste("  Standard 10-10 electrodes:", count_standard_10_10))
  if (biosemi_total > 0) {
    summary_lines <- c(summary_lines,
                       paste("  BioSemi A-series (A1-A32):", count_biosemi_a))
    summary_lines <- c(summary_lines,
                       paste("  BioSemi B-series (B1-B32):", count_biosemi_b))
    summary_lines <- c(summary_lines,
                       paste("  Total BioSemi electrodes:", biosemi_total))
  }
  summary_lines <- c(summary_lines,
                     paste("  External/Reference channels (EXG, EOG, etc.):", count_external))
  
  if (count_unrecognized > 0) {
    summary_lines <- c(summary_lines,
                       paste("  Unrecognized channels:", count_unrecognized))
    
    # List unrecognized channels
    unrecognized_chs <- channel_classification$channel_name[
      channel_classification$system_category == "Unknown"
    ]
    summary_lines <- c(summary_lines,
                       paste("    Channel names:", paste(unrecognized_chs, collapse = ", ")))
  }
  
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines, "RECOMMENDATIONS:")
  
  # Recommendations based on system detected
  if (naming_system == "10-20/10-10") {
    summary_lines <- c(summary_lines,
                       "✓ Standard international naming detected.")
    summary_lines <- c(summary_lines,
                       "✓ Compatible with most EEG analysis pipelines.")
    summary_lines <- c(summary_lines,
                       "✓ Ready for analysis without channel remapping.")
  } else if (naming_system == "BioSemi_AB") {
    summary_lines <- c(summary_lines,
                       "⚠ BioSemi internal naming detected (A1-A32 + B1-B32)")
    summary_lines <- c(summary_lines,
                       "⚠ Standard 10-20/10-10 equivalents provided in 'standard_equivalent' column")
    summary_lines <- c(summary_lines,
                       "⚠ Use channels_detected$standard_equivalent for standard names")
    summary_lines <- c(summary_lines,
                       "⚠ Export mapping: write.csv(result$channels_detected, 'mapping.csv')")
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "   BioSemi to Standard mapping examples:")
    summary_lines <- c(summary_lines, "   A1  -> Fp1  |  B1  -> Fpz")
    summary_lines <- c(summary_lines, "   A16 -> TP7  |  B16 -> Cz")
    summary_lines <- c(summary_lines, "   A32 -> CPz  |  B32 -> O2")
  } else if (naming_system == "Mixed") {
    summary_lines <- c(summary_lines,
                       "⚠ Mixed naming systems detected.")
    summary_lines <- c(summary_lines,
                       "⚠ Some channels use standard names, others use BioSemi names.")
    summary_lines <- c(summary_lines,
                       "⚠ Check 'standard_equivalent' column for mapped names.")
    summary_lines <- c(summary_lines,
                       "⚠ Manual inspection may be needed.")
  } else {
    summary_lines <- c(summary_lines,
                       "⚠ Unrecognized naming convention detected.")
    summary_lines <- c(summary_lines,
                       "⚠ Unable to automatically classify electrodes.")
    summary_lines <- c(summary_lines,
                       "⚠ Please verify electrode naming in your BDF file.")
  }
  
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines,
                     strrep("=", 70))
  
  # ========================================================================
  # BUILD DETAILED BREAKDOWN LIST
  # ========================================================================
  
  detailed_breakdown <- list(
    standard_10_20 = count_standard_10_20,
    standard_10_10 = count_standard_10_10,
    biosemi_a = count_biosemi_a,
    biosemi_b = count_biosemi_b,
    biosemi_total = biosemi_total,
    external = count_external,
    unrecognized = count_unrecognized,
    total_channels = n_channels
  )
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  return(list(
    naming_system = naming_system,
    confidence = confidence,
    channels_detected = channel_classification,
    summary = summary_lines,
    detailed_breakdown = detailed_breakdown
  ))
}
