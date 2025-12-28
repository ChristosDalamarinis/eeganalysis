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
#' @param data Either:
#'   - Character string: file path to BioSemi .bdf file
#'   - Character vector: channel names from loaded dataset
#'   - Data frame/matrix: dataset with channels as columns (uses colnames)
#'   - List: dataset with $channels component (extracts channel names)
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
#'   \item{summary}{Character vector with human-readable summary of findings}
#'   \item{detailed_breakdown}{List with counts by system type}
#' }
#'
#' @details
#' This function performs automatic detection of electrode naming conventions.
#' It intelligently handles multiple input types:
#'
#' **Usage Examples:**
#'
#' ```
#' # From file path
#' result <- detect_electrode_naming_system("data/subject_01.bdf")
#'
#' # From loaded dataset (list with $channels)
#' data2 <- read_biosemi("data/subject_01.bdf")
#' result <- detect_electrode_naming_system(data2$channels)
#'
#' # From data frame/matrix (uses column names)
#' result <- detect_electrode_naming_system(my_eeg_data)
#'
#' # From character vector of names
#' ch_names <- c("A1", "A2", "A3", ...)
#' result <- detect_electrode_naming_system(ch_names)
#' ```
#'
#' **BioSemi 64-Channel Layout:**
#' Maps A1-A32 + B1-B32 to standard 10-20/10-10 positions using official
#' BioSemi electrode cap layout.
#'
#' @examples
#' \dontrun{
#'   # Example 1: From file
#'   result <- detect_electrode_naming_system("data/subject_01.bdf")
#'
#'   # Example 2: From loaded dataset
#'   data2 <- read_biosemi("data/subject_01.bdf")
#'   result <- detect_electrode_naming_system(data2$channels)
#'
#'   # Example 3: See the mapping
#'   View(result$channels_detected[, c("channel_name", "standard_equivalent")])
#' }
#'
#' @keywords internal
#' @importFrom edfReader readEdfHeader
#' @export
#'
detect_electrode_naming_system <- function(data) {
  
  # ========================================================================
  # INPUT VALIDATION & TYPE DETECTION
  # ========================================================================
  
  channels_clean <- NULL
  source_info <- NULL
  
  # Determine input type and extract channel names accordingly
  if (is.character(data)) {
    # Case 1: Character vector
    if (length(data) == 1 && file.exists(data)) {
      # Single string that is a valid file path
      # MODE 1: Read from BDF file
      
      if (!grepl("\\.bdf$", data, ignore.case = TRUE)) {
        warning("WARNING: File does not have .bdf extension. ",
                "Expected BioSemi BDF format. Attempting to read anyway...")
      }
      
      header <- tryCatch({
        edfReader::readEdfHeader(data)
      }, error = function(e) {
        stop("ERROR reading BDF file header: ", e$message,
             "\nCheck that the file is a valid BioSemi BDF file.")
      })
      
      channels_raw <- header$sHeaders$label
      channels_clean <- trimws(channels_raw)
      source_info <- basename(data)
      
    } else {
      # MODE 2: Character vector of channel names
      channels_clean <- trimws(data)
      source_info <- "from provided channel names"
    }
    
  } else if (is.data.frame(data) || is.matrix(data)) {
    # MODE 3: Data frame or matrix - use column names
    channels_clean <- colnames(data)
    if (is.null(channels_clean)) {
      stop("ERROR: Data frame/matrix has no column names.\n",
           "Cannot extract channel names.")
    }
    channels_clean <- trimws(channels_clean)
    source_info <- "from data frame/matrix column names"
    
  } else if (is.list(data)) {
    # MODE 4: List - try to find channel names in common locations
    if (!is.null(data$channels)) {
      channels_clean <- trimws(data$channels)
      source_info <- "from data$channels"
    } else if (!is.null(data$channel_names)) {
      channels_clean <- trimws(data$channel_names)
      source_info <- "from data$channel_names"
    } else if (!is.null(data$channel_info) && !is.null(data$channel_info$label)) {
      channels_clean <- trimws(data$channel_info$label)
      source_info <- "from data$channel_info$label"
    } else {
      stop("ERROR: List provided but could not find channel names.\n",
           "Tried: data$channels, data$channel_names, data$channel_info$label\n",
           "Please provide channel names explicitly as a character vector.")
    }
    
  } else {
    # Unknown type
    stop("ERROR: Invalid input type.\n",
         "Expected one of:\n",
         "  - File path (character string)\n",
         "  - Channel names (character vector)\n",
         "  - Data frame/matrix (uses column names)\n",
         "  - List with $channels component\n",
         "Got: ", class(data)[1])
  }
  
  # Final validation
  if (is.null(channels_clean) || length(channels_clean) == 0) {
    stop("ERROR: No channel names found or extracted.")
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
  
  # External channel patterns
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
    standard_equivalent = character(n_channels),
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
      channel_classification$standard_equivalent[i] <- ch_name
      count_standard_10_20 <- count_standard_10_20 + 1
    }
    # Check if standard 10-10
    else if (ch_name %in% standard_10_10) {
      channel_classification$naming_type[i] <- "Standard 10-10"
      channel_classification$system_category[i] <- "10-20/10-10"
      channel_classification$standard_equivalent[i] <- ch_name
      count_standard_10_10 <- count_standard_10_10 + 1
    }
    # Check if BioSemi A-series
    else if (grepl(biosemi_a_pattern, ch_name)) {
      channel_classification$naming_type[i] <- "BioSemi A-series"
      channel_classification$system_category[i] <- "BioSemi_AB"
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
  
  standard_total <- count_standard_10_20 + count_standard_10_10
  biosemi_total <- count_biosemi_a + count_biosemi_b
  non_external <- n_channels - count_external
  non_external <- ifelse(non_external < 1, 1, non_external)
  
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
  summary_lines <- c(summary_lines,
                     paste("Electrode Naming System Detection for:", source_info))
  summary_lines <- c(summary_lines, strrep("=", 70))
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines,
                     paste("PRIMARY NAMING SYSTEM DETECTED:", naming_system))
  summary_lines <- c(summary_lines,
                     paste("Confidence Score:", round(confidence * 100, 1), "%"))
  summary_lines <- c(summary_lines, "")
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
                     paste("  External/Reference channels:", count_external))
  
  if (count_unrecognized > 0) {
    summary_lines <- c(summary_lines,
                       paste("  Unrecognized channels:", count_unrecognized))
    unrecognized_chs <- channel_classification$channel_name[
      channel_classification$system_category == "Unknown"
    ]
    summary_lines <- c(summary_lines,
                       paste("    Channel names:", paste(unrecognized_chs, collapse = ", ")))
  }
  
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines, "RECOMMENDATIONS:")
  
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
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "   BioSemi to Standard mapping examples:")
    summary_lines <- c(summary_lines, "   A1  -> Fp1  |  B1  -> Fpz")
    summary_lines <- c(summary_lines, "   A16 -> TP7  |  B16 -> Cz")
    summary_lines <- c(summary_lines, "   A32 -> CPz  |  B32 -> O2")
  } else if (naming_system == "Mixed") {
    summary_lines <- c(summary_lines,
                       "⚠ Mixed naming systems detected.")
    summary_lines <- c(summary_lines,
                       "⚠ Check 'standard_equivalent' column for mapped names.")
  } else {
    summary_lines <- c(summary_lines,
                       "⚠ Unrecognized naming convention detected.")
    summary_lines <- c(summary_lines,
                       "⚠ Please verify electrode naming.")
  }
  
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines, strrep("=", 70))
  
  # ========================================================================
  # BUILD DETAILED BREAKDOWN
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




#' ============================================================================
#' Plot 3D Electrode Layout on Scalp (BioSemi 64-Channel)
#'
#' Creates an interactive 3D visualization of electrode positions on a realistic
#' head model using OFFICIAL BioSemi 64-channel electrode coordinates.
#'
#' @param channels Character vector - channel names to plot. If NULL, plots all
#'   64 standard 10-20/10-10 electrodes. Can also accept BioSemi A/B naming.
#' @param highlight_channels Character vector - specific channels to highlight
#'   in a different color (optional)
#' @param show_labels Logical - whether to show electrode labels on the plot.
#'   Default is TRUE.
#' @param color_by Character - how to color electrodes. Options:
#'   - "region" = color by brain region (frontal, central, parietal, occipital, temporal)
#'   - "single" = all one color (default blue)
#'   Default is "region".
#' @param show_head Logical - whether to show the head outline. Default is TRUE.
#'
#' @return Returns a plotly object (interactive 3D plot viewable in RStudio)
#'
#' @details
#' This function visualizes electrode positions in 3D space using the official
#' BioSemi 64-channel electrode cap coordinates. All positions are exact Cartesian
#' coordinates (x, y, z) in millimeters from BioSemi's cap layout documentation.
#'
#' **Coordinate System:**
#' - Origin (0, 0, 0) at center of head
#' - X-axis: Left (-) to Right (+)
#' - Y-axis: Back (-) to Front (+)
#' - Z-axis: Bottom (-) to Top (+)
#' - Nose points in +Y direction
#'
#' @examples
#' \dontrun{
#'
#'   # Example 1: Plot all 64 electrodes
#'   plot_electrode_3d()
#'
#'   # Example 2: Plot specific channels
#'   plot_electrode_3d(channels = c("Cz", "Fz", "Pz", "Oz", "T7", "T8"))
#'
#'   # Example 3: Highlight midline electrodes
#'   plot_electrode_3d(highlight_channels = c("Fpz", "Fz", "Cz", "Pz", "Oz"))
#'
#'   # Example 4: No head outline, just electrodes
#'   plot_electrode_3d(show_head = FALSE)
#' }
#'
#' @keywords internal
#' @importFrom plotly plot_ly add_trace layout
#' @export
#'
plot_electrode_3d <- function(channels = NULL,
                              highlight_channels = NULL,
                              show_labels = TRUE,
                              color_by = "region",
                              show_head = TRUE) {
  
  # ========================================================================
  # LOAD REQUIRED PACKAGE
  # ========================================================================
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive 3D plots.\n",
         "Install it with: install.packages('plotly')")
  }
  
  # ========================================================================
  # OFFICIAL BIOSEMI 64-CHANNEL COORDINATES
  # ========================================================================
  
  electrode_coords <- data.frame(
    electrode = c(
      "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3",
      "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1",
      "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz",
      "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4",
      "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4",
      "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8",
      "P10", "PO8", "PO4", "O2"
    ),
    x = c(
      -27, -51, -36, -25, -48, -64, -71, -83, -78, -59,
      -33, -34, -63, -82, -87, -83, -78, -59, -33, -25,
      -48, -64, -71, -64, -51, -36, -27, 0, 0, 0,
      0, 0, 0, 27, 51, 36, 0, 0, 25, 48,
      64, 71, 83, 78, 59, 33, 0, 0, 34, 63,
      82, 87, 83, 78, 59, 33, 25, 48, 64, 71,
      64, 51, 36, 27
    ),
    y = c(
      83, 71, 76, 62, 59, 55, 51, 27, 30, 31,
      33, 0, 0, 0, 0, -27, -30, -31, -33, -62,
      -59, -55, -51, -47, -71, -76, -83, -79, -87, -82,
      -63, -34, 87, 83, 71, 76, 82, 63, 62, 59,
      55, 51, 27, 30, 31, 33, 34, 0, 0, 0,
      0, 0, -27, -30, -31, -33, -62, -59, -55, -51,
      -47, -71, -76, -83
    ),
    z = c(
      -3, -3, 24, 56, 44, 23, -3, -3, 27, 56,
      74, 81, 61, 31, -3, -3, 27, 56, 74, 56,
      44, 23, -3, -37, -3, 24, -3, -37, -3, 31,
      61, 81, -3, -3, -3, 24, 31, 61, 56, 44,
      23, -3, -3, 27, 56, 74, 81, 88, 81, 61,
      31, -3, -3, 27, 56, 74, 56, 44, 23, -3,
      -37, -3, 24, -3
    ),
    stringsAsFactors = FALSE
  )
  
  # ========================================================================
  # BIOSEMI A/B TO STANDARD MAPPING
  # ========================================================================
  
  biosemi_mapping <- data.frame(
    biosemi = c(
      "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
      "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20",
      "A21", "A22", "A23", "A24", "A25", "A26", "A27", "A28", "A29", "A30",
      "A31", "A32",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
      "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18", "B19", "B20",
      "B21", "B22", "B23", "B24", "B25", "B26", "B27", "B28", "B29", "B30",
      "B31", "B32"
    ),
    standard = c(
      "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3",
      "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1",
      "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz",
      "Pz", "CPz",
      "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8",
      "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8",
      "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8",
      "PO4", "O2"
    ),
    stringsAsFactors = FALSE
  )
  
  # ========================================================================
  # DETERMINE WHICH CHANNELS TO PLOT
  # ========================================================================
  
  if (is.null(channels)) {
    channels_to_plot <- electrode_coords$electrode
  } else {
    channels_to_plot <- channels
  }
  
  # Map BioSemi names to standard names if needed
  plot_data <- data.frame(
    original_name = channels_to_plot,
    standard_name = channels_to_plot,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(plot_data)) {
    ch <- plot_data$original_name[i]
    if (ch %in% biosemi_mapping$biosemi) {
      plot_data$standard_name[i] <- biosemi_mapping$standard[biosemi_mapping$biosemi == ch]
    }
  }
  
  # Merge with coordinates
  plot_data <- merge(plot_data, electrode_coords, 
                     by.x = "standard_name", by.y = "electrode",
                     all.x = TRUE)
  
  # Remove electrodes without coordinates
  missing_coords <- plot_data$original_name[is.na(plot_data$x)]
  if (length(missing_coords) > 0) {
    warning("No coordinates found for: ", paste(missing_coords, collapse = ", "))
    plot_data <- plot_data[!is.na(plot_data$x), ]
  }
  
  if (nrow(plot_data) == 0) {
    stop("No valid electrodes to plot.")
  }
  
  # ========================================================================
  # ASSIGN COLORS BY BRAIN REGION
  # ========================================================================
  
  if (color_by == "region") {
    plot_data$color_group <- "Other"
    plot_data$color_group[grepl("^Fp|^AF|^F[^CT]", plot_data$standard_name)] <- "Frontal"
    plot_data$color_group[grepl("^C|^FC", plot_data$standard_name)] <- "Central"
    plot_data$color_group[grepl("^P[^O]|^CP", plot_data$standard_name)] <- "Parietal"
    plot_data$color_group[grepl("^O|^PO|^Iz", plot_data$standard_name)] <- "Occipital"
    plot_data$color_group[grepl("^T[0-9]|^FT|^TP", plot_data$standard_name)] <- "Temporal"
    
    color_palette <- c(
      "Frontal" = "#e41a1c",
      "Central" = "#377eb8",
      "Parietal" = "#4daf4a",
      "Occipital" = "#984ea3",
      "Temporal" = "#ff7f00",
      "Other" = "#999999"
    )
    
  } else {
    plot_data$color_group <- "Electrode"
    color_palette <- c("Electrode" = "#377eb8")
  }
  
  # Highlight specific channels if requested
  if (!is.null(highlight_channels)) {
    plot_data$is_highlighted <- plot_data$original_name %in% highlight_channels |
      plot_data$standard_name %in% highlight_channels
  } else {
    plot_data$is_highlighted <- FALSE
  }
  
  # ========================================================================
  # CREATE REALISTIC HEAD OUTLINE (ELLIPSOID)
  # ========================================================================
  
  # Create ellipsoid head shape (more realistic than sphere)
  # Semi-axes: a (left-right), b (front-back), c (top-bottom)
  a <- 80   # Width (left-right)
  b <- 95   # Length (front-back)
  c <- 95   # Height (top-bottom, slightly elongated)
  
  # Parametric ellipsoid surface
  n_points <- 30
  u <- seq(0, 2*pi, length.out = n_points)
  v <- seq(0, pi, length.out = n_points)
  
  # Create wireframe circles (longitude and latitude lines)
  head_lines <- list()
  
  # Longitude lines (vertical circles from front to back)
  for (i in seq(1, length(u), by = 3)) {
    x_line <- a * sin(v) * cos(u[i])
    y_line <- b * sin(v) * sin(u[i])
    z_line <- c * cos(v)
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # Latitude lines (horizontal circles)
  for (j in seq(1, length(v), by = 3)) {
    x_line <- a * sin(v[j]) * cos(u)
    y_line <- b * sin(v[j]) * sin(u)
    z_line <- rep(c * cos(v[j]), length(u))
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # ========================================================================
  # CREATE NOSE INDICATOR
  # ========================================================================
  
  nose_x <- c(0, 8, -8, 0, 0)
  nose_y <- c(100, 88, 88, 100, 100)
  nose_z <- c(20, 12, 12, 20, 20)
  
  # ========================================================================
  # BUILD PLOTLY 3D FIGURE
  # ========================================================================
  
  p <- plotly::plot_ly()
  
  # Add head wireframe (if requested)
  if (show_head) {
    for (line in head_lines) {
      p <- plotly::add_trace(
        p,
        x = line$x, y = line$y, z = line$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = "lightgray", width = 1),
        hoverinfo = "none",
        showlegend = FALSE
      )
    }
  }
  
  # Add nose indicator
  p <- plotly::add_trace(
    p,
    x = nose_x, y = nose_y, z = nose_z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 5),
    hoverinfo = "none",
    showlegend = FALSE,
    name = "Nose"
  )
  
  # Add electrodes (by color group)
  for (group in unique(plot_data$color_group)) {
    group_data <- plot_data[plot_data$color_group == group, ]
    
    marker_sizes <- ifelse(group_data$is_highlighted, 14, 10)
    marker_borders <- ifelse(group_data$is_highlighted, "black", "white")
    border_widths <- ifelse(group_data$is_highlighted, 3, 2)
    
    p <- plotly::add_trace(
      p,
      data = group_data,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = ifelse(show_labels, "markers+text", "markers"),
      marker = list(
        size = marker_sizes,
        color = color_palette[group],
        line = list(color = marker_borders, width = border_widths)
      ),
      text = ~original_name,
      textposition = "top center",
      textfont = list(size = 10, color = "black", family = "Arial Black"),
      hovertext = ~paste0("<b>", original_name, "</b><br>",
                          "Standard: ", standard_name, "<br>",
                          "Region: ", color_group, "<br>",
                          "Position (mm):<br>",
                          "  X: ", round(x, 1), " (L/R)<br>",
                          "  Y: ", round(y, 1), " (Back/Front)<br>",
                          "  Z: ", round(z, 1), " (Bottom/Top)"),
      hoverinfo = "text",
      name = group,
      showlegend = TRUE
    )
  }
  
  # ========================================================================
  # CONFIGURE LAYOUT
  # ========================================================================
  
  p <- plotly::layout(
    p,
    title = list(
      text = "BioSemi 64-Channel Electrode Layout (3D)",
      font = list(size = 20, family = "Arial", color = "black")
    ),
    scene = list(
      xaxis = list(
        title = "← Left | Right →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-110, 110)
      ),
      yaxis = list(
        title = "← Back | Front (Nose) →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-110, 110)
      ),
      zaxis = list(
        title = "← Bottom | Top →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-60, 110)
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      ),
      aspectmode = "cube",
      bgcolor = "white"
    ),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    legend = list(
      x = 0.02,
      y = 0.98,
      bgcolor = "rgba(255, 255, 255, 0.9)",
      bordercolor = "black",
      borderwidth = 1,
      font = list(size = 12)
    )
  )
  
  return(p)
}








#' ============================================================================
#' Plot 3D Electrode Layout on Scalp (BioSemi 64-Channel - Spherical)
#'
#' Creates an interactive 3D visualization of electrode positions on a spherical
#' head model using OFFICIAL BioSemi 64-channel SPHERICAL coordinates.
#'
#' @param channels Character vector - channel names to plot. If NULL, plots all
#'   64 standard 10-20/10-10 electrodes. Can also accept BioSemi A/B naming.
#' @param highlight_channels Character vector - specific channels to highlight
#'   in a different color (optional)
#' @param show_labels Logical - whether to show electrode labels on the plot.
#'   Default is TRUE.
#' @param color_by Character - how to color electrodes. Options:
#'   - "region" = color by brain region (frontal, central, parietal, occipital, temporal)
#'   - "single" = all one color (default blue)
#'   Default is "region".
#' @param show_head Logical - whether to show the head outline. Default is TRUE.
#' @param head_radius Numeric - radius of head sphere in mm. Default is 87.54 mm
#'   (calculated from 55 cm circumference).
#'
#' @return Returns a plotly object (interactive 3D plot viewable in RStudio)
#'
#' @details
#' This function visualizes electrode positions in 3D space using the official
#' BioSemi 64-channel electrode cap SPHERICAL coordinates (θ, φ), which are then
#' converted to Cartesian (x, y, z) for plotting.
#'
#' **Spherical Coordinate System:**
#' - θ (inclination): range from +180° to -180° (at Cz, inclination = 0°)
#'   - Positive inclinations = right side of head
#'   - Negative inclinations = left side of head
#' - φ (azimuth): range from +90° to -90°
#'   - For positive inclinations: φ = angle from T8 (right ear)
#'   - For negative inclinations: φ = angle from T7 (left ear)
#'   - Positive = anti-clockwise, negative = clockwise
#' - r (radius): 87.54 mm (from 55 cm head circumference)
#'
#' **Conversion to Cartesian:**
#' - x = r * sin(θ) * cos(φ)
#' - y = r * sin(θ) * sin(φ)
#' - z = r * cos(θ)
#'
#' **Coordinate System:**
#' - Origin (0, 0, 0) at center of head
#' - X-axis: Left (-) to Right (+)
#' - Y-axis: Back (-) to Front (+)
#' - Z-axis: Bottom (-) to Top (+)
#' - Nose points in +Y direction
#'
#' @examples
#' \dontrun{
#'
#'   # Example 1: Plot all 64 electrodes
#'   plot_electrode_3d_spherical()
#'
#'   # Example 2: Plot specific channels
#'   plot_electrode_3d_spherical(channels = c("Cz", "Fz", "Pz", "Oz", "T7", "T8"))
#'
#'   # Example 3: Highlight midline electrodes
#'   plot_electrode_3d_spherical(highlight_channels = c("Fpz", "Fz", "Cz", "Pz", "Oz"))
#'
#'   # Example 4: Different head size
#'   plot_electrode_3d_spherical(head_radius = 90)
#' }
#'
#' @keywords internal
#' @importFrom plotly plot_ly add_trace layout
#' @export
#'
plot_electrode_3d_spherical <- function(channels = NULL,
                                        highlight_channels = NULL,
                                        show_labels = TRUE,
                                        color_by = "region",
                                        show_head = TRUE,
                                        head_radius = 87.54) {
  
  # ========================================================================
  # LOAD REQUIRED PACKAGE
  # ========================================================================
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive 3D plots.\n",
         "Install it with: install.packages('plotly')")
  }
  
  # ========================================================================
  # OFFICIAL BIOSEMI 64-CHANNEL SPHERICAL COORDINATES
  # ========================================================================
  # From BioSemi Cap_coords_all.xls file
  # θ (inclination) in degrees, φ (azimuth) in degrees
  
  electrode_coords_spherical <- data.frame(
    electrode = c(
      "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3",
      "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1",
      "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz",
      "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4",
      "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4",
      "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8",
      "P10", "PO8", "PO4", "O2"
    ),
    theta = c(
      -92, -92, -74, -50, -60, -75, -92, -92, -72, -50,
      -32, -23, -46, -69, -92, -92, -72, -50, -32, -50,
      -60, -75, -92, -115, -92, -74, -92, 115, 92, 69,
      46, 23, 92, 92, 92, 74, 69, 46, 50, 60,
      75, 92, 92, 72, 50, 32, 23, 0, 23, 46,
      69, 92, 92, 72, 50, 32, 50, 60, 75, 92,
      115, 92, 74, 92
    ),
    phi = c(
      -72, -54, -65, -68, -51, -41, -36, -18, -21, -28,
      -45, 0, 0, 0, 0, 18, 21, 28, 45, 68,
      51, 41, 36, 36, 54, 65, 72, -90, -90, -90,
      -90, -90, 90, 72, 54, 65, 90, 90, 68, 51,
      41, 36, 18, 21, 28, 45, 90, 0, 0, 0,
      0, 0, -18, -21, -28, -45, -68, -51, -41, -36,
      -36, -54, -65, -72
    ),
    stringsAsFactors = FALSE
  )
  
  # ========================================================================
  # CONVERT SPHERICAL TO CARTESIAN COORDINATES
  # ========================================================================
  
  # Convert degrees to radians
  deg_to_rad <- function(deg) {
    return(deg * pi / 180)
  }
  
  # Spherical to Cartesian conversion
  # x = r * sin(θ) * cos(φ)
  # y = r * sin(θ) * sin(φ)
  # z = r * cos(θ)
  
  electrode_coords_spherical$theta_rad <- deg_to_rad(electrode_coords_spherical$theta)
  electrode_coords_spherical$phi_rad <- deg_to_rad(electrode_coords_spherical$phi)
  
  electrode_coords_spherical$x <- head_radius * sin(electrode_coords_spherical$theta_rad) * 
    cos(electrode_coords_spherical$phi_rad)
  electrode_coords_spherical$y <- head_radius * sin(electrode_coords_spherical$theta_rad) * 
    sin(electrode_coords_spherical$phi_rad)
  electrode_coords_spherical$z <- head_radius * cos(electrode_coords_spherical$theta_rad)
  
  # Round to match BioSemi precision
  electrode_coords_spherical$x <- round(electrode_coords_spherical$x, 0)
  electrode_coords_spherical$y <- round(electrode_coords_spherical$y, 0)
  electrode_coords_spherical$z <- round(electrode_coords_spherical$z, 0)
  
  # ========================================================================
  # BIOSEMI A/B TO STANDARD MAPPING
  # ========================================================================
  
  biosemi_mapping <- data.frame(
    biosemi = c(
      "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
      "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20",
      "A21", "A22", "A23", "A24", "A25", "A26", "A27", "A28", "A29", "A30",
      "A31", "A32",
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
      "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18", "B19", "B20",
      "B21", "B22", "B23", "B24", "B25", "B26", "B27", "B28", "B29", "B30",
      "B31", "B32"
    ),
    standard = c(
      "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3",
      "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1",
      "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz",
      "Pz", "CPz",
      "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8",
      "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8",
      "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8",
      "PO4", "O2"
    ),
    stringsAsFactors = FALSE
  )
  
  # ========================================================================
  # DETERMINE WHICH CHANNELS TO PLOT
  # ========================================================================
  
  if (is.null(channels)) {
    channels_to_plot <- electrode_coords_spherical$electrode
  } else {
    channels_to_plot <- channels
  }
  
  # Map BioSemi names to standard names if needed
  plot_data <- data.frame(
    original_name = channels_to_plot,
    standard_name = channels_to_plot,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(plot_data)) {
    ch <- plot_data$original_name[i]
    if (ch %in% biosemi_mapping$biosemi) {
      plot_data$standard_name[i] <- biosemi_mapping$standard[biosemi_mapping$biosemi == ch]
    }
  }
  
  # Merge with coordinates
  plot_data <- merge(plot_data, electrode_coords_spherical, 
                     by.x = "standard_name", by.y = "electrode",
                     all.x = TRUE)
  
  # Remove electrodes without coordinates
  missing_coords <- plot_data$original_name[is.na(plot_data$x)]
  if (length(missing_coords) > 0) {
    warning("No coordinates found for: ", paste(missing_coords, collapse = ", "))
    plot_data <- plot_data[!is.na(plot_data$x), ]
  }
  
  if (nrow(plot_data) == 0) {
    stop("No valid electrodes to plot.")
  }
  
  # ========================================================================
  # ASSIGN COLORS BY BRAIN REGION
  # ========================================================================
  
  if (color_by == "region") {
    plot_data$color_group <- "Other"
    plot_data$color_group[grepl("^Fp|^AF|^F[^CT]", plot_data$standard_name)] <- "Frontal"
    plot_data$color_group[grepl("^C|^FC", plot_data$standard_name)] <- "Central"
    plot_data$color_group[grepl("^P[^O]|^CP", plot_data$standard_name)] <- "Parietal"
    plot_data$color_group[grepl("^O|^PO|^Iz", plot_data$standard_name)] <- "Occipital"
    plot_data$color_group[grepl("^T[0-9]|^FT|^TP", plot_data$standard_name)] <- "Temporal"
    
    color_palette <- c(
      "Frontal" = "#e41a1c",
      "Central" = "#377eb8",
      "Parietal" = "#4daf4a",
      "Occipital" = "#984ea3",
      "Temporal" = "#ff7f00",
      "Other" = "#999999"
    )
    
  } else {
    plot_data$color_group <- "Electrode"
    color_palette <- c("Electrode" = "#377eb8")
  }
  
  # Highlight specific channels if requested
  if (!is.null(highlight_channels)) {
    plot_data$is_highlighted <- plot_data$original_name %in% highlight_channels |
      plot_data$standard_name %in% highlight_channels
  } else {
    plot_data$is_highlighted <- FALSE
  }
  
  # ========================================================================
  # CREATE PERFECT SPHERE HEAD OUTLINE (from spherical coordinates)
  # ========================================================================
  
  # Create perfect sphere wireframe
  n_points <- 30
  u <- seq(0, 2*pi, length.out = n_points)
  v <- seq(0, pi, length.out = n_points)
  
  # Create wireframe circles (longitude and latitude lines)
  head_lines <- list()
  
  # Longitude lines (vertical circles from front to back)
  for (i in seq(1, length(u), by = 3)) {
    x_line <- head_radius * sin(v) * cos(u[i])
    y_line <- head_radius * sin(v) * sin(u[i])
    z_line <- head_radius * cos(v)
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # Latitude lines (horizontal circles)
  for (j in seq(1, length(v), by = 3)) {
    x_line <- head_radius * sin(v[j]) * cos(u)
    y_line <- head_radius * sin(v[j]) * sin(u)
    z_line <- rep(head_radius * cos(v[j]), length(u))
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # ========================================================================
  # CREATE NOSE INDICATOR
  # ========================================================================
  
  nose_length <- head_radius * 0.15  # 15% of radius
  nose_x <- c(0, head_radius * 0.09, -head_radius * 0.09, 0, 0)
  nose_y <- c(head_radius + nose_length, head_radius, head_radius, head_radius + nose_length, head_radius + nose_length)
  nose_z <- c(head_radius * 0.23, head_radius * 0.14, head_radius * 0.14, head_radius * 0.23, head_radius * 0.23)
  
  # ========================================================================
  # BUILD PLOTLY 3D FIGURE
  # ========================================================================
  
  p <- plotly::plot_ly()
  
  # Add head wireframe (if requested)
  if (show_head) {
    for (line in head_lines) {
      p <- plotly::add_trace(
        p,
        x = line$x, y = line$y, z = line$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = "lightgray", width = 1),
        hoverinfo = "none",
        showlegend = FALSE
      )
    }
  }
  
  # Add nose indicator
  p <- plotly::add_trace(
    p,
    x = nose_x, y = nose_y, z = nose_z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 5),
    hoverinfo = "none",
    showlegend = FALSE,
    name = "Nose"
  )
  
  # Add electrodes (by color group)
  for (group in unique(plot_data$color_group)) {
    group_data <- plot_data[plot_data$color_group == group, ]
    
    marker_sizes <- ifelse(group_data$is_highlighted, 14, 10)
    marker_borders <- ifelse(group_data$is_highlighted, "black", "white")
    border_widths <- ifelse(group_data$is_highlighted, 3, 2)
    
    p <- plotly::add_trace(
      p,
      data = group_data,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = ifelse(show_labels, "markers+text", "markers"),
      marker = list(
        size = marker_sizes,
        color = color_palette[group],
        line = list(color = marker_borders, width = border_widths)
      ),
      text = ~original_name,
      textposition = "top center",
      textfont = list(size = 10, color = "black", family = "Arial Black"),
      hovertext = ~paste0("<b>", original_name, "</b><br>",
                          "Standard: ", standard_name, "<br>",
                          "Region: ", color_group, "<br>",
                          "Spherical (degrees):<br>",
                          "  θ (inclination): ", theta, "°<br>",
                          "  φ (azimuth): ", phi, "°<br>",
                          "Cartesian (mm):<br>",
                          "  X: ", round(x, 1), " (L/R)<br>",
                          "  Y: ", round(y, 1), " (Back/Front)<br>",
                          "  Z: ", round(z, 1), " (Bottom/Top)"),
      hoverinfo = "text",
      name = group,
      showlegend = TRUE
    )
  }
  
  # ========================================================================
  # CONFIGURE LAYOUT
  # ========================================================================
  
  axis_range <- head_radius * 1.3
  
  p <- plotly::layout(
    p,
    title = list(
      text = "BioSemi 64-Channel Electrode Layout (3D - Spherical)",
      font = list(size = 20, family = "Arial", color = "black")
    ),
    scene = list(
      xaxis = list(
        title = "← Left | Right →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range, axis_range)
      ),
      yaxis = list(
        title = "← Back | Front (Nose) →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range, axis_range)
      ),
      zaxis = list(
        title = "← Bottom | Top →",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range * 0.7, axis_range)
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      ),
      aspectmode = "cube",
      bgcolor = "white"
    ),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    legend = list(
      x = 0.02,
      y = 0.98,
      bgcolor = "rgba(255, 255, 255, 0.9)",
      bordercolor = "black",
      borderwidth = 1,
      font = list(size = 12)
    )
  )
  
  return(p)
}

