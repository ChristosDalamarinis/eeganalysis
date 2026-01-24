#' ==============================================================================
#' Extract EVENT codes from BDF STATUS channel
#' ==============================================================================
#'
#' Function to read BioSemi BDF files and extract trigger/event codes 
#' from the STATUS channel using native R implementation (no edfReader dependency)
#'
#' Author: Christos Dalamarinis
#' Date: Jan - s2026
#' 
#' BioSemi BDF STATUS Channel Structure (per BioSemi format spec):
#' - Bits 1-16: Trigger codes (up to 65535)
#' - Bits 16-23: System status (CMS, battery, etc.)
#' - Each sample = 24-bit value (3 bytes, little-endian)
#' 
#' Reference: https://www.biosemi.com/faq/file_format.htm
#' ==============================================================================
#'
#' Extract EVENT codes from BDF STATUS channel (native)
#'
#' Convenience wrapper around \code{read_bdf_native()} and
#' \code{extract_events_native()} to return only the events table.
#'
#' Requires that \code{read_bdf_native()} and \code{extract_events_native()}
#' from \code{read_bdf_native.R} are available in the package.
#'
#' @param bdf_file Path to BDF file (character string).
#' @param verbose Logical. If TRUE (default), print detailed progress messages.
#'
#' @return Data frame with columns:
#'   - onset: sample index of the event
#'   - onset_time: time in seconds
#'   - type: trigger code as character
#'   - description: text label (e.g. "Trigger: 255")
#'
#' The returned data frame has attributes:
#'   - bdf_file: file path
#'   - sample_rate: sampling rate in Hz
#'   - n_samples: total number of samples
#'   - duration_sec: total duration in seconds
#'
#' @examples
#' \dontrun{
#'   events <- extract_bdf_events("sub-01_task_eeg.bdf")
#'   head(events)
#'   table(events$type)
#' }
#'
#' @export
extract_bdf_events <- function(bdf_file, verbose = TRUE) {
  
  # ========== STEP 1: Validate File ==========
  if (verbose) cat("\n========== BDF Event Extraction ==========\n")
  if (verbose) cat("[1/5] Validating BDF file...\n")
  
  if (!file.exists(bdf_file)) {
    stop("BDF file not found: ", bdf_file, call. = FALSE)
  }
  
  if (!grepl("\\.bdf$", bdf_file, ignore.case = TRUE)) {
    warning("File does not have .bdf extension.")
  }
  
  file_size_mb <- file.info(bdf_file)$size / (1024^2)
  if (verbose) {
    cat("  ✓ File found: ", basename(bdf_file), "\n", sep = "")
    cat("  ✓ File size: ", round(file_size_mb, 1), " MB\n", sep = "")
  }
  
  # ========== STEP 2: Read BDF File ==========
  if (verbose) cat("\n[2/5] Reading BDF file (native reader)...\n")
  
  start_time <- Sys.time()
  
  # Use native reader (from read_bdf_native.R)
  # Note: read_bdf_native has its own verbose output, so we pass verbose through
  eeg_obj <- read_bdf_native(bdf_file, verbose = verbose)
  
  read_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (verbose && !verbose) {  # Only print this if read_bdf_native wasn't verbose
    cat("  ✓ File read successfully in ", round(read_time, 2), " seconds\n", sep = "")
  }
  
  # ========== STEP 3: Extract Events ==========
  if (verbose) cat("\n[3/5] Extracting events from EEG object...\n")
  
  events <- eeg_obj$events
  
  if (is.null(events)) {
    warning("No events structure found in EEG object for file: ", bdf_file)
    events <- data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  if (verbose) {
    if (nrow(events) == 0) {
      cat("  ⚠ No events found in STATUS channel\n")
    } else {
      cat("  ✓ Events extracted: ", nrow(events), "\n", sep = "")
    }
  }
  
  # ========== STEP 4: Add Metadata Attributes ==========
  if (verbose) cat("\n[4/5] Adding metadata attributes...\n")
  
  attr(events, "bdf_file")    <- bdf_file
  attr(events, "sample_rate") <- eeg_obj$sampling_rate
  attr(events, "n_samples")   <- length(eeg_obj$times)
  attr(events, "duration_sec")<- max(eeg_obj$times)
  
  if (verbose) {
    cat("  ✓ Sample rate: ", eeg_obj$sampling_rate, " Hz\n", sep = "")
    cat("  ✓ Total samples: ", length(eeg_obj$times), "\n", sep = "")
    cat("  ✓ Duration: ", round(max(eeg_obj$times), 2), " seconds\n", sep = "")
  }
  
  # ========== STEP 5: Summary Statistics ==========
  if (verbose && nrow(events) > 0) {
    cat("\n[5/5] Event extraction summary:\n")
    
    unique_codes <- unique(events$type)
    cat("  • Total events: ", nrow(events), "\n", sep = "")
    cat("  • Unique trigger codes: ", length(unique_codes), "\n", sep = "")
    
    if (length(unique_codes) <= 10) {
      cat("  • Trigger codes: ", paste(sort(as.numeric(unique_codes)), collapse = ", "), "\n", sep = "")
    } else {
      cat("  • Trigger code range: ", min(as.numeric(unique_codes)), 
          " to ", max(as.numeric(unique_codes)), "\n", sep = "")
    }
    
    cat("  • First event at: ", round(min(events$onset_time), 3), " sec\n", sep = "")
    cat("  • Last event at: ", round(max(events$onset_time), 3), " sec\n", sep = "")
    
    if (nrow(events) > 1) {
      mean_interval <- mean(diff(events$onset_time))
      cat("  • Mean inter-event interval: ", round(mean_interval, 3), " sec\n", sep = "")
      cat("  • Approximate event rate: ", round(1/mean_interval, 2), " events/sec\n", sep = "")
    }
    
    cat("\n✓ Event extraction complete!\n")
    cat("==========================================\n\n")
  } else if (verbose) {
    cat("\n[5/5] Complete (no events found)\n")
    cat("==========================================\n\n")
  }
  
  return(events)
}


#' Summary of extracted BDF events
#'
#' @param events Data frame returned from \code{extract_bdf_events()}
#' @export
summary_bdf_events <- function(events) {
  
  cat("====== BDF Event Summary ======\n")
  cat("File: ", attr(events, "bdf_file"), "\n", sep = "")
  cat("Sample rate: ", attr(events, "sample_rate"), " Hz\n", sep = "")
  cat("Recording duration: ",
      round(attr(events, "duration_sec"), 2), " seconds\n", sep = "")
  cat("Total samples: ", attr(events, "n_samples"), "\n", sep = "")
  cat("\n")
  
  cat("Event Summary:\n")
  cat("  Total events: ", nrow(events), "\n", sep = "")
  cat("  Unique trigger codes: ", length(unique(events$type)), "\n", sep = "")
  cat("\n")
  
  cat("Trigger Code Distribution:\n")
  print(table(events$type))
  
  if (nrow(events) > 0) {
    cat("\nTiming:\n")
    cat("  First event at: ",
        round(min(events$onset_time), 3), " seconds\n", sep = "")
    cat("  Last event at: ",
        round(max(events$onset_time), 3), " seconds\n", sep = "")
    if (nrow(events) > 1) {
      cat("  Mean interval: ",
          round(mean(diff(events$onset_time)), 3), " seconds\n", sep = "")
    }
  }
}


#' Plot BDF event timeline
#'
#' @param events Data frame returned from \code{extract_bdf_events()}
#' @param show_codes Logical. If TRUE, label each event with its trigger code.
#' @export
plot_bdf_events <- function(events, show_codes = FALSE) {
  
  if (is.null(events) || !nrow(events)) {
    warning("No events to plot.")
    return(invisible(NULL))
  }
  
  plot(events$onset_time,
       as.numeric(events$type),
       type = "h",
       xlab = "Time (seconds)",
       yl
       
#' Extract EVENT codes from BDF STATUS channel (native)
#'
#' Convenience wrapper around \code{read_bdf_native()} and
#' \code{extract_events_native()} to return only the events table.
#'
#' Requires that \code{read_bdf_native()} and \code{extract_events_native()}
#' from \code{read_bdf_native.R} are available in the package.
#'
#' @param bdf_file Path to BDF file (character string).
#' @param verbose Logical. If TRUE (default), print detailed progress messages.
#'
#' @return Data frame with columns:
#'   - onset: sample index of the event
#'   - onset_time: time in seconds
#'   - type: trigger code as character
#'   - description: text label (e.g. "Trigger: 255")
#'
#' The returned data frame has attributes:
#'   - bdf_file: file path
#'   - sample_rate: sampling rate in Hz
#'   - n_samples: total number of samples
#'   - duration_sec: total duration in seconds
#'
#' @examples
#' \dontrun{
#'   events <- extract_bdf_events("sub-01_task_eeg.bdf")
#'   head(events)
#'   table(events$type)
#' }
#'
#' @export
extract_bdf_events <- function(bdf_file, verbose = TRUE) {
  
  # ========== STEP 1: Validate File ==========
  if (verbose) cat("\n========== BDF Event Extraction ==========\n")
  if (verbose) cat("[1/5] Validating BDF file...\n")
  
  if (!file.exists(bdf_file)) {
    stop("BDF file not found: ", bdf_file, call. = FALSE)
  }
  
  if (!grepl("\\.bdf$", bdf_file, ignore.case = TRUE)) {
    warning("File does not have .bdf extension.")
  }
  
  file_size_mb <- file.info(bdf_file)$size / (1024^2)
  if (verbose) {
    cat("  ✓ File found: ", basename(bdf_file), "\n", sep = "")
    cat("  ✓ File size: ", round(file_size_mb, 1), " MB\n", sep = "")
  }
  
  # ========== STEP 2: Read BDF File ==========
  if (verbose) cat("\n[2/5] Reading BDF file (native reader)...\n")
  
  start_time <- Sys.time()
  
  # Use native reader (from read_bdf_native.R)
  # Note: read_bdf_native has its own verbose output, so we pass verbose through
  eeg_obj <- read_bdf_native(bdf_file, verbose = verbose)
  
  read_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (verbose && !verbose) {  # Only print this if read_bdf_native wasn't verbose
    cat("  ✓ File read successfully in ", round(read_time, 2), " seconds\n", sep = "")
  }
  
  # ========== STEP 3: Extract Events ==========
  if (verbose) cat("\n[3/5] Extracting events from EEG object...\n")
  
  events <- eeg_obj$events
  
  if (is.null(events)) {
    warning("No events structure found in EEG object for file: ", bdf_file)
    events <- data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  if (verbose) {
    if (nrow(events) == 0) {
      cat("  ⚠ No events found in STATUS channel\n")
    } else {
      cat("  ✓ Events extracted: ", nrow(events), "\n", sep = "")
    }
  }
  
  # ========== STEP 4: Add Metadata Attributes ==========
  if (verbose) cat("\n[4/5] Adding metadata attributes...\n")
  
  attr(events, "bdf_file")    <- bdf_file
  attr(events, "sample_rate") <- eeg_obj$sampling_rate
  attr(events, "n_samples")   <- length(eeg_obj$times)
  attr(events, "duration_sec")<- max(eeg_obj$times)
  
  if (verbose) {
    cat("  ✓ Sample rate: ", eeg_obj$sampling_rate, " Hz\n", sep = "")
    cat("  ✓ Total samples: ", length(eeg_obj$times), "\n", sep = "")
    cat("  ✓ Duration: ", round(max(eeg_obj$times), 2), " seconds\n", sep = "")
  }
  
  # ========== STEP 5: Summary Statistics ==========
  if (verbose && nrow(events) > 0) {
    cat("\n[5/5] Event extraction summary:\n")
    
    unique_codes <- unique(events$type)
    cat("  • Total events: ", nrow(events), "\n", sep = "")
    cat("  • Unique trigger codes: ", length(unique_codes), "\n", sep = "")
    
    if (length(unique_codes) <= 10) {
      cat("  • Trigger codes: ", paste(sort(as.numeric(unique_codes)), collapse = ", "), "\n", sep = "")
    } else {
      cat("  • Trigger code range: ", min(as.numeric(unique_codes)), 
          " to ", max(as.numeric(unique_codes)), "\n", sep = "")
    }
    
    cat("  • First event at: ", round(min(events$onset_time), 3), " sec\n", sep = "")
    cat("  • Last event at: ", round(max(events$onset_time), 3), " sec\n", sep = "")
    
    if (nrow(events) > 1) {
      mean_interval <- mean(diff(events$onset_time))
      cat("  • Mean inter-event interval: ", round(mean_interval, 3), " sec\n", sep = "")
      cat("  • Approximate event rate: ", round(1/mean_interval, 2), " events/sec\n", sep = "")
    }
    
    cat("\n✓ Event extraction complete!\n")
    cat("==========================================\n\n")
  } else if (verbose) {
    cat("\n[5/5] Complete (no events found)\n")
    cat("==========================================\n\n")
  }
  
  return(events)
}


#' Summary of extracted BDF events
#'
#' @param events Data frame returned from \code{extract_bdf_events()}
#' @export
summary_bdf_events <- function(events) {
  
  cat("====== BDF Event Summary ======\n")
  cat("File: ", attr(events, "bdf_file"), "\n", sep = "")
  cat("Sample rate: ", attr(events, "sample_rate"), " Hz\n", sep = "")
  cat("Recording duration: ",
      round(attr(events, "duration_sec"), 2), " seconds\n", sep = "")
  cat("Total samples: ", attr(events, "n_samples"), "\n", sep = "")
  cat("\n")
  
  cat("Event Summary:\n")
  cat("  Total events: ", nrow(events), "\n", sep = "")
  cat("  Unique trigger codes: ", length(unique(events$type)), "\n", sep = "")
  cat("\n")
  
  cat("Trigger Code Distribution:\n")
  print(table(events$type))
  
  if (nrow(events) > 0) {
    cat("\nTiming:\n")
    cat("  First event at: ",
        round(min(events$onset_time), 3), " seconds\n", sep = "")
    cat("  Last event at: ",
        round(max(events$onset_time), 3), " seconds\n", sep = "")
    if (nrow(events) > 1) {
      cat("  Mean interval: ",
          round(mean(diff(events$onset_time)), 3), " seconds\n", sep = "")
    }
  }
}

#' Plot BDF event timeline
#'
#' @param events Data frame returned from \code{extract_bdf_events()}
#' @param show_codes Logical. If TRUE, label each event with its trigger code.
#' @export
plot_bdf_events <- function(events, show_codes = FALSE) {
  
  if (is.null(events) || !nrow(events)) {
    warning("No events to plot.")
    return(invisible(NULL))
  }
  
  plot(events$onset_time,
       as.numeric(events$type),
       type = "h",
       xlab = "Time (seconds)",
       ylab = "Trigger code",
       main = paste("BDF Events from",
                    basename(attr(events, "bdf_file"))),
       col = "darkblue",
       lwd = 2)
  
  points(events$onset_time,
         as.numeric(events$type),
         pch = 19,
         col = "red",
         cex = 0.8)
  
  if (show_codes && nrow(events) < 50) {
    text(events$onset_time,
         as.numeric(events$type),
         labels = events$type,
         pos = 3,
         cex = 0.7)
  }
  
  grid()
}
