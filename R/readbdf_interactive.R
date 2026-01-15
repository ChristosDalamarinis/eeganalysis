#' ============================================================================
#' BioSemi BDF File Reader Functions (with Smart Downsampling & Chunked Backend)
#' ============================================================================
#'
#' This module provides functions to read and import EEG data recorded from
#' BioSemi systems (BDF format). BDF (BioSemi Data Format) is a binary format
#' for storing EEG data with metadata and event markers.
#'
#' BioSemi systems typically use:
#' - CMS/DRL reference scheme (Common Mode Sense / Driven Right Leg)
#' - High sampling rates (256-2048 Hz typical)
#' - 64, 128, or 256 electrode channels
#' - Status channel for event markers/triggers
#'
#' This version adds:
#' - Smart interactive downsampling decisions
#' - A fast, chunked low-level reader for huge/high-rate files
#' - Seamless fallback to edfReader::readEdfSignals() for standard files
#'
#' Author: Christos Dalamarinis
#' Date: Jan - 2026
#' ============================================================================

#' Internal Helper: Downsample During Import
#'
#' Simplified downsampling optimized for import process.
#' Uses decimate-like method with anti-aliasing filter to prevent spectral aliasing.
#'
#' @param eeg_obj EEG object to downsample
#' @param target_rate Target sampling rate in Hz
#' @param verbose Print progress messages
#' @return Downsampled eeg object
#' @keywords internal
downsample_during_import <- function(eeg_obj, target_rate, verbose = TRUE) {
  
  # Calculate downsampling parameters
  downsample_factor <- round(eeg_obj$sampling_rate / target_rate)
  actual_rate <- eeg_obj$sampling_rate / downsample_factor
  
  if (verbose) {
    cat("\n")
    cat(strrep("=", 75), "\n")
    cat("DOWNSAMPLING\n")
    cat(strrep("=", 75), "\n")
  }
  
  # Anti-aliasing filter parameters
  new_nyquist <- actual_rate / 2
  filter_cutoff <- new_nyquist * 0.9
  
  if (verbose) {
    cat(" Original rate:", eeg_obj$sampling_rate, "Hz\n")
    cat(" Target rate:", target_rate, "Hz\n")
    cat(" Actual rate:", actual_rate, "Hz\n")
    cat(" Downsample factor:", downsample_factor, "\n")
    cat(" Filter cutoff:", round(filter_cutoff, 1), "Hz\n\n")
  }
  
  # Require signal package
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop(
      "'signal' package required for downsampling.\n",
      "Install with: install.packages('signal')",
      call. = FALSE
    )
  }
  
  # Design Butterworth filter
  nyquist_freq <- eeg_obj$sampling_rate / 2
  normalized_cutoff <- filter_cutoff / nyquist_freq
  if (normalized_cutoff >= 1) {
    stop(
      "Filter cutoff too high. This should not happen - please report as bug.",
      call. = FALSE
    )
  }
  
  if (verbose) {
    cat(" [1/3] Applying anti-aliasing filter...\n")
  }
  
  butter_filter <- signal::butter(n = 8, W = normalized_cutoff, type = "low")
  
  # Apply filter and decimate
  n_channels <- nrow(eeg_obj$data)
  n_timepoints <- ncol(eeg_obj$data)
  filtered_data <- matrix(0, nrow = n_channels, ncol = n_timepoints)
  
  for (ch in 1:n_channels) {
    filtered_data[ch, ] <- signal::filtfilt(butter_filter, eeg_obj$data[ch, ])
  }
  
  if (verbose) {
    cat(" [2/3] Decimating signal...\n")
  }
  
  # Decimate
  decimate_indices <- seq(1, n_timepoints, by = downsample_factor)
  downsampled_data <- filtered_data[, decimate_indices, drop = FALSE]
  new_times <- eeg_obj$times[decimate_indices]
  
  if (verbose) {
    cat(" [3/3] Updating events...\n")
  }
  
  # Update events
  new_events <- eeg_obj$events
  if (nrow(new_events) > 0) {
    # Adjust onset indices
    new_events$onset <- round(new_events$onset / downsample_factor)
    # Clamp to valid range
    new_events$onset <- pmax(1, pmin(new_events$onset, length(new_times)))
    # Update onset_time if present
    if ("onset_time" %in% names(new_events)) {
      new_events$onset_time <- new_times[new_events$onset]
    }
  }
  
  # Update preprocessing history
  preprocessing_step <- paste0(
    "Downsampled during import from ", eeg_obj$sampling_rate, " Hz to ",
    actual_rate, " Hz (factor: ", downsample_factor, ")"
  )
  
  if (verbose) {
    cat(
      "\n Samples reduced:", n_timepoints, "→", length(new_times),
      "(", round((1 - length(new_times) / n_timepoints) * 100, 1),
      "% reduction)\n"
    )
    cat(" Events adjusted:", nrow(new_events), "events\n")
    cat(strrep("=", 75), "\n\n")
  }
  
  # Create new object
  new_eeg(
    data = downsampled_data,
    channels = eeg_obj$channels,
    sampling_rate = actual_rate,
    times = new_times,
    events = new_events,
    metadata = eeg_obj$metadata,
    reference = eeg_obj$reference,
    preprocessing_history = c(eeg_obj$preprocessing_history, list(preprocessing_step))
  )
}

#' Read BioSemi BDF File (Interactive + Chunked Backend)
#'
#' This is the main function to import EEG data from a BioSemi .bdf file.
#' It reads the binary file, extracts all EEG channels and metadata,
#' parses event triggers, and returns a properly structured eeg object
#' ready for analysis.
#'
#' It can:
#' - Use edfReader::readEdfSignals() (standard backend)
#' - Use a chunked low-level reader for huge/high-rate files (fast path)
#' - Prompt intelligently about downsampling based on sampling rate and file size
#'
#' @param file_path Character string with full path to the .bdf file
#'   Example: "/Users/user/data/subject_01.bdf"
#' @param interactive Logical. If TRUE (default), prompts user about downsampling
#'   when high sampling rates (>512 Hz) are detected. Set FALSE for batch processing.
#' @param downsample_to Numeric. Target sampling rate in Hz for automatic downsampling.
#'   If NULL (default), uses interactive prompts. If specified, automatically
#'   downsamples without prompting.
#'
#' @details
#' The function performs the following steps:
#' 1. Validates the file exists and has .bdf extension
#' 2. Reads the BDF header using edfReader::readEdfHeader()
#' 3. Chooses between edfReader backend and chunked low-level reader
#' 4. Extracts channel names, sampling rate, and recording metadata
#' 5. Separates EEG channels from status/trigger channel (if using edfReader path)
#' 6. Extracts and parses event markers from status channel
#' 7. (Optional) Intelligently prompts for downsampling if rate > 512 Hz
#' 8. Returns an eeg object with all data and metadata
#'
#' The status channel is typically:
#' - Named "Status" in BioSemi recordings
#' - Contains trigger codes corresponding to experimental events
#' - Bit values encode different triggers (0, 1, 2, 4, 8, etc.)
#' - Used to epoch data around events of interest
#'
#' @return An object of class 'eeg' containing:
#' \describe{
#'   \item{data}{Matrix of EEG signal values (channels × time points)}
#'   \item{channels}{Character vector of EEG channel names}
#'   \item{sampling_rate}{Numeric - sampling rate in Hz}
#'   \item{times}{Numeric vector - time points in seconds}
#'   \item{events}{Data frame - event markers (onset, type, description)}
#'   \item{metadata}{List - recording information and file metadata}
#'   \item{reference}{Character - reference scheme (Biosemi CMS/DRL)}
#'   \item{preprocessing_history}{List - processing steps applied}
#' }
#'
#' @examples
#' \dontrun{
#' # Load with interactive downsampling prompt
#' eeg <- read_biosemi_interactive("~/data/experiment/subject_01.bdf")
#'
#' # Load and downsample automatically to 512 Hz
#' eeg <- read_biosemi_interactive("~/data/subject_01.bdf", downsample_to = 512)
#'
#' # Batch processing without prompts
#' files <- list.files(pattern = "\\.bdf$")
#' data_list <- lapply(files, function(f) {
#'   read_biosemi_interactive(f, interactive = FALSE, downsample_to = 256)
#' })
#'
#' # View summary information
#' print(eeg)
#'
#' # Access components
#' head(eeg$data)        # First rows of EEG data
#' eeg$channels          # Channel names
#' eeg$sampling_rate     # Sampling rate
#' eeg$events            # Event table
#' eeg$metadata          # Recording metadata
#' }
#'
#' @importFrom edfReader readEdfHeader readEdfSignals
#' @importFrom signal butter filtfilt
#' @export
read_biosemi_interactive <- function(file_path,
                                     interactive = TRUE,
                                     downsample_to = NULL) {
  
  # ========== INPUT VALIDATION ==========
  if (!file.exists(file_path)) {
    stop("ERROR: File not found at path: ", file_path, call. = FALSE)
  }
  
  if (!grepl("\\.bdf$", file_path, ignore.case = TRUE)) {
    warning(
      "WARNING: File does not have .bdf extension. ",
      "Expected BioSemi BDF format. Attempting to read anyway..."
    )
  }
  
  # ========== READ BDF FILE HEADER ==========
  cat("[1/5] Reading BDF header...")
  header <- tryCatch({
    edfReader::readEdfHeader(file_path)
  }, error = function(e) {
    stop(
      "ERROR reading BDF file header: ", e$message,
      "\nCheck that the file is a valid BioSemi BDF file.",
      call. = FALSE
    )
  })
  cat(" Done!\n")
  
  # ========== EXTRACT BASIC INFO FROM HEADER ==========
  channels_raw   <- header$sHeaders$label
  channels_clean <- trimws(channels_raw)
  sampling_rate  <- header$sHeaders$sRate[1]
  n_channels     <- length(channels_clean)
  
  # Estimate file size and duration
  n_samples_per_channel <- header$sHeaders$sLength[1]
  duration_sec      <- n_samples_per_channel / sampling_rate
  estimated_size_mb <- (n_channels * n_samples_per_channel * 8) / (1024^2)
  
  cat(" Channels: ", n_channels, "\n", sep = "")
  cat(" Sampling rate: ", sampling_rate, " Hz\n", sep = "")
  cat(" Duration: ~", round(duration_sec, 1), " seconds\n", sep = "")
  cat(" Estimated size: ~", round(estimated_size_mb, 1), " MB\n", sep = "")
  
  # ========== EARLY DOWNSAMPLING PROMPT ==========
  # Determine target rate BEFORE reading the large data file
  target_downsample_rate <- downsample_to # Start with parameter value
  
  if (interactive && is.null(downsample_to) && sampling_rate > 512) {
    cat("\n")
    cat(strrep("=", 75), "\n")
    cat("HIGH SAMPLING RATE DETECTED\n")
    cat(strrep("=", 75), "\n")
    cat(" Current sampling rate:", sampling_rate, "Hz\n")
    cat(" Estimated data size: ~", round(estimated_size_mb, 1), "MB\n")
    cat(" Duration: ~", round(duration_sec, 1), "seconds\n\n")
    cat("⚠ NOTICE: Reading this file at full resolution may take several minutes.\n\n")
    
    # Frequency preservation guidance
    cat("FREQUENCY PRESERVATION GUIDE:\n")
    cat(" • 512 Hz → preserves up to ~230 Hz (high gamma, RIFT/SSVEP)\n")
    cat(" • 256 Hz → preserves up to ~115 Hz (standard EEG + some high gamma)\n")
    cat(" • 128 Hz → preserves up to ~58 Hz (standard EEG bands only)\n\n")
    
    cat("⚠ IMPORTANT WARNINGS:\n")
    cat(" • DO NOT downsample if analyzing RIFT/SSVEP (60-120 Hz)\n")
    cat(" • DO NOT downsample if analyzing high gamma (>60 Hz)\n")
    cat(" • DO NOT downsample if analyzing ripples (80-250 Hz)\n\n")
    
    cat("RECOMMENDATIONS:\n")
    cat(" • For standard ERP analysis (0.5-50 Hz): 256 Hz is sufficient\n")
    cat(" • For time-frequency analysis up to 80 Hz: 512 Hz minimum\n")
    cat(" • For high-frequency research: keep original sampling rate\n\n")
    
    cat("Would you like to downsample during import for faster loading? (y/n): ")
    user_response <- tolower(trimws(readline()))
    
    if (user_response %in% c("y", "yes")) {
      cat("\nEnter target sampling rate in Hz (e.g., 256, 512): ")
      target_input <- trimws(readline())
      target_rate  <- suppressWarnings(as.numeric(target_input))
      if (!is.na(target_rate) && target_rate > 0 && target_rate < sampling_rate) {
        target_downsample_rate <- target_rate
        cat("\n✓ Will downsample to", target_rate, "Hz during import\n")
      } else {
        cat("\n✗ Invalid target rate. Will load at original sampling rate.\n")
      }
    } else {
      cat("\n✓ Will load at original sampling rate:", sampling_rate, "Hz\n")
    }
    
    cat(strrep("=", 75), "\n\n")
  }
  
  # ========== CHOOSE IMPORT BACKEND ==========
  use_chunked <- FALSE
  if (sampling_rate >= 2048 && estimated_size_mb > 1500) {
    if (interactive && is.null(downsample_to)) {
      cat("\n")
      cat(strrep("=", 75), "\n")
      cat("LARGE HIGH-RATE FILE DETECTED\n")
      cat(strrep("=", 75), "\n")
      cat(" Estimated size: ~", round(estimated_size_mb, 1), " MB\n", sep = "")
      cat(" Sampling rate :", sampling_rate, " Hz\n")
      cat("This may load much faster with the chunked reader.\n")
      cat("Use chunked low-level reader instead of edfReader? (y/n): ")
      resp <- tolower(trimws(readline()))
      use_chunked <- resp %in% c("y", "yes")
    } else {
      # Non-interactive: default to chunked for huge files
      use_chunked <- TRUE
    }
  }
  
  # ========== READ BDF SIGNAL DATA ==========
  if (!use_chunked) {
    # ---------- ORIGINAL PATH (edfReader) ----------
    cat("[2/5] Reading EEG signals...")
    signals <- tryCatch({
      edfReader::readEdfSignals(header)
    }, error = function(e) {
      stop("ERROR reading BDF signals: ", e$message, call. = FALSE)
    })
    cat(" Done!\n")
    
    # ========== CONVERT SIGNALS TO MATRIX ==========
    cat("[3/5] Converting to matrix...")
    n_channels <- length(signals)
    sample_lengths <- sapply(1:n_channels, function(i) length(signals[[i]]$signal))
    if (length(unique(sample_lengths)) > 1) {
      stop("ERROR: Channels have different numbers of samples!", call. = FALSE)
    }
    n_samples <- sample_lengths[1]
    eeg_matrix <- do.call(rbind, lapply(signals, function(sig) sig$signal))
    cat(" Done!\n")
    
    # ========== IDENTIFY AND SEPARATE STATUS CHANNEL ==========
    cat("[4/5] Processing channels...")
    status_channel_idx <- which(tolower(channels_clean) == "status")
    if (length(status_channel_idx) == 0) {
      status_channel_idx <- n_channels
    }
    
    status_signal <- as.numeric(eeg_matrix[status_channel_idx[1], ])
    eeg_channels_idx <- setdiff(1:n_channels, status_channel_idx)
    eeg_data <- eeg_matrix[eeg_channels_idx, , drop = FALSE]
    channels_final <- channels_clean[eeg_channels_idx]
    cat(" Done!\n")
    
    # ========== CREATE TIME VECTOR ==========
    times <- (0:(n_samples - 1)) / sampling_rate
    
    # ========== EXTRACT EVENTS ==========
    cat("[5/5] Extracting events...")
    events_df <- extract_biosemi_events(status_signal, sampling_rate)
    cat(" Done!\n")
    cat(" Events found: ", nrow(events_df), "\n\n", sep = "")
    
  } else {
    # ---------- CHUNKED FAST PATH ----------
    cat("[2/5] Reading EEG signals with chunked reader...\n")
    chunk_res <- read_bdf_chunked_lowlevel(
      file_path   = file_path,
      header      = header,
      target_rate = target_downsample_rate,
      chunk_mb    = 200
    )
    
    eeg_data      <- chunk_res$eeg_data      # channels x samples, EEG channels only
    channels_final <- chunk_res$channels
    sampling_rate <- chunk_res$srate
    status_signal <- chunk_res$status
    times         <- chunk_res$times
    events_df     <- chunk_res$events
    n_channels    <- nrow(eeg_data)
    n_samples     <- ncol(eeg_data)
    
    cat("[3/5] Converting to matrix... Done! (chunked backend)\n")
    cat("[4/5] Processing channels... Done! (status already separated)\n")
    cat("[5/5] Extracting events... Done! (from status in chunked backend)\n")
    cat(" Events found: ", nrow(events_df), "\n\n", sep = "")
  }
  
  # ========== PREPARE METADATA ==========
  metadata <- list(
    file_path = file_path,
    file_name = basename(file_path),
    date_recorded = tryCatch(header$startDate, error = function(e) NA),
    time_recorded = tryCatch(header$startTime, error = function(e) NA),
    device = "BioSemi",
    n_channels_recorded = length(channels_raw),
    n_eeg_channels = n_channels,
    status_channel_name = NA,  # status already handled
    reference_scheme = "Biosemi CMS/DRL (Common Mode Sense / Driven Right Leg)",
    duration_seconds = max(times),
    n_samples = n_samples,
    n_events = nrow(events_df)
  )
  
  # ========== CREATE EEG OBJECT ==========
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
  
  cat("✓ Successfully imported BioSemi BDF file!\n")
  
  # ========== APPLY DOWNSAMPLING IF REQUESTED (edfReader path) ==========
  if (!use_chunked &&
      !is.null(target_downsample_rate) &&
      target_downsample_rate < sampling_rate) {
    cat("\nApplying downsampling...\n")
    eeg_obj <- downsample_during_import(
      eeg_obj,
      target_rate = target_downsample_rate,
      verbose = TRUE
    )
  }
  
  return(eeg_obj)
}

#' Internal: unpack 24-bit signed integers from raw vector
#' @keywords internal
unpack_int24 <- function(raw_vec) {
  n_ints <- length(raw_vec) / 3L
  if (n_ints == 0L) return(integer(0))
  m <- matrix(as.integer(raw_vec), nrow = 3L)
  vals <- m[1, ] + bitwShiftL(m[2, ], 8L) + bitwShiftL(m[3, ], 16L)
  # convert from unsigned 24-bit to signed
  neg_mask <- vals >= 2^23
  vals[neg_mask] <- vals[neg_mask] - 2^24
  vals
}

#' Low-level chunked BDF reader (full-rate events, optional downsampling)
#'
#' Reads BDF data in chunks to avoid loading the full file into memory.
#' 1) Builds full-rate EEG and Status.
#' 2) Detects events on the full-rate Status.
#' 3) If requested, downsamples EEG (and times) only.
#' 4) Rescales event onsets to the new sampling rate.
#'
#' @param file_path Path to BDF file
#' @param header Edf header object from edfReader::readEdfHeader()
#' @param target_rate Optional target sampling rate (Hz). If NULL, keep full rate.
#' @param chunk_mb Chunk size in megabytes for reading
#' @keywords internal
read_bdf_chunked_lowlevel <- function(file_path, header,
                                      target_rate = NULL,
                                      chunk_mb = 200) {
  
  con <- file(file_path, "rb")
  on.exit(close(con), add = TRUE)
  
  # ----- basic metadata from header -----
  labels <- trimws(header$sHeaders$label)
  srate_full <- header$sHeaders$sRate[1]
  nsamp_full <- header$sHeaders$sLength[1]
  nchan      <- length(labels)
  
  # Status channel index
  status_idx <- which(tolower(labels) == "status")
  if (length(status_idx) == 0) status_idx <- nchan
  
  # BDF data start
  data_offset <- header$headLen
  seek(con, where = data_offset, origin = "start")
  
  # 24-bit ints: total samples (all channels × time)
  total_samples <- nsamp_full * nchan
  bytes_per_int <- 3L
  
  chunk_bytes    <- chunk_mb * 1024^2
  ints_per_chunk <- floor(chunk_bytes / bytes_per_int)
  
  # ----- storage for full-rate EEG + status -----
  eeg_list_full   <- list()   # EEG channels only, full rate
  status_full_vec <- numeric(0)
  read_so_far     <- 0L
  
  # ----- Step 1: read full-rate EEG + status in chunks -----
  while (read_so_far < total_samples) {
    to_read <- min(ints_per_chunk, total_samples - read_so_far)
    raw_chunk <- readBin(con, what = "raw", n = to_read * bytes_per_int)
    if (!length(raw_chunk)) break
    
    vals <- unpack_int24(raw_chunk)
    # ensure full columns only
    n_full <- floor(length(vals) / nchan) * nchan
    if (n_full == 0L) break
    vals <- vals[seq_len(n_full)]
    
    mat <- matrix(vals, nrow = nchan, byrow = TRUE)
    
    status_full_vec <- c(status_full_vec, mat[status_idx, ])
    eeg_list_full[[length(eeg_list_full) + 1L]] <- mat[-status_idx, , drop = FALSE]
    
    read_so_far <- read_so_far + (n_full)
  }
  
  eegdata_full <- do.call(cbind, eeg_list_full)  # EEG channels x timepoints
  channels_eeg <- labels[-status_idx]
  
  # full-rate time vector
  times_full <- (0:(ncol(eegdata_full) - 1)) / srate_full
  
  # ----- Step 2: detect events on full-rate status -----
  events_full <- extract_biosemi_events(status_full_vec, srate_full)
  
  # ----- Step 3: if no downsampling requested -----
  if (is.null(target_rate) || target_rate >= srate_full) {
    return(list(
      eeg_data = eegdata_full,
      channels = channels_eeg,
      srate    = srate_full,
      status   = status_full_vec,
      times    = times_full,
      events   = events_full
    ))
  }
  
  # ----- Step 4: if downsampling is requested -----
  
  # integer factor (assumes target_rate divides srate_full)
  factor <- round(srate_full / target_rate)
  if (factor < 1) factor <- 1
  srate_ds <- srate_full / factor
  
  # 4.1 Downsample EEG + times only (keep status_full at full rate)
  idx       <- seq(1, ncol(eegdata_full), by = factor)
  eegdata_ds <- eegdata_full[, idx, drop = FALSE]
  times_ds   <- times_full[idx]
  
  # 4.2 Rescale event onsets to new rate
  events_ds <- events_full
  if (nrow(events_ds) > 0) {
    events_ds$onset <- round(events_full$onset / factor)
    events_ds$onset <- pmax(1, events_ds$onset)
    events_ds$onset <- pmin(events_ds$onset, length(times_ds))
    if ("onset_time" %in% names(events_ds)) {
      events_ds$onset_time <- times_ds[events_ds$onset]
    }
  }
  
  # Return downsampled set
  list(
    eeg_data = eegdata_ds,
    channels = channels_eeg,
    srate    = srate_ds,
    status   = status_full_vec,  # still full-rate, kept for completeness
    times    = times_ds,
    events   = events_ds
  )
}


#' Extract Events from BioSemi Status Channel
#'
#' This internal function parses event markers from the BioSemi status channel.
#' BioSemi systems encode experiment triggers as numeric codes in the status
#' channel. This function identifies when these values change (indicating events)
#' and creates an event table.
#'
#' @param status_signal Numeric vector containing the status channel values
#'   at each time point
#' @param sampling_rate Numeric - sampling rate in Hz
#'   Used to convert sample indices to time (seconds)
#'
#' @return Data frame with event information:
#' \describe{
#'   \item{onset}{Integer - event onset in samples (1-indexed)}
#'   \item{onset_time}{Numeric - event onset in seconds}
#'   \item{type}{Character - trigger code as string (e.g., "255")}
#'   \item{description}{Character - human-readable description}
#' }
#'
#' @keywords internal
extract_biosemi_events <- function(status_signal, sampling_rate) {
  
  # ========== HANDLE EMPTY STATUS ==========
  if (is.null(status_signal) || length(status_signal) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # ========== DETECT EVENT TRANSITIONS ==========
  status_numeric <- as.numeric(status_signal)
  
  # Find where status changes and is non-zero
  status_changes <- c(TRUE, diff(status_numeric) != 0)
  event_idx_unique <- which(status_changes & status_numeric != 0)
  
  # ========== HANDLE CASE WITH NO EVENTS ==========
  if (length(event_idx_unique) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # ========== CREATE EVENTS DATA FRAME ==========
  events_df <- data.frame(
    onset = event_idx_unique,
    onset_time = (event_idx_unique - 1) / sampling_rate,
    type = as.character(status_numeric[event_idx_unique]),
    description = paste0(
      "Trigger code: ",
      as.character(status_numeric[event_idx_unique])
    ),
    stringsAsFactors = FALSE
  )
  
  events_df
}

#' Generate Quality Summary of BioSemi Import
#'
#' After importing a BioSemi file, this function provides a comprehensive
#' quality report including data statistics, channel information, event counts,
#' and any issues detected during import.
#'
#' @param eeg An eeg object imported via read_biosemi_interactive()
#'
#' @return Invisibly returns a list with quality metrics
#'
#' @export
summarize_biosemi_import_interactive <- function(eeg) {
  
  # ========== HEADER ==========
  cat("\n")
  cat(strrep("=", 75), "\n")
  cat("BioSemi Import Quality Summary\n")
  cat(strrep("=", 75), "\n\n")
  
  # ========== DATA STATISTICS ==========
  cat("DATA STATISTICS:\n")
  cat(" Channels: ", length(eeg$channels), "\n", sep = "")
  cat(" Time points: ", ncol(eeg$data), "\n", sep = "")
  cat(" Duration: ", round(max(eeg$times), 3), " seconds\n", sep = "")
  cat(" Sampling rate: ", eeg$sampling_rate, " Hz\n", sep = "")
  
  # ========== CHANNEL INFORMATION ==========
  cat("\nCHANNEL INFORMATION:\n")
  channels_display <- paste(
    eeg$channels[1:min(8, length(eeg$channels))],
    collapse = ", "
  )
  if (length(eeg$channels) > 8) {
    channels_display <- paste0(
      channels_display, " ... (+",
      length(eeg$channels) - 8, " more)"
    )
  }
  cat(" Channels: ", channels_display, "\n", sep = "")
  
  # ========== DATA QUALITY ANALYSIS ==========
  cat("\nDATA QUALITY METRICS:\n")
  channel_vars <- apply(eeg$data, 1, var, na.rm = TRUE)
  channel_means <- apply(eeg$data, 1, mean, na.rm = TRUE)
  cat(
    " Amplitude range: [",
    round(min(eeg$data), 2), ", ",
    round(max(eeg$data), 2), "] µV\n", sep = ""
  )
  cat(
    " Mean across channels: ",
    round(mean(channel_means), 2), " µV\n", sep = ""
  )
  
  # Identify problematic channels
  flat_channels <- which(channel_vars < 0.01)
  if (length(flat_channels) > 0) {
    cat(
      " WARNING: ", length(flat_channels),
      " potentially flat channels: ",
      paste(eeg$channels[flat_channels], collapse = ", "),
      "\n", sep = ""
    )
  }
  
  # ========== EVENT INFORMATION ==========
  cat("\nEVENT INFORMATION:\n")
  cat(" Total events: ", nrow(eeg$events), "\n", sep = "")
  if (nrow(eeg$events) > 0) {
    event_types <- unique(eeg$events$type)
    cat(" Event types: ", paste(event_types, collapse = ", "), "\n", sep = "")
    cat(" Event counts:\n")
    for (et in event_types) {
      count <- sum(eeg$events$type == et)
      cat(" Trigger ", et, ": ", count, " occurrences\n", sep = "")
    }
  } else {
    cat(" WARNING: No events detected in status channel!\n")
  }
  
  # ========== REFERENCE INFORMATION ==========
  cat("\nREFERENCE INFORMATION:\n")
  cat(" Reference: ", eeg$reference, "\n", sep = "")
  
  # ========== FOOTER ==========
  cat(strrep("=", 75), "\n\n")
  
  # ========== RETURN METRICS INVISIBLY ==========
  quality_metrics <- list(
    n_channels = length(eeg$channels),
    n_samples = ncol(eeg$data),
    duration_sec = max(eeg$times),
    sampling_rate = eeg$sampling_rate,
    channel_vars = channel_vars,
    n_events = nrow(eeg$events),
    event_types = if (nrow(eeg$events) > 0) unique(eeg$events$type) else character(0),
    flat_channels = if (length(flat_channels) > 0) eeg$channels[flat_channels] else character(0)
  )
  
  invisible(quality_metrics)
}

# End of read_biosemi_interactive.R
