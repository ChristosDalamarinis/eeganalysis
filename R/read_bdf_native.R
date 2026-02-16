#' ============================================================================
#'                          Native BDF File Reader
#' ============================================================================
#'
#' This module provides a native R implementation for reading BioSemi BDF
#' files without external dependencies like edfReader.
#'
#' Depends on: eeg_class.R (for new_eeg() constructor)
#'
#' Author: Christos Dalamarinis
#' Date: Jan 2026
#' ============================================================================
#'
#' Read BDF File Header
#'
#' @param file_path Character string with path to .bdf file
#' @return List containing header information
#' @keywords internal
read_bdf_header_native <- function(file_path) {
  
  con <- file(file_path, "rb")
  on.exit(close(con), add = TRUE)
  
  # Read fixed header (256 bytes)
  version <- readBin(con, "raw", n = 8)
  patient_id <- rawToChar(readBin(con, "raw", n = 80))
  recording_id <- rawToChar(readBin(con, "raw", n = 80))
  start_date <- rawToChar(readBin(con, "raw", n = 8))
  start_time <- rawToChar(readBin(con, "raw", n = 8))
  header_bytes <- as.integer(rawToChar(readBin(con, "raw", n = 8)))
  reserved <- rawToChar(readBin(con, "raw", n = 44))
  n_data_records <- as.integer(rawToChar(readBin(con, "raw", n = 8)))
  record_duration <- as.numeric(rawToChar(readBin(con, "raw", n = 8)))
  n_channels <- as.integer(rawToChar(readBin(con, "raw", n = 4)))
  
  # Read variable header (per-channel information)
  labels <- sapply(1:n_channels, function(i) {
    trimws(rawToChar(readBin(con, "raw", n = 16)))
  })
  
  transducer <- sapply(1:n_channels, function(i) {
    trimws(rawToChar(readBin(con, "raw", n = 80)))
  })
  
  physical_dim <- sapply(1:n_channels, function(i) {
    trimws(rawToChar(readBin(con, "raw", n = 8)))
  })
  
  physical_min <- sapply(1:n_channels, function(i) {
    as.numeric(trimws(rawToChar(readBin(con, "raw", n = 8))))
  })
  
  physical_max <- sapply(1:n_channels, function(i) {
    as.numeric(trimws(rawToChar(readBin(con, "raw", n = 8))))
  })
  
  digital_min <- sapply(1:n_channels, function(i) {
    as.numeric(trimws(rawToChar(readBin(con, "raw", n = 8))))
  })
  
  digital_max <- sapply(1:n_channels, function(i) {
    as.numeric(trimws(rawToChar(readBin(con, "raw", n = 8))))
  })
  
  prefilter <- sapply(1:n_channels, function(i) {
    trimws(rawToChar(readBin(con, "raw", n = 80)))
  })
  
  n_samples_per_record <- sapply(1:n_channels, function(i) {
    as.integer(trimws(rawToChar(readBin(con, "raw", n = 8))))
  })
  
  reserved_per_channel <- sapply(1:n_channels, function(i) {
    readBin(con, "raw", n = 32)
  })
  
  # Calculate sampling rates
  s_rate <- n_samples_per_record / record_duration
  
  # Calculate total samples per channel
  total_samples <- n_samples_per_record * n_data_records
  
  list(
    nChannels = n_channels,
    labels = labels,
    sRate = s_rate,
    nSamplesPerRecord = n_samples_per_record,
    nDataRecords = n_data_records,
    recordDuration = record_duration,
    digitalMin = digital_min,
    digitalMax = digital_max,
    physicalMin = physical_min,
    physicalMax = physical_max,
    prefilter = prefilter,
    transducer = transducer,
    physicalDim = physical_dim,
    headerBytes = header_bytes,
    totalSamples = total_samples,
    startDate = trimws(start_date),
    startTime = trimws(start_time),
    patientID = trimws(patient_id),
    recordingID = trimws(recording_id)
  )
}

#' Unpack 24-bit Signed Integers
#'
#' @param raw_vec Raw vector of bytes (length must be multiple of 3)
#' @return Integer vector of unpacked values
#' @keywords internal
unpack_int24_native <- function(raw_vec) {
  
  n_ints <- length(raw_vec) %/% 3L
  if (n_ints == 0L) return(integer(0))
  
  # Ensure we have complete triplets
  raw_vec <- raw_vec[1:(n_ints * 3L)]
  
  # Reshape into 3 rows (one per byte)
  m <- matrix(as.integer(raw_vec), nrow = 3L)
  
  # Combine bytes: byte1 + byte2*256 + byte3*65536
  vals <- m[1, ] + bitwShiftL(m[2, ], 8L) + bitwShiftL(m[3, ], 16L)
  
  # Convert from unsigned 24-bit to signed (two's complement)
  neg_mask <- vals >= 8388608L  # 2^23
  vals[neg_mask] <- vals[neg_mask] - 16777216L  # 2^24
  
  vals
}

#' Convert Digital Values to Physical Units
#'
#' @param digital Integer vector of digital values
#' @param digital_min Digital minimum from header
#' @param digital_max Digital maximum from header
#' @param physical_min Physical minimum from header (microV)
#' @param physical_max Physical maximum from header (microV)
#' @return Numeric vector of physical values
#' @keywords internal
digital_to_physical <- function(digital, digital_min, digital_max, 
                                physical_min, physical_max) {
  
  # Force symmetric range 16/02/2026
  digital_max_adjusted <- -digital_min
  
  gain <- (physical_max - physical_min) / (digital_max_adjusted - digital_min)
  offset <- physical_min - gain * digital_min
  
  gain * digital + offset
}

#' Extract Events from STATUS Channel
#'
#' @param status_signal Numeric vector of STATUS channel values
#' @param sampling_rate Sampling rate in Hz
#' @return Data frame with columns: onset, onset_time, type, description
#' @keywords internal
extract_events_native <- function(status_signal, sampling_rate) {
  
  if (is.null(status_signal) || length(status_signal) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = integer(0),
      description = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # CRITICAL: Mask to lower 16 bits only
  status_masked <- bitwAnd(as.integer(status_signal), 65535L)
  
  # Detect changes in trigger values
  status_diff <- c(0, diff(status_masked))
  event_indices <- which(status_diff != 0 & status_masked != 0)
  
  # Adjust indices
  # event_indices <- event_indices + 1L
  
  if (length(event_indices) == 0) {
    return(data.frame(
      onset = integer(0),
      onset_time = numeric(0),
      type = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create events data frame
  trigger_values <- status_masked[event_indices]
  
  data.frame(
    onset = event_indices,
    onset_time = (event_indices - 1) / sampling_rate,
    type = as.integer(trigger_values),
    description = paste0("Trigger: ", trigger_values),
    stringsAsFactors = FALSE
  )
}

#' Read BDF File with Native R Implementation (Memory-Efficient Chunked Version)
#'
#' @param file_path Character string with path to .bdf file
#' @param verbose Logical, print progress messages
#' @param chunk_records Number of data records to read per chunk (default: 100)
#' @return An object of class 'eeg' created with new_eeg()
#' @export
read_bdf_native <- function(file_path, verbose = TRUE, chunk_records = 100) {
  
  # ========== VALIDATION ==========
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path, call. = FALSE)
  }
  
  if (!grepl("\\.bdf$", file_path, ignore.case = TRUE)) {
    warning("File does not have .bdf extension. Attempting to read stopped.")
    return(invisible(NULL))
  }
  
  # ========== READ HEADER ==========
  if (verbose) cat("[1/5] Reading BDF header...\n")
  
  header <- tryCatch({
    read_bdf_header_native(file_path)
  }, error = function(e) {
    stop("Failed to read BDF header: ", e$message, call. = FALSE)
  })
  
  if (verbose) {
    cat("  Channels: ", header$nChannels, "\n", sep = "")
    cat("  Sampling rate: ", header$sRate[1], " Hz\n", sep = "")
    cat("  Duration: ~", round(header$totalSamples[1] / header$sRate[1], 1), 
        " seconds\n", sep = "")
    cat("  Data records: ", header$nDataRecords, "\n", sep = "")
  }
  
  # ========== IDENTIFY STATUS CHANNEL ==========
  status_idx <- which(tolower(header$labels) == "status")
  if (length(status_idx) == 0) {
    status_idx <- header$nChannels
    if (verbose) cat("  WARNING: No 'Status' channel found, using last channel\n")
  }
  
  eeg_idx <- setdiff(1:header$nChannels, status_idx)
  
  # ========== PREPARE OUTPUT STORAGE ==========
  n_samples_final <- header$totalSamples[1]
  n_eeg_channels <- length(eeg_idx)
  
  if (verbose) cat("[2/5] Allocating memory for ", n_eeg_channels, " EEG channels...\n", sep = "")
  
  # Pre-allocate output matrices
  eeg_data_physical <- matrix(0, nrow = n_eeg_channels, ncol = n_samples_final)
  status_data <- numeric(n_samples_final)
  
  # ========== READ DATA IN CHUNKS ==========
  if (verbose) cat("[3/5] Reading data in chunks...\n")
  
  con <- file(file_path, "rb")
  on.exit(close(con), add = TRUE)
  
  seek(con, where = header$headerBytes, origin = "start")
  
  # Calculate chunk parameters
  samples_per_record <- sum(header$nSamplesPerRecord)
  bytes_per_record <- samples_per_record * 3L
  n_chunks <- ceiling(header$nDataRecords / chunk_records)
  
  output_col_idx <- 1L
  
  for (chunk_idx in 1:n_chunks) {
    
    # Calculate records to read in this chunk
    records_remaining <- header$nDataRecords - (chunk_idx - 1) * chunk_records
    records_this_chunk <- min(chunk_records, records_remaining)
    
    if (verbose && n_chunks > 1) {
      cat("  Chunk ", chunk_idx, "/", n_chunks, " (", records_this_chunk, " records)...\r", sep = "")
    }
    
    # Read chunk of raw data
    chunk_bytes <- bytes_per_record * records_this_chunk
    raw_chunk <- readBin(con, "raw", n = chunk_bytes)
    
    if (length(raw_chunk) == 0) break
    
    # Unpack 24-bit integers
    chunk_values <- unpack_int24_native(raw_chunk)
    
    # Process each record in this chunk
    value_idx <- 1L
    
    for (rec in 1:records_this_chunk) {
      
      # Process each channel in this record
      for (ch in 1:header$nChannels) {
        n_samp <- header$nSamplesPerRecord[ch]
        
        if (value_idx + n_samp - 1L > length(chunk_values)) break
        
        ch_data <- chunk_values[value_idx:(value_idx + n_samp - 1L)]
        value_idx <- value_idx + n_samp
        
        # Store in appropriate location
        col_range <- output_col_idx:(output_col_idx + n_samp - 1L)
        
        if (ch %in% eeg_idx) {
          eeg_row <- which(eeg_idx == ch)
          
          # Convert to physical units immediately to save memory
          eeg_data_physical[eeg_row, col_range] <- digital_to_physical(
            ch_data,
            header$digitalMin[ch],
            header$digitalMax[ch],
            header$physicalMin[ch],
            header$physicalMax[ch]
          )
        } else {
          # STATUS channel - keep as digital values
          status_data[col_range] <- ch_data
        }
      }
      
      output_col_idx <- output_col_idx + header$nSamplesPerRecord[1]
    }
  }
  
  if (verbose) cat("\n")
  
  # ========== EXTRACT EVENTS ==========
  if (verbose) cat("[4/5] Extracting events from STATUS channel...\n")
  
  sampling_rate <- header$sRate[1]
  times <- (0:(n_samples_final - 1)) / sampling_rate
  events_df <- extract_events_native(status_data, sampling_rate)
  
  if (verbose) {
    cat("  Events found: ", nrow(events_df), "\n", sep = "")
  }
  
  # ========== PREPARE METADATA ==========
  metadata <- list(
    file_path = file_path,
    file_name = basename(file_path),
    date_recorded = header$startDate,
    time_recorded = header$startTime,
    device = "BioSemi",
    n_channels_recorded = header$nChannels,
    n_eeg_channels = n_eeg_channels,
    status_channel_name = header$labels[status_idx],
    reference_scheme = "Biosemi CMS/DRL (Common Mode Sense / Driven Right Leg)",
    duration_seconds = max(times),
    n_samples = n_samples_final,
    n_events = nrow(events_df),
    patient_id = header$patientID,
    recording_id = header$recordingID
  )
  
  # ========== CREATE EEG OBJECT ==========
  if (verbose) cat("[5/5] Creating EEG object...\n")
  
  eeg_obj <- new_eeg(
    data = eeg_data_physical,
    channels = header$labels[eeg_idx],
    sampling_rate = sampling_rate,
    times = times,
    events = events_df,
    metadata = metadata,
    reference = "Biosemi CMS/DRL",
    preprocessing_history = list(
      paste0("Imported from BDF file using native chunked reader: ", basename(file_path))
    )
  )
  
  if (verbose) cat("\n[OK] Successfully imported BDF file!\n\n")
  
  return(eeg_obj)
}

