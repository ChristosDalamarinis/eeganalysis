#' ============================================================================
#'                          EEG Downsampling Functions
#' ============================================================================
#'
#' This module provides functions for downsampling EEG data to reduce sampling
#' rate while preserving signal integrity. Downsampling is essential for
#' reducing computational overhead, file sizes, and analysis time without
#' losing meaningful neural information.
#'
#' The module implements anti-aliasing filters to prevent spectral aliasing
#' artifacts when reducing sampling rate. Most EEG analyses work well with
#' 250-512 Hz sampling, even when data is recorded at higher rates like
#' 2048 Hz.
#'
#' Author: Christos Dalamarinis
#' Date: 2026
#' ============================================================================
#'
#' Downsample EEG Data
#'
#' Reduces the sampling rate of EEG data by applying an anti-aliasing low-pass
#' filter followed by decimation. This function preserves signal integrity while
#' reducing computational overhead and file size.
#'
#' @param eeg_obj An object of class 'eeg' containing EEG data
#' 
#' @param target_rate Numeric value specifying the desired sampling rate in Hz.
#'                    Must be lower than current sampling rate and satisfy
#'                    Nyquist theorem (at least 2 times highest frequency of interest).
#'                    Common values: 128, 256, 512 Hz
#' 
#' @param method Character string specifying downsampling method:
#'               \describe{
#'                 \item{"decimate"}{(default) Applies anti-aliasing filter then decimates.
#'                                   Recommended for all research applications.}
#'                 \item{"simple"}{Direct downsampling without filtering.
#'                                 Not recommended - may cause spectral aliasing.}
#'               }
#' 
#' @param filter_order Integer specifying the order of the anti-aliasing Butterworth filter.
#'                     Default: 8. Higher values create sharper cutoff (more computation).
#'                     Typical range: 4-12.
#' 
#' @param cutoff_ratio Numeric value between 0 and 1 specifying the filter cutoff
#'                     as a proportion of the new Nyquist frequency.
#'                     Default: 0.9 (90% of Nyquist).
#'                     Formula: cutoff_freq = (target_rate / 2) times cutoff_ratio
#'                     Lower values (0.8) = more conservative, higher values (0.95) = more aggressive.
#' 
#' @param preserve_bands Numeric vector of length 2: c(min_freq, max_freq) in Hz.
#'                       If provided, validates that target sampling rate can preserve
#'                       these frequencies. Throws error if impossible.
#'                       Example: c(0.5, 50) for standard EEG bands.
#'                       Default: NULL (no validation).
#' 
#' @param event_strategy Character string specifying how event timing is adjusted:
#'                       \describe{
#'                         \item{"round"}{(default) Round to nearest sample}
#'                         \item{"nearest"}{Find temporally nearest sample}
#'                         \item{"floor"}{Always round down (event at or before actual time)}
#'                         \item{"ceiling"}{Always round up (event at or after actual time)}
#'                       }
#' 
#' @param remove_outbound_events Logical. If TRUE (default), removes events that fall
#'                               outside the downsampled time range with a warning.
#'                               If FALSE, clamps event indices to valid range.
#' 
#' @param verbose Logical indicating whether to print progress messages and summary.
#'                Default: TRUE. Set to FALSE for batch processing.
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input parameters and data quality
#' 2. Calculates downsampling factor (must result in integer ratio)
#' 3. Applies anti-aliasing low-pass Butterworth filter at cutoff frequency
#' 4. Decimates the signal by selecting every nth sample
#' 5. Updates time vector and event onset times accordingly
#' 6. Validates and adjusts events for new sampling rate
#'
#' \strong{Anti-aliasing filtering} prevents spectral aliasing by removing frequencies
#' above the new Nyquist frequency before decimation. The Butterworth filter is
#' applied with zero-phase filtering (forward-backward) using filtfilt().
#'
#' \strong{Frequency Preservation:}
#' The cutoff frequency determines what frequencies are preserved:
#' \itemize{
#'   \item{Target 512 Hz -> cutoff ~230 Hz -> suitable for high gamma}
#'   \item{Target 256 Hz -> cutoff ~115 Hz -> suitable for RIFT/SSVEP studies}
#'   \item{Target 128 Hz -> cutoff ~58 Hz -> suitable only for standard EEG bands}
#' }
#'
#' \strong{DO NOT downsample if you need to preserve:}
#' \itemize{
#'   \item{RIFT/SSVEP signals (60-120 Hz) -> keep at least 256 Hz sampling}
#'   \item{High gamma (60-100 Hz) -> keep at least 200 Hz sampling}
#'   \item{Ripples (80-250 Hz) -> keep  at least 500 Hz sampling}
#' }
#'
#' @return An 'eeg' object with downsampled data and updated parameters:
#'   \describe{
#'     \item{data}{Downsampled signal matrix (channels times new_timepoints)}
#'     \item{sampling_rate}{Updated to actual achieved sampling rate}
#'     \item{times}{Adjusted time vector for new sampling rate}
#'     \item{events}{Event data frame with recalculated onset times}
#'     \item{channels}{Unchanged channel names}
#'     \item{metadata}{Unchanged metadata}
#'     \item{reference}{Unchanged reference scheme}
#'     \item{preprocessing_history}{Updated with downsampling record}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Standard ERP analysis (preserve 0.5-50 Hz)
#'   eeg_down <- downsample_eeg(eeg, target_rate = 256,
#'                              preserve_bands = c(0.5, 50))
#'
#'   # RIFT study (preserve 60-100 Hz)
#'   eeg_down <- downsample_eeg(eeg, target_rate = 512,
#'                              preserve_bands = c(60, 100),
#'                              filter_order = 10)
#'   
#'   # Conservative downsampling
#'   eeg_down <- downsample_eeg(eeg, target_rate = 256,
#'                              cutoff_ratio = 0.8,
#'                              event_strategy = "floor")
#'   
#'   # Batch processing (silent)
#'   eeg_down <- downsample_eeg(eeg, target_rate = 256, verbose = FALSE)
#' }
#'
#' @seealso \code{\link{new_eeg}} for EEG object structure
#'
#' @importFrom signal butter filtfilt
#' @export
downsample <- function(eeg_obj, 
                           target_rate, 
                           method = c("decimate", "simple"),
                           filter_order = 8,
                           cutoff_ratio = 0.9,
                           preserve_bands = NULL,
                           event_strategy = c("round", "nearest", "floor", "ceiling"),
                           remove_outbound_events = TRUE,
                           verbose = TRUE) {
  
  # ========== MATCH ARGUMENTS ==========
  
  method <- match.arg(method)
  event_strategy <- match.arg(event_strategy)
  
  # ========== INPUT VALIDATION ==========
  
  # Validate eeg_obj
  if (!inherits(eeg_obj, "eeg")) {
    stop("Input must be an object of class 'eeg'", call. = FALSE)
  }
  
  # Validate target_rate
  if (!is.numeric(target_rate) || length(target_rate) != 1 || target_rate <= 0) {
    stop("target_rate must be a single positive numeric value", call. = FALSE)
  }
  
  if (target_rate >= eeg_obj$sampling_rate) {
    stop("target_rate (", target_rate, " Hz) must be less than current sampling rate (", 
         eeg_obj$sampling_rate, " Hz)", call. = FALSE)
  }
  
  # Validate filter_order
  if (!is.numeric(filter_order) || length(filter_order) != 1 || 
      filter_order <= 0 || filter_order %% 1 != 0) {
    stop("filter_order must be a positive integer", call. = FALSE)
  }
  
  # Validate cutoff_ratio
  if (!is.numeric(cutoff_ratio) || length(cutoff_ratio) != 1 || 
      cutoff_ratio <= 0 || cutoff_ratio > 1) {
    stop("cutoff_ratio must be a numeric value between 0 and 1", call. = FALSE)
  }
  
  # Validate preserve_bands
  if (!is.null(preserve_bands)) {
    if (!is.numeric(preserve_bands) || length(preserve_bands) != 2) {
      stop("preserve_bands must be a numeric vector of length 2: c(min_freq, max_freq)", 
           call. = FALSE)
    }
    if (preserve_bands[1] >= preserve_bands[2]) {
      stop("preserve_bands[1] must be less than preserve_bands[2]", call. = FALSE)
    }
    if (any(preserve_bands < 0)) {
      stop("preserve_bands values must be positive", call. = FALSE)
    }
  }
  
  # Validate remove_outbound_events
  if (!is.logical(remove_outbound_events) || length(remove_outbound_events) != 1) {
    stop("remove_outbound_events must be a single logical value (TRUE or FALSE)", 
         call. = FALSE)
  }
  
  # Validate verbose
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }
  
  # Check for NA/Inf values in data
  if (any(is.na(eeg_obj$data))) {
    stop("EEG data contains NA values. Please clean data before downsampling.", 
         call. = FALSE)
  }
  
  if (any(!is.finite(eeg_obj$data))) {
    stop("EEG data contains non-finite values (Inf or -Inf). Please clean data before downsampling.", 
         call. = FALSE)
  }
  
  # ========== CALCULATE DOWNSAMPLING FACTOR ==========
  
  downsample_factor <- eeg_obj$sampling_rate / target_rate
  
  # Check if factor is too small
  if (downsample_factor < 1.5) {
    stop("Downsampling factor too small (", round(downsample_factor, 2), 
         "). Consider a lower target_rate for meaningful downsampling.", call. = FALSE)
  }
  
  # Check if factor is close to an integer
  if (abs(downsample_factor - round(downsample_factor)) > 0.01) {
    warning("Downsampling factor (", round(downsample_factor, 2), 
            ") is not an integer. Rounding to nearest integer for optimal results.", 
            call. = FALSE, immediate. = TRUE)
  }
  
  downsample_factor <- round(downsample_factor)
  actual_target_rate <- eeg_obj$sampling_rate / downsample_factor
  
  if (abs(actual_target_rate - target_rate) > 1 && verbose) {
    message("NOTE: Actual achieved sampling rate will be ", actual_target_rate, 
            " Hz (factor: ", downsample_factor, ")")
  }
  
  # ========== CALCULATE FILTER PARAMETERS ==========
  
  new_nyquist <- actual_target_rate / 2
  filter_cutoff <- new_nyquist * cutoff_ratio
  
  # ========== VALIDATE FREQUENCY PRESERVATION ==========
  
  if (!is.null(preserve_bands)) {
    max_freq_needed <- max(preserve_bands)
    min_freq_needed <- min(preserve_bands)
    
    if (max_freq_needed > filter_cutoff) {
      stop("Cannot preserve frequencies up to ", max_freq_needed, " Hz with target rate of ", 
           actual_target_rate, " Hz (cutoff: ", round(filter_cutoff, 1), " Hz).\n",
           "Minimum required sampling rate: ", ceiling(max_freq_needed * 2.2), " Hz\n",
           "Recommendation: Use target_rate >= ", ceiling(max_freq_needed * 2.2), " Hz",
           call. = FALSE)
    }
    
    if (verbose) {
      message("[OK] Frequency preservation validated: ", min_freq_needed, "-", 
              max_freq_needed, " Hz can be preserved at ", actual_target_rate, " Hz")
    }
  }
  
  # Warn about high-frequency loss
  if (filter_cutoff < 60 && verbose) {
    message("\n WARNING: High-frequency content will be removed!")
    message("  Filter cutoff: ", round(filter_cutoff, 1), " Hz")
    message("  Frequencies above ~", round(filter_cutoff, 1), " Hz will be lost.")
    message("  This is NOT suitable for:")
    message("    - RIFT/SSVEP studies (typically 60-120 Hz)")
    message("    - High gamma analysis (>60 Hz)")
    message("    - Ripple detection (80-250 Hz)")
    message("  Consider keeping original sampling rate for high-frequency research.\n")
  }
  
  # ========== EXTRACT DATA DIMENSIONS ==========
  
  n_channels <- nrow(eeg_obj$data)
  n_timepoints <- ncol(eeg_obj$data)
  
  # Check for very short recordings
  if (n_timepoints < downsample_factor * 10) {
    warning("Very short recording (", n_timepoints, " samples). ",
            "Downsampling may produce unreliable results.", 
            call. = FALSE, immediate. = TRUE)
  }
  
  # ========== DOWNSAMPLING ==========
  
  if (method == "decimate") {
    # Apply anti-aliasing filter before decimation
    
    if (verbose) {
      cat("Applying anti-aliasing filter (order ", filter_order, ", cutoff ", 
          round(filter_cutoff, 1), " Hz)...\n", sep = "")
    }
    
    # Design low-pass Butterworth filter
    nyquist_freq <- eeg_obj$sampling_rate / 2
    normalized_cutoff <- filter_cutoff / nyquist_freq
    
    # Ensure normalized cutoff is valid (< 1)
    if (normalized_cutoff >= 1) {
      stop("Filter cutoff frequency is too high for the current sampling rate. ",
           "This should not happen - please report this as a bug.", call. = FALSE)
    }
    
    # Ensure signal package is available
    if (!requireNamespace("signal", quietly = TRUE)) {
      stop("'signal' package is required for decimation method.\n",
           "Install with: install.packages('signal')", call. = FALSE)
    }
    
    # Create filter
    butter_filter <- signal::butter(n = filter_order, 
                                    W = normalized_cutoff, 
                                    type = "low")
    
    # Apply filter to each channel
    filtered_data <- matrix(0, nrow = n_channels, ncol = n_timepoints)
    for (ch in 1:n_channels) {
      filtered_data[ch, ] <- signal::filtfilt(butter_filter, eeg_obj$data[ch, ])
    }
    
    if (verbose) {
      cat("Decimating signal (factor: ", downsample_factor, ")...\n", sep = "")
    }
    
    # Decimate: select every nth sample
    decimate_indices <- seq(1, n_timepoints, by = downsample_factor)
    downsampled_data <- filtered_data[, decimate_indices, drop = FALSE]
    
  } else if (method == "simple") {
    # Direct downsampling without filtering (not recommended)
    warning("Simple downsampling without anti-aliasing filter may cause spectral aliasing. ",
            "Consider using method = 'decimate' for research applications.", 
            call. = FALSE, immediate. = TRUE)
    
    if (verbose) {
      cat("Performing simple downsampling (no anti-aliasing filter)...\n")
    }
    
    decimate_indices <- seq(1, n_timepoints, by = downsample_factor)
    downsampled_data <- eeg_obj$data[, decimate_indices, drop = FALSE]
  }
  
  # ========== UPDATE TIME VECTOR ==========
  
  new_times <- eeg_obj$times[decimate_indices]
  
  # ========== UPDATE EVENTS ==========
  
  # Adjust event onset indices for new sampling rate
  if (nrow(eeg_obj$events) > 0 && ncol(eeg_obj$events) > 0) {
    new_events <- eeg_obj$events
    
    # Save original onset for bounds checking
    original_onset <- new_events$onset
    
    # Apply event strategy to recalculate onset sample indices
    if (event_strategy == "round") {
      new_events$onset <- round(original_onset / downsample_factor)
    } else if (event_strategy == "floor") {
      new_events$onset <- floor(original_onset / downsample_factor)
    } else if (event_strategy == "ceiling") {
      new_events$onset <- ceiling(original_onset / downsample_factor)
    } else if (event_strategy == "nearest") {
      # Find nearest sample in downsampled data
      new_events$onset <- round(original_onset / downsample_factor)
    }
    
    # Check for events outside the downsampled time range
    out_of_bounds <- new_events$onset < 1 | new_events$onset > length(new_times)
    
    if (any(out_of_bounds)) {
      n_removed <- sum(out_of_bounds)
      
      if (remove_outbound_events) {
        warning("Removed ", n_removed, " event(s) that fell outside the downsampled time range", 
                call. = FALSE, immediate. = TRUE)
        new_events <- new_events[!out_of_bounds, , drop = FALSE]
      } else {
        warning(n_removed, " event(s) fell outside the downsampled time range. ",
                "Clamping to valid range (may cause event clustering at boundaries).",
                call. = FALSE, immediate. = TRUE)
        # Clamp to valid range
        new_events$onset <- pmax(1, pmin(new_events$onset, length(new_times)))
      }
    }
    
    # Update onset_time based on new indices (if any events remain)
    if (nrow(new_events) > 0) {
      # Ensure onset indices are within valid range
      new_events$onset <- pmax(1, pmin(new_events$onset, length(new_times)))
      
      # Update onset_time based on new time vector
      if ("onset_time" %in% names(new_events)) {
        new_events$onset_time <- new_times[new_events$onset]
      }
    }
  } else {
    new_events <- eeg_obj$events
  }
  
  # ========== UPDATE PREPROCESSING HISTORY ==========
  
  preprocessing_step <- paste0(
    "Downsampled from ", eeg_obj$sampling_rate, " Hz to ", 
    actual_target_rate, " Hz using '", method, "' method (factor: ", 
    downsample_factor, ", cutoff: ", round(filter_cutoff, 1), " Hz)"
  )
  
  new_history <- c(eeg_obj$preprocessing_history, list(preprocessing_step))
  
  # ========== CREATE NEW EEG OBJECT ==========
  
  downsampled_eeg <- new_eeg(
    data = downsampled_data,
    channels = eeg_obj$channels,
    sampling_rate = actual_target_rate,
    times = new_times,
    events = new_events,
    metadata = eeg_obj$metadata,
    reference = eeg_obj$reference,
    preprocessing_history = new_history
  )
  
  # ========== SUMMARY MESSAGE ==========
  
  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("Downsampling Complete\n")
    cat(strrep("=", 70), "\n")
    cat("  Original sampling rate:  ", eeg_obj$sampling_rate, " Hz\n")
    cat("  Original samples:        ", n_timepoints, "\n")
    cat("  New sampling rate:       ", actual_target_rate, " Hz\n")
    cat("  New samples:             ", length(new_times), "\n")
    cat("  Reduction:               ", round((1 - length(new_times)/n_timepoints) * 100, 1), "%\n")
    cat("  Downsampling factor:     ", downsample_factor, "\n")
    cat("  Method:                  ", method, "\n")
    if (method == "decimate") {
      cat("  Filter order:            ", filter_order, "\n")
      cat("  Filter cutoff:           ", round(filter_cutoff, 1), " Hz\n")
    }
    cat("  Event strategy:          ", event_strategy, "\n")
    cat("  Events retained:         ", nrow(new_events), " / ", nrow(eeg_obj$events), "\n")
    cat(strrep("=", 70), "\n\n")
  }
  
  return(downsampled_eeg)
}

# End of downsample.R
