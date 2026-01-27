#' ============================================================================
#'            EEG Epoching Functions - Enhanced Trigger Detection
#' ============================================================================
#'
#' This module provides functions for epoching (segmenting) continuous EEG data
#' into time-locked segments around experimental events. This enhanced version
#' implements BioSemi-aware trigger detection that exploits the 24-bit file
#' structure to automatically separate experimental triggers (bits 0-15) from
#' system status codes (bits 16-23), ensuring robust epoching regardless of
#' hardware noise in the status channel.
#'
#' The module extracts epochs based on trigger codes/event markers, applies
#' baseline correction, and optionally rejects epochs with artifacts or missing
#' data. It supports flexible epoch windows, multiple event types, and various
#' baseline correction methods.
#'
#' Key Enhancement: Bitwise masking for BioSemi trigger isolation using the
#' 0xFFFF mask to filter experimental triggers from system triggers, following
#' industry standards. This approach eliminates the need for dataset-specific 
#' tuning and works universally across all BioSemi recordings.
#'
#' Author: Christos Dalamarinis
#' Date: Jan - 2026 NOT TESTED
#' ============================================================================
#'
#' ======================= INSPECT TRIGGERS FUNCTION ==========================
#' Inspect Event Triggers in EEG Data
#'
#' Provides a concise diagnostic summary of event triggers in EEG recordings.
#' Automatically distinguishes experimental triggers (bits 0-15) from system
#' status triggers (bits 16-23) using BioSemi's 24-bit file structure.
#'
#' @param eeg_obj An object of class 'eeg' containing EEG data with events
#' @param mode Character string: "auto" (default), "raw", or "manual"
#' @param trigger_threshold Numeric. Minimum trigger duration in seconds (default: 0.002)
#' @param min_iei Numeric. Minimum inter-event interval in seconds (default: 0.01)
#' @param exclude_types Character vector of event types to exclude (default: NULL)
#' @param plot Logical. Generate diagnostic plots (default: TRUE)
#' @param export_csv Character. Path to export cleaned events (default: NULL)
#'
#' @return List with diagnostic information (invisible)
#'
#' @export
inspect_triggers <- function(eeg_obj,
                             mode = c("auto", "raw", "manual"),
                             trigger_threshold = 0.002,
                             min_iei = 0.01,
                             exclude_types = NULL,
                             plot = TRUE,
                             export_csv = NULL) {
  
  mode <- match.arg(mode)
  
  # ========== VALIDATION ==========
  if (!inherits(eeg_obj, "eeg")) {
    stop("Input must be an object of class 'eeg'", call. = FALSE)
  }
  
  if (is.null(eeg_obj$events) || nrow(eeg_obj$events) == 0) {
    cat("\n⚠ No events found in EEG recording!\n")
    cat("Cannot epoch data without triggers.\n\n")
    return(invisible(list(n_events_raw = 0, n_events_cleaned = 0)))
  }
  
  # ========== FILTER EXPERIMENTAL TRIGGERS (BITS 0-15) ==========
  events_raw <- eeg_obj$events
  n_events_raw <- nrow(events_raw)
  recording_duration <- max(eeg_obj$times)
  
  # Create bit mask for experimental triggers (bits 0-15)
  # Mask: 0xFFFF = 65535 (binary: 0000000011111111 11111111)
  EXPERIMENTAL_TRIGGER_MASK <- 0xFFFF  # 2^16 - 1
  
  # Filter to keep only experimental triggers
  # Apply bitwise AND to extract only bits 0-15
  experimental_triggers <- bitwAnd(as.integer(events_raw$type), EXPERIMENTAL_TRIGGER_MASK)
  
  # Keep only events where experimental trigger bits are non-zero
  events_cleaned <- events_raw[experimental_triggers != 0, ]
  
  # Update trigger codes to reflect only experimental bits
  events_cleaned$type <- experimental_triggers[experimental_triggers != 0]
  
  # Track how many system triggers were removed
  n_system_triggers_removed <- n_events_raw - nrow(events_cleaned)
  
  filtering_notes <- character(0)
  if (n_system_triggers_removed > 0) {
    filtering_notes <- c(filtering_notes,
                         paste0("Removed ", format(n_system_triggers_removed, big.mark = ","),
                                " system status triggers (bits 16-23)"))
  }
  
  # ========== ADDITIONAL FILTERING (IF NOT RAW MODE) ==========
  if (mode != "raw" && nrow(events_cleaned) > 0) {
    
    # Step 2: Exclude user-specified types
    if (!is.null(exclude_types)) {
      n_before_exclude <- nrow(events_cleaned)
      events_cleaned <- events_cleaned[!events_cleaned$type %in% exclude_types, ]
      n_excluded <- n_before_exclude - nrow(events_cleaned)
      if (n_excluded > 0) {
        filtering_notes <- c(filtering_notes,
                             paste0("Excluded ", n_excluded, " events: ",
                                    paste(exclude_types, collapse = ", ")))
      }
    }
    
    if (nrow(events_cleaned) > 1) {
      # Step 3: Merge consecutive identical events
      n_before_merge <- nrow(events_cleaned)
      type_changes <- c(TRUE, events_cleaned$type[-1] != events_cleaned$type[-nrow(events_cleaned)])
      events_cleaned <- events_cleaned[type_changes, ]
      n_merged <- n_before_merge - nrow(events_cleaned)
      if (n_merged > 0) {
        filtering_notes <- c(filtering_notes,
                             paste0("Merged ", format(n_merged, big.mark = ","),
                                    " consecutive identical events"))
      }
      
      # Step 4: Remove very brief events (likely noise spikes)
      if (mode == "auto") {
        n_before_brief <- nrow(events_cleaned)
        durations <- c(diff(events_cleaned$onset_time), Inf)
        valid <- durations >= trigger_threshold
        events_cleaned <- events_cleaned[valid, ]
        n_brief <- n_before_brief - nrow(events_cleaned)
        if (n_brief > 0) {
          filtering_notes <- c(filtering_notes,
                               paste0("Removed ", n_brief, " brief events (< ",
                                      round(trigger_threshold * 1000, 1), " ms)"))
        }
      }
      
      # Step 5: Filter by minimum inter-event interval
      if (mode == "auto") {
        n_before_iei <- nrow(events_cleaned)
        iei <- diff(events_cleaned$onset_time)
        keep <- c(TRUE, iei >= min_iei)
        events_cleaned <- events_cleaned[keep, ]
        n_removed_iei <- n_before_iei - nrow(events_cleaned)
        if (n_removed_iei > 0) {
          filtering_notes <- c(filtering_notes,
                               paste0("Removed ", n_removed_iei, " events with IEI < ",
                                      round(min_iei * 1000, 1), " ms"))
        }
      }
    }
  }
  
  n_events_cleaned <- nrow(events_cleaned)
  
  # ========== CALCULATE STATISTICS ==========
  if (n_events_cleaned > 0) {
    event_types <- unique(events_cleaned$type)
    event_counts <- table(events_cleaned$type)
    
    if (n_events_cleaned > 1) {
      iei <- diff(events_cleaned$onset_time)
      mean_iei <- mean(iei)
      median_iei <- median(iei)
      event_rate <- (n_events_cleaned / recording_duration) * 60
    } else {
      mean_iei <- NA
      median_iei <- NA
      event_rate <- NA
    }
  } else {
    event_types <- character(0)
    event_counts <- integer(0)
    mean_iei <- NA
    median_iei <- NA
    event_rate <- NA
  }
  
  # ========== COMPACT OUTPUT ==========
  cat("\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  TRIGGER INSPECTION (BioSemi Experimental Triggers: Bits 0-15)\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  
  cat("Recording: ", round(recording_duration/60, 1), " min | ",
      eeg_obj$sampling_rate, " Hz | Mode: ", toupper(mode), "\n", sep = "")
  
  if (mode == "raw") {
    cat("Events: ", format(n_events_raw, big.mark = ","), " total → ",
        format(n_events_cleaned, big.mark = ","), " experimental triggers\n\n", sep = "")
  } else {
    reduction_pct <- round((1 - n_events_cleaned/n_events_raw)*100, 1)
    cat("Events: ", format(n_events_raw, big.mark = ","), " total → ",
        format(n_events_cleaned, big.mark = ","), " triggers (",
        reduction_pct, "% filtered)\n\n", sep = "")
  }
  
  # Show filtering steps
  if (length(filtering_notes) > 0) {
    cat("Filtering:\n")
    for (note in filtering_notes) {
      cat("  • ", note, "\n", sep = "")
    }
    cat("\n")
  }
  
  if (n_events_cleaned == 0) {
    cat("⚠ No experimental triggers found!\n\n")
    cat("Possible issues:\n")
    cat("  • All events were system status triggers (bits 16-23)\n")
    cat("  • No experimental triggers sent during recording\n")
    cat("  • Try: inspect_triggers(eeg, mode=\"raw\") to see all events\n\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
    return(invisible(list(
      n_events_raw = n_events_raw,
      n_events_cleaned = 0,
      n_system_removed = n_system_triggers_removed,
      events_cleaned = data.frame()
    )))
  }
  
  # Trigger summary table
  cat("EXPERIMENTAL TRIGGERS:\n")
  cat(sprintf("  %-10s %8s %7s %9s %9s\n", "Code", "Count", "%", "First(s)", "Last(s)"))
  cat("  ─────────────────────────────────────────────────────\n")
  
  for (event_type in names(sort(event_counts, decreasing = TRUE))) {
    count <- event_counts[event_type]
    pct <- (count / n_events_cleaned) * 100
    type_events <- events_cleaned[events_cleaned$type == event_type, ]
    cat(sprintf("  %-10s %8s %6.1f%% %9.2f %9.2f\n",
                event_type,
                format(count, big.mark = ","),
                pct,
                min(type_events$onset_time),
                max(type_events$onset_time)))
  }
  
  cat("\n")
  
  # Timing info
  if (!is.na(mean_iei)) {
    cat("Timing: Mean ITI = ", sprintf("%.3f", mean_iei), " s | ",
        "Median ITI = ", sprintf("%.3f", median_iei), " s | ",
        "Rate = ", round(event_rate, 1), "/min\n", sep = "")
    cat("\n")
  }
  
  # Epoching recommendation
  cat("Recommended epoch parameters:\n")
  suggested_tmax <- if (!is.na(mean_iei) && mean_iei > 1.0) {
    min(0.8, mean_iei * 0.6)
  } else if (!is.na(mean_iei) && mean_iei > 0.5) {
    0.5
  } else {
    0.3
  }
  
  if (length(event_types) == 1) {
    cat("  epoch_eeg(eeg, events = \"all\", tmin = -0.2, tmax = ",
        round(suggested_tmax, 2), ")\n", sep = "")
  } else if (length(event_types) <= 3) {
    cat("  epoch_eeg(eeg, events = c(\"", paste(event_types, collapse = "\", \""),
        "\"), tmin = -0.2, tmax = ", round(suggested_tmax, 2), ")\n", sep = "")
  } else {
    cat("  epoch_eeg(eeg, events = c(\"", event_types[1], "\", \"", event_types[2],
        "\", ...), tmin = -0.2, tmax = ", round(suggested_tmax, 2), ")\n", sep = "")
  }
  
  cat("\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  
  # ========== PLOTS ==========
  if (plot && n_events_cleaned > 0) {
    tryCatch({
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      
      # Plot 1: Timeline
      if (n_events_cleaned <= 5000) {
        plot(events_cleaned$onset_time, rep(1, n_events_cleaned),
             type = "h", lwd = 2, col = "steelblue",
             xlab = "Time (s)", ylab = "", yaxt = "n",
             main = paste0("Trigger Timeline (n=", format(n_events_cleaned, big.mark = ","), ")"),
             ylim = c(0, 1.2))
        abline(h = 0, col = "gray70")
      } else {
        bins <- hist(events_cleaned$onset_time, breaks = 100, plot = FALSE)
        plot(bins$mids, bins$counts, type = "h", lwd = 2, col = "steelblue",
             xlab = "Time (s)", ylab = "Count per bin",
             main = paste0("Trigger Timeline (n=", format(n_events_cleaned, big.mark = ","), ")"))
      }
      
      # Plot 2: Frequencies
      counts_sorted <- sort(event_counts, decreasing = TRUE)
      if (length(counts_sorted) > 15) counts_sorted <- counts_sorted[1:15]
      
      barplot(counts_sorted,
              col = rainbow(length(counts_sorted), alpha = 0.7),
              xlab = if (length(event_counts) > 15) "Trigger Code (Top 15)" else "Trigger Code",
              ylab = "Count",
              main = "Trigger Frequencies",
              las = 2, cex.names = 0.8)
      
      par(mfrow = c(1, 1))
    }, error = function(e) {
      warning("Could not generate plots: ", e$message, call. = FALSE)
      par(mfrow = c(1, 1))
    })
  }
  
  # ========== EXPORT ==========
  if (!is.null(export_csv)) {
    tryCatch({
      write.csv(events_cleaned, file = export_csv, row.names = FALSE)
      cat("✓ Exported to:", export_csv, "\n\n")
    }, error = function(e) {
      warning("Could not export: ", e$message, call. = FALSE)
    })
  }
  
  # ========== RETURN ==========
  invisible(list(
    n_events_raw = n_events_raw,
    n_events_cleaned = n_events_cleaned,
    n_system_removed = n_system_triggers_removed,
    event_types = event_types,
    event_counts = event_counts,
    mean_iei = mean_iei,
    median_iei = median_iei,
    event_rate = event_rate,
    events_cleaned = events_cleaned
  ))
}




#' ======================= EPOCH EEG FUNCTION ==========================
#' Epoch EEG Data Around Events
#'
#' Segments continuous EEG data into time-locked epochs (trials) centered on
#' experimental events or triggers. Uses BioSemi-aware trigger detection to
#' automatically filter experimental triggers from system status codes.
#'
#' @param eeg_obj An object of class 'eeg' containing continuous EEG data
#' @param events Character vector or numeric vector specifying which event types
#'   to epoch around. Can be trigger codes, event labels, or "all"
#' @param tmin Numeric value specifying epoch start time relative to event onset
#'   in seconds (default: -0.2)
#' @param tmax Numeric value specifying epoch end time relative to event onset
#'   in seconds (default: 0.8)
#' @param baseline Numeric vector c(start, end) specifying baseline correction
#'   window in seconds. Set to NULL to skip baseline correction (default: c(-0.2, 0))
#' @param baseline_method Character: "mean", "median", or "none" (default: "mean")
#' @param reject_threshold Numeric. Amplitude threshold in µV for epoch rejection.
#'   Set to NULL to disable (default: NULL)
#' @param preload Logical. If TRUE, extract all epoch data into memory (default: TRUE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A list of class 'eeg_epochs' containing:
#'   \describe{
#'     \item{data}{3D array (channels × timepoints × trials) of epoched data}
#'     \item{channels}{Character vector of channel names}
#'     \item{times}{Numeric vector of time points relative to event (seconds)}
#'     \item{events}{Data frame with event information for each epoch}
#'     \item{sampling_rate}{Numeric sampling rate in Hz}
#'     \item{tmin}{Start time of epoch window}
#'     \item{tmax}{End time of epoch window}
#'     \item{baseline}{Baseline correction window used}
#'     \item{n_epochs}{Total number of epochs extracted}
#'     \item{rejected}{Logical vector indicating rejected epochs}
#'     \item{rejection_log}{Data frame with rejection reasons}
#'   }
#'
#' @export
epoch_eeg <- function(eeg_obj,
                      events = "all",
                      tmin = -0.2,
                      tmax = 0.8,
                      baseline = c(-0.2, 0),
                      baseline_method = c("mean", "median", "none"),
                      reject_threshold = NULL,
                      preload = TRUE,
                      verbose = TRUE) {
  
  baseline_method <- match.arg(baseline_method)
  
  # ========== INPUT VALIDATION ==========
  if (!inherits(eeg_obj, "eeg")) {
    stop("Input must be an object of class 'eeg'", call. = FALSE)
  }
  
  if (is.null(eeg_obj$events) || nrow(eeg_obj$events) == 0) {
    stop("No events found in EEG object. Run inspect_triggers() first.", call. = FALSE)
  }
  
  if (tmin >= tmax) {
    stop("tmin (", tmin, ") must be less than tmax (", tmax, ")", call. = FALSE)
  }
  
  if (!is.null(baseline)) {
    if (length(baseline) != 2 || baseline[1] >= baseline[2]) {
      stop("baseline must be NULL or c(start, end) with start < end", call. = FALSE)
    }
    if (baseline[1] < tmin || baseline[2] > tmax) {
      stop("Baseline window must be within epoch window [", tmin, ", ", tmax, "]", call. = FALSE)
    }
  }
  
  if (baseline_method == "none") {
    baseline <- NULL
  }
  
  # ========== FILTER EXPERIMENTAL TRIGGERS ==========
  if (verbose) cat("Filtering experimental triggers (bits 0-15)...\n")
  
  EXPERIMENTAL_TRIGGER_MASK <- 0xFFFF
  experimental_triggers <- bitwAnd(as.integer(eeg_obj$events$type), EXPERIMENTAL_TRIGGER_MASK)
  valid_events <- eeg_obj$events[experimental_triggers != 0, ]
  valid_events$type <- experimental_triggers[experimental_triggers != 0]
  
  if (nrow(valid_events) == 0) {
    stop("No experimental triggers found. All events were system triggers.", call. = FALSE)
  }
  
  # ========== SELECT EVENTS TO EPOCH ==========
  if (is.character(events) && length(events) == 1 && events == "all") {
    selected_events <- valid_events
  } else {
    selected_events <- valid_events[valid_events$type %in% events, ]
  }
  
  if (nrow(selected_events) == 0) {
    stop("No events matched the selection criteria: ", paste(events, collapse = ", "), call. = FALSE)
  }
  
  if (verbose) {
    cat("Selected", nrow(selected_events), "events for epoching\n")
  }
  
  # ========== CALCULATE EPOCH PARAMETERS ==========
  smin <- round(tmin * eeg_obj$sampling_rate)
  smax <- round(tmax * eeg_obj$sampling_rate)
  n_samples_epoch <- smax - smin + 1
  
  epoch_times <- seq(tmin, tmax, length.out = n_samples_epoch)
  
  n_channels <- length(eeg_obj$channels)
  n_timepoints <- ncol(eeg_obj$data)
  n_trials <- nrow(selected_events)
  
  # ========== CALCULATE BASELINE INDICES ==========
  if (!is.null(baseline)) {
    baseline_smin <- which.min(abs(epoch_times - baseline[1]))
    baseline_smax <- which.min(abs(epoch_times - baseline[2]))
    baseline_indices <- baseline_smin:baseline_smax
  } else {
    baseline_indices <- NULL
  }
  
  # ========== EXTRACT EPOCHS ==========
  if (verbose) cat("Extracting epochs...\n")
  
  if (preload) {
    epoch_array <- array(NA, dim = c(n_channels, n_samples_epoch, n_trials))
    dimnames(epoch_array)[[1]] <- eeg_obj$channels
  } else {
    epoch_array <- NULL
  }
  
  valid_epochs <- rep(TRUE, n_trials)
  rejection_log <- data.frame(
    epoch_id = integer(0),
    event_type = character(0),
    event_time = numeric(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n_trials)) {
    event_onset <- selected_events$onset[i]
    epoch_start <- event_onset + smin
    epoch_end <- event_onset + smax
    
    # Check boundaries
    if (epoch_start < 1 || epoch_end > n_timepoints) {
      valid_epochs[i] <- FALSE
      rejection_log <- rbind(rejection_log, data.frame(
        epoch_id = i,
        event_type = selected_events$type[i],
        event_time = selected_events$onset_time[i],
        reason = "Extends beyond data boundaries",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Extract epoch data
    if (preload) {
      epoch_data <- eeg_obj$data[, epoch_start:epoch_end, drop = FALSE]
      
      # Apply baseline correction
      if (!is.null(baseline)) {
        for (ch in 1:n_channels) {
          if (baseline_method == "mean") {
            baseline_value <- mean(epoch_data[ch, baseline_indices], na.rm = TRUE)
          } else if (baseline_method == "median") {
            baseline_value <- median(epoch_data[ch, baseline_indices], na.rm = TRUE)
          }
          epoch_data[ch, ] <- epoch_data[ch, ] - baseline_value
        }
      }
      
      # Check rejection threshold
      if (!is.null(reject_threshold)) {
        epoch_max <- max(abs(epoch_data), na.rm = TRUE)
        if (epoch_max > reject_threshold) {
          valid_epochs[i] <- FALSE
          rejection_log <- rbind(rejection_log, data.frame(
            epoch_id = i,
            event_type = selected_events$type[i],
            event_time = selected_events$onset_time[i],
            reason = paste0("Amplitude exceeds threshold (", round(epoch_max, 1), " µV)"),
            stringsAsFactors = FALSE
          ))
          next
        }
      }
      
      epoch_array[, , i] <- epoch_data
    }
  }
  
  # ========== FILTER VALID EPOCHS ==========
  n_rejected <- sum(!valid_epochs)
  n_accepted <- sum(valid_epochs)
  
  if (n_accepted == 0) {
    stop("All epochs were rejected. Check data quality and rejection criteria.", call. = FALSE)
  }
  
  if (preload) {
    epoch_array <- epoch_array[, , valid_epochs, drop = FALSE]
  }
  
  selected_events <- selected_events[valid_epochs, ]
  selected_events$epoch_id <- seq_len(n_accepted)
  
  # ========== CREATE EPOCHS OBJECT ==========
  epochs_obj <- structure(
    list(
      data = epoch_array,
      channels = eeg_obj$channels,
      times = epoch_times,
      events = selected_events,
      sampling_rate = eeg_obj$sampling_rate,
      tmin = tmin,
      tmax = tmax,
      baseline = baseline,
      baseline_method = if (!is.null(baseline)) baseline_method else "none",
      n_epochs = n_accepted,
      rejected = !valid_epochs,
      rejection_log = rejection_log,
      metadata = eeg_obj$metadata
    ),
    class = "eeg_epochs"
  )
  
  # ========== SUMMARY ==========
  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("Epoching Complete\n")
    cat(strrep("=", 70), "\n")
    cat("  Epoch window:     [", round(tmin, 3), ", ", round(tmax, 3), "] s\n", sep = "")
    cat("  Channels:         ", n_channels, "\n")
    cat("  Epochs extracted: ", n_accepted, "\n")
    cat("  Epochs rejected:  ", n_rejected, "\n")
    if (!is.null(baseline)) {
      cat("  Baseline:         ", baseline_method, " [", 
          round(baseline[1], 3), ", ", round(baseline[2], 3), "] s\n", sep = "")
    }
    cat(strrep("=", 70), "\n\n")
  }
  
  return(epochs_obj)
}


#' ======================= AVERAGE EPOCHS FUNCTION ==========================
#' Average Epochs to Create ERPs
#'
#' Computes event-related potentials (ERPs) by averaging epoched data across
#' trials. Can average by condition or create grand average across all epochs.
#'
#' @param epochs_obj An object of class 'eeg_epochs' from epoch_eeg()
#' @param by Character. Group epochs by "event_type" or "all" for grand average
#'   (default: "event_type")
#' @param method Character. "mean" or "median" averaging (default: "mean")
#' @param return_se Logical. Return standard error of the mean (default: TRUE)
#' @param verbose Logical. Print summary (default: TRUE)
#'
#' @return A list of class 'eeg_evoked' containing:
#'   \describe{
#'     \item{evoked}{List of averaged data (one per condition)}
#'     \item{data}{3D array or matrix of averaged data}
#'     \item{se}{Standard error (if return_se = TRUE)}
#'     \item{n_trials}{Number of trials averaged per condition}
#'     \item{channels}{Channel names}
#'     \item{times}{Time vector}
#'     \item{conditions}{Condition names}
#'     \item{sampling_rate}{Sampling rate}
#'   }
#'
#' @export
average_epochs <- function(epochs_obj,
                           by = c("event_type", "all"),
                           method = c("mean", "median"),
                           return_se = TRUE,
                           verbose = TRUE) {
  
  by <- match.arg(by)
  method <- match.arg(method)
  
  # ========== VALIDATION ==========
  if (!inherits(epochs_obj, "eeg_epochs")) {
    stop("Input must be an object of class 'eeg_epochs'", call. = FALSE)
  }
  
  if (is.null(epochs_obj$data)) {
    stop("Epoch data not loaded. Set preload = TRUE in epoch_eeg()", call. = FALSE)
  }
  
  # ========== DETERMINE GROUPING ==========
  if (by == "all") {
    conditions <- list(all = seq_len(epochs_obj$n_epochs))
    condition_names <- "all"
  } else {
    event_types <- unique(epochs_obj$events$type)
    conditions <- lapply(event_types, function(et) {
      which(epochs_obj$events$type == et)
    })
    names(conditions) <- as.character(event_types)
    condition_names <- names(conditions)
  }
  
  n_channels <- length(epochs_obj$channels)
  n_times <- length(epochs_obj$times)
  n_conditions <- length(conditions)
  
  # ========== COMPUTE AVERAGES ==========
  if (verbose) cat("Computing ERPs for", n_conditions, "condition(s)...\n")
  
  evoked_list <- list()
  se_list <- list()
  n_trials_list <- integer(n_conditions)
  
  for (i in seq_along(conditions)) {
    cond_name <- names(conditions)[i]
    trial_indices <- conditions[[i]]
    n_trials <- length(trial_indices)
    n_trials_list[i] <- n_trials
    
    if (n_trials == 0) {
      warning("No trials found for condition: ", cond_name, call. = FALSE)
      next
    }
    
    # Extract trials for this condition
    cond_data <- epochs_obj$data[, , trial_indices, drop = FALSE]
    
    # Compute average
    if (method == "mean") {
      evoked_data <- apply(cond_data, c(1, 2), mean, na.rm = TRUE)
    } else if (method == "median") {
      evoked_data <- apply(cond_data, c(1, 2), median, na.rm = TRUE)
    }
    
    evoked_list[[cond_name]] <- evoked_data
    
    # Compute standard error
    if (return_se) {
      se_data <- apply(cond_data, c(1, 2), function(x) {
        sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      })
      se_list[[cond_name]] <- se_data
    }
    
    if (verbose) {
      cat("  ", cond_name, ": ", n_trials, " trials averaged\n", sep = "")
    }
  }
  
  # ========== CREATE EVOKED OBJECT ==========
  evoked_obj <- structure(
    list(
      evoked = evoked_list,
      se = if (return_se) se_list else NULL,
      n_trials = n_trials_list,
      conditions = condition_names,
      channels = epochs_obj$channels,
      times = epochs_obj$times,
      sampling_rate = epochs_obj$sampling_rate,
      tmin = epochs_obj$tmin,
      tmax = epochs_obj$tmax,
      baseline = epochs_obj$baseline,
      method = method
    ),
    class = "eeg_evoked"
  )
  
  if (verbose) {
    cat("\nERP averaging complete.\n")
    cat("Total trials across conditions:", sum(n_trials_list), "\n\n")
  }
  
  return(evoked_obj)
}


#' ======================= PLOT EPOCHS FUNCTION (OPTIMIZED) ==========================
#' Plot Epoched EEG Data
#'
#' Visualize epoched data for quality control and inspection. Optimized for
#' large datasets with many trials and channels.
#'
#' @param epochs_obj An object of class 'eeg_epochs' from epoch_eeg()
#' @param plot_type Character. Type of plot: "butterfly", "image", "trials", or "evoked"
#'   (default: "butterfly")
#' @param channels Character vector of channel names to plot. NULL = all channels
#'   (default: NULL)
#' @param events Numeric vector of event types to plot. NULL = all events
#'   (default: NULL)
#' @param max_trials Numeric. Maximum number of trials to plot (default: 50)
#' @param n_trials Numeric. For "trials" plot only (default: 20)
#' @param color_by Character. Color by "channel" or "event_type" (default: "channel")
#' @param add_legend Logical. Add legend (default: TRUE)
#' @param alpha Numeric. Transparency for individual trials (0-1, default: 0.1)
#'
#' @return None (generates plot)
#'
#' @export
plot_epochs <- function(epochs_obj,
                        plot_type = c("butterfly", "image", "trials", "evoked"),
                        channels = NULL,
                        events = NULL,
                        max_trials = 50,
                        n_trials = 20,
                        color_by = c("channel", "event_type"),
                        add_legend = TRUE,
                        alpha = 0.1) {
  
  plot_type <- match.arg(plot_type)
  color_by <- match.arg(color_by)
  
  # ========== VALIDATION ==========
  if (!inherits(epochs_obj, "eeg_epochs")) {
    stop("Input must be an object of class 'eeg_epochs'", call. = FALSE)
  }
  
  if (is.null(epochs_obj$data)) {
    stop("Epoch data not loaded. Set preload = TRUE in epoch_eeg()", call. = FALSE)
  }
  
  # ========== SELECT CHANNELS ==========
  if (is.null(channels)) {
    channel_indices <- seq_along(epochs_obj$channels)
    channel_names <- epochs_obj$channels
  } else {
    channel_indices <- which(epochs_obj$channels %in% channels)
    if (length(channel_indices) == 0) {
      stop("No matching channels found", call. = FALSE)
    }
    channel_names <- epochs_obj$channels[channel_indices]
  }
  
  # ========== SELECT EVENTS ==========
  if (is.null(events)) {
    trial_indices <- seq_len(epochs_obj$n_epochs)
  } else {
    trial_indices <- which(epochs_obj$events$type %in% events)
    if (length(trial_indices) == 0) {
      stop("No matching events found", call. = FALSE)
    }
  }
  
  # LIMIT TRIALS FOR BUTTERFLY PLOT
  if (plot_type == "butterfly" && length(trial_indices) > max_trials) {
    cat("Note: Plotting", max_trials, "of", length(trial_indices), 
        "trials for speed. Use max_trials to adjust.\n")
    trial_indices <- sample(trial_indices, max_trials)
  }
  
  plot_data <- epochs_obj$data[channel_indices, , trial_indices, drop = FALSE]
  
  # ========== BUTTERFLY PLOT (OPTIMIZED) ==========
  if (plot_type == "butterfly") {
    
    # Compute average first (faster)
    avg_data <- apply(plot_data, c(1, 2), mean, na.rm = TRUE)
    
    plot(range(epochs_obj$times), range(plot_data, na.rm = TRUE),
         type = "n",
         xlab = "Time (s)",
         ylab = "Amplitude (µV)",
         main = paste("Butterfly Plot -", length(trial_indices), "trials,",
                      length(channel_indices), "channels"))
    abline(h = 0, col = "gray70", lty = 2)
    abline(v = 0, col = "red", lty = 2, lwd = 2)
    
    # Plot individual trials (semi-transparent)
    trial_colors <- adjustcolor("gray30", alpha.f = alpha)
    for (ch_idx in seq_along(channel_indices)) {
      for (trial in seq_len(dim(plot_data)[3])) {
        lines(epochs_obj$times, plot_data[ch_idx, , trial],
              col = trial_colors, lwd = 0.3)
      }
    }
    
    # Plot channel averages (bold, colored)
    colors <- rainbow(length(channel_indices))
    for (ch_idx in seq_along(channel_indices)) {
      lines(epochs_obj$times, avg_data[ch_idx, ],
            col = colors[ch_idx], lwd = 2)
    }
    
    if (add_legend && length(channel_indices) <= 15) {
      legend("topright", legend = channel_names,
             col = colors, lty = 1, lwd = 2,
             bg = "white", cex = 0.7, ncol = 2)
    }
  }
  
  # ========== IMAGE PLOT ==========
  else if (plot_type == "image") {
    # Average across channels
    trial_data <- apply(plot_data, c(2, 3), mean, na.rm = TRUE)
    
    image(epochs_obj$times,
          seq_len(ncol(trial_data)),
          trial_data,
          xlab = "Time (s)",
          ylab = "Trial Number",
          main = paste("Epoch Image -", length(channel_indices), "channels averaged"),
          col = colorRampPalette(c("blue", "white", "red"))(100))
    abline(v = 0, col = "black", lwd = 2)
    contour(epochs_obj$times, seq_len(ncol(trial_data)), trial_data,
            add = TRUE, col = "black", lwd = 0.5, labcex = 0)
  }
  
  # ========== TRIAL-BY-TRIAL PLOT ==========
  else if (plot_type == "trials") {
    n_plot <- min(n_trials, length(trial_indices))
    n_rows <- ceiling(sqrt(n_plot))
    n_cols <- ceiling(n_plot / n_rows)
    
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 2, 1))
    
    for (i in seq_len(n_plot)) {
      trial_idx <- trial_indices[i]
      trial_data <- plot_data[, , i]
      
      plot(range(epochs_obj$times), range(trial_data, na.rm = TRUE),
           type = "n",
           xlab = "", ylab = "",
           main = paste("Trial", trial_idx))
      abline(h = 0, col = "gray70", lty = 2)
      abline(v = 0, col = "red", lty = 2)
      
      colors <- rainbow(length(channel_indices), alpha = 0.5)
      for (ch in seq_along(channel_indices)) {
        lines(epochs_obj$times, trial_data[ch, ],
              col = colors[ch], lwd = 0.8)
      }
    }
    
    par(mfrow = c(1, 1))
  }
  
  # ========== EVOKED PLOT ==========
  else if (plot_type == "evoked") {
    # Compute average and SE
    avg_data <- apply(plot_data, c(1, 2), mean, na.rm = TRUE)
    se_data <- apply(plot_data, c(1, 2), function(x) {
      sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    })
    
    plot(range(epochs_obj$times), 
         range(c(avg_data - 1.96*se_data, avg_data + 1.96*se_data), na.rm = TRUE),
         type = "n",
         xlab = "Time (s)",
         ylab = "Amplitude (µV)",
         main = paste("Evoked Response -", length(trial_indices), "trials averaged"))
    abline(h = 0, col = "gray70", lty = 2)
    abline(v = 0, col = "red", lty = 2, lwd = 2)
    
    colors <- rainbow(length(channel_indices))
    for (ch_idx in seq_along(channel_indices)) {
      # Plot 95% confidence interval
      polygon(c(epochs_obj$times, rev(epochs_obj$times)),
              c(avg_data[ch_idx, ] + 1.96*se_data[ch_idx, ],
                rev(avg_data[ch_idx, ] - 1.96*se_data[ch_idx, ])),
              col = adjustcolor(colors[ch_idx], alpha.f = 0.2),
              border = NA)
      
      # Plot mean
      lines(epochs_obj$times, avg_data[ch_idx, ],
            col = colors[ch_idx], lwd = 2)
    }
    
    if (add_legend && length(channel_indices) <= 10) {
      legend("topright", legend = channel_names,
             col = colors, lty = 1, lwd = 2,
             bg = "white", cex = 0.8)
    }
  }
}





