#' ============================================================================
#'                            ERP Analysis Functions
#' ============================================================================
#'
#' This module provides functions for Event-Related Potential (ERP) analysis
#' of epoched EEG data. It includes tools for computing ERPs, extracting
#' ERP components (peaks, latencies, amplitudes), measuring ERP features,
#' and conducting statistical analyses on ERP waveforms.
#'
#' The module supports common ERP analysis workflows including:
#' - Grand average computation across subjects/conditions
#' - Peak detection and latency measurement
#' - Mean amplitude extraction in predefined time windows
#' - Component identification (P1, N1, P2, N2, P3, etc.)
#' - Statistical comparison between conditions
#' - Topographic mapping of ERP features
#'
#' Author: Christos Dalamarinis
#' Date: February 2026
#' ============================================================================
#' 
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