#' ============================================================================
#'                         EEG Object Summary Function
#' ============================================================================
#'
#' This module provides a dedicated inspection function for eeg objects.
#' Unlike print.eeg(), which gives a quick overview, eeg_summary() produces
#' a full diagnostic report with:
#'   - Stats computed separately for EEG and EXG channels
#'   - Suspicious channel flagging (flat, noisy, outlier amplitude)
#'
#' Channel classification relies on detect_external_channels() from
#' setexchannels.R. The Status channel is automatically excluded from
#' all statistical computations.
#'
#' Author: Christos Dalamarinis
#' Date: March 2026
#' Status: Not tested with a testfile - Pending 21/03/2026
#' ============================================================================
#'
#' Full Diagnostic Summary of an EEG Object
#'
#' Produces a detailed report of an eeg object, split into a clean EEG
#' channel stats block, an EXG channel stats block, and a suspicious
#' channel flagging section. Unlike print.eeg(), stats are never polluted
#' by the Status channel or EXG signals.
#'
#' @param eeg_obj An object of class 'eeg'.
#'
#' @param flag_flat_threshold Numeric. Channels with std below this value
#'   (in µV) are flagged as potentially flat / disconnected. Default: 0.5.
#'
#' @param flag_amplitude_threshold Numeric. Channels with any absolute
#'   amplitude exceeding this value (in µV) are flagged as potentially
#'   noisy or saturated. Default: 500.
#'
#' @param flag_outlier_sd_multiplier Numeric. A channel whose std exceeds
#'   this multiple of the median std across all EEG channels is flagged
#'   as an outlier. Default: 3.
#'
#' @param show_exg Logical. If TRUE (default), prints the EXG channel
#'   stats block. Set to FALSE to suppress it when EXG channels are not
#'   of interest.
#'
#' @return Invisibly returns a list with three elements:
#'   \item{eeg_stats}{Data frame of per-channel stats for EEG channels}
#'   \item{exg_stats}{Data frame of per-channel stats for EXG channels
#'     (NULL if no EXG channels present)}
#'   \item{flags}{Data frame of flagged channels with reason column
#'     (NULL if no flags raised)}
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   eeg_summary(sub1)
#'
#'   # Stricter flat-channel threshold
#'   eeg_summary(sub1, flag_flat_threshold = 1.0)
#'
#'   # Suppress EXG block
#'   eeg_summary(sub1, show_exg = FALSE)
#'
#'   # Capture the return value for downstream use
#'   report <- eeg_summary(sub1)
#'   report$flags   # inspect flagged channels programmatically
#' }
#'
#' @export
eeg_summary <- function(eeg_obj,
                        flag_flat_threshold        = 0.5,
                        flag_amplitude_threshold   = 500,
                        flag_outlier_sd_multiplier = 3,
                        show_exg                   = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("eeg_summary() requires an object of class 'eeg'.\n",
         "Got: ", class(eeg_obj)[1])
  }
  
  if (!is.matrix(eeg_obj$data)) {
    stop("eeg_obj$data must be a numeric matrix (channels x timepoints).")
  }
  
  # ========== CLASSIFY CHANNELS ==========
  
  all_channels <- eeg_obj$channels
  
  # Status channel — always excluded from stats
  status_idx <- which(grepl("^status$", tolower(all_channels)))
  
  # EXG channels — two-pass detection:
  #
  # Pass 1: detect_external_channels() catches original names (e.g. "EXG1")
  # Pass 2: regex fallback catches renamed channels where the original name
  #         is preserved in parentheses by apply_external_labels(keep_original=TRUE)
  #         e.g. "MASTOID LEFT (EXG5)", "EOG LEFT (EXG1)", "ECG (EXG3)"
  #
  # Known external name patterns to look for inside parentheses:
  #   EXG1-EXG8, GSR1-GSR2, Plet, Temp, Resp, Erg1-Erg2
  
  exg_pattern <- paste0(
    "\\(",                          # opening parenthesis
    "(",                            # start capture group
    "EXG[1-8]",                     # BioSemi external electrodes
    "|GSR[12]",                     # Galvanic skin response
    "|Plet",                        # Plethysmograph
    "|Temp",                        # Temperature
    "|Resp",                        # Respiration
    "|Erg[12]",                     # Ergo/AUX
    ")",                            # end capture group
    "\\)"                           # closing parenthesis
  )
  
  # Pass 1 - original names
  exg_names_pass1 <- detect_external_channels(all_channels)
  
  # Pass 2 — renamed channels (keep_original = TRUE pattern)
  exg_names_pass2 <- all_channels[grepl(exg_pattern, all_channels,
                                        ignore.case = TRUE)]
  
  # Combine both passes, remove duplicates
  exg_names <- unique(c(exg_names_pass1, exg_names_pass2))
  exg_idx   <- which(all_channels %in% exg_names)
  
  # EEG channels — everything that is not Status and not EXG
  eeg_idx <- setdiff(seq_along(all_channels), c(status_idx, exg_idx))
  
  if (length(eeg_idx) == 0) {
    stop("No EEG channels found after excluding EXG and Status channels.\n",
         "Check that channel names follow standard 10-20/10-10 conventions.")
  }
  
  # ========== COMPUTE PER-CHANNEL STATS HELPER ==========
  
  .channel_stats <- function(data_matrix, idx, channel_names) {
    if (length(idx) == 0) return(NULL)
    
    result <- data.frame(
      channel  = channel_names[idx],
      min_uv   = numeric(length(idx)),
      max_uv   = numeric(length(idx)),
      mean_uv  = numeric(length(idx)),
      std_uv   = numeric(length(idx)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(idx)) {
      ch_data           <- data_matrix[idx[i], ]
      result$min_uv[i]  <- round(min(ch_data),  2)
      result$max_uv[i]  <- round(max(ch_data),  2)
      result$mean_uv[i] <- round(mean(ch_data), 2)
      result$std_uv[i]  <- round(sd(ch_data),   2)
    }
    
    return(result)
  }
  
  eeg_stats <- .channel_stats(eeg_obj$data, eeg_idx, all_channels)
  exg_stats <- .channel_stats(eeg_obj$data, exg_idx, all_channels)
  
  # ========== FLAG SUSPICIOUS EEG CHANNELS ==========
  
  flags      <- data.frame(
    channel = character(0),
    reason  = character(0),
    value   = character(0),
    stringsAsFactors = FALSE
  )
  
  median_std <- median(eeg_stats$std_uv, na.rm = TRUE)
  
  for (i in seq_len(nrow(eeg_stats))) {
    
    ch   <- eeg_stats$channel[i]
    std  <- eeg_stats$std_uv[i]
    amax <- max(abs(eeg_stats$min_uv[i]), abs(eeg_stats$max_uv[i]))
    
    # Flat / disconnected
    if (std < flag_flat_threshold) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Flat / possibly disconnected",
        value   = paste0("std = ", std, " µV  (threshold: < ",
                         flag_flat_threshold, " µV)"),
        stringsAsFactors = FALSE
      ))
    }
    
    # Saturated / excessive amplitude
    if (amax > flag_amplitude_threshold) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Excessive amplitude",
        value   = paste0("peak |amplitude| = ", round(amax, 2),
                         " µV  (threshold: > ", flag_amplitude_threshold, " µV)"),
        stringsAsFactors = FALSE
      ))
    }
    
    # Outlier std relative to the rest of the array
    if (std > flag_outlier_sd_multiplier * median_std) {
      flags <- rbind(flags, data.frame(
        channel = ch,
        reason  = "Outlier noise level",
        value   = paste0("std = ", std, " µV  (", flag_outlier_sd_multiplier,
                         " x median std = ",
                         round(flag_outlier_sd_multiplier * median_std, 2), " µV)"),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(flags) == 0) flags <- NULL
  
  # ========== PRINT REPORT ==========
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("EEG Summary Report\n")
  cat(strrep("=", 70), "\n\n")
  
  # ---- Recording info ----
  cat("RECORDING INFORMATION:\n")
  cat("  Sampling rate  : ", eeg_obj$sampling_rate, " Hz\n", sep = "")
  cat("  Duration       : ", sprintf("%.2f", ncol(eeg_obj$data) /
                                       eeg_obj$sampling_rate), " seconds\n", sep = "")
  cat("  Total channels : ", length(all_channels), "\n", sep = "")
  cat("    EEG channels : ", length(eeg_idx), "\n", sep = "")
  cat("    EXG channels : ", length(exg_idx), "\n", sep = "")
  if (length(status_idx) > 0) {
    cat("    Status channel : present (excluded from all stats)\n")
  }
  cat("  Reference      : ", eeg_obj$reference, "\n", sep = "")
  
  # ---- EEG global stats ----
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("EEG CHANNEL STATISTICS  (n = ", length(eeg_idx), " channels)\n", sep = "")
  cat(strrep("-", 70), "\n")
  
  cat("  Amplitude range : ", round(min(eeg_stats$min_uv), 2),
      "  to  ", round(max(eeg_stats$max_uv), 2), "  µV\n", sep = "")
  cat("  Mean amplitude  : ", round(mean(eeg_stats$mean_uv), 2), "  µV\n", sep = "")
  cat("  Median std      : ", round(median(eeg_stats$std_uv), 2), "  µV\n", sep = "")
  cat("  Std range       : ", round(min(eeg_stats$std_uv), 2),
      "  to  ", round(max(eeg_stats$std_uv), 2), "  µV\n", sep = "")
  
  # ---- EEG per-channel table ----
  cat("\n  Per-channel breakdown:\n")
  cat(sprintf("  %-12s  %10s  %10s  %10s  %10s\n",
              "Channel", "Min (µV)", "Max (µV)", "Mean (µV)", "Std (µV)"))
  cat("  ", strrep("-", 58), "\n", sep = "")
  for (i in seq_len(nrow(eeg_stats))) {
    cat(sprintf("  %-12s  %10.2f  %10.2f  %10.2f  %10.2f\n",
                eeg_stats$channel[i],
                eeg_stats$min_uv[i],
                eeg_stats$max_uv[i],
                eeg_stats$mean_uv[i],
                eeg_stats$std_uv[i]))
  }
  
  # ---- EXG stats ----
  if (show_exg && !is.null(exg_stats)) {
    cat("\n")
    cat(strrep("-", 70), "\n")
    cat("EXG CHANNEL STATISTICS  (n = ", length(exg_idx), " channels)\n", sep = "")
    cat("  Note: stats reported for reference only — EXG channels are not\n")
    cat("  included in EEG stats above and are not flagged below.\n")
    cat(strrep("-", 70), "\n")
    
    cat(sprintf("  %-30s  %10s  %10s  %10s  %10s\n",
                "Channel", "Min (µV)", "Max (µV)", "Mean (µV)", "Std (µV)"))
    cat("  ", strrep("-", 74), "\n", sep = "")
    for (i in seq_len(nrow(exg_stats))) {
      cat(sprintf("  %-30s  %10.2f  %10.2f  %10.2f  %10.2f\n",
                  exg_stats$channel[i],
                  exg_stats$min_uv[i],
                  exg_stats$max_uv[i],
                  exg_stats$mean_uv[i],
                  exg_stats$std_uv[i]))
    }
  } else if (show_exg && length(exg_idx) == 0) {
    cat("\n  No EXG channels detected.\n")
  }
  
  # ---- Suspicious channels ----
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("SUSPICIOUS CHANNEL FLAGS\n")
  cat(strrep("-", 70), "\n")
  cat("  Thresholds used:\n")
  cat("    Flat            : std < ", flag_flat_threshold, " µV\n", sep = "")
  cat("    Excessive amp.  : peak |amplitude| > ", flag_amplitude_threshold, " µV\n", sep = "")
  cat("    Outlier noise   : std > ", flag_outlier_sd_multiplier,
      " x median std (", round(flag_outlier_sd_multiplier * median_std, 2), " µV)\n", sep = "")
  cat("\n")
  
  if (is.null(flags)) {
    cat("  [OK] No suspicious channels detected.\n")
  } else {
    cat("  [!] ", nrow(flags), " flag(s) raised:\n\n", sep = "")
    for (i in seq_len(nrow(flags))) {
      cat("  Channel : ", flags$channel[i], "\n", sep = "")
      cat("  Reason  : ", flags$reason[i],  "\n", sep = "")
      cat("  Detail  : ", flags$value[i],   "\n", sep = "")
      if (i < nrow(flags)) cat("\n")
    }
  }
  
  cat("\n")
  cat(strrep("=", 70), "\n\n")
  
  # ========== RETURN INVISIBLY ==========
  
  invisible(list(
    eeg_stats = eeg_stats,
    exg_stats = exg_stats,
    flags     = flags
  ))
}