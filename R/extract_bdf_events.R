#' ==============================================================================
#'              Extract EVENT codes from BDF STATUS channel
#' ==============================================================================
#'
#' Function to read BioSemi BDF files and extract trigger/event codes 
#' from the STATUS channel using native R implementation (no edfReader dependency)
#'
#' Author: Christos Dalamarinis
#' Date: Jan - 2026
#' 
#' BioSemi BDF STATUS Channel Structure (per BioSemi format spec):
#' - Bits 1-16: Trigger codes (up to 65535)
#' - Bits 16-23: System status (CMS, battery, etc.)
#' - Each sample = 24-bit value (3 bytes, little-endian)
#' 
#' Reference: https://www.biosemi.com/faq/file_format.htm
#' 
#' 
#' 
#' ==============================================================================
#'
#' Extract EVENT codes from BDF file or EEG object
#'
#' Extracts trigger/event codes from either a BDF file path or an already-loaded
#' EEG object. When a file path is provided, it uses \code{read_bdf_native()}
#' internally. When an EEG object is provided, it extracts events directly.
#'
#' @param data Either:
#'   \itemize{
#'     \item A character string with path to .bdf file, OR
#'     \item An \code{eeg} object (from \code{read_bdf_native()})
#'   }
#' @param verbose Logical. If TRUE (default), print detailed progress messages.
#'   Only applies when reading from file.
#'
#' @return Data frame with columns:
#'   - onset: sample index of the event
#'   - onset_time: time in seconds
#'   - type: trigger code as character
#'   - description: text label (e.g. "Trigger: 255")
#'
#' The returned data frame has attributes:
#'   - bdf_file: file path (if applicable)
#'   - sample_rate: sampling rate in Hz
#'   - n_samples: total number of samples
#'   - duration_sec: total duration in seconds
#'
#' @examples
#' \dontrun{
#'   # Option 1: Extract from file directly
#'   events <- extract_bdf_events("sub-01_task_eeg.bdf")
#'   
#'   # Option 2: Extract from already-loaded EEG object
#'   eeg_data <- read_bdf_native("sub-01_task_eeg.bdf")
#'   events <- extract_bdf_events(eeg_data, verbose = FALSE)
#'   
#'   # View results
#'   head(events)
#'   table(events$type)
#' }
#'
#' @export
extract_bdf_events <- function(data, verbose = TRUE) {
  
  # ========== STEP 0: Determine input type ==========
  is_file <- is.character(data) && length(data) == 1
  is_eeg_obj <- inherits(data, "eeg") || (is.list(data) && all(c("events", "sampling_rate") %in% names(data)))
  
  if (!is_file && !is_eeg_obj) {
    stop("'data' must be either a file path (character string) or an eeg object", 
         call. = FALSE)
  }
  
  # ========== BRANCH 1: File path provided ==========
  if (is_file) {
    bdf_file <- data
    
    # Validate file
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
      cat("  [OK] File found: ", basename(bdf_file), "\n", sep = "")
      cat("  [OK] File size: ", round(file_size_mb, 1), " MB\n", sep = "")
    }
    
    # Read BDF file
    if (verbose) cat("\n[2/5] Reading BDF file (native reader)...\n")
    
    eeg_obj <- read_bdf_native(bdf_file, verbose = verbose)
    
    if (verbose) cat("\n[3/5] Extracting events from EEG object...\n")
    
  } else {
    # ========== BRANCH 2: EEG object provided ==========
    eeg_obj <- data
    bdf_file <- eeg_obj$metadata$file_path %||% "unknown"
    
    if (verbose) {
      cat("\n========== BDF Event Extraction (from EEG object) ==========\n")
      cat("[1/2] Extracting events from loaded EEG object...\n")
    }
  }
  
  # ========== COMMON: Extract Events ==========
  events <- eeg_obj$events
  
  if (is.null(events)) {
    warning("No events structure found in EEG object")
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
      cat("  [WARN] No events found in STATUS channel\n")
    } else {
      cat("  [OK] Events extracted: ", nrow(events), "\n", sep = "")
    }
  }
  
  # ========== Add Metadata Attributes ==========
  step_label <- if (is_file) "[4/5]" else "[2/2]"
  if (verbose) cat("\n", step_label, " Adding metadata attributes...\n", sep = "")
  
  attr(events, "bdf_file")    <- bdf_file
  attr(events, "sample_rate") <- eeg_obj$sampling_rate
  attr(events, "n_samples")   <- length(eeg_obj$times)
  attr(events, "duration_sec")<- max(eeg_obj$times)
  
  if (verbose) {
    cat("  [OK] Sample rate: ", eeg_obj$sampling_rate, " Hz\n", sep = "")
    cat("  [OK] Total samples: ", length(eeg_obj$times), "\n", sep = "")
    cat("  [OK] Duration: ", round(max(eeg_obj$times), 2), " seconds\n", sep = "")
  }
  
  # ========== Summary Statistics ==========
  if (verbose && nrow(events) > 0) {
    step_label <- if (is_file) "[5/5]" else "[Summary]"
    cat("\n", step_label, " Event extraction summary:\n", sep = "")
    
    unique_codes <- unique(events$type)
    cat("  - Total events: ", nrow(events), "\n", sep = "")
    cat("  - Unique trigger codes: ", length(unique_codes), "\n", sep = "")
    
    if (length(unique_codes) <= 10) {
      cat("  - Trigger codes: ", paste(sort(as.numeric(unique_codes)), collapse = ", "), "\n", sep = "")
    } else {
      cat("  - Trigger code range: ", min(as.numeric(unique_codes)), 
          " to ", max(as.numeric(unique_codes)), "\n", sep = "")
    }
    
    cat("  - First event at: ", round(min(events$onset_time), 3), " sec\n", sep = "")
    cat("  - Last event at: ", round(max(events$onset_time), 3), " sec\n", sep = "")
    
    if (nrow(events) > 1) {
      mean_interval <- mean(diff(events$onset_time))
      cat("  - Mean inter-event interval: ", round(mean_interval, 3), " sec\n", sep = "")
      cat("  - Approximate event rate: ", round(1/mean_interval, 2), " events/sec\n", sep = "")
    }
    
    cat("\n [OK] Event extraction complete!\n")
    cat("==========================================\n\n")
  } else if (verbose) {
    step_label <- if (is_file) "[5/5]" else ""
    cat("\n", step_label, " Complete (no events found)\n", sep = "")
    cat("==========================================\n\n")
  }
  
  return(events)
}

# Helper for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x



#' Summary of extracted BDF events
#'
#' @param events Data frame returned from \code{extract_bdf_events()}
#' 
#' @export
summary_bdf_events <- function(events) {
  
  cat("====== BDF Event Summary ======\n")
  cat("File: ", attr(events, "bdf_file"), "\n", sep = "")
  cat("Sample rate: ", attr(events, "sample_rate"), " Hz\n", sep = "")
  cat("Recording duration: ",
      round(attr(events, "duration_sec") %||% 0, 2), " seconds\n", sep = "")
  cat("Total samples: ", attr(events, "n_samples") %||% "N/A", "\n", sep = "")
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


#' Validate and summarize BDF events
#'
#' Provides a concise quality-control summary of extracted events.
#' Only creates plots if issues are detected or explicitly requested.
#'
#' @param data Either a file path, eeg object, or events data frame
#' @param verbose Logical. Print detailed summary
#' @param plot Logical or character. Options:
#'   \itemize{
#'     \item FALSE - no plot (default, just print summary)
#'     \item TRUE - auto-decide based on detected issues
#'     \item "always" - always create diagnostic plot
#'     \item "timeline" - show event timeline
#'     \item "sequence" - show event sequence patterns
#'   }
#'
#' @return Invisibly returns a list with:
#'   - summary: data frame with per-trigger statistics
#'   - issues: character vector of detected problems
#'   - recommendation: suggested next steps
#'
#' @export
validate_bdf_events <- function(data, verbose = TRUE, plot = FALSE) {
  
  # Extract events (reuse existing logic)
  events <- .extract_events_from_input(data)
  
  if (nrow(events) == 0) {
    message(" [WARN] No events found")
    return(invisible(list(summary = NULL, issues = "no_events")))
  }
  
  # ========== ANALYSIS ==========
  n_events <- nrow(events)
  unique_triggers <- sort(as.numeric(unique(events$type)))
  duration_sec <- attr(events, "duration_sec") %||% max(events$onset_time)
  sample_rate <- attr(events, "sample_rate") %||% NA
  
  # Per-trigger statistics
  trigger_summary <- as.data.frame(table(events$type))
  names(trigger_summary) <- c("trigger", "count")
  trigger_summary$trigger <- as.numeric(as.character(trigger_summary$trigger))
  trigger_summary$percentage <- round(trigger_summary$count / n_events * 100, 1)
  trigger_summary <- trigger_summary[order(-trigger_summary$count), ]
  
  # Timing analysis
  if (n_events > 1) {
    intervals_ms <- diff(events$onset_time) * 1000
    mean_isi <- mean(intervals_ms)
    sd_isi <- sd(intervals_ms)
    min_isi <- min(intervals_ms)
    max_isi <- max(intervals_ms)
    cv_isi <- sd_isi / mean_isi  # Coefficient of variation
  }
  
  # ========== ISSUE DETECTION ==========
  issues <- character(0)
  
  # Check 1: Very few events
  if (n_events < 10) {
    issues <- c(issues, "very_few_events")
  }
  
  # Check 2: Imbalanced conditions (if >1 trigger type)
  if (length(unique_triggers) > 1) {
    max_count <- max(trigger_summary$count)
    min_count <- min(trigger_summary$count)
    if (max_count / min_count > 5) {
      issues <- c(issues, "imbalanced_conditions")
    }
  }
  
  # Check 3: Suspiciously short intervals (<50ms)
  # Use isTRUE() to safely handle NA values
  if (n_events > 1 && exists("min_isi") && isTRUE(min_isi < 50)) {
    issues <- c(issues, "double_triggers")
  }
  
  # Check 4: High variability in timing
  # Use isTRUE() to safely handle NA/NaN/Inf values
  if (n_events > 1 && exists("cv_isi") && isTRUE(cv_isi > 0.5)) {
    issues <- c(issues, "irregular_timing")
  }
  
  # Check 5: Unexpected trigger codes (e.g., >255 for 8-bit systems)
  if (any(unique_triggers > 255)) {
    issues <- c(issues, "unusual_trigger_values")
  }
  
  # ========== PRINT SUMMARY ==========
  if (verbose) {
    cat("\n")
    cat("========================================\n")
    cat("           EVENT VALIDATION            \n")
    cat("========================================\n\n")
    
    cat("File:", basename(attr(events, "bdf_file") %||% "Unknown"), "\n")
    cat("Duration:", round(duration_sec, 1), "sec\n")
    cat("Total events:", n_events, "\n")
    cat("Event rate:", round(n_events / duration_sec, 2), "events/sec\n\n")
    
    cat("TRIGGER SUMMARY:\n")
    print(trigger_summary, row.names = FALSE)
    cat("\n")
    
    if (n_events > 1) {
      cat("TIMING STATISTICS:\n")
      cat("  Mean ISI:", round(mean_isi, 1), "ms\n")
      cat("  SD:", round(sd_isi, 1), "ms\n")
      cat("  Range:", round(min_isi, 1), "-", round(max_isi, 1), "ms\n")
      cat("  Regularity:", ifelse(cv_isi < 0.2, "Regular", 
                                  ifelse(cv_isi < 0.5, "Moderate", "Irregular")), 
          "(CV =", round(cv_isi, 2), ")\n\n")
    }
    
    # Report issues
    if (length(issues) > 0) {
      cat("[WARN] POTENTIAL ISSUES DETECTED:\n")
      if ("very_few_events" %in% issues) {
        cat("  - Very few events (<10) - check recording\n")
      }
      if ("imbalanced_conditions" %in% issues) {
        cat("  - Imbalanced conditions - check experimental design\n")
      }
      if ("double_triggers" %in% issues) {
        cat("  - Very short intervals (<50ms) - possible double-triggers\n")
        short_intervals <- which(intervals_ms < 50)
        cat("    Occurs at events:", paste(head(short_intervals, 5), collapse = ", "))
        if (length(short_intervals) > 5) cat(" ...")
        cat("\n")
      }
      if ("irregular_timing" %in% issues) {
        cat("  - Highly variable timing - check experimental control\n")
      }
      if ("unusual_trigger_values" %in% issues) {
        cat("  - Trigger values >255 detected:", paste(unique_triggers[unique_triggers > 255], collapse = ", "), "\n")
      }
      cat("\n")
    } else {
      cat("[OK] No issues detected\n\n")
    }
    
    cat("========================================\n\n")
  }
  
  # ========== DECIDE ON PLOTTING ==========
  should_plot <- FALSE
  if (is.logical(plot)) {
    if (plot == TRUE && length(issues) > 0) should_plot <- TRUE
    if (plot == TRUE && length(issues) == 0 && verbose) {
      cat("No issues detected - skipping plot. Use plot='always' to force.\n\n")
    }
  } else if (is.character(plot)) {
    should_plot <- TRUE
  }
  
  # ========== CREATE DIAGNOSTIC PLOT (ONLY IF NEEDED) ==========
  if (should_plot) {
    
    plot_type <- if (is.character(plot)) plot else "timeline"
    
    if (plot_type == "timeline" || plot_type == "always") {
      # Simple, informative timeline
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
      
      # Panel 1: Event sequence
      plot(events$onset_time, as.numeric(events$type),
           pch = "|", cex = 2, col = "darkblue",
           xlab = "Time (sec)", ylab = "Trigger",
           main = "Event Timeline")
      grid()
      
      # Panel 2: Inter-event intervals (if >1 event)
      if (n_events > 1) {
        plot(events$onset_time[-1], intervals_ms,
             type = "h", col = ifelse(intervals_ms < 50, "red", "darkgreen"),
             xlab = "Time (sec)", ylab = "Interval (ms)",
             main = "Inter-Event Intervals (red = <50ms)")
        abline(h = mean_isi, col = "blue", lty = 2, lwd = 2)
        grid()
      }
      
      par(mfrow = c(1, 1))
    }
    
    if (plot_type == "sequence") {
      # Show trigger sequence pattern
      plot(seq_along(events$type), as.numeric(events$type),
           type = "b", pch = 19, cex = 0.8,
           xlab = "Event number", ylab = "Trigger code",
           main = "Trigger Sequence Pattern",
           col = as.numeric(as.factor(events$type)))
      grid()
    }
  }
  
  # ========== RETURN ==========
  result <- list(
    summary = trigger_summary,
    issues = if (length(issues) > 0) issues else NULL,
    events = events,
    timing = if (n_events > 1) list(
      mean_isi_ms = mean_isi,
      sd_isi_ms = sd_isi,
      cv = cv_isi
    ) else NULL
  )
  
  invisible(result)
}


# Helper to extract events from various input types
.extract_events_from_input <- function(data) {
  if (is.character(data) && length(data) == 1) {
    return(extract_bdf_events(data, verbose = FALSE))
  } else if (inherits(data, "eeg") || (is.list(data) && "events" %in% names(data))) {
    return(data$events)
  } else if (is.data.frame(data) && all(c("onset", "onset_time", "type") %in% names(data))) {
    return(data)
  } else {
    stop("Invalid input type")
  }
}


# ==============================================================================
# label_bdf_events.R
# ==============================================================================
#
# Interactive labeling of BDF trigger codes
#
# Author: Christos Dalamarinis
# Date: March 2026
# ==============================================================================


#' Interactively label BDF trigger codes
#'
#' Displays a summary table of all unique trigger codes found in an events
#' data frame, then prompts the user to provide a plain-text label for each
#' one. Pressing Enter without typing a label assigns \code{"unknown"} to that
#' code. The resulting labels are stored in a new \code{label} column and the
#' full code-to-label mapping is attached as the \code{trigger_labels}
#' attribute so it can be inspected or reused programmatically.
#'
#' Author: Christos Dalamarinis
#' Date: March 2026
#'
#' @param events  Data frame returned by \code{extract_bdf_events()}.
#'   Must contain at minimum the columns \code{type} and \code{onset_time}.
#' @param verbose Logical. If \code{TRUE} (default), print progress messages
#'   and the final mapping summary.
#'
#' @return The input data frame with one additional column:
#'   \describe{
#'     \item{\code{label}}{Character. User-supplied label for each event, or
#'       \code{"unknown"} when the user skipped that trigger code.}
#'   }
#'   The returned data frame also carries a new attribute:
#'   \describe{
#'     \item{\code{trigger_labels}}{Named character vector mapping each trigger
#'       code (name) to its label (value), e.g.
#'       \code{c("1" = "standard", "2" = "target", "128" = "unknown")}.}
#'   }
#'
#' @details
#' This function requires an interactive R session because it uses
#' \code{\link[base]{readline}} to collect input. It will stop with an
#' informative error if called in a non-interactive context (e.g. during
#' \code{R CMD check} or automated testing). Use the
#' \code{apply_trigger_labels()} helper to apply a pre-built mapping without
#' interactivity.
#'
#' @examples
#' \dontrun{
#'   # Typical workflow
#'   events <- extract_bdf_events("sub-01_task_eeg.bdf")
#'   events <- label_bdf_events(events)
#'
#'   # Inspect the mapping that was created
#'   attr(events, "trigger_labels")
#'
#'   # Apply the same mapping to another participant's events
#'   mapping <- attr(events, "trigger_labels")
#'   events_p2 <- extract_bdf_events("sub-02_task_eeg.bdf")
#'   events_p2 <- apply_trigger_labels(events_p2, mapping)
#' }
#'
#' @seealso \code{\link{extract_bdf_events}}, \code{\link{apply_trigger_labels}}
#'
#' @export
label_bdf_events <- function(events, verbose = TRUE) {

  # ========== STEP 1: Validate input ==========

  if (!is.data.frame(events)) {
    stop("'events' must be a data frame returned by extract_bdf_events().",
         call. = FALSE)
  }

  required_cols <- c("type", "onset_time")
  missing_cols  <- setdiff(required_cols, names(events))
  if (length(missing_cols) > 0) {
    stop("'events' is missing required column(s): ",
         paste(missing_cols, collapse = ", "),
         ". Was it created by extract_bdf_events()?",
         call. = FALSE)
  }

  if (nrow(events) == 0) {
    warning("'events' contains zero rows. Returning unchanged.", call. = FALSE)
    events$label <- character(0)
    return(invisible(events))
  }

  # ========== STEP 2: Check interactive context ==========

  if (!interactive()) {
    stop(
      "label_bdf_events() requires an interactive R session.\n",
      "To apply labels programmatically use apply_trigger_labels().",
      call. = FALSE
    )
  }

  # ========== STEP 3: Build summary table of unique trigger codes ==========

  unique_codes <- sort(unique(events$type),
                       method = "radix",
                       decreasing = FALSE)
  n_codes      <- length(unique_codes)

  # Summary statistics per code
  summary_tbl <- do.call(rbind, lapply(unique_codes, function(code) {
    idx <- events$type == code
    data.frame(
      Code       = code,
      Count      = sum(idx),
      First_sec  = round(min(events$onset_time[idx]), 3),
      Last_sec   = round(max(events$onset_time[idx]), 3),
      stringsAsFactors = FALSE
    )
  }))

  # ========== STEP 4: Print header + summary table ==========

  if (verbose) {
    cat("\n")
    cat("========== BDF Trigger Labeling ==========\n")
    cat("  Unique trigger codes : ", n_codes,  "\n", sep = "")
    cat("  Total events         : ", nrow(events), "\n", sep = "")
    cat("\n")
    cat("  Trigger summary:\n")
    cat("  ", strrep("-", 52), "\n", sep = "")
    cat(sprintf("  %-10s  %8s  %12s  %12s\n",
                "Code", "Count", "First (sec)", "Last (sec)"))
    cat("  ", strrep("-", 52), "\n", sep = "")
    for (i in seq_len(nrow(summary_tbl))) {
      cat(sprintf("  %-10s  %8d  %12.3f  %12.3f\n",
                  summary_tbl$Code[i],
                  summary_tbl$Count[i],
                  summary_tbl$First_sec[i],
                  summary_tbl$Last_sec[i]))
    }
    cat("  ", strrep("-", 52), "\n", sep = "")
    cat("\n")
    cat("  For each code, type a label and press Enter.\n")
    cat("  Press Enter with no text to assign 'unknown'.\n")
    cat("\n")
  }

  # ========== STEP 5: Interactive prompt — one code at a time ==========

  mapping <- setNames(character(n_codes), unique_codes)  # named character vector

  for (i in seq_len(n_codes)) {
    code  <- unique_codes[i]
    count <- summary_tbl$Count[i]

    prompt_str <- sprintf(
      "  [%d/%d]  Code: %-8s  (%d events)  Label: ",
      i, n_codes, code, count
    )

    raw_input <- readline(prompt = prompt_str)
    raw_input <- trimws(raw_input)

    # Empty input → auto-assign "unknown"
    if (nchar(raw_input) == 0) {
      mapping[code] <- "unknown"
      if (verbose) cat("         → assigned: 'unknown'\n")
    } else {
      mapping[code] <- raw_input
    }
  }

  # ========== STEP 6: Apply mapping — add 'label' column ==========

  events$label <- mapping[events$type]

  # ========== STEP 7: Attach mapping as attribute ==========

  attr(events, "trigger_labels") <- mapping

  # ========== STEP 8: Print final mapping summary ==========

  if (verbose) {
    cat("\n")
    cat("  Labeling complete. Final mapping:\n")
    cat("  ", strrep("-", 36), "\n", sep = "")
    cat(sprintf("  %-12s  %-20s\n", "Code", "Label"))
    cat("  ", strrep("-", 36), "\n", sep = "")
    for (code in names(mapping)) {
      cat(sprintf("  %-12s  %-20s\n", code, mapping[code]))
    }
    cat("  ", strrep("-", 36), "\n", sep = "")
    n_unknown <- sum(mapping == "unknown")
    if (n_unknown > 0) {
      cat("\n  [NOTE]", n_unknown, "code(s) left as 'unknown'.\n")
    }
    cat("==========================================\n\n")
  }

  invisible(events)
}


# ==============================================================================
#' Apply a pre-built trigger label mapping to an events data frame
#'
#' Non-interactive companion to \code{\link{label_bdf_events}}. Applies an
#' existing code-to-label mapping (typically the \code{trigger_labels}
#' attribute from a previously labeled events data frame) to a new events
#' data frame. Useful for applying the same experimental mapping across
#' multiple participants without re-entering labels each time.
#'
#' @param events  Data frame returned by \code{extract_bdf_events()}.
#' @param mapping Named character vector of the form
#'   \code{c("code" = "label", ...)}, as returned by
#'   \code{attr(labeled_events, "trigger_labels")}.
#' @param unmapped_label Character scalar. Label to assign to any trigger code
#'   in \code{events} that is absent from \code{mapping}.
#'   Default: \code{"unknown"}.
#' @param verbose Logical. If \code{TRUE} (default), report any unmapped codes.
#'
#' @return The input data frame with a \code{label} column added and the
#'   \code{trigger_labels} attribute attached (identical to \code{mapping},
#'   plus any auto-filled unmapped codes).
#'
#' @examples
#' \dontrun{
#'   # Build mapping interactively for participant 1
#'   events_p1 <- extract_bdf_events("sub-01_task_eeg.bdf")
#'   events_p1 <- label_bdf_events(events_p1)
#'   mapping   <- attr(events_p1, "trigger_labels")
#'
#'   # Re-use for participant 2 — no prompts needed
#'   events_p2 <- extract_bdf_events("sub-02_task_eeg.bdf")
#'   events_p2 <- apply_trigger_labels(events_p2, mapping)
#' }
#'
#' @seealso \code{\link{label_bdf_events}}
#'
#' @export
apply_trigger_labels <- function(events,
                                 mapping,
                                 unmapped_label = "unknown",
                                 verbose        = TRUE) {

  # ========== Validate inputs ==========

  if (!is.data.frame(events)) {
    stop("'events' must be a data frame.", call. = FALSE)
  }
  if (!"type" %in% names(events)) {
    stop("'events' must contain a 'type' column.", call. = FALSE)
  }
  if (!is.character(mapping) || is.null(names(mapping))) {
    stop("'mapping' must be a named character vector.", call. = FALSE)
  }

  if (nrow(events) == 0) {
    warning("'events' contains zero rows. Returning unchanged.", call. = FALSE)
    events$label <- character(0)
    return(invisible(events))
  }

  # ========== Detect any codes not covered by the mapping ==========

  present_codes <- unique(events$type)
  unmapped      <- setdiff(present_codes, names(mapping))

  if (length(unmapped) > 0) {
    if (verbose) {
      cat("[apply_trigger_labels] WARNING: ",
          length(unmapped),
          " code(s) not found in mapping — assigned '",
          unmapped_label, "':\n", sep = "")
      cat("  ", paste(unmapped, collapse = ", "), "\n\n", sep = "")
    }
    extra          <- setNames(rep(unmapped_label, length(unmapped)), unmapped)
    full_mapping   <- c(mapping, extra)
  } else {
    full_mapping <- mapping
  }

  # ========== Apply ==========

  events$label <- full_mapping[events$type]

  attr(events, "trigger_labels") <- full_mapping

  invisible(events)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
#
#  # --- Interactive labeling (requires console) ---
#  events        <- extract_bdf_events("sub-01_task_eeg.bdf")
#  events        <- label_bdf_events(events)
#
#  # Inspect result
#  head(events[, c("onset_time", "type", "label")])
#  attr(events, "trigger_labels")
#
#  # --- Re-use mapping for a second participant ---
#  mapping       <- attr(events, "trigger_labels")
#  events_p2     <- extract_bdf_events("sub-02_task_eeg.bdf")
#  events_p2     <- apply_trigger_labels(events_p2, mapping)
#
#  # --- Save mapping to disk for future sessions ---
#  saveRDS(mapping, "task_trigger_mapping.rds")
#  mapping       <- readRDS("task_trigger_mapping.rds")
# ==============================================================================




# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
#
# # Extract events from a BDF file
# events <- extract_bdf_events("sub-01_task_eeg.bdf")
#
# # View the first few events
# head(events)
#
# # Get summary statistics
# summary_bdf_events(events)
#
# # Filter specific trigger codes (e.g., only code 255)
# target_events <- events[events$type == "255", ]
# nrow(target_events)
#
# # Export to CSV for further analysis
# write.csv(events, "extracted_events.csv", row.names = FALSE)
#
# ==============================================================================
