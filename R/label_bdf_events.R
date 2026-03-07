# ==============================================================================
#                             Label EEG Events
# ==============================================================================
#
# Interactive labeling of BDF trigger codes
#
# Author: Christos Dalamarinis
# Date: March 2026
# ==============================================================================
#'
#' Displays a summary table of all unique trigger codes found in an events
#' data frame, then prompts the user to provide a plain-text label for each
#' one. Pressing Enter without typing a label assigns \code{"unknown"} to that
#' code. The resulting labels are stored in a new \code{label} column and the
#' full code-to-label mapping is attached as the \code{trigger_labels}
#' attribute so it can be inspected or reused programmatically.
#' 
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


#' ==============================================================================
#'        Apply a pre-built trigger label mapping to an events data frame
#' ==============================================================================
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