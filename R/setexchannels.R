#' Interactively Identify and Label External Channels (Database-Driven)
#'
#' @description
#' This function detects external channels in imported EEG data using the 
#' official electrode database and allows the user to interactively label 
#' each one (e.g., "EMG", "ECG", "EOG").
#'
#' @param data A data frame or matrix containing EEG data with channel names,
#'   or a character vector of channel names
#' @param channel_col Character string specifying the column name containing 
#'   channel labels. If NULL, uses column names. Default is NULL.
#'
#' @return A named list with two elements:
#'   \item{labels}{Named vector mapping external channel names to their labels}
#'   \item{summary}{Data frame summarizing the labeling}
#'
#' @details
#' This function relies on get_electrode_database() from channel_info2.R to 
#' accurately identify which channels are external electrodes. It recognizes:
#' - EXG1-EXG8 (BioSemi external electrodes)
#' - GSR1/GSR2 (Galvanic skin response)
#' - Erg1/Erg2 (Ergo/AUX inputs)
#' - Resp (Respiration belt)
#' - Plet (Plethysmograph)
#' - Temp (Temperature)
#'
#' @examples
#' \dontrun{
#' # After importing BioSemi data
#' eeg_data <- import_biosemi("myfile.bdf")
#' external_labels <- identify_external_channels(eeg_data)
#' }
#'
#' @export
identify_external_channels <- function(data, channel_col = NULL) {
  
  # ========================================================================
  # EXTRACT CHANNEL NAMES
  # ========================================================================
  
  if (is.character(data)) {
    # Character vector of channel names
    channel_names <- trimws(data)
  } else if (is.data.frame(data) || is.matrix(data)) {
    # Data frame or matrix - use column names
    if (is.null(channel_col)) {
      channel_names <- colnames(data)
    } else {
      if (!channel_col %in% colnames(data)) {
        stop("Specified channel_col '", channel_col, "' not found in data.")
      }
      channel_names <- unique(data[[channel_col]])
    }
  } else if (is.list(data)) {
    # List - try common locations
    if (!is.null(data$channels)) {
      channel_names <- trimws(data$channels)
    } else if (!is.null(data$channel_names)) {
      channel_names <- trimws(data$channel_names)
    } else {
      stop("Cannot find channel names in list. Please provide channel names directly.")
    }
  } else {
    stop("Invalid input type. Expected data frame, matrix, list, or character vector.")
  }
  
  if (is.null(channel_names) || length(channel_names) == 0) {
    stop("No channel names found or extracted.")
  }
  
  # ========================================================================
  # GET ELECTRODE DATABASE
  # ========================================================================
  
  electrode_db <- get_electrode_database()
  
  # ========================================================================
  # IDENTIFY EXTERNAL CHANNELS USING DATABASE
  # ========================================================================
  
  external_channels <- character()
  
  for (ch_name in channel_names) {
    # Case-insensitive lookup
    ch_lower <- tolower(ch_name)
    
    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]
      
      # Check if it's an external channel
      # External channels have position_type = "External" or similar categories
      if (electrode_info$position_type %in% c("External", "GSR", "Ergo/AUX", 
                                              "Respiration", "Plethysmograph", 
                                              "Temperature")) {
        external_channels <- c(external_channels, ch_name)
      }
    }
  }
  
  # ========================================================================
  # REPORT FINDINGS
  # ========================================================================
  
  if (length(external_channels) == 0) {
    cat("\n✓ No external channels detected in the dataset.\n\n")
    return(list(labels = character(0), summary = data.frame()))
  }
  
  # Report number of external channels found
  cat("\n========================================\n")
  cat("  EXTERNAL CHANNELS DETECTED\n")
  cat("========================================\n\n")
  cat("Found", length(external_channels), "external channel(s):\n")
  
  # Show channel names with their database descriptions
  for (ch_name in external_channels) {
    ch_lower <- tolower(ch_name)
    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]
      cat("  - ", ch_name, " (", electrode_info$position_name, ")\n", sep = "")
    } else {
      cat("  - ", ch_name, "\n", sep = "")
    }
  }
  cat("\n")
  
  # ========================================================================
  # INTERACTIVE LABELING
  # ========================================================================
  
  cat("Please label each external channel.\n")
  cat("Common labels: EMG (muscle), EOG (eye), ECG (heart), GSR (skin conductance)\n")
  cat("Press ENTER to skip a channel.\n\n")
  
  labels <- character(length(external_channels))
  names(labels) <- external_channels
  
  for (i in seq_along(external_channels)) {
    channel <- external_channels[i]
    
    # Get database info for context
    ch_lower <- tolower(channel)
    context_info <- ""
    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]
      context_info <- paste0(" [", electrode_info$position_type, "]")
    }
    
    repeat {
      user_input <- readline(prompt = sprintf("[%d/%d] %s%s = ", 
                                              i, 
                                              length(external_channels), 
                                              channel,
                                              context_info))
      
      # Trim whitespace
      user_input <- trimws(user_input)
      
      # If empty, mark as unlabeled
      if (user_input == "") {
        labels[channel] <- "Unlabeled"
        cat("  → Skipped (marked as Unlabeled)\n\n")
        break
      }
      
      # Validate input (only letters, numbers, underscores, hyphens, spaces)
      if (grepl("^[A-Za-z0-9_ -]+$", user_input)) {
        labels[channel] <- user_input
        cat("  → Labeled as:", user_input, "\n\n")
        break
      } else {
        cat("  ✗ Invalid input. Please use only letters, numbers, hyphens, underscores, or spaces.\n")
      }
    }
  }
  
  # ========================================================================
  # CREATE SUMMARY
  # ========================================================================
  
  summary_df <- data.frame(
    Channel = names(labels),
    Label = unname(labels),
    Type = character(length(labels)),
    Description = character(length(labels)),
    stringsAsFactors = FALSE
  )
  
  # Add database info to summary
  for (i in 1:nrow(summary_df)) {
    ch_name <- summary_df$Channel[i]
    ch_lower <- tolower(ch_name)
    
    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]
      summary_df$Type[i] <- electrode_info$position_type
      summary_df$Description[i] <- electrode_info$position_name
    } else {
      summary_df$Type[i] <- "Unknown"
      summary_df$Description[i] <- "Not in database"
    }
  }
  
  # ========================================================================
  # DISPLAY SUMMARY
  # ========================================================================
  
  cat("========================================\n")
  cat("  LABELING SUMMARY\n")
  cat("========================================\n\n")
  print(summary_df, row.names = FALSE)
  cat("\n")
  
  # Count labeled vs unlabeled
  n_labeled <- sum(labels != "Unlabeled")
  n_unlabeled <- sum(labels == "Unlabeled")
  
  cat("Total labeled: ", n_labeled, " / ", length(external_channels), "\n", sep = "")
  if (n_unlabeled > 0) {
    cat("⚠ ", n_unlabeled, " channel(s) remain unlabeled.\n", sep = "")
  }
  cat("\n")
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  invisible(list(
    labels = labels,
    summary = summary_df
  ))
}


#' Detect External Channels (Non-Interactive, Database-Driven)
#'
#' @description
#' Detects external channels using the electrode database without user interaction.
#' Useful for scripting and automation.
#'
#' @param data A data frame, matrix, list, or character vector containing channel names
#' @param channel_col Character string specifying column with channel names (for data frames)
#'
#' @return Character vector of detected external channel names
#'
#' @export
detect_external_channels <- function(data, channel_col = NULL) {
  
  # Extract channel names
  if (is.character(data)) {
    channel_names <- trimws(data)
  } else if (is.data.frame(data) || is.matrix(data)) {
    if (is.null(channel_col)) {
      channel_names <- colnames(data)
    } else {
      if (!channel_col %in% colnames(data)) {
        stop("Specified channel_col '", channel_col, "' not found in data.")
      }
      channel_names <- unique(data[[channel_col]])
    }
  } else if (is.list(data)) {
    if (!is.null(data$channels)) {
      channel_names <- trimws(data$channels)
    } else if (!is.null(data$channel_names)) {
      channel_names <- trimws(data$channel_names)
    } else {
      stop("Cannot find channel names in list.")
    }
  } else {
    stop("Invalid input type.")
  }
  
  # Get electrode database
  electrode_db <- get_electrode_database()
  
  # Identify external channels
  external_channels <- character()
  
  for (ch_name in channel_names) {
    ch_lower <- tolower(ch_name)
    
    if (ch_lower %in% names(electrode_db)) {
      electrode_info <- electrode_db[[ch_lower]]
      
      if (electrode_info$position_type %in% c("External", "GSR", "Ergo/AUX", 
                                              "Respiration", "Plethysmograph", 
                                              "Temperature")) {
        external_channels <- c(external_channels, ch_name)
      }
    }
  }
  
  return(external_channels)
}


#' Apply External Channel Labels to Data
#'
#' @description
#' Renames external channels in the dataset based on user-provided labels.
#' Handles eeg class objects, data frames/matrices, and list structures.
#'
#' @param data An eeg object, data frame, matrix, or list containing EEG data
#' @param labels Named vector of labels (output from identify_external_channels)
#' @param keep_original Logical; if TRUE, keeps original names in parentheses
#'
#' @return Modified data with renamed channels
#'
#' @examples
#' \dontrun{
#' external_info <- identify_external_channels(eeg_data)
#' eeg_data_labeled <- apply_external_labels(eeg_data, external_info$labels)
#' }
#'
#' @export
apply_external_labels <- function(data, labels, keep_original = TRUE) {
  
  if (length(labels) == 0) {
    message("No labels to apply.")
    return(data)
  }
  
  # Handle different data structures
  if (is.list(data) && !is.data.frame(data)) {
    # List structure (eeg object or other list)
    
    # Find where channel names are stored
    if ("channels" %in% names(data)) {
      current_names <- data$channels
      channel_location <- "channels"
    } else if ("channel_names" %in% names(data)) {
      current_names <- data$channel_names
      channel_location <- "channel_names"
    } else {
      stop("Cannot find channel names in data structure.\n",
           "Expected data$channels or data$channel_names\n",
           "Available names: ", paste(names(data), collapse = ", "))
    }
    
    # Create new names
    new_names <- current_names
    n_changed <- 0
    
    for (i in seq_along(current_names)) {
      if (current_names[i] %in% names(labels)) {
        label <- labels[current_names[i]]
        
        if (label != "Unlabeled") {
          if (keep_original) {
            new_names[i] <- paste0(label, " (", current_names[i], ")")
          } else {
            new_names[i] <- label
          }
          n_changed <- n_changed + 1
        }
      }
    }
    
    # Apply new names back to the data structure
    data[[channel_location]] <- new_names
    
    cat("✓ Applied", n_changed, "external channel label(s).\n")
    
    return(data)
    
  } else if (is.data.frame(data) || is.matrix(data)) {
    # Data frame or matrix
    
    # Get current column names
    current_names <- colnames(data)
    
    if (is.null(current_names)) {
      stop("Data frame/matrix has no column names.")
    }
    
    # Create new names
    new_names <- current_names
    n_changed <- 0
    
    for (i in seq_along(current_names)) {
      if (current_names[i] %in% names(labels)) {
        label <- labels[current_names[i]]
        
        if (label != "Unlabeled") {
          if (keep_original) {
            new_names[i] <- paste0(label, " (", current_names[i], ")")
          } else {
            new_names[i] <- label
          }
          n_changed <- n_changed + 1
        }
      }
    }
    
    # Apply new names
    colnames(data) <- new_names
    
    cat("✓ Applied", n_changed, "external channel label(s).\n")
    
    return(data)
    
  } else {
    stop("Invalid data type. Expected eeg object, data frame, matrix, or list.\n",
         "Got: ", class(data)[1])
  }
}

