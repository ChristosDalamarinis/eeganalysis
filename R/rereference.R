#' ============================================================
#'                  EEG Re-referencing Utilities
#' ============================================================
#'
#' This module provides functions to re-reference EEG recordings to
#' different reference schemes (e.g., common average or selected channels)
#' within the eeganalysis framework for downstream preprocessing and analysis.
#'
#' Author: Christos Dalamarinis
#' Date: Jan - 2026
#' ============================================================
#'
#' Re-reference EEG data
#'
#' Change the reference of EEG signals to a new set of reference channels
#' or to the common average across channels.
#'
#' @param eeg An EEG object as used in this package.
#' @param ref Character or integer vector specifying the new reference
#'   channel(s). Use \code{"average"} to apply a common average reference.
#' @param exclude Optional character or integer vector of channels to
#'   exclude from computing the reference (e.g., EOG, EMG).
#' @param copy Logical. If \code{TRUE}, return a new EEG object and leave
#'   \code{eeg} unchanged. If \code{FALSE}, re-reference \code{eeg} in place
#'   (if your class supports this).
#'
#' @return An EEG object with updated reference.
#' @details
#' EEG data are measured as voltage differences between each electrode
#' and a reference electrode. Re-referencing applies a linear transformation
#' so that all channels are expressed relative to the new reference defined.
#'
#' If \code{ref = "average"}, the reference signal is computed as the mean
#' across all non-excluded EEG channels at each time point (common average
#' reference).
#' When one or several channels are specified in \code{ref}, the reference
#' signal is the mean of those channels at each time point.
#'
#' @examples
#' \dontrun{
#' # Common average reference
#' eeg_avg <- eeg_rereference(eeg, ref = "average")
#'
#' # Re-reference to linked mastoids M1/M2
#' eeg_mastoids <- eeg_rereference(eeg, ref = c("M1", "M2"))
#' }
#'
#' @export
eeg_rereference <- function(eeg,
                            ref = "average",
                            exclude = NULL,
                            copy = TRUE) {
  
  # make a copy if requested
  if (copy) {
    eeg_out <- eeg
  } else {
    eeg_out <- eeg
  }
  
  # ---- get data matrix and channel names ----
  signals   <- eeg_out$data       # channels x time
  chan_names <- eeg_out$channels
  
  # ---- validate dimensions ----
  if (nrow(signals) != length(chan_names)) {
    stop("Number of columns in signals does not match number of channel names.")
  }
  
  # ---- handle exclude ----
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      excl_idx <- exclude
    } else {
      excl_idx <- match(exclude, chan_names)
    }
    excl_idx <- excl_idx[!is.na(excl_idx)]
  } else {
    excl_idx <- integer(0)
  }
  
  # channels available for referencing (contributors)
  contrib_idx <- setdiff(seq_along(chan_names), excl_idx)
  
  if (length(contrib_idx) == 0) {
    stop("No channels left to compute reference after applying 'exclude'.")
  }
  
  # ---- determine reference channels ----
  if (identical(ref, "average")) {
    ref_idx <- contrib_idx
  } else {
    if (is.numeric(ref)) {
      ref_idx <- ref
    } else {
      ref_idx <- match(ref, chan_names)
    }
    if (any(is.na(ref_idx))) {
      stop("Some reference channels specified in 'ref' were not found.")
    }
    # also ensure ref channels are part of contributors
    ref_idx <- intersect(ref_idx, contrib_idx)
    if (length(ref_idx) == 0) {
      stop("Reference channels are all excluded by 'exclude'.")
    }
  }
  
  # ---- compute reference signal (1 value per time point) ----
  ref_signal <- colMeans(signals[ref_idx, , drop = FALSE], na.rm = TRUE)
  
  # ---- subtract reference from all non-excluded channels ----
  apply_idx <- contrib_idx
  
  signals[apply_idx, ] <- sweep(signals[apply_idx, , drop = FALSE],
                                2,
                                ref_signal,
                                FUN = "-")
  
  # ---- update signals in object ----
  eeg_out$data <- signals
  
  # ---- build reference label ----
  ref_label <- if (identical(ref, "average")) {
    "Common Average"
  } else {
    paste(chan_names[ref_idx], collapse = "+")
  }
  
  # ---- update $reference field (read by print.eeg) ----
  eeg_out$reference <- ref_label
  
  # ---- update metadata$reference_scheme if metadata exists ----
  if (!is.null(eeg_out$metadata)) {
    eeg_out$metadata$reference_scheme <- ref_label
  }
  
  # ---- append to preprocessing history ----
  eeg_out$preprocessing_history <- c(
    eeg_out$preprocessing_history,
    list(paste0("Re-referenced to: ", ref_label))
  )
  
  eeg_out
}