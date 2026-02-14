#' ============================================================
#' EEG Re-referencing Utilities  
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
  
  # ---- get data matrix and channel names (adapt to your class) ----
  # Example assumptions:
  #   - signal matrix: eeg_out$signals  (channels x time or time x channels; adapt)
  #   - channel names: eeg_out$chan_info$name
  #
  # Replace these with the actual field names of your package.
  signals <- eeg_out$data # assuming channels x time from the "eeg_class.R"
  chan_names <- eeg_out$channels
  
  # check dimensions and orientation
  # Here assume signals is samples x channels
  if (nrow(signals) != length(chan_names)) {
    stop("Number of columns in signals does not match number of channel names.")
  }
  
  # determine which channels to use as reference and as contributors
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
  
  # ---- compute reference signal (samples x 1) ----
  ref_signal <- colMeans(signals[ref_idx, , drop = FALSE], na.rm = TRUE)
  
  # ---- subtract reference from all channels except excluded? ----
  # Common practice: subtract from all EEG channels (including ref channels),
  # but typically still exclude non-EEG channels (EOG, triggers) from being modified.
  # Here we only modify contrib_idx and ref_idx, not excluded channels.
  # If you want to modify all, replace 'contrib_idx' with 'seq_along(chan_names)'.
  apply_idx <- contrib_idx
  
  signals[apply_idx, ] <- sweep(signals[apply_idx, , drop = FALSE],
                                1,
                                ref_signal,
                                FUN = "-")
  
  # update signals in object
  eeg_out$data <- signals
  
  # optionally update metadata about reference
  if (!is.null(eeg_out$meta)) {
    eeg_out$meta$reference <- if (identical(ref, "average")) {
      "average"
    } else {
      paste(chan_names[ref_idx], collapse = "+")
    }
  }
  
  eeg_out
}
