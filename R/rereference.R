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
#'   exclude from computing the reference (e.g., EOG, EMG). Excluded
#'   channels are left untouched entirely - they are not part of the EEG
#'   reference scheme, so the reference is not subtracted from them either.
#'   Channels listed in \code{eeg$bads} are handled differently - see
#'   Details.
#' @param copy Logical. If \code{TRUE}, return a new EEG object and leave
#'   \code{eeg} unchanged. If \code{FALSE}, re-reference \code{eeg} in place
#'   (if your class supports this).
#' @param drop_ref Logical. If \code{TRUE}, remove the channel(s) used as
#'   \code{ref} from the returned data/channels after re-referencing. Has no
#'   effect when \code{ref = "average"} (a warning is issued instead, since
#'   the reference channels are then the entire contributing set). Default
#'   \code{FALSE}, matching the behavior of MNE-Python's
#'   \code{set_eeg_reference()}, which also keeps reference channels in the
#'   data unless dropped explicitly.
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
#' \strong{Bad channels (\code{eeg$bads}):} channels marked bad are excluded
#' from computing the reference signal (so noisy/broken data does not
#' contaminate the average), but the reference is still subtracted from
#' them, same as any other EEG channel. This differs from \code{exclude},
#' which is meant for channels that were never on the EEG reference scheme
#' (EOG/EMG/etc.) and are left untouched by this function entirely. Keeping
#' bad channels in the same reference frame means they stay consistent with
#' the rest of the data once they are later fixed (e.g., interpolated). A
#' channel requested directly via \code{ref} that is also listed in
#' \code{eeg$bads} is treated as unavailable, the same as an excluded one.
#'
#' After re-referencing to specific channels (e.g., linked mastoids), those
#' channels remain in the data by default, transformed into mirror-image
#' signals of each other (e.g., for \code{ref = c("M1", "M2")}, the new
#' \code{M1 == -M2}). This makes them perfectly (anti)correlated, reducing
#' the effective rank of the data by one — relevant before rank-sensitive
#' steps such as ICA or PCA/whitening. Set \code{drop_ref = TRUE} to remove
#' them from the returned object instead.
#'
#' @examples
#' \dontrun{
#' # Common average reference
#' eeg_avg <- eeg_rereference(eeg, ref = "average")
#'
#' # Re-reference to linked mastoids M1/M2, keeping them in the data
#' eeg_mastoids <- eeg_rereference(eeg, ref = c("M1", "M2"))
#'
#' # Re-reference to linked mastoids and drop them afterward
#' eeg_mastoids <- eeg_rereference(eeg, ref = c("M1", "M2"), drop_ref = TRUE)
#' }
#'
#' @export
eeg_rereference <- function(eeg,
                            ref = "average",
                            exclude = NULL,
                            copy = TRUE,
                            drop_ref = FALSE) {
  
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
  
  # ---- handle exclude (channels never on the EEG reference scheme) ----
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

  # ---- handle bads (shared eeg$bads list: noisy/broken EEG channels) ----
  # Bad channels must not pollute the reference computation, but the
  # reference is still subtracted from them (see @details) - so they get
  # their own index set instead of being folded into 'exclude'.
  if (!is.null(eeg_out$bads) && length(eeg_out$bads) > 0) {
    bads_idx <- match(eeg_out$bads, chan_names)
    bads_idx <- bads_idx[!is.na(bads_idx)]
  } else {
    bads_idx <- integer(0)
  }

  # channels that receive the reference subtraction: everything except
  # 'exclude' (bad channels ARE still re-referenced)
  apply_idx <- setdiff(seq_along(chan_names), excl_idx)

  # channels available for computing the reference: also drops bads
  contrib_idx <- setdiff(apply_idx, bads_idx)

  if (length(contrib_idx) == 0) {
    stop("No channels left to compute reference after applying 'exclude' and 'bads'.")
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
    # also ensure ref channels are part of contributors (not excluded, not bad)
    ref_idx <- intersect(ref_idx, contrib_idx)
    if (length(ref_idx) == 0) {
      stop("Reference channels are all excluded by 'exclude' or marked bad in eeg$bads.")
    }
  }

  # ---- compute reference signal (1 value per time point) ----
  ref_signal <- colMeans(signals[ref_idx, , drop = FALSE], na.rm = TRUE)

  # ---- subtract reference from all non-excluded channels (bads included) ----
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

  # ---- optionally drop the reference channel(s) from the output ----
  if (drop_ref) {
    if (identical(ref, "average")) {
      warning("drop_ref is ignored when ref = 'average' (the reference ",
              "channels are the entire contributing set).")
    } else {
      dropped_names <- chan_names[ref_idx]

      if (length(ref_idx) >= nrow(eeg_out$data)) {
        stop("Dropping the reference channel(s) would remove all channels ",
             "from the data.")
      }

      eeg_out$data     <- eeg_out$data[-ref_idx, , drop = FALSE]
      eeg_out$channels <- eeg_out$channels[-ref_idx]

      eeg_out$preprocessing_history <- c(
        eeg_out$preprocessing_history,
        list(paste0("Dropped reference channel(s): ",
                    paste(dropped_names, collapse = ", ")))
      )
    }
  }

  eeg_out
}