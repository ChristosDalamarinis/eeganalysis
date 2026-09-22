# ============================================================================
#              Bipolar Referencing: Derive a Channel by Subtraction
# ============================================================================
#
# One helper, one job: given two raw electrode channels, build a new channel
# that is their difference (anode - cathode) and add it to the recording -
# the standard way to turn a pair of single-ended EXG inputs (BioSemi's
# extra electrode ports) into a clean eye-movement signal. Two electrodes
# placed above and below the eye pick up nearly the same noise and the same
# stray brain activity; subtracting one from the other cancels that shared
# part and leaves mostly the eye's own up/down movement.
#
# Sits one step earlier in the pipeline than regression.R, the same way
# eeg_rereference() does: hand it a recording with raw EXG channels, get
# back an ordinary eeg object with one derived channel in their place. Name
# that channel with "EOG" in it (e.g. "VEOG") and fit_eog_regression() finds
# it automatically, the same lookup find_bads_eog() uses
# (.resolve_reference_channels(), R/ica_detect.R) - nothing in regression.R
# needs to know this channel was built by subtraction rather than recorded
# directly.
#
# Reference-invariant by construction: anode and cathode are both expressed
# relative to whatever reference is currently active, so subtracting them
# cancels that reference term out algebraically. Unlike fit_eog_regression(),
# this function does not need a particular reference scheme and can run
# before or after eeg_rereference() with the same result.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# Status: Built.
# Tested: see tests/testthat/test-bipolar.R
# ============================================================================

# ----------------------------------------------------------------------------
# set_bipolar_reference() - derive a new channel as anode minus cathode
# ----------------------------------------------------------------------------
#' Derive a Bipolar Channel (Anode Minus Cathode)
#'
#' Subtracts one raw channel from another and adds the result to the
#' recording as a new channel - the standard way to turn a pair of
#' single-ended electrodes (e.g. BioSemi's EXG inputs) into a clean
#' eye-movement signal: place one electrode above the eye (\code{anode}) and
#' one below it (\code{cathode}), and \code{anode - cathode} cancels the
#' noise and stray brain activity the two pick up in common, leaving mostly
#' the eye's own movement. A port of the idea behind MNE-Python's
#' \code{set_bipolar_reference()}.
#'
#' @param eeg An object of class \code{eeg} (see \code{\link{new_eeg}}) -
#'   continuous data only.
#' @param anode Character scalar: the channel name to keep the sign of (e.g.
#'   the electrode above the eye). Must exist in \code{eeg$channels}.
#' @param cathode Character scalar: the channel name to subtract (e.g. the
#'   electrode below the eye). Must exist in \code{eeg$channels} and differ
#'   from \code{anode}.
#' @param ch_name Character scalar: the name for the new derived channel
#'   (e.g. \code{"VEOG"}). Required - there is no default, so every call
#'   picks a deliberate, meaningful name. Include \code{"EOG"} in it if you
#'   want \code{\link{find_bads_eog}}/\code{\link{fit_eog_regression}} to
#'   find it automatically later (see Details). Must not collide with a
#'   channel name that will still exist once this call finishes.
#' @param drop Logical, default \code{TRUE}. Removes \code{anode} and
#'   \code{cathode} from the returned object once the new channel is built,
#'   since keeping the two raw electrodes around alongside their own
#'   combination is normally just clutter (and would confuse later
#'   automatic channel lookups). Set \code{FALSE} to keep all three.
#'
#' @return A new object of class \code{eeg}, the same shape as \code{eeg}
#'   plus one new channel named \code{ch_name} (and minus \code{anode}/
#'   \code{cathode} if \code{drop = TRUE}), with a note appended to
#'   \code{preprocessing_history}. Since R does not change arguments in
#'   place, reassign the result (\code{eeg <- set_bipolar_reference(eeg,
#'   ...)}); the input object itself is left untouched.
#'
#' @details
#' \strong{Channel type.} The new channel is recorded as \code{"external"} in
#' \code{eeg$channel_types}, the same category every other EXG/EOG/ECG
#' channel already has (set directly, not by re-running the name-based
#' classifier - see \code{\link{new_eeg}}).
#'
#' \strong{Reference.} \code{anode} and \code{cathode} are both expressed
#' relative to whatever reference is currently active, so the subtraction
#' cancels that reference term out algebraically. Unlike
#' \code{\link{fit_eog_regression}}, this function does not need - and does
#' not check for - a particular reference scheme; it can run before or after
#' \code{\link{eeg_rereference}} with the same result.
#'
#' \strong{Bad channels.} If \code{anode} or \code{cathode} is listed in
#' \code{eeg$bads}, a warning names which one and the new channel is added
#' to \code{eeg$bads} in the result - it was built from a known-noisy
#' electrode, so later steps that automatically skip bad channels (e.g.
#' \code{\link{fit_eog_regression}}'s auto-detection) skip this one too,
#' instead of silently trusting it. If a dropped source channel was itself
#' in \code{eeg$bads}, its now-nonexistent name is removed from the result.
#'
#' @examples
#' \dontrun{
#'   # electrodes above/below the left eye, on a BioSemi recording
#'   eeg <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
#'                                ch_name = "VEOG")
#'
#'   # fit_eog_regression() finds "VEOG" automatically, same as a recorded one
#'   model <- fit_eog_regression(eeg)
#'
#'   # keep the two raw electrodes instead of dropping them
#'   eeg <- set_bipolar_reference(eeg, anode = "EXG1", cathode = "EXG2",
#'                                ch_name = "VEOG", drop = FALSE)
#' }
#'
#' @seealso \code{\link{fit_eog_regression}}, \code{\link{find_bads_eog}},
#'   \code{\link{eeg_rereference}}
#'
#' @export
set_bipolar_reference <- function(eeg, anode, cathode, ch_name, drop = TRUE) {

  # ========== VALIDATE eeg ==========

  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  # ========== VALIDATE REQUIRED ARGUMENTS WERE SUPPLIED ==========

  missing_args <- c("anode", "cathode", "ch_name")[
    c(missing(anode), missing(cathode), missing(ch_name))
  ]
  if (length(missing_args) > 0) {
    stop("ERROR: missing required argument(s): ",
         paste(missing_args, collapse = ", "), " - there are no defaults.",
         call. = FALSE)
  }

  # ========== VALIDATE anode / cathode ==========

  if (!is.character(anode) || length(anode) != 1 || is.na(anode)) {
    stop("ERROR: 'anode' must be a single, non-missing character string ",
         "naming one channel.", call. = FALSE)
  }
  if (!is.character(cathode) || length(cathode) != 1 || is.na(cathode)) {
    stop("ERROR: 'cathode' must be a single, non-missing character string ",
         "naming one channel.", call. = FALSE)
  }
  if (identical(anode, cathode)) {
    stop("ERROR: 'anode' and 'cathode' are both '", anode, "' - a channel ",
         "cannot be derived from itself.", call. = FALSE)
  }

  anode_idx <- match(anode, eeg$channels)
  if (is.na(anode_idx)) {
    stop("ERROR: 'anode' channel '", anode, "' not found in eeg$channels.",
         call. = FALSE)
  }
  cathode_idx <- match(cathode, eeg$channels)
  if (is.na(cathode_idx)) {
    stop("ERROR: 'cathode' channel '", cathode, "' not found in ",
         "eeg$channels.", call. = FALSE)
  }

  # ========== VALIDATE ch_name ==========

  if (!is.character(ch_name) || length(ch_name) != 1 || is.na(ch_name)) {
    stop("ERROR: 'ch_name' must be a single, non-missing character string ",
         "naming the new channel (e.g. \"VEOG\").", call. = FALSE)
  }

  surviving <- if (isTRUE(drop)) {
    setdiff(eeg$channels, c(anode, cathode))
  } else {
    eeg$channels
  }
  if (ch_name %in% surviving) {
    stop("ERROR: 'ch_name' (\"", ch_name, "\") already names a channel ",
         "that would still exist after this call. Choose a different name.",
         call. = FALSE)
  }

  # ========== BUILD THE NEW CHANNEL ==========

  out <- eeg
  new_row <- eeg$data[anode_idx, ] - eeg$data[cathode_idx, ]

  out$data            <- rbind(out$data, new_row)
  out$channels        <- c(out$channels, ch_name)
  out$channel_types   <- c(out$channel_types, "external")
  rownames(out$data)  <- NULL     # rbind() would otherwise label the new row
                                   # "new_row" (deparse.level = 1) while every
                                   # other row stays "" - channels are looked
                                   # up positionally via eeg$channels anyway.

  # ========== OPTIONALLY DROP THE SOURCE CHANNELS ==========

  if (isTRUE(drop)) {
    drop_idx <- match(c(anode, cathode), out$channels)
    out$data          <- out$data[-drop_idx, , drop = FALSE]
    out$channels      <- out$channels[-drop_idx]
    out$channel_types <- out$channel_types[-drop_idx]
    out$bads          <- setdiff(out$bads, c(anode, cathode))
  }

  # ========== BAD-SOURCE HANDLING ==========

  bad_sources <- intersect(c(anode, cathode), eeg$bads)
  if (length(bad_sources) > 0) {
    warning("Source channel(s) marked bad in eeg$bads: ",
            paste(bad_sources, collapse = ", "), ". '", ch_name,
            "' was built from them and has been marked bad too.",
            call. = FALSE)
    out$bads <- union(out$bads, ch_name)
  }

  # ========== PREPROCESSING HISTORY ==========

  out$preprocessing_history <- c(
    out$preprocessing_history,
    list(paste0("Bipolar channel '", ch_name, "' derived as '", anode,
                "' minus '", cathode, "'",
                if (isTRUE(drop)) {
                  paste0(" (sources '", anode, "' and '", cathode,
                         "' dropped)")
                } else {
                  " (sources kept)"
                }))
  )

  out
}
