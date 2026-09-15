#' ============================================================================
#'                       Bad EEG Channel Interpolation
#' ============================================================================
#'
#' This module rebuilds the signal of channels listed in eeg_obj$bads (see
#' R/eeg_class.R and R/bad_channels.R) using spherical spline interpolation
#' (Perrin et al., 1989) - the "fixer" counterpart to find_bad_channels()'s
#' "detector". A bad channel's signal is reconstructed as a linear
#' combination of the good channels, weighted by a Legendre-polynomial
#' function of the angular distance between electrodes on the scalp sphere
#' (using positions from eeg_obj$montage, see R/montage.R) - not raw
#' Euclidean distance. Ported directly from MNE-Python's
#' mne/channels/interpolation.py (_calc_g / _make_interpolation_matrix /
#' _interpolate_bads_eeg).
#'
#' By default, successfully interpolated channels are cleared back out of
#' eeg_obj$bads, since they are no longer bad.
#'
#' Author: Christos Dalamarinis
#' Date: Sep - 2026
#' Status: Ready
#' ============================================================================
#'
#' Interpolate Bad EEG Channels
#'
#' Rebuilds the signal of every channel listed in \code{eeg_obj$bads} (by
#' hand or by \code{\link{find_bad_channels}}) as a linear combination of
#' the good channels, using spherical spline interpolation (Perrin et al.,
#' 1989) - the same method as MNE-Python's \code{raw.interpolate_bads()}.
#' By default, interpolated channels are cleared back out of \code{bads}
#' once fixed, since they are no longer bad.
#'
#' @param eeg_obj An object of class 'eeg', with a montage attached (see
#'   \code{\link{set_montage}}).
#'
#' @param exclude Character vector of channel names to leave out of
#'   interpolation entirely - neither fixed if bad, nor used as a donor for
#'   other channels. Default \code{NULL} (none).
#'
#' @param origin Numeric vector of length 3 (x, y, z, in the same units as
#'   \code{eeg_obj$montage$positions}) giving the center of the sphere
#'   channel positions are measured against. Default \code{c(0, 0, 0)},
#'   correct for any montage built with \code{\link{create_montage}} - see
#'   Details.
#'
#' @param reset_bads Logical. If \code{TRUE} (default), successfully
#'   interpolated channels are removed from \code{eeg_obj$bads} (they are
#'   no longer bad). If \code{FALSE}, their data is still fixed but they
#'   stay listed in \code{bads}.
#'
#' @return The input \code{eeg_obj} with interpolated channels' data
#'   replaced, \code{$bads} updated (see \code{reset_bads}), and a step
#'   appended to \code{$preprocessing_history}. If there were no eligible
#'   bad channels to interpolate, \code{eeg_obj} is returned unchanged
#'   (with a message).
#'
#' @details
#' Only channels classified as \code{"eeg"} in \code{eeg_obj$channel_types}
#' are ever touched - mirrors MNE's own EEG-only interpolation path
#' (\code{mne.channels.interpolation._interpolate_bads_eeg}); external
#' (EOG/ECG/EMG/GSR/etc.) and status channels are never interpolated even
#' if listed in \code{bads}. A bad channel with no position in
#' \code{eeg_obj$montage} cannot be interpolated and is skipped with a
#' warning naming it.
#'
#' \code{origin} is the center of the sphere channel positions are assumed
#' to lie on - needed because the interpolation weights are a function of
#' the \emph{angle} between two electrodes as seen from that center, not
#' their straight-line distance. MNE has to estimate this from real
#' (imperfect) head-digitization data (its \code{origin = "auto"}); this
#' package's electrode database (see \code{\link{get_electrode_database}})
#' instead gives every channel a fixed \code{radius} from a single center,
#' so any montage built with \code{\link{create_montage}} is already
#' exactly spherical around \code{c(0, 0, 0)} - nothing to estimate. Only
#' override this for a custom, non-template montage. A sanity-check warning
#' fires if the supplied channel positions turn out not to be close to
#' spherical around \code{origin}.
#'
#' @examples
#' \dontrun{
#'   eeg <- set_montage(eeg, create_montage())
#'   eeg <- find_bad_channels(eeg)
#'   eeg <- interpolate_bads(eeg)
#'   eeg$bads  # interpolated channels are gone
#'
#'   # Fix the data but keep the channels flagged
#'   eeg <- interpolate_bads(eeg, reset_bads = FALSE)
#' }
#'
#' @seealso \code{\link{find_bad_channels}}, \code{\link{create_montage}},
#'   \code{\link{set_montage}}, \code{\link{eeg_rereference}}
#'
#' @export
interpolate_bads <- function(eeg_obj,
                              exclude = NULL,
                              origin = c(0, 0, 0),
                              reset_bads = TRUE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!is.matrix(eeg_obj$data)) {
    stop("ERROR: eeg_obj$data must be a numeric matrix (channels x timepoints).",
         call. = FALSE)
  }
  if (is.null(eeg_obj$montage) || !inherits(eeg_obj$montage, "montage")) {
    stop("ERROR: No montage attached - interpolation needs channel scalp ",
         "positions. Attach one with set_montage().", call. = FALSE)
  }
  if (!is.numeric(origin) || length(origin) != 3) {
    stop("ERROR: 'origin' must be a numeric vector of length 3 (x, y, z).",
         call. = FALSE)
  }

  # ========== DETERMINE PICKS / BAD / GOOD CHANNELS ==========

  picks <- eeg_obj$channels[eeg_obj$channel_types == "eeg"]
  if (!is.null(exclude)) {
    picks <- setdiff(picks, exclude)
  }

  bads_idx <- picks %in% eeg_obj$bads

  if (length(picks) == 0 || !any(bads_idx)) {
    message("interpolate_bads(): no bad EEG channel(s) to interpolate - ",
            "returning eeg_obj unchanged.")
    return(eeg_obj)
  }

  has_position <- picks %in% eeg_obj$montage$channels

  bad_no_position <- picks[bads_idx & !has_position]
  if (length(bad_no_position) > 0) {
    warning("interpolate_bads(): ", length(bad_no_position),
            " bad channel(s) have no position in the attached montage and ",
            "cannot be interpolated: ", paste(bad_no_position, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  bad_channels  <- picks[bads_idx & has_position]
  good_channels <- picks[!bads_idx & has_position]

  if (length(bad_channels) == 0) {
    message("interpolate_bads(): no bad EEG channel(s) with a montage ",
            "position left to interpolate - returning eeg_obj unchanged.")
    return(eeg_obj)
  }
  if (length(good_channels) == 0) {
    stop("ERROR: No good EEG channels with a montage position are available ",
         "to interpolate from.", call. = FALSE)
  }

  # ========== CHANNEL POSITIONS, CENTERED ON origin ==========

  mont_pos <- eeg_obj$montage$positions
  pos_good <- as.matrix(mont_pos[match(good_channels, mont_pos$channel),
                                  c("x", "y", "z")])
  pos_bad  <- as.matrix(mont_pos[match(bad_channels, mont_pos$channel),
                                  c("x", "y", "z")])

  origin   <- as.numeric(origin)
  pos_good <- sweep(pos_good, 2, origin, "-")
  pos_bad  <- sweep(pos_bad, 2, origin, "-")

  # Sanity check mirroring MNE's own _interpolate_bads_eeg: if the channel
  # positions are not (roughly) equidistant from 'origin', they are not
  # close to spherical and the spline fit below will be poor.
  all_dist   <- sqrt(rowSums(rbind(pos_good, pos_bad)^2))
  dist_ratio <- mean(all_dist / mean(all_dist))
  if (abs(1 - dist_ratio) > 0.1) {
    warning("interpolate_bads(): channel positions are not close to ",
            "spherical around 'origin' - interpolation results may be ",
            "inaccurate.", call. = FALSE, immediate. = TRUE)
  }

  # ========== BUILD INTERPOLATION MATRIX AND APPLY ==========

  interp_matrix <- make_interpolation_matrix(pos_good, pos_bad)

  good_idx <- match(good_channels, eeg_obj$channels)
  bad_idx  <- match(bad_channels, eeg_obj$channels)

  eeg_obj$data[bad_idx, ] <- interp_matrix %*% eeg_obj$data[good_idx, , drop = FALSE]

  # ========== RESET bads / HISTORY ==========

  if (isTRUE(reset_bads)) {
    eeg_obj$bads <- setdiff(eeg_obj$bads, bad_channels)
  }

  history_entry <- paste0(
    "interpolate_bads(): interpolated ", length(bad_channels),
    " channel(s) (", paste(bad_channels, collapse = ", "), ") from ",
    length(good_channels), " good channel(s)",
    if (isTRUE(reset_bads)) {
      " - cleared from bads."
    } else {
      " - kept in bads (reset_bads = FALSE)."
    })

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history,
                                      list(history_entry))

  # ========== RETURN ==========

  eeg_obj
}

#' Build a Spherical-Spline Interpolation Matrix
#'
#' Internal helper computing the linear map from a set of "good" channel
#' positions onto a set of target positions (e.g. bad channels), using
#' regularized spherical spline interpolation (Perrin et al., 1989). Direct
#' port of MNE-Python's
#' \code{mne.channels.interpolation._make_interpolation_matrix}.
#'
#' @param pos_from Numeric matrix, \code{n_from x 3} (x/y/z), positions to
#'   interpolate from. Should already be centered on the sphere's origin
#'   (see \code{\link{interpolate_bads}}) - normalized onto the unit sphere
#'   internally.
#' @param pos_to Numeric matrix, \code{n_to x 3}, positions to interpolate
#'   to. Same centering requirement as \code{pos_from}.
#' @param stiffness Numeric, spline stiffness (\code{m} in Perrin et al.),
#'   passed to \code{\link{calc_g}}. Default 4, matching MNE.
#' @param n_legendre_terms Integer, number of Legendre polynomial terms,
#'   passed to \code{\link{calc_g}}. Default 50, matching MNE.
#' @param alpha Numeric, regularization added to the diagonal of the
#'   good-good kernel matrix before inversion. Default \code{1e-5},
#'   matching MNE.
#'
#' @return Numeric matrix, \code{n_to x n_from} - the weights mapping good
#'   channel signals onto the \code{pos_to} locations.
#'
#' @seealso \code{\link{calc_g}}, \code{\link{interpolate_bads}}
#' @keywords internal
make_interpolation_matrix <- function(pos_from, pos_to,
                                       stiffness = 4,
                                       n_legendre_terms = 50,
                                       alpha = 1e-5) {

  n_from <- nrow(pos_from)
  n_to   <- nrow(pos_to)

  # Position dimnames (e.g. channel names inherited from montage$positions)
  # carry no meaning for this purely positional matrix algebra - drop them
  # up front so they can't ride along into G_from/G_to_from/C_inv below.
  dimnames(pos_from) <- NULL
  dimnames(pos_to)   <- NULL

  # Normalize each row (channel position) onto the unit sphere.
  pos_from <- sweep(pos_from, 1, sqrt(rowSums(pos_from^2)), "/")
  pos_to   <- sweep(pos_to, 1, sqrt(rowSums(pos_to^2)), "/")

  cosang_from    <- pos_from %*% t(pos_from)
  cosang_to_from <- pos_to %*% t(pos_from)

  G_from    <- calc_g(cosang_from, stiffness = stiffness,
                       n_legendre_terms = n_legendre_terms)
  G_to_from <- calc_g(cosang_to_from, stiffness = stiffness,
                       n_legendre_terms = n_legendre_terms)

  diag(G_from) <- diag(G_from) + alpha

  # Augmented system (Perrin et al., 1989): the extra row/column encodes
  # the constraint sum(spline coefficients) = 0 needed for a unique fit.
  # A provable consequence is that each row of the returned matrix sums to
  # 1, so a spatially constant input field is reproduced exactly.
  C <- rbind(
    cbind(G_from, rep(1, n_from)),
    c(rep(1, n_from), 0)
  )
  C_inv <- MASS::ginv(C)

  cbind(G_to_from, rep(1, n_to)) %*% C_inv[, -ncol(C_inv), drop = FALSE]
}

#' Spherical Spline Green's Function
#'
#' Internal helper evaluating the spherical-spline "g" function (Perrin et
#' al., 1989) at cosine-of-angle values between pairs of points on a
#' sphere - the Legendre-series kernel \code{\link{make_interpolation_matrix}}
#' is built from. Direct port of MNE-Python's
#' \code{mne.channels.interpolation._calc_g}.
#'
#' @param cosang Numeric vector or matrix of cosine-of-angle values (i.e.
#'   the dot product of unit vectors) between pairs of points on a sphere.
#' @param stiffness Numeric, spline stiffness (\code{m} in Perrin et al.).
#'   Default 4, matching MNE.
#' @param n_legendre_terms Integer, number of Legendre polynomial terms to
#'   sum. Default 50, matching MNE.
#'
#' @return A numeric vector or matrix the same shape as \code{cosang}.
#'
#' @seealso \code{\link{make_interpolation_matrix}}, \code{\link{interpolate_bads}}
#' @keywords internal
calc_g <- function(cosang, stiffness = 4, n_legendre_terms = 50) {
  n_seq   <- seq_len(n_legendre_terms)
  factors <- (2 * n_seq + 1) /
    (n_seq^stiffness * (n_seq + 1)^stiffness * 4 * pi)
  .legendre_series_eval(cosang, c(0, factors))
}

#' Evaluate a Legendre Series (internal)
#'
#' Evaluates \code{sum(coeffs[n+1] * P_n(x))} for \code{n = 0..length(coeffs)-1},
#' where \code{P_n} is the Legendre polynomial of degree \code{n}, via the
#' standard 3-term recurrence (\code{P_0 = 1}, \code{P_1 = x},
#' \code{n*P_n = (2n-1)*x*P_(n-1) - (n-1)*P_(n-2)}). Equivalent to numpy's
#' \code{numpy.polynomial.legendre.legval(x, coeffs)}, used by
#' \code{\link{calc_g}}. Works elementwise whether \code{x} is a plain
#' vector or a matrix - shape is preserved throughout.
#'
#' @param x Numeric vector or matrix, points to evaluate at (typically in
#'   \code{[-1, 1]}, cosine-of-angle values).
#' @param coeffs Numeric vector of coefficients, \code{coeffs[1]} is the
#'   \code{P_0} term, \code{coeffs[2]} the \code{P_1} term, etc.
#'
#' @return Numeric vector or matrix, same shape as \code{x}.
#' @keywords internal
.legendre_series_eval <- function(x, coeffs) {

  n_terms <- length(coeffs)

  ones  <- x
  ones[] <- 1
  total <- coeffs[1] * ones

  if (n_terms >= 2) {

    p_prev <- ones
    p_curr <- x
    total  <- total + coeffs[2] * p_curr

    if (n_terms >= 3) {
      for (n in 2:(n_terms - 1)) {
        p_next <- ((2 * n - 1) * x * p_curr - (n - 1) * p_prev) / n
        total  <- total + coeffs[n + 1] * p_next
        p_prev <- p_curr
        p_curr <- p_next
      }
    }
  }

  total
}
