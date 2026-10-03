# ============================================================================
#           Current Source Density (Surface Laplacian) Transform
# ============================================================================
#
# This module replaces the scalp voltage of every EEG channel with its
# current source density (CSD): the spherical-spline surface Laplacian of the
# scalp potential (Perrin et al., 1987, 1989; Kayser & Tenke, 2015). CSD is
# reference-free - a constant added to every channel, which is all a change of
# reference does, disappears - and it sharpens broad, smeared-out voltage maps
# into more focal ones. It is built from three steps: the weight table
# (.csd_matrix), its application to the data (.csd_apply_matrix) and the
# user-facing wrapper (compute_csd); it reuses the spline kernels calc_g() and
# calc_h() from R/interpolate.R.
#
# The transform is linear: one fixed channels-by-channels weight table that
# depends only on the electrode positions (from a montage, see R/montage.R),
# applied unchanged to every time point and every epoch.
#
# Author: Christos Dalamarinis
# Date: Oct - 2026
# ============================================================================

#' Compute Current Source Density (Surface Laplacian)
#'
#' Replaces the scalp voltage of every EEG channel with its current source
#' density (CSD): the spherical-spline surface Laplacian of the scalp
#' potential (Perrin et al., 1987, 1989; Kayser & Tenke, 2015). CSD is
#' reference-free and sharpens broad, smeared-out voltage maps into more focal
#' ones, so it is a standard companion to ERP and time-frequency topographies.
#'
#' @param x An object of class \code{'eeg'} (continuous data) or
#'   \code{'eeg_epochs'} (epoched data, see \code{\link{epoch_eeg}}). Only
#'   channels classified as \code{"eeg"} in \code{x$channel_types} are
#'   transformed; EOG, status and other channels are left exactly as they are.
#'
#' @param montage An object of class \code{'montage'} giving the scalp position
#'   of every EEG channel (see \code{\link{create_montage}}). Default
#'   \code{NULL}, which uses \code{x$montage}. Epoched objects do not carry a
#'   montage, so pass one explicitly for those (for example
#'   \code{montage = eeg$montage}).
#'
#' @param head_radius Numeric, radius of the head sphere in the same units as
#'   the montage positions (millimetres for \code{\link{create_montage}}). It
#'   only scales the result, by 1 / radius squared, and so sets the unit of the
#'   output. Default \code{NULL}: the mean distance of the EEG electrodes from
#'   \code{origin} (about 87.5 mm for the biosemi64 template).
#'
#' @param origin Numeric vector of length 3 (x, y, z, in the same units as the
#'   montage positions) giving the centre of the sphere the electrodes lie on.
#'   Default \code{c(0, 0, 0)}, correct for any montage built with
#'   \code{\link{create_montage}} - see \code{\link{interpolate_bads}}.
#'
#' @param lambda2 Numeric in [0, 1), regularisation of the spline fit. Larger
#'   values give a smoother result. Default \code{1e-5}.
#'
#' @param stiffness Numeric, non-negative, stiffness of the spline (\code{m} in
#'   Perrin et al.). Default 4.
#'
#' @param n_legendre_terms Integer, number of Legendre terms to sum. Default 50.
#'
#' @param verbose Logical. If \code{TRUE} (default), print a one-line summary.
#'
#' @return An object of the same class as \code{x}, with the EEG channels'
#'   data replaced by their current source density, \code{$reference} set to
#'   \code{"CSD"} (and \code{$metadata$reference_scheme} too), and a step
#'   appended to \code{$preprocessing_history}. Everything else (other
#'   channels, \code{bads}, \code{montage}, \code{annotations},
#'   \code{channel_types}, events, times) is unchanged. The input object itself
#'   is left untouched.
#'
#' @details
#' \strong{What it does.} Each EEG channel's voltage is replaced by how sharply
#' the scalp potential is peaked at that electrode compared with its
#' surroundings (the second spatial derivative over the scalp surface). Broad,
#' smooth patterns come out near zero; sharp local peaks come out large.
#' Because a constant added to every channel - which is all a change of
#' reference does - has no peaks at all, CSD is reference-free.
#'
#' \strong{How.} The head is treated as a sphere and each electrode is placed
#' on it from its montage position. A spherical spline (Perrin et al., 1989) is
#' fitted through the voltages and its surface Laplacian is evaluated at every
#' electrode. All of this is linear, so it collapses into one fixed weight
#' table (channels x channels) that depends only on the electrode positions
#' and is applied unchanged to every time point and every epoch.
#'
#' \strong{Unit.} The output is in microvolts per square metre when the input
#' is in microvolts (the unit the package assumes, e.g. for BioSemi data read
#' with \code{\link{read_bdf_native}}): the table is divided by the square of
#' the head radius, converted to metres. Divide by 10,000 for microvolts per
#' square centimetre. A local voltage maximum gives a positive value.
#'
#' \strong{Requirements.} Every EEG channel needs a position in the montage
#' (positions are assumed to be in millimetres, as in the template). There must
#' be no bad EEG channels: a bad channel would corrupt the values of its
#' neighbours, so this function refuses to run on them. Use
#' \code{\link{interpolate_bads}} first. The EEG data must be finite.
#'
#' \strong{Order of steps.} Do CSD last: after filtering, bad-channel handling,
#' ICA and re-referencing. Afterwards the data are no longer in microvolts, so
#' amplitude thresholds (for example the \code{reject_threshold} of
#' \code{\link{epoch_eeg}}) no longer apply, and re-referencing the data makes
#' no sense. CSD cannot be applied twice (\code{x$reference} is
#' \code{"CSD"} afterwards).
#'
#' \strong{Epochs and averages.} CSD is linear, so applying it to the epochs
#' and then averaging gives the same result as averaging first and then
#' applying CSD (for mean averaging). \code{eeg_evoked} objects are not
#' supported directly because they carry neither the montage nor the channel
#' types: apply CSD to the epochs first.
#'
#' \strong{Noise.} CSD emphasises local detail, including spatial noise, so a
#' noisy channel stands out more than in the voltage map. Spatial resolution
#' also depends on electrode density: a dense montage (such as 64 channels) is
#' what the method is meant for.
#'
#' @examples
#' \dontrun{
#'   # Continuous data
#'   eeg <- set_montage(eeg, create_montage())
#'   eeg <- find_bad_channels(eeg)
#'   eeg <- interpolate_bads(eeg)          # CSD refuses bad channels
#'   eeg_csd <- compute_csd(eeg)
#'
#'   # Epoched data: CSD first, then average (same result as the other order)
#'   epochs     <- epoch_eeg(eeg, tmin = -0.2, tmax = 0.8)
#'   epochs_csd <- compute_csd(epochs, montage = eeg$montage)
#'   erp_csd    <- average_epochs(epochs_csd)
#' }
#'
#' @seealso \code{\link{interpolate_bads}}, \code{\link{create_montage}},
#'   \code{\link{set_montage}}, \code{\link{epoch_eeg}},
#'   \code{\link{average_epochs}}, \code{\link{plot_topography}}
#'
#' @export
compute_csd <- function(x,
                        montage = NULL,
                        head_radius = NULL,
                        origin = c(0, 0, 0),
                        lambda2 = 1e-5,
                        stiffness = 4,
                        n_legendre_terms = 50,
                        verbose = TRUE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(x, "eeg") && !inherits(x, "eeg_epochs")) {
    stop("ERROR: 'x' must be an object of class 'eeg' or 'eeg_epochs'.",
         call. = FALSE)
  }

  is_epochs <- inherits(x, "eeg_epochs")

  if (is_epochs) {
    if (is.null(x$data)) {
      stop("ERROR: epoch data not loaded - re-run epoch_eeg() with ",
           "preload = TRUE.", call. = FALSE)
    }
    if (length(dim(x$data)) != 3) {
      stop("ERROR: x$data must be a 3D array (channels x timepoints x ",
           "epochs).", call. = FALSE)
    }
  } else {
    if (!is.matrix(x$data)) {
      stop("ERROR: x$data must be a numeric matrix (channels x timepoints).",
           call. = FALSE)
    }
  }

  if (!is.numeric(lambda2) || length(lambda2) != 1 || is.na(lambda2) ||
      lambda2 < 0 || lambda2 >= 1) {
    stop("ERROR: 'lambda2' must be a single number from 0 up to (but not ",
         "including) 1.", call. = FALSE)
  }
  if (!is.numeric(stiffness) || length(stiffness) != 1 ||
      is.na(stiffness) || stiffness < 0) {
    stop("ERROR: 'stiffness' must be a single non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(n_legendre_terms) || length(n_legendre_terms) != 1 ||
      is.na(n_legendre_terms) || n_legendre_terms < 1 ||
      n_legendre_terms != round(n_legendre_terms)) {
    stop("ERROR: 'n_legendre_terms' must be a single whole number of at ",
         "least 1.", call. = FALSE)
  }
  if (!is.numeric(origin) || length(origin) != 3 || any(!is.finite(origin))) {
    stop("ERROR: 'origin' must be a numeric vector of length 3 (x, y, z).",
         call. = FALSE)
  }
  if (!is.null(head_radius) &&
      (!is.numeric(head_radius) || length(head_radius) != 1 ||
       !is.finite(head_radius) || head_radius <= 0)) {
    stop("ERROR: 'head_radius' must be NULL or a single positive number.",
         call. = FALSE)
  }
  if (!isTRUE(verbose) && !isFALSE(verbose)) {
    stop("ERROR: 'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(montage)) montage <- x$montage
  if (is.null(montage) || !inherits(montage, "montage")) {
    stop("ERROR: No montage available - CSD needs channel scalp positions. ",
         "Attach one with set_montage(), or pass one explicitly via the ",
         "'montage' argument.", call. = FALSE)
  }

  if (identical(x$reference, "CSD")) {
    stop("ERROR: CSD has already been applied to this object (x$reference ",
         "is 'CSD') - it must not be applied twice.", call. = FALSE)
  }

  # ========== PICK EEG CHANNELS ==========

  picks <- which(x$channel_types == "eeg")
  if (length(picks) == 0) {
    stop("ERROR: No EEG channels found (x$channel_types has no \"eeg\" ",
         "entries).", call. = FALSE)
  }
  pick_names <- x$channels[picks]

  bad_eeg <- intersect(pick_names, x$bads)
  if (length(bad_eeg) > 0) {
    stop("ERROR: CSD cannot be computed with bad EEG channels: ",
         paste(bad_eeg, collapse = ", "),
         ". Interpolate them first with interpolate_bads().", call. = FALSE)
  }

  mont_pos <- montage$positions
  rows     <- match(pick_names, mont_pos$channel)
  if (anyNA(rows)) {
    stop("ERROR: EEG channel(s) with no position in the montage: ",
         paste(pick_names[is.na(rows)], collapse = ", "),
         ". CSD needs a position for every EEG channel.", call. = FALSE)
  }

  eeg_data <- if (is_epochs) x$data[picks, , ] else x$data[picks, ]
  if (any(!is.finite(eeg_data))) {
    stop("ERROR: EEG data contains NA or non-finite values. Clean the data ",
         "before computing CSD.", call. = FALSE)
  }

  # ========== CHANNEL POSITIONS, CENTERED ON origin ==========

  pos  <- as.matrix(mont_pos[rows, c("x", "y", "z")])
  pos  <- sweep(pos, 2, as.numeric(origin), "-")
  dist <- sqrt(rowSums(pos^2))

  bad_pos <- !is.finite(dist) | dist == 0
  if (any(bad_pos)) {
    stop("ERROR: Zero or non-finite channel position found for: ",
         paste(pick_names[bad_pos], collapse = ", "), ".", call. = FALSE)
  }

  unit_pos <- round(sweep(pos, 1, dist, "/"), 6)
  dup <- duplicated(unit_pos) | duplicated(unit_pos, fromLast = TRUE)
  if (any(dup)) {
    stop("ERROR: These EEG channels share the same position, which CSD ",
         "cannot handle: ", paste(pick_names[dup], collapse = ", "), ".",
         call. = FALSE)
  }

  # Sanity check: electrodes should lie (roughly) on a sphere around 'origin'
  # (any electrode more than 10 % from the mean distance triggers it).
  if (max(abs(dist / mean(dist) - 1)) > 0.1) {
    warning("compute_csd(): channel positions are not close to spherical ",
            "around 'origin' - CSD values may be inaccurate.",
            call. = FALSE, immediate. = TRUE)
  }

  # ========== HEAD RADIUS ==========

  radius_source <- if (is.null(head_radius)) "estimated from the montage" else
    "user supplied"
  if (is.null(head_radius)) head_radius <- mean(dist)
  radius_m <- head_radius / 1000

  # ========== BUILD WEIGHT TABLE AND APPLY ==========

  X   <- .csd_matrix(pos, lambda2, stiffness, n_legendre_terms, radius_m)
  out <- .csd_apply_matrix(x, picks, X)

  # ========== MARK AS CSD ==========

  out$reference <- "CSD"
  if (!is.null(out$metadata)) out$metadata$reference_scheme <- "CSD"

  history_entry <- paste0(
    "compute_csd(): current source density (spherical-spline surface Laplacian) ",
    "applied to ", length(picks), " EEG channel(s); head radius ",
    format(round(head_radius, 2), nsmall = 2), " mm (", radius_source,
    "), lambda2 = ", format(lambda2), ", stiffness = ", format(stiffness),
    ", n_legendre_terms = ", n_legendre_terms,
    "; output unit is uV/m^2 for input in uV.")
  out$preprocessing_history <- c(out$preprocessing_history,
                                 list(history_entry))

  if (verbose) {
    message("compute_csd(): CSD applied to ", length(picks),
            " EEG channel(s); head radius ",
            format(round(head_radius, 2), nsmall = 2), " mm (", radius_source,
            "). Output unit: uV/m^2 for input in uV.")
  }

  out
}

#' Build the CSD Weight Table (internal)
#'
#' Internal helper computing the square matrix that turns the voltages of a
#' set of EEG channels into their current source density: row i holds the
#' weights that make up channel i's CSD value. The G kernel
#' (\code{\link{calc_g}}) is regularised and inverted to fit the spline, and
#' the H kernel (\code{\link{calc_h}}) reads off its curvature at each
#' electrode.
#'
#' @param pos Numeric matrix, \code{n x 3} (x/y/z), one row per channel,
#'   already centred on the sphere's origin (see \code{\link{compute_csd}}) -
#'   normalised onto the unit sphere internally.
#' @param lambda2 Numeric, regularisation added to the diagonal of the G
#'   kernel matrix before inversion. Default \code{1e-5}.
#' @param stiffness Numeric, spline stiffness, passed to the kernels. Default
#'   4.
#' @param n_legendre_terms Integer, number of Legendre terms, passed to the
#'   kernels. Default 50.
#' @param radius Numeric, head radius in metres. The whole table is divided
#'   by its square. Default 1.
#'
#' @return Numeric matrix, \code{n x n}. Every row sums to (about) zero, which
#'   is why a constant added to every channel has no effect.
#'
#' @seealso \code{\link{compute_csd}}, \code{\link{calc_g}},
#'   \code{\link{calc_h}}
#' @keywords internal
.csd_matrix <- function(pos, lambda2 = 1e-5, stiffness = 4,
                        n_legendre_terms = 50, radius = 1) {

  dimnames(pos) <- NULL

  # Positions onto the unit sphere, then cosines of the angles between all pairs
  pos    <- sweep(pos, 1, sqrt(rowSums(pos^2)), "/")
  cosang <- pmin(pmax(pos %*% t(pos), -1), 1)

  G <- calc_g(cosang, stiffness = stiffness, n_legendre_terms = n_legendre_terms)
  H <- calc_h(cosang, stiffness = stiffness, n_legendre_terms = n_legendre_terms)
  diag(G) <- diag(G) + lambda2

  Gi <- tryCatch(solve(G), error = function(e) {
    stop("ERROR: the electrode kernel matrix could not be inverted - check for ",
         "two electrodes at the same position, or raise lambda2.", call. = FALSE)
  })

  total_col <- colSums(Gi)
  total_all <- sum(total_col)
  n         <- nrow(H)
  Cp2       <- Gi %*% (diag(n) - 1 / n)
  c02       <- colSums(Cp2) / total_all
  C2        <- Cp2 - total_col %o% c02

  (H %*% C2) / radius^2
}

#' Apply the CSD Weight Table to an eeg or eeg_epochs Object (internal)
#'
#' Multiplies the EEG channels of \code{x} by the weight table \code{X}.
#' For an \code{eeg} object that is one matrix product on the channels x time
#' matrix. For an \code{eeg_epochs} object all epochs go through in a single
#' product (the channels x time x epochs array is laid out as a channels x
#' (time * epochs) matrix and folded back), which gives exactly the same
#' result as applying the table epoch by epoch. Channels outside
#' \code{picks} are not touched.
#'
#' @param x An object of class \code{'eeg'} or \code{'eeg_epochs'}.
#' @param picks Integer vector, row indices of the channels to transform (in
#'   the same order as the rows and columns of \code{X}).
#' @param X Numeric matrix, \code{length(picks) x length(picks)}, from
#'   \code{.csd_matrix}.
#'
#' @return \code{x} with \code{x$data[picks, ]} replaced.
#'
#' @seealso \code{\link{compute_csd}}
#' @keywords internal
.csd_apply_matrix <- function(x, picks, X) {

  if (inherits(x, "eeg_epochs")) {
    n_time   <- dim(x$data)[2]
    n_epochs <- dim(x$data)[3]
    stacked  <- matrix(x$data[picks, , , drop = FALSE], nrow = length(picks))
    x$data[picks, , ] <- array(X %*% stacked,
                               dim = c(length(picks), n_time, n_epochs))
  } else {
    x$data[picks, ] <- X %*% x$data[picks, , drop = FALSE]
  }

  x
}
