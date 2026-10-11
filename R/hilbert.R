# ============================================================================
#              Hilbert Transform: Envelope and Phase of Each Channel
# ============================================================================
#
# The Hilbert transform turns a signal into its "analytic signal": for every
# sample, the voltage together with a partner signal that is the same wave
# shifted by a quarter of a cycle. Voltage and partner are like the two
# coordinates of a point going round a circle: the distance of the point from
# the centre is the instantaneous amplitude of the wave (its envelope) and its
# angle is the instantaneous phase. Envelope and phase are the basis of
# band-limited power over time, of phase-amplitude coupling and of phase-based
# connectivity.
#
# apply_hilbert() replaces the chosen channels of an 'eeg' or 'eeg_epochs'
# object by their envelope (in the unit of the data) or their phase (in
# radians). The analytic signal is computed through the FFT (Marple, 1999): the
# signal is padded with zeros to a fast FFT length, the negative frequencies
# are dropped and the positive ones doubled, and the result is cut back to the
# original length. The complex analytic signal itself is never stored in an
# 'eeg' object: complex numbers would break other functions of the package or
# silently lose their imaginary part, and the envelope and the phase cover what
# the analytic signal is used for. .analytic_signal() is also the core of the
# older helpers .hilbert_envelope() (R/annotations.R) and .hilbert_phase()
# (R/ica_detect.R).
#
# Author: Christos Dalamarinis
# Date: Oct - 2026
# ============================================================================

#' Apply the Hilbert Transform: Envelope or Phase of Each Channel
#'
#' Replaces the chosen channels by their instantaneous amplitude (the
#' \emph{envelope}: how big a rhythm is at each moment) or by their
#' instantaneous \emph{phase} (where in its up-and-down cycle the rhythm is at
#' each moment), both taken from the analytic signal that the Hilbert transform
#' builds. Filter the data to the band of interest first (see
#' \code{\link{eeg_bandpass}}): the envelope and the phase of broadband EEG are
#' not the power and the phase of a rhythm. The squared envelope is the power of
#' the rhythm over time; the phase is what phase-amplitude coupling and
#' phase-based connectivity measures are built on.
#'
#' @param eeg_obj An object of class \code{'eeg'} (continuous data) or
#'   \code{'eeg_epochs'} (epoched data, see \code{\link{epoch_eeg}}).
#'
#' @param output \code{"envelope"} (default) or \code{"phase"}. See Details.
#'
#' @param picks \code{NULL} (default), a character vector of channel names or a
#'   numeric vector of channel numbers: the channels to transform. \code{NULL}
#'   means every channel typed \code{"eeg"} in \code{eeg_obj$channel_types}, bad
#'   channels included (every channel is transformed on its own, so a bad
#'   channel does not affect the others). Channels that are not picked are left
#'   exactly as they are.
#'
#' @param n_fft \code{"auto"} (default), \code{NULL} or a single whole number:
#'   the length the signal is zero-padded to before the FFT (the padding is cut
#'   off afterwards). \code{"auto"} uses the next length whose only prime factors
#'   are 2, 3 and 5, which the FFT handles fastest (\code{\link[stats]{nextn}});
#'   \code{NULL} does not pad; a number must be at least the number of time
#'   points. Padding changes the result a little, mostly near the ends of the
#'   data.
#'
#' @param verbose Logical. If \code{TRUE} (default), print a one-line summary.
#'
#' @return An object of the same class as \code{eeg_obj}, with the picked
#'   channels replaced by their envelope or phase and a step appended to
#'   \code{$preprocessing_history}. Everything else (other channels,
#'   \code{bads}, \code{montage}, \code{annotations}, \code{channel_types},
#'   events, times, \code{reference}) is unchanged. The input object itself is
#'   left untouched.
#'
#' @details
#' \strong{Envelope.} The magnitude of the analytic signal, in the unit of the
#' data (microvolts for data read with \code{\link{read_bdf_native}}). It is
#' never negative. For a rhythm that swells and shrinks, such as alpha, it is the
#' curve that follows the peaks of the wave; for a steady sine wave of amplitude
#' 5 it is 5 at every sample.
#'
#' \strong{Phase.} The angle of the analytic signal, in radians, from just above
#' \eqn{-\pi} to \eqn{\pi}. It is 0 at a positive peak of the wave,
#' \eqn{\pi/2} at a falling zero crossing, \eqn{\pi} (or \eqn{-\pi}) at a trough and
#' \eqn{-\pi/2} at a rising zero crossing. Phases are circular: to compare two
#' of them, take the angle of the difference of their complex exponentials, not
#' the plain difference (\eqn{-\pi + 0.1} and \eqn{\pi - 0.1} are close
#' together). Do not baseline-correct or average phases like voltages.
#'
#' \strong{How.} The analytic signal is \eqn{x + iy}, where \eqn{y} is the
#' Hilbert transform of \eqn{x}: the same wave shifted by a quarter of a cycle.
#' It is computed through the FFT: the negative frequencies are dropped, the
#' positive ones doubled, and the result is transformed back (Marple, 1999).
#'
#' \strong{Edges.} The transform treats the signal as if it repeated around its
#' ends, so the first and last second or so are less reliable (zero padding
#' moves the problem rather than removing it). Apply it to the \emph{continuous}
#' recording, after filtering, and epoch the result: the envelope and phase
#' objects can be passed to \code{\link{epoch_eeg}} (for a phase object use
#' \code{baseline = NULL}, and no amplitude thresholds). Applying it to short
#' epochs gives edge effects in every epoch.
#'
#' \strong{No complex values.} The complex analytic signal itself is not
#' returned: an \code{'eeg'} object holds real numbers only. Envelope and phase
#' are what band power, phase-amplitude coupling and phase-based connectivity
#' need.
#'
#' Every picked channel must be finite; locate gaps with
#' \code{\link{annotate_nan}} and deal with them first. Applying the function
#' twice to the same object gives a warning, because the envelope or the phase
#' of an envelope or a phase rarely means anything.
#'
#' @references Marple, S. L. (1999). Computing the discrete-time "analytic"
#'   signal via FFT. \emph{IEEE Transactions on Signal Processing}, 47(9),
#'   2600-2603.
#'
#' @examples
#' \dontrun{
#'   # Alpha power over time
#'   alpha     <- eeg_bandpass(eeg, l_freq = 8, h_freq = 12)
#'   alpha_env <- apply_hilbert(alpha, output = "envelope")
#'   alpha_pow <- alpha_env$data^2
#'
#'   # Phase of a slow rhythm, cut into epochs (no baseline for a phase)
#'   theta_phase <- apply_hilbert(eeg_bandpass(eeg, l_freq = 4, h_freq = 8),
#'                                output = "phase")
#'   epochs <- epoch_eeg(theta_phase, events = "all", tmin = -0.2, tmax = 0.8,
#'                       baseline = NULL, baseline_method = "none")
#' }
#'
#' @seealso \code{\link{eeg_bandpass}}, \code{\link{epoch_eeg}},
#'   \code{\link{eeg_band_power}}, \code{\link{annotate_nan}}
#'
#' @export
apply_hilbert <- function(eeg_obj,
                          output  = c("envelope", "phase"),
                          picks   = NULL,
                          n_fft   = "auto",
                          verbose = TRUE) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg") && !inherits(eeg_obj, "eeg_epochs")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg' or 'eeg_epochs'.",
         call. = FALSE)
  }
  is_epochs <- inherits(eeg_obj, "eeg_epochs")

  if (!is_epochs) {
    if (!is.matrix(eeg_obj$data)) {
      stop("ERROR: eeg_obj$data must be a numeric matrix (channels x timepoints).",
           call. = FALSE)
    }
  } else {
    if (is.null(eeg_obj$data)) {
      stop("ERROR: epoch data not loaded - re-run epoch_eeg() with preload = TRUE.",
           call. = FALSE)
    }
    if (length(dim(eeg_obj$data)) != 3) {
      stop("ERROR: eeg_obj$data must be a 3D array (channels x timepoints x epochs).",
           call. = FALSE)
    }
  }

  if (identical(output, c("envelope", "phase"))) output <- "envelope"
  if (!is.character(output) || length(output) != 1 ||
      !(output %in% c("envelope", "phase"))) {
    stop("ERROR: 'output' must be \"envelope\" or \"phase\".", call. = FALSE)
  }

  n_times <- if (is_epochs) dim(eeg_obj$data)[2] else ncol(eeg_obj$data)
  if (is.null(n_fft)) {
    n_fft <- n_times
  } else if (identical(n_fft, "auto")) {
    n_fft <- stats::nextn(n_times)
  } else if (is.numeric(n_fft) && length(n_fft) == 1 && is.finite(n_fft) &&
             n_fft == round(n_fft) && n_fft >= 1) {
    n_fft <- as.integer(n_fft)
  } else {
    stop("ERROR: 'n_fft' must be \"auto\", NULL or a single whole number.",
         call. = FALSE)
  }
  if (n_fft < n_times) {
    stop("ERROR: 'n_fft' (", n_fft, ") must be at least the number of time points (",
         n_times, ").", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("ERROR: 'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  # ========== PICK THE CHANNELS ==========

  n_ch <- length(eeg_obj$channels)
  if (is.null(picks)) {
    idx <- which(eeg_obj$channel_types == "eeg")
  } else if (is.character(picks)) {
    unknown <- setdiff(picks, eeg_obj$channels)
    if (length(unknown) > 0) {
      stop("ERROR: These channels in 'picks' are not in eeg_obj: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    idx <- sort(unique(match(picks, eeg_obj$channels)))
  } else if (is.numeric(picks) && all(is.finite(picks)) && all(picks == round(picks))) {
    if (any(picks < 1 | picks > n_ch)) {
      stop("ERROR: 'picks' numbers must be between 1 and ", n_ch, ".", call. = FALSE)
    }
    idx <- sort(unique(as.integer(picks)))
  } else {
    stop("ERROR: 'picks' must be NULL, channel names or channel numbers.",
         call. = FALSE)
  }

  if (length(idx) == 0) {
    stop("ERROR: No channels to transform (no channel is typed \"eeg\", or 'picks' is empty).",
         call. = FALSE)
  }

  channel_data <- function(i) if (is_epochs) eeg_obj$data[i, , ] else eeg_obj$data[i, ]

  not_finite <- vapply(idx, function(i) !all(is.finite(channel_data(i))), logical(1))
  if (any(not_finite)) {
    stop("ERROR: Channel(s) with NA or non-finite values cannot be transformed: ",
         paste(eeg_obj$channels[idx][not_finite], collapse = ", "),
         ". Clean the data first (see annotate_nan()).", call. = FALSE)
  }

  # as.character() because startsWith() fails on the NULL of an empty history
  if (any(startsWith(as.character(unlist(eeg_obj$preprocessing_history)),
                     "apply_hilbert():"))) {
    warning("apply_hilbert(): this object already went through apply_hilbert() (see ",
            "preprocessing_history); transforming it again rarely makes sense.",
            call. = FALSE, immediate. = TRUE)
  }

  # ========== TRANSFORM, CHANNEL BY CHANNEL ==========

  for (i in idx) {
    z   <- .analytic_signal(channel_data(i), n_fft)
    res <- if (output == "envelope") Mod(z) else Arg(z)
    if (is_epochs) eeg_obj$data[i, , ] <- res else eeg_obj$data[i, ] <- res
  }

  # ========== HISTORY AND MESSAGE ==========

  what <- if (output == "envelope") "amplitude envelope" else "instantaneous phase (radians)"
  pad  <- if (n_fft > n_times) {
    paste0("zero-padded from ", n_times, " to ", n_fft, " samples")
  } else {
    "no zero padding"
  }

  eeg_obj$preprocessing_history <- c(eeg_obj$preprocessing_history, list(paste0(
    "apply_hilbert(): replaced ", length(idx), " channel(s) by their ", what,
    " (", pad, ").")))

  if (verbose) {
    message("apply_hilbert(): ", length(idx), " channel(s) replaced by their ", what,
            " (", pad, ").")
  }

  eeg_obj
}


#' Analytic Signal via the FFT (internal)
#'
#' Internal helper computing the analytic signal of a vector, or of every
#' column of a matrix, through the FFT: the signal is zero-padded to
#' \code{n_fft} samples, its negative frequencies are dropped and its positive
#' ones doubled (the zero frequency and, for an even length, the Nyquist
#' frequency are kept once), and the result is transformed back and cut to the
#' original length (Marple, 1999). It is the shared core of
#' \code{\link{apply_hilbert}}, \code{\link{.hilbert_envelope}} and
#' \code{\link{.hilbert_phase}}.
#'
#' @param x Numeric vector, or numeric matrix with one signal per column (time
#'   along the rows).
#' @param n_fft Integer, the FFT length, at least the number of time points.
#'   Default \code{NULL}: no padding (the length of the signal).
#'
#' @return Complex vector (for a vector) or complex matrix (for a matrix), the
#'   same shape as \code{x}. Its real part is \code{x} itself.
#'
#' @seealso \code{\link{apply_hilbert}}
#' @keywords internal
.analytic_signal <- function(x, n_fft = NULL) {
  is_vec <- is.null(dim(x))
  x      <- as.matrix(x)                          # time x signals
  n      <- nrow(x)
  if (is.null(n_fft)) n_fft <- n

  if (n_fft > n) x <- rbind(x, matrix(0, n_fft - n, ncol(x)))

  # Weights: keep the zero frequency (and the Nyquist frequency of an even
  # length) once, double the positive frequencies, drop the negative ones.
  h <- numeric(n_fft)
  if (n_fft %% 2 == 0) {
    h[1] <- 1
    h[n_fft / 2 + 1] <- 1
    if (n_fft > 2) h[2:(n_fft / 2)] <- 2
  } else {
    h[1] <- 1
    if (n_fft > 1) h[2:((n_fft + 1) / 2)] <- 2
  }

  z <- stats::mvfft(stats::mvfft(x) * h, inverse = TRUE)[seq_len(n), , drop = FALSE] / n_fft
  if (is_vec) z[, 1] else z
}
