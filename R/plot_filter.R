# ============================================================================
#                    Looking at a filter: plot_filter()
# ============================================================================
#
# A filter is easiest to understand from three pictures:
#
#   1. Its impulse response (the "spike test"): what comes out when one single
#      spike goes in. It shows how long the filter rings and whether it rings
#      before the spike (a zero-phase filter does, a causal one does not).
#   2. Its magnitude response: how much of every frequency is let through, in
#      dB (0 dB = untouched, -6 dB = half, -20 dB = a tenth).
#   3. Its group delay: how late every frequency comes out.
#
# plot_filter() does not describe the filter with a formula. It sends one
# spike through eeg_bandpass() or eeg_notch() themselves and measures what
# comes out, so the pictures always show the filter exactly as it is applied
# to data: both passes of a zero-phase IIR filter, the shift of a zero-phase
# FIR filter, the padding, everything. The other two pictures follow from the
# impulse response h[m] (m = the lag in samples, 0 = the spike) by Fourier
# analysis:
#
#   H(w)   = sum_m h[m] exp(-i w m)                       frequency response
#   tau(w) = Re( sum_m m h[m] exp(-i w m) / H(w) )        group delay, samples
#
# with w = 2 pi f / sampling rate.
#
# Autor: Christos Dalamarinis
# Date: Oct - 2026

# ----------------------------------------------------------------------------
# .measure_impulse_response() - send one spike through the real filter
# ----------------------------------------------------------------------------

#' Measure the impulse response of a filter function
#'
#' Builds a recording of zeros with one single spike (value 1) in the middle,
#' runs it through \code{fun}, and returns what comes out. The window must be
#' long enough to hold the whole response, so the two outer fifths of the
#' output are checked: if either holds more than \code{tol} of the peak, the
#' window is doubled and the measurement is repeated, up to \code{max_window}.
#' Whatever comes out is then cut to the stretch that holds the response:
#' symmetric around the spike for a zero-phase filter, from the spike onwards
#' for a causal one.
#'
#' @param sfreq      Numeric. Sampling rate in Hz.
#' @param fun        Function taking one \code{'eeg'} object and returning the
#'   filtered \code{'eeg'} object (e.g. a call to \code{\link{eeg_bandpass}}).
#' @param causal     Logical. TRUE: the response starts at the spike. FALSE
#'   (zero phase): the response is symmetric around the spike.
#' @param window     Numeric. First window length to try, in seconds.
#' @param max_window Numeric. Longest window to try, in seconds.
#' @param tol        Numeric. How much of the peak the two outer fifths of the
#'   window may hold.
#' @param trim       Numeric. Samples smaller than this fraction of the peak
#'   are cut off the two ends of the returned response.
#' @return A list with \code{h} (the response), \code{lags} (its time axis in
#'   samples, 0 = the spike), \code{window} (seconds, the window that was
#'   used) and \code{contained} (FALSE if the response did not fit into
#'   \code{max_window}).
#' @keywords internal
.measure_impulse_response <- function(sfreq, fun, causal, window = 32,
                                      max_window = 1024, tol = 1e-9,
                                      trim = 1e-10) {

  dur   <- window
  first <- TRUE
  repeat {
    n  <- 2L * as.integer(ceiling(dur * sfreq / 2)) + 1L   # odd: spike in the middle
    c0 <- (n + 1L) %/% 2L
    x  <- numeric(n)
    x[c0] <- 1
    eeg <- new_eeg(data = matrix(x, nrow = 1L), channels = "Cz",
                   sampling_rate = sfreq)

    # a warning of the real filter (e.g. an unused argument) is shown once
    y <- if (first) fun(eeg)$data[1, ] else suppressWarnings(fun(eeg)$data[1, ])
    first <- FALSE

    peak  <- max(abs(y))
    if (!is.finite(peak) || peak == 0) {
      stop("The filter returned no response to a single spike.", call. = FALSE)
    }
    fifth <- n %/% 5L
    outer <- max(abs(y[c(seq_len(fifth), (n - fifth + 1L):n)])) / peak
    contained <- outer <= tol
    if (contained || dur * 2 > max_window) break
    dur <- dur * 2
  }

  big <- which(abs(y) > trim * peak)
  lo  <- min(big)
  hi  <- max(big)
  if (causal) {
    lo <- c0                                   # nothing can come before the spike
  } else {
    reach <- max(c0 - lo, hi - c0)             # symmetric around the spike
    lo <- c0 - reach
    hi <- c0 + reach
  }
  hi <- max(hi, lo + 2L)                       # at least three samples

  list(h = y[lo:hi], lags = (lo:hi) - c0, window = dur, contained = contained)
}

# ----------------------------------------------------------------------------
# .frequency_response() - magnitude and group delay of an impulse response
# ----------------------------------------------------------------------------

#' Frequency response and group delay of an impulse response
#'
#' Evaluates \eqn{H(w) = \sum_m h[m] e^{-i w m}} and the group delay
#' \eqn{\tau(w) = Re(\sum_m m h[m] e^{-i w m} / H(w))} at the requested
#' frequencies. For a response of moderate length the sums are evaluated
#' exactly at every requested frequency (the powers of \eqn{e^{-i w}} are built
#' by repeated multiplication). A very long response (more than 20000
#' samples, for example a 0.1 Hz IIR high-pass at 2048 Hz) would take
#' seconds that way, so there the two sums are computed with an FFT that is
#' zero-padded to many times the length of the response (the frequency spacing
#' of the FFT is then far finer than anything the response varies on) and read
#' off at the requested frequencies by linear interpolation of the complex
#' values; the difference from the exact sums is far below what a plot can
#' show.
#'
#' @param h     Numeric vector. The impulse response.
#' @param lags  Integer vector, same length as \code{h}: consecutive lags in
#'   samples (0 = the spike).
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param freq  Numeric vector. Frequencies in Hz, between 0 and the Nyquist
#'   frequency.
#' @param exact Logical. TRUE: exact sums; FALSE: the FFT version. Default:
#'   exact for responses of up to 20000 samples.
#' @return A list with \code{H} (complex response), \code{magnitude_db}
#'   (\eqn{10 \log_{10}|H|^2}, floored at -200 dB) and \code{delay} (group
#'   delay in seconds, \code{NA} where the magnitude is below -100 dB: there
#'   almost nothing is left and the delay is dominated by rounding).
#' @keywords internal
.frequency_response <- function(h, lags, sfreq, freq, exact = length(h) <= 20000L) {

  n_h <- length(h)
  lh  <- lags * h

  if (exact) {
    w <- 2 * pi * freq / sfreq
    H <- complex(length(w))
    D <- complex(length(w))
    for (k in seq_along(w)) {
      z  <- exp(-1i * w[k])
      p  <- exp(-1i * w[k] * lags[1]) * cumprod(c(1+0i, rep(z, n_h - 1L)))
      H[k] <- sum(h * p)
      D[k] <- sum(lh * p)
    }
  } else {
    n_fft <- 2^min(22, ceiling(log2(max(64 * n_h, 2^17))))
    slot  <- (lags %% n_fft) + 1L              # lag m sits at place m (mod n_fft)
    zh <- numeric(n_fft); zh[slot] <- h
    zd <- numeric(n_fft); zd[slot] <- lh
    keep   <- seq_len(n_fft %/% 2L + 1L)
    f_bins <- (keep - 1L) * sfreq / n_fft
    Hf <- fft(zh)[keep]
    Df <- fft(zd)[keep]
    interp <- function(v) {                    # rule = 2: a frequency one rounding
      stats::approx(f_bins, Re(v), xout = freq, rule = 2)$y +    # step above the
        1i * stats::approx(f_bins, Im(v), xout = freq, rule = 2)$y   # Nyquist is kept
    }
    H <- interp(Hf)
    D <- interp(Df)
  }

  mag   <- 10 * log10(pmax(Mod(H)^2, 1e-20))
  delay <- Re(D / H) / sfreq
  delay[mag < -100] <- NA_real_                # not meaningful, see above

  list(H = H, magnitude_db = mag, delay = delay)
}

# ----------------------------------------------------------------------------
# small helpers for the drawing
# ----------------------------------------------------------------------------

#' Frequency tick marks for a plot
#'
#' @param flim   Numeric length 2. Frequency limits in Hz.
#' @param fscale "log" or "linear".
#' @return Numeric vector of tick positions inside \code{flim}.
#' @keywords internal
.freq_ticks <- function(flim, fscale) {
  if (fscale == "linear") return(pretty(flim))
  decades <- 10^seq(floor(log10(flim[1])), ceiling(log10(flim[2])))
  ticks   <- sort(c(decades, 2 * decades, 5 * decades))
  ticks[ticks >= flim[1] * (1 - 1e-9) & ticks <= flim[2] * (1 + 1e-9)]
}

#' Title of a filter plot
#'
#' @param kind   "bandpass" or "notch".
#' @param method "fir" or "iir".
#' @param phase  The phase word.
#' @param l_freq,h_freq,freqs The frequencies given by the user.
#' @param order  Integer. Butterworth order (used for IIR).
#' @return A character string.
#' @keywords internal
.filter_title <- function(kind, method, phase, l_freq, h_freq, freqs, order) {

  what <- if (kind == "notch") {
    paste0("notch ", paste(freqs, collapse = ", "), " Hz")
  } else if (!is.null(l_freq) && !is.null(h_freq)) {
    paste0("band-pass ", l_freq, "-", h_freq, " Hz")
  } else if (!is.null(l_freq)) {
    paste0("high-pass ", l_freq, " Hz")
  } else {
    paste0("low-pass ", h_freq, " Hz")
  }

  family <- if (method == "iir") {
    paste0("IIR (Butterworth, order ", order, ") ")
  } else {
    "FIR "
  }

  how <- switch(phase,
                zero    = if (method == "iir") "zero phase (forward and backward)"
                          else "zero phase",
                forward = "causal (one pass)",
                minimum = "causal (minimum phase)",
                phase)

  paste0(family, what, ", ", how)
}

# ----------------------------------------------------------------------------
# plot_filter()
# ----------------------------------------------------------------------------

#' Plot What a Filter Does (Impulse, Magnitude and Delay)
#'
#' Shows a filter of \code{\link{eeg_bandpass}} or \code{\link{eeg_notch}} in
#' three pictures, one under the other: its impulse response (what comes out
#' of the filter when one single spike goes in), its magnitude response (how
#' much of every frequency is let through, in dB) and its group delay (how
#' late every frequency comes out). Use it to check a filter before it is
#' applied to data, or to compare FIR and IIR filters and zero-phase and
#' causal ones.
#'
#' The filter is shown exactly as it is applied to data: the function sends a
#' single spike through \code{eeg_bandpass()} (when \code{l_freq} and/or
#' \code{h_freq} are given) or through \code{eeg_notch()} (when \code{freqs}
#' is given) and measures what comes out. Everything the real filter does is
#' therefore in the pictures: for \code{phase = "zero"} both passes of an IIR
#' filter (so its cutoffs are 6 dB down instead of 3 dB) and the delay
#' correction of an FIR filter (so there is no delay); for a causal filter its
#' true delay.
#'
#' \strong{Impulse response.} For a zero-phase filter the spike is at time 0
#' and the response is symmetric: it rings before the spike as well as after
#' it. For a causal filter (\code{phase = "minimum"} for FIR, \code{"forward"}
#' for IIR) the response starts at the spike and nothing comes before it.
#'
#' \strong{Magnitude.} \eqn{10 \log_{10} |H|^2} of the frequency response
#' \eqn{H}, so 0 dB means a frequency passes untouched, -3 dB about 70%,
#' -6 dB half, -20 dB a tenth, -40 dB a hundredth of the amplitude.
#'
#' \strong{Delay.} The group delay in seconds. A zero-phase filter has no
#' delay (a flat line at 0); a causal filter delays different frequencies by
#' different amounts, usually most near its cutoffs. Where the magnitude is
#' below -40 dB (grey bands) there is almost no signal left and the delay is
#' not meaningful.
#'
#' The measurement uses a window of zeros around the spike that is made longer
#' until the response fits into it (at least 32 seconds, more for a filter
#' with a low high-pass cutoff). A very long causal FIR filter, such as a 0.1
#' Hz high-pass, is as slow to build here as it is when filtering.
#'
#' @param x An object of class \code{'eeg'} (only its sampling rate is used),
#'   or a single number: the sampling rate in Hz.
#' @param l_freq Numeric or \code{NULL}. High-pass cutoff in Hz, as in
#'   \code{\link{eeg_bandpass}}. With \code{h_freq} it makes a band-pass.
#' @param h_freq Numeric or \code{NULL}. Low-pass cutoff in Hz, as in
#'   \code{\link{eeg_bandpass}}.
#' @param freqs Numeric or \code{NULL}. Notch frequencies in Hz, as in
#'   \code{\link{eeg_notch}}. Give either \code{freqs} or \code{l_freq} and/or
#'   \code{h_freq}, not both.
#' @param method Character. \code{"fir"} (default) or \code{"iir"}.
#' @param phase Character. \code{"zero"} (default); \code{"minimum"} for a
#'   causal FIR filter, \code{"forward"} for a causal IIR filter.
#' @param ... Further arguments for the filter function, such as
#'   \code{iir_params}, \code{l_trans_bandwidth}, \code{h_trans_bandwidth},
#'   \code{filter_length}, \code{fir_window}, \code{pad} (for
#'   \code{eeg_bandpass()}) or \code{notch_widths}, \code{trans_bandwidth}
#'   (for \code{eeg_notch()}). They mean exactly what they mean there.
#' @param panels Character vector. Which pictures to draw, from
#'   \code{"time"} (impulse response), \code{"magnitude"} and \code{"delay"},
#'   from top to bottom. Default: all three.
#' @param fscale Character. Frequency axis: \code{"log"} (default) or
#'   \code{"linear"}.
#' @param flim Numeric vector of length 2 or \code{NULL}. Frequency limits in
#'   Hz. Default: 0.1 Hz to half the sampling rate (0 Hz to half the sampling
#'   rate for \code{fscale = "linear"}). A log axis needs a start above 0; for
#'   a high-pass below 0.5 Hz start lower, e.g. \code{flim = c(0.01, 256)}.
#' @param alim Numeric vector of length 2. Magnitude limits in dB. Default
#'   \code{c(-80, 10)}.
#' @param tlim Numeric vector of length 2 or \code{NULL}. Time limits in
#'   seconds for the impulse response. Default: the stretch where the response
#'   is larger than 0.1% of its peak, with a little margin.
#' @param dlim Numeric vector of length 2 or \code{NULL}. Delay limits in
#'   seconds. Default: chosen from the delay where the magnitude is above
#'   -40 dB.
#' @param main Character or \code{NULL}. Title above the first picture.
#'   Default: a description of the filter.
#' @param col Colour of the curves.
#'
#' @return Invisibly, a list with:
#'  \describe{
#'    \item{time}{Numeric vector, time axis of the impulse response in seconds
#'      (0 = the spike)}
#'    \item{impulse}{Numeric vector, the impulse response}
#'    \item{freq}{Numeric vector, the 1000 frequencies in Hz at which the
#'      response was evaluated}
#'    \item{magnitude_db}{Numeric vector, magnitude in dB at \code{freq}}
#'    \item{delay}{Numeric vector, group delay in seconds at \code{freq}
#'      (\code{NA} where the magnitude is below -100 dB: there almost nothing
#'      is left and the delay is dominated by rounding)}
#'    \item{sfreq}{Numeric, the sampling rate in Hz}
#'    \item{window}{Numeric, length in seconds of the window of zeros that was
#'      needed}
#'  }
#'
#' @examples
#' \dontrun{
#'   # Which one lets less of the 50 Hz hum through, and which one rings less?
#'   plot_filter(eeg, h_freq = 40)                       # FIR low-pass
#'   plot_filter(eeg, h_freq = 40, method = "iir")       # IIR low-pass
#'
#'   # A causal filter: no ringing before the spike, but a delay
#'   plot_filter(eeg, l_freq = 1, h_freq = 40, method = "iir", phase = "forward")
#'
#'   # A notch, with the sampling rate given as a plain number
#'   plot_filter(512, freqs = 50, method = "iir")
#'
#'   # Only the magnitude, on a linear frequency axis
#'   plot_filter(eeg, l_freq = 1, h_freq = 40, panels = "magnitude",
#'               fscale = "linear", flim = c(0, 60))
#' }
#'
#' @references
#' Oppenheim, A. V., & Schafer, R. W. (2010). Discrete-Time Signal Processing
#' (3rd ed.). Pearson.
#'
#' @seealso \code{\link{eeg_bandpass}}, \code{\link{eeg_notch}}
#'
#' @importFrom graphics plot abline lines par grid polygon axis
#' @importFrom grDevices adjustcolor
#'
#' @export
plot_filter <- function(x,
                        l_freq = NULL,
                        h_freq = NULL,
                        freqs  = NULL,
                        method = c("fir", "iir"),
                        phase  = "zero",
                        ...,
                        panels = c("time", "magnitude", "delay"),
                        fscale = c("log", "linear"),
                        flim   = NULL,
                        alim   = c(-80, 10),
                        tlim   = NULL,
                        dlim   = NULL,
                        main   = NULL,
                        col    = "#1f77b4") {

  # ========== INPUT VALIDATION ==========

  if (inherits(x, "eeg")) {
    sfreq <- x$sampling_rate
  } else if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
    sfreq <- x
  } else {
    stop("'x' must be an object of class 'eeg' or a sampling rate in Hz.",
         call. = FALSE)
  }
  if (!is.numeric(sfreq) || length(sfreq) != 1L || !is.finite(sfreq) ||
      sfreq <= 0) {
    stop("The sampling rate must be a single positive number.", call. = FALSE)
  }

  use_notch <- !is.null(freqs)
  use_band  <- !is.null(l_freq) || !is.null(h_freq)
  if (use_notch && use_band) {
    stop("Give either 'freqs' (notch) or 'l_freq' and/or 'h_freq' ",
         "(band-pass, high-pass, low-pass), not both.", call. = FALSE)
  }
  if (!use_notch && !use_band) {
    stop("Give 'l_freq' and/or 'h_freq' (band-pass, high-pass, low-pass) ",
         "or 'freqs' (notch).", call. = FALSE)
  }

  method <- match.arg(method)
  phase  <- .check_phase(phase, method)
  fscale <- match.arg(fscale)

  all_panels <- c("time", "magnitude", "delay")
  if (!is.character(panels) || length(panels) == 0L ||
      !all(panels %in% all_panels)) {
    stop("'panels' must be one or more of \"time\", \"magnitude\", \"delay\".",
         call. = FALSE)
  }
  panels <- unique(panels)

  nyquist <- sfreq / 2
  if (is.null(flim)) {
    flim <- c(if (fscale == "log") 0.1 else 0, nyquist)
  }
  if (!is.numeric(flim) || length(flim) != 2L || anyNA(flim) ||
      flim[1] >= flim[2]) {
    stop("'flim' must be two increasing frequencies in Hz.", call. = FALSE)
  }
  if (fscale == "log" && flim[1] <= 0) {
    stop("'flim[1]' must be above 0 Hz for a log frequency axis.", call. = FALSE)
  }
  if (flim[1] < 0 || flim[2] > nyquist) {
    stop("'flim' must lie between 0 Hz and the Nyquist frequency (", nyquist,
         " Hz).", call. = FALSE)
  }
  if (!is.numeric(alim) || length(alim) != 2L || anyNA(alim) ||
      alim[1] >= alim[2]) {
    stop("'alim' must be two increasing values in dB.", call. = FALSE)
  }
  for (lim in list(tlim = tlim, dlim = dlim)) {
    if (!is.null(lim) && (!is.numeric(lim) || length(lim) != 2L ||
                          anyNA(lim) || lim[1] >= lim[2])) {
      stop("'tlim' and 'dlim' must be NULL or two increasing values in seconds.",
           call. = FALSE)
    }
  }

  # ========== MEASURE THE FILTER ==========

  kind <- if (use_notch) "notch" else "bandpass"
  fun  <- if (use_notch) {
    function(e) eeg_notch(e, freqs = freqs, method = method, phase = phase,
                          ..., verbose = FALSE)
  } else {
    function(e) eeg_bandpass(e, l_freq = l_freq, h_freq = h_freq,
                             method = method, phase = phase, ...,
                             verbose = FALSE)
  }

  # a slow high-pass rings for a long time: start with a longer window
  window0 <- if (!is.null(l_freq) && is.numeric(l_freq) && l_freq > 0) {
    min(600, max(32, 20 / l_freq))
  } else {
    32
  }

  resp <- .measure_impulse_response(sfreq, fun, causal = (phase != "zero"),
                                    window = window0)
  if (!resp$contained) {
    warning("The response of this filter is longer than the longest window ",
            "tried (", resp$window, " s); the pictures show it cut off.",
            call. = FALSE)
  }

  # ========== FREQUENCY RESPONSE AND GROUP DELAY ==========

  f_grid <- if (fscale == "log") {
    10^seq(log10(flim[1]), log10(flim[2]), length.out = 1000L)
  } else {
    seq(flim[1], flim[2], length.out = 1000L)
  }
  f_grid[c(1L, 1000L)] <- flim       # the ends exactly as asked: 10^log10(x) is not always x
  fr <- .frequency_response(resp$h, resp$lags, sfreq, f_grid)
  t_axis <- resp$lags / sfreq

  # ========== DRAW ==========

  order <- list(...)$iir_params$order
  if (is.null(order)) order <- 4L
  if (is.null(main)) {
    main <- .filter_title(kind, method, phase, l_freq, h_freq, freqs, order)
  }
  marks <- if (use_notch) freqs else c(l_freq, h_freq)

  # time window: where the response is larger than 0.1 % of its peak
  if (is.null(tlim)) {
    big  <- which(abs(resp$h) > 1e-3 * max(abs(resp$h)))
    tlim <- t_axis[range(big)]
    pad  <- 0.1 * diff(tlim)
    if (pad == 0) pad <- 1 / sfreq
    tlim <- tlim + c(-pad, pad)
  }

  # delay window: from the delay where something is left of the signal
  if (is.null(dlim)) {
    ok  <- fr$magnitude_db > -40 & !is.na(fr$delay)
    rng <- if (any(ok)) range(fr$delay[ok]) else c(0, 0)
    if (diff(rng) < 1e-4) {
      dlim <- mean(rng) + c(-0.05, 0.05)
    } else {
      pad  <- 0.1 * diff(rng)
      dlim <- rng + c(-pad, pad)
      if (dlim[1] > 0 && dlim[1] < 0.5 * dlim[2]) dlim[1] <- 0
    }
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(length(panels), 1L), mar = c(4.2, 4.6, 2.4, 1),
      mgp = c(2.7, 0.7, 0))

  log_arg  <- if (fscale == "log") "x" else ""
  ticks    <- .freq_ticks(flim, fscale)
  nulled   <- fr$magnitude_db <= -39.9
  first    <- TRUE

  for (panel in panels) {
    if (panel == "time") {
      plot(t_axis, resp$h, type = "n", xlim = tlim,
           xlab = "Time (s)", ylab = "Amplitude", main = if (first) main else "")
      grid(col = "grey88")
      abline(h = 0, col = "grey70")
      abline(v = 0, col = "grey70", lty = 2)
      lines(t_axis, resp$h, col = col, lwd = 1.4)
    } else if (panel == "magnitude") {
      plot(f_grid, fr$magnitude_db, type = "n", log = log_arg, xlim = flim,
           ylim = alim, xaxt = "n", xlab = "Frequency (Hz)",
           ylab = "Magnitude (dB)", main = if (first) main else "")
      axis(1, at = ticks, labels = ticks)
      abline(v = ticks, h = pretty(alim), col = "grey88")
      abline(v = marks, col = "grey55", lty = 3)
      lines(f_grid, fr$magnitude_db, col = col, lwd = 1.4)
    } else {
      plot(f_grid, fr$delay, type = "n", log = log_arg, xlim = flim,
           ylim = dlim, xaxt = "n", xlab = "Frequency (Hz)",
           ylab = "Group delay (s)", main = if (first) main else "")
      axis(1, at = ticks, labels = ticks)
      abline(v = ticks, h = pretty(dlim), col = "grey88")
      abline(v = marks, col = "grey55", lty = 3)
      runs <- rle(nulled)
      ends <- cumsum(runs$lengths)
      begs <- ends - runs$lengths + 1L
      for (r in which(runs$values)) {
        xs <- f_grid[c(begs[r], ends[r])]
        polygon(c(xs, rev(xs)), rep(dlim, each = 2L),
                col = adjustcolor("black", 0.08), border = NA)
      }
      abline(h = 0, col = "grey70")
      lines(f_grid, fr$delay, col = col, lwd = 1.4)
    }
    first <- FALSE
  }

  invisible(list(time         = t_axis,
                 impulse      = resp$h,
                 freq         = f_grid,
                 magnitude_db = fr$magnitude_db,
                 delay        = fr$delay,
                 sfreq        = sfreq,
                 window       = resp$window))
}
