#' ============================================================================
#'                    Fourier / Spectral Analysis Functions
#' ============================================================================
#'
#' This module provides functions for frequency-domain analysis of EEG data,
#' implementing four complementary spectral estimation approaches:
#'
#'   eeg_fft()        - Raw FFT: amplitude, power, and phase spectra using the
#'                      Fast Fourier Transform with optional windowing.
#'                      Best for short, stationary segments or ERPs.
#'
#'   eeg_psd_welch()  - Welch averaged periodogram (PSD in uV^2/Hz). Reduces
#'                      spectral variance by averaging overlapping windowed
#'                      segments. The standard method for resting-state EEG.
#'
#'   eeg_multitaper() - Multitaper spectral estimation using Discrete Prolate
#'                      Spheroidal Sequences (DPSS / Slepian tapers). Provides
#'                      an optimal bias-variance trade-off via multiple
#'                      orthogonal tapers. Output in uV^2/Hz.
#'
#'   eeg_band_power() - Extracts total or relative power in the canonical EEG
#'                      frequency bands (delta, theta, alpha, beta, gamma).
#'                      Accepts raw EEG objects or a pre-computed eeg_spectrum.
#'
#' All four functions accept objects of class 'eeg', 'eeg_epochs', and
#' 'eeg_evoked', enabling seamless integration with the eeganalysis pipeline.
#' All return objects of class 'eeg_spectrum' (except eeg_band_power which
#' returns a data frame or list of data frames for multi-condition evoked data).
#'
#' Return object structure (eeg_spectrum):
#'   $power         - matrix [channels x frequencies] or named list for evoked
#'   $amplitude     - sqrt(power)
#'   $phase         - phase angles in radians (eeg_fft only; NULL otherwise)
#'   $frequencies   - numeric vector of frequencies in Hz
#'   $channels      - character vector of channel names
#'   $sampling_rate - recording sampling rate in Hz
#'   $method        - "fft" | "welch" | "multitaper"
#'   $input_class   - "eeg" | "eeg_epochs" | "eeg_evoked"
#'   $n_fft         - FFT length used
#'   $window_type   - window function applied ("hann", "hamming", etc.)
#'   $per_epoch     - logical; TRUE if per-epoch spectra were computed
#'   $epoch_power   - 3D array [channels x freqs x epochs] (per_epoch = TRUE)
#'   $n_epochs      - number of epochs (eeg_epochs input only)
#'   $conditions    - condition names (eeg_evoked input only)
#'   + method-specific fields (window_length, overlap, time_bandwidth, etc.)
#'
#  Author: Christos Dalamarinis
#  Date: March 2026
#  Status: In development
#' ============================================================================
#'
#'
#' ============================================================================
#'                          INTERNAL HELPER FUNCTIONS
#' ============================================================================
#'
#' Generate a spectral window vector
#'
#' @param n       Integer. Window length in samples.
#' @param type    Character. One of "hann", "hamming", "blackman", "none".
#' @return Numeric vector of length n.
#' @keywords internal
.make_window <- function(n, type = c("hann", "hamming", "blackman", "none")) {
  type <- match.arg(type)
  idx  <- 0:(n - 1)
  switch(type,
         hann     = {
           w <- 0.5 * (1 - cos(2 * pi * idx / (n - 1)))
           w / max(w)                          # normalise to unit peak
         },
         hamming  = 0.54 - 0.46 * cos(2 * pi * idx / (n - 1)),
         blackman = {
           w <- 0.42 - 0.50 * cos(2 * pi * idx / (n - 1)) +
             0.08 * cos(4 * pi * idx / (n - 1))
           pmax(w, 0)
         },
         none     = rep(1.0, n)
  )
}


#' Compute one-sided FFT amplitude, power, and phase for a single channel
#'
#' The amplitude spectrum is corrected for the one-sided representation
#' (all bins except DC and Nyquist are doubled) and normalised by n_fft,
#' yielding peak amplitudes in the same units as the input signal (uV).
#' Power = amplitude^2.  Phase = Arg(FFT).
#'
#' @param x     Numeric vector. Signal (already windowed if desired).
#' @param sr    Numeric. Sampling rate in Hz.
#' @param n_fft Integer. FFT length (zero-pads or truncates x to this length).
#' @return Named list: power, amplitude, phase (all length floor(n_fft/2)+1),
#'         and freqs (frequency vector in Hz).
#' @keywords internal
.fft_one_sided <- function(x, sr, n_fft) {
  n <- length(x)
  
  # Zero-pad or truncate to n_fft
  if (n_fft > n) {
    x <- c(x, rep(0.0, n_fft - n))
  } else if (n_fft < n) {
    x <- x[seq_len(n_fft)]
  }
  
  X      <- fft(x)
  n_half <- floor(n_fft / 2L) + 1L
  X_os   <- X[seq_len(n_half)]
  
  # Amplitude: normalise by n_fft, then double interior bins (one-sided correction)
  amp      <- Mod(X_os) / n_fft
  amp[seq(2L, n_half - 1L)] <- 2.0 * amp[seq(2L, n_half - 1L)]
  
  freqs <- seq(0, sr / 2, length.out = n_half)
  
  list(
    power     = amp^2,
    amplitude = amp,
    phase     = Arg(X_os),
    freqs     = freqs
  )
}


#' Compute Discrete Prolate Spheroidal Sequences (DPSS / Slepian tapers)
#'
#' Uses the symmetric tridiagonal commuting matrix approach (Slepian 1978,
#' Thomson 1982). The eigenvectors of this tridiagonal matrix are the DPSS.
#' This is numerically stable for typical EEG epoch lengths (n <= ~4096).
#' For very long signals (n > 4096) consider segmenting first.
#'
#' @param n  Integer. Signal length.
#' @param nw Numeric. Time-bandwidth product (e.g. 4 gives NW=4).
#' @param k  Integer. Number of tapers to return (k <= 2*nw - 1 recommended).
#' @return Numeric matrix of dimensions [n x k]. Columns are the tapers,
#'         normalised to unit energy. Sign convention: first element positive.
#' @keywords internal
.compute_dpss <- function(n, nw, k) {
  W <- nw / n           # normalised half-bandwidth
  
  # Diagonal and off-diagonal elements of the tridiagonal commuting matrix
  i_main  <- 0L:(n - 1L)
  d_main  <- ((n - 1L - 2L * i_main) / 2)^2 * cos(2 * pi * W)
  i_off   <- seq_len(n - 1L)
  d_off   <- i_off * (n - i_off) / 2
  
  # Build full symmetric tridiagonal matrix
  T_mat <- diag(d_main)
  for (i in i_off) {
    T_mat[i, i + 1L] <- d_off[i]
    T_mat[i + 1L, i] <- d_off[i]
  }
  
  # Eigen-decomposition (symmetric); take the k vectors with largest |eigenvalue|
  eig    <- eigen(T_mat, symmetric = TRUE)
  ord    <- order(abs(eig$values), decreasing = TRUE)
  tapers <- eig$vectors[, ord[seq_len(k)], drop = FALSE]
  
  # Normalise each taper: unit energy, positive first element (sign convention)
  for (j in seq_len(k)) {
    if (tapers[1L, j] < 0) tapers[, j] <- -tapers[, j]
    tapers[, j] <- tapers[, j] / sqrt(sum(tapers[, j]^2))
  }
  
  tapers   # [n x k]
}


#' Apply pre-computed DPSS tapers to a signal and return multitaper PSD
#'
#' @param x      Numeric vector of length n (must match nrow of tapers).
#' @param sr     Numeric. Sampling rate in Hz.
#' @param tapers Numeric matrix [n x k] from .compute_dpss().
#' @return Numeric vector of PSD estimates (uV^2/Hz), length floor(n/2)+1.
#' @keywords internal
.mtaper_psd <- function(x, sr, tapers) {
  n      <- length(x)
  k      <- ncol(tapers)
  n_half <- floor(n / 2L) + 1L
  acc    <- numeric(n_half)
  
  for (j in seq_len(k)) {
    X   <- fft(x * tapers[, j])[seq_len(n_half)]
    p   <- (Mod(X)^2) * 2 / (n * sr)
    # DC and Nyquist bins must not be doubled
    p[c(1L, n_half)] <- p[c(1L, n_half)] / 2
    acc <- acc + p
  }
  
  acc / k    # average across tapers
}


# ============================================================================
#                         RAW FFT SPECTRUM  (eeg_fft)
# ============================================================================

#' Compute FFT Amplitude, Power, and Phase Spectra of EEG Data
#'
#' Applies the Fast Fourier Transform to each channel of an EEG object,
#' returning amplitude, power, and phase spectra. Supports optional spectral
#' windowing to control spectral leakage. Accepts 'eeg', 'eeg_epochs', and
#' 'eeg_evoked' objects.
#'
#' For 'eeg_epochs' input the behaviour is controlled by \code{per_epoch}:
#' \itemize{
#'   \item FALSE (default): epochs are averaged in the time domain first
#'         (equivalent to FFT of the ERP), suitable for phase-locked activity.
#'   \item TRUE: FFT is computed per epoch; the mean power across epochs is
#'         stored in \code{$power} and individual epoch power in
#'         \code{$epoch_power} (3-D array).
#' }
#' For 'eeg_evoked' input one spectrum is computed per condition.
#'
#' @param input_obj An object of class 'eeg', 'eeg_epochs', or 'eeg_evoked'.
#' @param n_fft     Integer or NULL. FFT length. NULL (default) uses the signal
#'                  length. Values > signal length zero-pad (interpolates the
#'                  spectrum). Values < signal length truncate.
#' @param window    Character. Spectral window: "hann" (default), "hamming",
#'                  "blackman", or "none". Windowing reduces spectral leakage
#'                  at the cost of frequency resolution. Use "none" only when
#'                  the signal length is an exact multiple of the frequency of
#'                  interest.
#' @param per_epoch Logical. Only relevant for 'eeg_epochs' input. If TRUE,
#'                  compute FFT per epoch and store results in $epoch_power.
#'                  Default: FALSE.
#' @param verbose   Logical. Print a processing summary. Default: TRUE.
#'
#' @return An object of class 'eeg_spectrum'. See module header for full field
#'   descriptions. The \code{$phase} field is populated for raw FFT output;
#'   it is NULL for Welch and multitaper outputs.
#'
#' @examples
#' \dontrun{
#'   # Continuous EEG
#'   spec <- eeg_fft(eeg, window = "hann")
#'
#'   # Epoched data - per-epoch power (e.g. for induced oscillations)
#'   spec <- eeg_fft(epochs, per_epoch = TRUE)
#'
#'   # Evoked (condition-wise spectra)
#'   spec <- eeg_fft(evoked_obj)
#'   spec$power[["cond1"]]   # [channels x frequencies]
#' }
#'
#' @seealso \code{\link{eeg_psd_welch}}, \code{\link{eeg_multitaper}},
#'   \code{\link{eeg_band_power}}
#' @export
eeg_fft <- function(input_obj,
                    n_fft     = NULL,
                    window    = c("hann", "hamming", "blackman", "none"),
                    per_epoch = FALSE,
                    verbose   = TRUE) {
  
  window      <- match.arg(window)
  valid_cls   <- c("eeg", "eeg_epochs", "eeg_evoked")
  input_class <- class(input_obj)[1L]
  
  # ========== INPUT VALIDATION ==========
  if (!input_class %in% valid_cls) {
    stop("Input must be one of: ", paste(valid_cls, collapse = ", "),
         ". Got: '", input_class, "'.", call. = FALSE)
  }
  
  if (!is.null(n_fft) && (!is.numeric(n_fft) || n_fft < 1)) {
    stop("n_fft must be a positive integer or NULL.", call. = FALSE)
  }
  
  # ========== DISPATCH: eeg ==========
  if (input_class == "eeg") {
    
    sr         <- input_obj$sampling_rate
    data_mat   <- input_obj$data           # [channels x timepoints]
    chan_names <- input_obj$channels
    n_ch       <- nrow(data_mat)
    n_tp       <- ncol(data_mat)
    n_fft_use  <- if (is.null(n_fft)) n_tp else as.integer(n_fft)
    n_half     <- floor(n_fft_use / 2L) + 1L
    freqs      <- seq(0, sr / 2, length.out = n_half)
    win_vec    <- .make_window(n_tp, window)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("FFT Spectral Analysis  [eeg]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Signal length:     ", n_tp, "samples (", round(n_tp / sr, 2), "s)\n")
      cat("  FFT length:        ", n_fft_use, "\n")
      cat("  Window:            ", window, "\n")
      cat("  Freq. resolution:  ", round(freqs[2L], 4), "Hz\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    pwr_mat   <- matrix(0, nrow = n_ch, ncol = n_half,
                        dimnames = list(chan_names, NULL))
    amp_mat   <- matrix(0, nrow = n_ch, ncol = n_half,
                        dimnames = list(chan_names, NULL))
    phase_mat <- matrix(0, nrow = n_ch, ncol = n_half,
                        dimnames = list(chan_names, NULL))
    
    for (i in seq_len(n_ch)) {
      res              <- .fft_one_sided(data_mat[i, ] * win_vec, sr, n_fft_use)
      pwr_mat[i, ]   <- res$power
      amp_mat[i, ]   <- res$amplitude
      phase_mat[i, ] <- res$phase
    }
    
    return(structure(
      list(
        power         = pwr_mat,
        amplitude     = amp_mat,
        phase         = phase_mat,
        frequencies   = freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "fft",
        input_class   = input_class,
        n_fft         = n_fft_use,
        window_type   = window,
        per_epoch     = FALSE,
        epoch_power   = NULL,
        conditions    = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_epochs ==========
  if (input_class == "eeg_epochs") {
    
    if (is.null(input_obj$data)) {
      stop("Epoch data not loaded. Set preload = TRUE in epoch_eeg().", call. = FALSE)
    }
    
    sr         <- input_obj$sampling_rate
    data_arr   <- input_obj$data      # [channels x samples x epochs]
    chan_names <- input_obj$channels
    n_ch       <- dim(data_arr)[1L]
    n_samp     <- dim(data_arr)[2L]
    n_ep       <- dim(data_arr)[3L]
    n_fft_use  <- if (is.null(n_fft)) n_samp else as.integer(n_fft)
    n_half     <- floor(n_fft_use / 2L) + 1L
    freqs      <- seq(0, sr / 2, length.out = n_half)
    win_vec    <- .make_window(n_samp, window)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("FFT Spectral Analysis  [eeg_epochs]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Samples per epoch: ", n_samp, "\n")
      cat("  Epochs:            ", n_ep, "\n")
      cat("  FFT length:        ", n_fft_use, "\n")
      cat("  Window:            ", window, "\n")
      cat("  Per epoch:         ", per_epoch, "\n")
      cat("  Freq. resolution:  ", round(freqs[2L], 4), "Hz\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    if (per_epoch) {
      # Compute FFT for every epoch; store 3-D array [channels x freqs x epochs]
      epoch_pwr   <- array(0, dim = c(n_ch, n_half, n_ep))
      epoch_amp   <- array(0, dim = c(n_ch, n_half, n_ep))
      epoch_phase <- array(0, dim = c(n_ch, n_half, n_ep))
      dimnames(epoch_pwr)[[1L]] <- chan_names
      
      for (ep in seq_len(n_ep)) {
        for (i in seq_len(n_ch)) {
          res                    <- .fft_one_sided(data_arr[i, , ep] * win_vec,
                                                   sr, n_fft_use)
          epoch_pwr[i, , ep]   <- res$power
          epoch_amp[i, , ep]   <- res$amplitude
          epoch_phase[i, , ep] <- res$phase
        }
      }
      
      # Mean power and amplitude across epochs
      pwr_mat   <- apply(epoch_pwr,   c(1L, 2L), mean)
      amp_mat   <- apply(epoch_amp,   c(1L, 2L), mean)
      rownames(pwr_mat) <- rownames(amp_mat) <- chan_names
      
    } else {
      # Average time domain first, then FFT
      avg_data  <- apply(data_arr, c(1L, 2L), mean)   # [channels x samples]
      pwr_mat   <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      amp_mat   <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      phase_mat <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      
      for (i in seq_len(n_ch)) {
        res              <- .fft_one_sided(avg_data[i, ] * win_vec, sr, n_fft_use)
        pwr_mat[i, ]   <- res$power
        amp_mat[i, ]   <- res$amplitude
        phase_mat[i, ] <- res$phase
      }
      epoch_pwr <- NULL
    }
    
    return(structure(
      list(
        power         = pwr_mat,
        amplitude     = amp_mat,
        phase         = if (per_epoch) NULL else phase_mat,
        frequencies   = freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "fft",
        input_class   = input_class,
        n_fft         = n_fft_use,
        window_type   = window,
        per_epoch     = per_epoch,
        epoch_power   = if (per_epoch) epoch_pwr else NULL,
        n_epochs      = n_ep,
        conditions    = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_evoked ==========
  if (input_class == "eeg_evoked") {
    
    conditions <- input_obj$conditions
    sr         <- input_obj$sampling_rate
    chan_names <- input_obj$channels
    n_ch       <- length(chan_names)
    n_samp     <- ncol(input_obj$evoked[[1L]])
    n_fft_use  <- if (is.null(n_fft)) n_samp else as.integer(n_fft)
    n_half     <- floor(n_fft_use / 2L) + 1L
    freqs      <- seq(0, sr / 2, length.out = n_half)
    win_vec    <- .make_window(n_samp, window)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("FFT Spectral Analysis  [eeg_evoked]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Conditions:        ", paste(conditions, collapse = ", "), "\n")
      cat("  FFT length:        ", n_fft_use, "\n")
      cat("  Window:            ", window, "\n")
      cat("  Freq. resolution:  ", round(freqs[2L], 4), "Hz\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    pwr_list   <- list()
    amp_list   <- list()
    phase_list <- list()
    
    for (cond in conditions) {
      cond_data <- input_obj$evoked[[cond]]   # [channels x times]
      pwr_c     <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      amp_c     <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      phase_c   <- matrix(0, nrow = n_ch, ncol = n_half,
                          dimnames = list(chan_names, NULL))
      
      for (i in seq_len(n_ch)) {
        res            <- .fft_one_sided(cond_data[i, ] * win_vec, sr, n_fft_use)
        pwr_c[i, ]   <- res$power
        amp_c[i, ]   <- res$amplitude
        phase_c[i, ] <- res$phase
      }
      
      pwr_list[[cond]]   <- pwr_c
      amp_list[[cond]]   <- amp_c
      phase_list[[cond]] <- phase_c
    }
    
    return(structure(
      list(
        power         = pwr_list,
        amplitude     = amp_list,
        phase         = phase_list,
        frequencies   = freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "fft",
        input_class   = input_class,
        n_fft         = n_fft_use,
        window_type   = window,
        per_epoch     = FALSE,
        epoch_power   = NULL,
        conditions    = conditions
      ),
      class = "eeg_spectrum"
    ))
  }
}


# ============================================================================
#                    WELCH POWER SPECTRAL DENSITY  (eeg_psd_welch)
# ============================================================================

#' Compute Welch Power Spectral Density of EEG Data
#'
#' Estimates the Power Spectral Density (PSD) using the Welch averaged
#' periodogram method. The signal is divided into overlapping segments;
#' each segment is windowed and FFT'd, and the resulting power spectra are
#' averaged. This substantially reduces spectral variance compared to a
#' single FFT, making it the recommended method for resting-state EEG.
#'
#' Output is in uV^2/Hz, calibrated by the window power and sampling rate.
#'
#' For 'eeg_epochs' input the \code{per_epoch} flag controls whether the
#' Welch PSD is computed independently for each epoch (TRUE) or on the
#' trial-averaged waveform (FALSE, default).
#'
#' @param input_obj     An object of class 'eeg', 'eeg_epochs', or 'eeg_evoked'.
#' @param window_length Numeric. Segment length in seconds (default: 2.0).
#'   Determines frequency resolution: delta_f = 1 / window_length.
#'   Longer windows give finer resolution but fewer segments (more variance).
#'   If window_length exceeds epoch length, the whole epoch is used as one
#'   segment and a warning is issued.
#' @param overlap  Numeric in [0, 1). Fractional overlap between consecutive
#'   segments (default: 0.5 = 50%). Higher overlap smooths the PSD estimate
#'   but does not increase effective degrees of freedom proportionally.
#' @param window   Character. Spectral window: "hann" (default), "hamming",
#'   "blackman", or "none". The window is applied per-segment before FFT.
#' @param n_fft    Integer or NULL. FFT length per segment. NULL (default) uses
#'   the segment length. A value larger than the segment length zero-pads each
#'   segment, interpolating the spectrum (does not increase resolution).
#' @param per_epoch Logical. Only for 'eeg_epochs'. Compute Welch PSD per epoch
#'   then average, rather than averaging epochs in the time domain first.
#'   Default: FALSE.
#' @param verbose   Logical. Print processing summary. Default: TRUE.
#'
#' @return An object of class 'eeg_spectrum'. \code{$power} contains the PSD
#'   in uV^2/Hz. The \code{$phase} field is NULL (not meaningful for Welch).
#'
#' @examples
#' \dontrun{
#'   # Default: 2 s windows, 50% overlap, Hann window
#'   psd <- eeg_psd_welch(eeg)
#'
#'   # Finer resolution with 4 s windows
#'   psd <- eeg_psd_welch(eeg, window_length = 4, overlap = 0.75)
#'
#'   # Per-epoch PSD (induced oscillations)
#'   psd <- eeg_psd_welch(epochs, per_epoch = TRUE)
#' }
#'
#' @seealso \code{\link{eeg_fft}}, \code{\link{eeg_multitaper}},
#'   \code{\link{eeg_band_power}}
#' @export
eeg_psd_welch <- function(input_obj,
                          window_length = 2.0,
                          overlap       = 0.5,
                          window        = c("hann", "hamming", "blackman", "none"),
                          n_fft         = NULL,
                          per_epoch     = FALSE,
                          verbose       = TRUE) {
  
  window      <- match.arg(window)
  valid_cls   <- c("eeg", "eeg_epochs", "eeg_evoked")
  input_class <- class(input_obj)[1L]
  
  # ========== INPUT VALIDATION ==========
  if (!input_class %in% valid_cls) {
    stop("Input must be one of: ", paste(valid_cls, collapse = ", "),
         ". Got: '", input_class, "'.", call. = FALSE)
  }
  
  if (!is.numeric(window_length) || length(window_length) != 1 || window_length <= 0) {
    stop("window_length must be a single positive number (seconds).", call. = FALSE)
  }
  
  if (!is.numeric(overlap) || length(overlap) != 1 || overlap < 0 || overlap >= 1) {
    stop("overlap must be a single number in [0, 1).", call. = FALSE)
  }
  
  # ========== INTERNAL WELCH WORKER ==========
  # Computes Welch PSD for a single channel time series x.
  # Returns numeric vector of length n_half = floor(n_fft_w/2) + 1.
  .welch_channel <- function(x, sr, win_samp, step, win_vec, n_fft_w) {
    n        <- length(x)
    n_half   <- floor(n_fft_w / 2L) + 1L
    win_norm <- sum(win_vec^2) * sr   # = Fs * sum(w^2), normalises to uV^2/Hz
    
    starts <- seq(1L, n - win_samp + 1L, by = step)
    if (length(starts) < 1L) {
      stop("Signal too short for window_length = ", window_length,
           " s. Reduce window_length or shorten overlap.", call. = FALSE)
    }
    
    acc <- numeric(n_half)
    
    for (s in starts) {
      seg <- x[s:(s + win_samp - 1L)] * win_vec
      if (n_fft_w > win_samp) seg <- c(seg, rep(0.0, n_fft_w - win_samp))
      
      X   <- fft(seg)[seq_len(n_half)]
      p   <- (Mod(X)^2) / win_norm
      # Double one-sided bins; leave DC and Nyquist as-is
      p[seq(2L, n_half - 1L)] <- 2.0 * p[seq(2L, n_half - 1L)]
      acc <- acc + p
    }
    
    acc / length(starts)
  }
  
  # ========== SHARED SETUP HELPER ==========
  .welch_setup <- function(sr, n_samp) {
    win_samp <- round(window_length * sr)
    if (win_samp > n_samp) {
      warning("window_length (", window_length, " s = ", win_samp, " samples) ",
              "exceeds signal length (", n_samp, " samples). ",
              "Using full signal as a single segment.", call. = FALSE)
      win_samp <- n_samp
    }
    step     <- max(1L, round(win_samp * (1 - overlap)))
    n_fft_w  <- if (is.null(n_fft)) win_samp else max(as.integer(n_fft), win_samp)
    win_vec  <- .make_window(win_samp, window)
    n_half   <- floor(n_fft_w / 2L) + 1L
    freqs    <- seq(0, sr / 2, length.out = n_half)
    list(win_samp = win_samp, step = step, n_fft_w = n_fft_w,
         win_vec = win_vec, n_half = n_half, freqs = freqs)
  }
  
  # ========== DISPATCH: eeg ==========
  if (input_class == "eeg") {
    
    sr         <- input_obj$sampling_rate
    data_mat   <- input_obj$data
    chan_names <- input_obj$channels
    n_ch       <- nrow(data_mat)
    n_tp       <- ncol(data_mat)
    s          <- .welch_setup(sr, n_tp)
    n_seg      <- floor((n_tp - s$win_samp) / s$step) + 1L
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Welch PSD  [eeg]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Signal length:     ", n_tp, "samples (", round(n_tp / sr, 2), "s)\n")
      cat("  Window length:     ", window_length, "s (", s$win_samp, "samples)\n")
      cat("  Overlap:           ", overlap * 100, "%\n")
      cat("  Segments:          ", n_seg, "\n")
      cat("  FFT length:        ", s$n_fft_w, "\n")
      cat("  Freq. resolution:  ", round(sr / s$n_fft_w, 4), "Hz\n")
      cat("  Window:            ", window, "\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    psd_mat <- matrix(0, nrow = n_ch, ncol = s$n_half,
                      dimnames = list(chan_names, NULL))
    for (i in seq_len(n_ch)) {
      psd_mat[i, ] <- .welch_channel(data_mat[i, ], sr,
                                     s$win_samp, s$step, s$win_vec, s$n_fft_w)
    }
    
    return(structure(
      list(
        power         = psd_mat,
        amplitude     = sqrt(psd_mat),
        phase         = NULL,
        frequencies   = s$freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "welch",
        input_class   = input_class,
        n_fft         = s$n_fft_w,
        window_type   = window,
        window_length = window_length,
        overlap       = overlap,
        n_segments    = n_seg,
        per_epoch     = FALSE,
        epoch_power   = NULL,
        conditions    = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_epochs ==========
  if (input_class == "eeg_epochs") {
    
    if (is.null(input_obj$data)) {
      stop("Epoch data not loaded. Set preload = TRUE in epoch_eeg().", call. = FALSE)
    }
    
    sr         <- input_obj$sampling_rate
    data_arr   <- input_obj$data    # [channels x samples x epochs]
    chan_names <- input_obj$channels
    n_ch       <- dim(data_arr)[1L]
    n_samp     <- dim(data_arr)[2L]
    n_ep       <- dim(data_arr)[3L]
    s          <- .welch_setup(sr, n_samp)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Welch PSD  [eeg_epochs]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Epochs:            ", n_ep, "\n")
      cat("  Samples per epoch: ", n_samp, "\n")
      cat("  Window length:     ", window_length, "s (", s$win_samp, "samples)\n")
      cat("  Overlap:           ", overlap * 100, "%\n")
      cat("  Per epoch:         ", per_epoch, "\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    if (per_epoch) {
      epoch_pwr <- array(0, dim = c(n_ch, s$n_half, n_ep))
      dimnames(epoch_pwr)[[1L]] <- chan_names
      
      for (ep in seq_len(n_ep)) {
        for (i in seq_len(n_ch)) {
          epoch_pwr[i, , ep] <- .welch_channel(data_arr[i, , ep], sr,
                                               s$win_samp, s$step,
                                               s$win_vec, s$n_fft_w)
        }
      }
      
      psd_mat <- apply(epoch_pwr, c(1L, 2L), mean)
      rownames(psd_mat) <- chan_names
      
    } else {
      avg_data <- apply(data_arr, c(1L, 2L), mean)   # [channels x samples]
      psd_mat  <- matrix(0, nrow = n_ch, ncol = s$n_half,
                         dimnames = list(chan_names, NULL))
      for (i in seq_len(n_ch)) {
        psd_mat[i, ] <- .welch_channel(avg_data[i, ], sr,
                                       s$win_samp, s$step, s$win_vec, s$n_fft_w)
      }
      epoch_pwr <- NULL
    }
    
    return(structure(
      list(
        power         = psd_mat,
        amplitude     = sqrt(psd_mat),
        phase         = NULL,
        frequencies   = s$freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "welch",
        input_class   = input_class,
        n_fft         = s$n_fft_w,
        window_type   = window,
        window_length = window_length,
        overlap       = overlap,
        n_segments    = floor((n_samp - s$win_samp) / s$step) + 1L,
        per_epoch     = per_epoch,
        epoch_power   = if (per_epoch) epoch_pwr else NULL,
        n_epochs      = n_ep,
        conditions    = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_evoked ==========
  if (input_class == "eeg_evoked") {
    
    conditions <- input_obj$conditions
    sr         <- input_obj$sampling_rate
    chan_names <- input_obj$channels
    n_ch       <- length(chan_names)
    n_samp     <- ncol(input_obj$evoked[[1L]])
    s          <- .welch_setup(sr, n_samp)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Welch PSD  [eeg_evoked]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:    ", n_ch, "\n")
      cat("  Conditions:  ", paste(conditions, collapse = ", "), "\n")
      cat("  Window:      ", window_length, "s, overlap", overlap * 100, "%\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    psd_list <- list()
    for (cond in conditions) {
      cond_data   <- input_obj$evoked[[cond]]   # [channels x times]
      psd_c       <- matrix(0, nrow = n_ch, ncol = s$n_half,
                            dimnames = list(chan_names, NULL))
      for (i in seq_len(n_ch)) {
        psd_c[i, ] <- .welch_channel(cond_data[i, ], sr,
                                     s$win_samp, s$step, s$win_vec, s$n_fft_w)
      }
      psd_list[[cond]] <- psd_c
    }
    
    return(structure(
      list(
        power         = psd_list,
        amplitude     = lapply(psd_list, sqrt),
        phase         = NULL,
        frequencies   = s$freqs,
        channels      = chan_names,
        sampling_rate = sr,
        method        = "welch",
        input_class   = input_class,
        n_fft         = s$n_fft_w,
        window_type   = window,
        window_length = window_length,
        overlap       = overlap,
        n_segments    = floor((n_samp - s$win_samp) / s$step) + 1L,
        per_epoch     = FALSE,
        epoch_power   = NULL,
        conditions    = conditions
      ),
      class = "eeg_spectrum"
    ))
  }
}


# ============================================================================
#                   MULTITAPER PSD  (eeg_multitaper)
# ============================================================================

#' Compute Multitaper Power Spectral Density of EEG Data
#'
#' Estimates the PSD using the multitaper method (Thomson, 1982). Multiple
#' orthogonal Discrete Prolate Spheroidal Sequences (DPSS / Slepian tapers)
#' are applied to the signal; the resulting periodograms are averaged. This
#' approach simultaneously controls spectral leakage (bias) and reduces
#' variance without the frequency-resolution loss of simple averaging.
#'
#' The key parameter is \code{time_bandwidth} (NW): larger values give more
#' spectral smoothing but coarser frequency resolution. For EEG, NW = 4
#' (with k = 7 tapers) is a common default. The spectral half-bandwidth is
#' W = NW / (n * Fs) Hz, controlling the minimum resolvable frequency
#' separation.
#'
#' For 'eeg' input the signal is divided into segments of length
#' \code{segment_length} seconds (with \code{overlap}) and multitaper PSD
#' is computed per segment then averaged. For 'eeg_epochs' input multitaper
#' PSD is computed on the full epoch (no segmentation).
#'
#' @param input_obj      An object of class 'eeg', 'eeg_epochs', or 'eeg_evoked'.
#' @param time_bandwidth Numeric. Time-bandwidth product NW (default: 4).
#'   Controls bias-variance trade-off. Larger values = more smoothing.
#'   Typical EEG range: 2-8.
#' @param n_tapers       Integer or NULL. Number of tapers. NULL (default)
#'   uses 2 * time_bandwidth - 1, the maximum that maintains good spectral
#'   concentration. Must be >= 1 and <= 2 * time_bandwidth - 1.
#' @param segment_length Numeric. Segment length in seconds for 'eeg' input
#'   (default: 2.0). Ignored for 'eeg_epochs' and 'eeg_evoked'.
#' @param overlap        Numeric in [0, 1). Fractional overlap between segments
#'   for 'eeg' input (default: 0.5). Ignored for epochs/evoked.
#' @param per_epoch      Logical. For 'eeg_epochs': compute multitaper PSD per
#'   epoch then average (TRUE), or average epochs first (FALSE). Default: FALSE.
#' @param verbose        Logical. Print processing summary. Default: TRUE.
#'
#' @return An object of class 'eeg_spectrum'. Output units are uV^2/Hz.
#'   \code{$phase} is NULL. Additional fields: \code{$time_bandwidth},
#'   \code{$n_tapers}, \code{$segment_length} (for 'eeg' input).
#'
#' @references
#'   Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
#'   \emph{Proc. IEEE}, 70, 1055--1096.
#'
#'   Slepian, D. (1978). Prolate spheroidal wave functions, Fourier analysis,
#'   and uncertainty V. \emph{Bell Syst. Tech. J.}, 57, 1371--1430.
#'
#' @examples
#' \dontrun{
#'   # Default: NW = 4, k = 7 tapers
#'   psd <- eeg_multitaper(eeg)
#'
#'   # More smoothing (wider half-bandwidth)
#'   psd <- eeg_multitaper(eeg, time_bandwidth = 6)
#'
#'   # Epochs: per-epoch multitaper PSD
#'   psd <- eeg_multitaper(epochs, per_epoch = TRUE)
#' }
#'
#' @seealso \code{\link{eeg_fft}}, \code{\link{eeg_psd_welch}},
#'   \code{\link{eeg_band_power}}
#' @export
eeg_multitaper <- function(input_obj,
                           time_bandwidth = 4,
                           n_tapers       = NULL,
                           segment_length = 2.0,
                           overlap        = 0.5,
                           per_epoch      = FALSE,
                           verbose        = TRUE) {
  
  valid_cls   <- c("eeg", "eeg_epochs", "eeg_evoked")
  input_class <- class(input_obj)[1L]
  
  # ========== INPUT VALIDATION ==========
  if (!input_class %in% valid_cls) {
    stop("Input must be one of: ", paste(valid_cls, collapse = ", "),
         ". Got: '", input_class, "'.", call. = FALSE)
  }
  
  if (!is.numeric(time_bandwidth) || time_bandwidth <= 0) {
    stop("time_bandwidth must be a positive number (e.g. 4).", call. = FALSE)
  }
  
  if (is.null(n_tapers)) {
    n_tapers <- max(1L, floor(2 * time_bandwidth - 1))
  }
  
  n_tapers <- as.integer(n_tapers)
  if (n_tapers < 1L) {
    stop("n_tapers must be a positive integer.", call. = FALSE)
  }
  
  max_tapers <- floor(2 * time_bandwidth - 1)
  if (n_tapers > max_tapers) {
    warning("n_tapers (", n_tapers, ") exceeds the recommended maximum ",
            "(2 * time_bandwidth - 1 = ", max_tapers, "). ",
            "Excess tapers have poor spectral concentration.", call. = FALSE)
  }
  
  if (!is.numeric(overlap) || overlap < 0 || overlap >= 1) {
    stop("overlap must be a number in [0, 1).", call. = FALSE)
  }
  
  # ========== DISPATCH: eeg ==========
  if (input_class == "eeg") {
    
    sr         <- input_obj$sampling_rate
    data_mat   <- input_obj$data
    chan_names <- input_obj$channels
    n_ch       <- nrow(data_mat)
    n_tp       <- ncol(data_mat)
    
    win_samp <- min(round(segment_length * sr), n_tp)
    step     <- max(1L, round(win_samp * (1 - overlap)))
    starts   <- seq(1L, n_tp - win_samp + 1L, by = step)
    n_half   <- floor(win_samp / 2L) + 1L
    freqs    <- seq(0, sr / 2, length.out = n_half)
    tapers   <- .compute_dpss(win_samp, time_bandwidth, n_tapers)
    h_bw     <- round(time_bandwidth / win_samp * sr, 3)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Multitaper PSD  [eeg]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Signal length:     ", n_tp, "samples (", round(n_tp / sr, 2), "s)\n")
      cat("  Time-bandwidth NW: ", time_bandwidth, "\n")
      cat("  Tapers (K):        ", n_tapers, "\n")
      cat("  Spectral half-bw:  +-", h_bw, "Hz\n")
      cat("  Segment length:    ", segment_length, "s (", win_samp, "samples)\n")
      cat("  Segments:          ", length(starts), "\n")
      cat("  Freq. resolution:  ", round(sr / win_samp, 4), "Hz\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    psd_mat <- matrix(0, nrow = n_ch, ncol = n_half,
                      dimnames = list(chan_names, NULL))
    
    for (i in seq_len(n_ch)) {
      seg_acc <- numeric(n_half)
      for (s_start in starts) {
        seg      <- data_mat[i, s_start:(s_start + win_samp - 1L)]
        seg_acc  <- seg_acc + .mtaper_psd(seg, sr, tapers)
      }
      psd_mat[i, ] <- seg_acc / length(starts)
    }
    
    return(structure(
      list(
        power          = psd_mat,
        amplitude      = sqrt(psd_mat),
        phase          = NULL,
        frequencies    = freqs,
        channels       = chan_names,
        sampling_rate  = sr,
        method         = "multitaper",
        input_class    = input_class,
        n_fft          = win_samp,
        window_type    = "dpss",
        time_bandwidth = time_bandwidth,
        n_tapers       = n_tapers,
        segment_length = segment_length,
        overlap        = overlap,
        n_segments     = length(starts),
        per_epoch      = FALSE,
        epoch_power    = NULL,
        conditions     = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_epochs ==========
  if (input_class == "eeg_epochs") {
    
    if (is.null(input_obj$data)) {
      stop("Epoch data not loaded. Set preload = TRUE in epoch_eeg().", call. = FALSE)
    }
    
    sr         <- input_obj$sampling_rate
    data_arr   <- input_obj$data
    chan_names <- input_obj$channels
    n_ch       <- dim(data_arr)[1L]
    n_samp     <- dim(data_arr)[2L]
    n_ep       <- dim(data_arr)[3L]
    n_half     <- floor(n_samp / 2L) + 1L
    freqs      <- seq(0, sr / 2, length.out = n_half)
    tapers     <- .compute_dpss(n_samp, time_bandwidth, n_tapers)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Multitaper PSD  [eeg_epochs]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Epochs:            ", n_ep, "\n")
      cat("  Samples per epoch: ", n_samp, "\n")
      cat("  Time-bandwidth NW: ", time_bandwidth, "\n")
      cat("  Tapers (K):        ", n_tapers, "\n")
      cat("  Per epoch:         ", per_epoch, "\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    if (per_epoch) {
      epoch_pwr <- array(0, dim = c(n_ch, n_half, n_ep))
      dimnames(epoch_pwr)[[1L]] <- chan_names
      
      for (ep in seq_len(n_ep)) {
        for (i in seq_len(n_ch)) {
          epoch_pwr[i, , ep] <- .mtaper_psd(data_arr[i, , ep], sr, tapers)
        }
      }
      
      psd_mat <- apply(epoch_pwr, c(1L, 2L), mean)
      rownames(psd_mat) <- chan_names
      
    } else {
      avg_data <- apply(data_arr, c(1L, 2L), mean)
      psd_mat  <- matrix(0, nrow = n_ch, ncol = n_half,
                         dimnames = list(chan_names, NULL))
      for (i in seq_len(n_ch)) {
        psd_mat[i, ] <- .mtaper_psd(avg_data[i, ], sr, tapers)
      }
      epoch_pwr <- NULL
    }
    
    return(structure(
      list(
        power          = psd_mat,
        amplitude      = sqrt(psd_mat),
        phase          = NULL,
        frequencies    = freqs,
        channels       = chan_names,
        sampling_rate  = sr,
        method         = "multitaper",
        input_class    = input_class,
        n_fft          = n_samp,
        window_type    = "dpss",
        time_bandwidth = time_bandwidth,
        n_tapers       = n_tapers,
        per_epoch      = per_epoch,
        epoch_power    = if (per_epoch) epoch_pwr else NULL,
        n_epochs       = n_ep,
        conditions     = NULL
      ),
      class = "eeg_spectrum"
    ))
  }
  
  # ========== DISPATCH: eeg_evoked ==========
  if (input_class == "eeg_evoked") {
    
    conditions <- input_obj$conditions
    sr         <- input_obj$sampling_rate
    chan_names <- input_obj$channels
    n_ch       <- length(chan_names)
    n_samp     <- ncol(input_obj$evoked[[1L]])
    n_half     <- floor(n_samp / 2L) + 1L
    freqs      <- seq(0, sr / 2, length.out = n_half)
    tapers     <- .compute_dpss(n_samp, time_bandwidth, n_tapers)
    
    if (verbose) {
      cat("\n", strrep("=", 60), "\n", sep = "")
      cat("Multitaper PSD  [eeg_evoked]\n")
      cat(strrep("=", 60), "\n")
      cat("  Channels:          ", n_ch, "\n")
      cat("  Conditions:        ", paste(conditions, collapse = ", "), "\n")
      cat("  Time-bandwidth NW: ", time_bandwidth, "\n")
      cat("  Tapers (K):        ", n_tapers, "\n")
      cat(strrep("=", 60), "\n\n")
    }
    
    psd_list <- list()
    for (cond in conditions) {
      cond_data   <- input_obj$evoked[[cond]]
      psd_c       <- matrix(0, nrow = n_ch, ncol = n_half,
                            dimnames = list(chan_names, NULL))
      for (i in seq_len(n_ch)) {
        psd_c[i, ] <- .mtaper_psd(cond_data[i, ], sr, tapers)
      }
      psd_list[[cond]] <- psd_c
    }
    
    return(structure(
      list(
        power          = psd_list,
        amplitude      = lapply(psd_list, sqrt),
        phase          = NULL,
        frequencies    = freqs,
        channels       = chan_names,
        sampling_rate  = sr,
        method         = "multitaper",
        input_class    = input_class,
        n_fft          = n_samp,
        window_type    = "dpss",
        time_bandwidth = time_bandwidth,
        n_tapers       = n_tapers,
        per_epoch      = FALSE,
        epoch_power    = NULL,
        conditions     = conditions
      ),
      class = "eeg_spectrum"
    ))
  }
}


# ============================================================================
#                        BAND POWER  (eeg_band_power)
# ============================================================================

#' Extract EEG Frequency Band Power
#'
#' Computes the power in canonical EEG frequency bands (delta, theta, alpha,
#' beta, gamma) using trapezoidal integration over the PSD. Accepts either
#' a pre-computed \code{eeg_spectrum} object or a raw EEG object (in which
#' case a Welch PSD is computed internally).
#'
#' Default frequency bands follow standard EEG conventions:
#' \itemize{
#'   \item delta: 0.5-4 Hz
#'   \item theta: 4-8 Hz
#'   \item alpha: 8-13 Hz
#'   \item beta:  13-30 Hz
#'   \item gamma: 30-100 Hz
#' }
#' Custom bands can be supplied as a named list of \code{c(low_hz, high_hz)} pairs.
#'
#' @param input_obj An \code{eeg_spectrum} object (output of \code{eeg_fft()},
#'   \code{eeg_psd_welch()}, or \code{eeg_multitaper()}), or a raw 'eeg',
#'   'eeg_epochs', or 'eeg_evoked' object.
#' @param bands A named list where each element is \code{c(low_hz, high_hz)}.
#'   Defaults to the five canonical EEG bands. Custom example:
#'   \code{list(low_alpha = c(8, 10), high_alpha = c(10, 13))}.
#' @param method  Character. Method used when \code{input_obj} is a raw EEG
#'   object: "welch" (default), "fft", or "multitaper". Ignored when
#'   \code{input_obj} is already an \code{eeg_spectrum}.
#' @param relative Logical. If TRUE, return relative band power (each band
#'   divided by the sum of all requested bands). Default: FALSE (absolute
#'   power in uV^2).
#' @param verbose  Logical. Print a band power summary table. Default: TRUE.
#' @param ...      Additional arguments passed to the underlying spectral
#'   function when \code{input_obj} is a raw EEG object (e.g. window_length).
#'
#' @return For 'eeg' and 'eeg_epochs' input (or an \code{eeg_spectrum} derived
#'   from them): a data frame with one row per channel and one column per band.
#'   For 'eeg_evoked' input: a named list of such data frames, one per condition.
#'
#' @examples
#' \dontrun{
#'   # From raw EEG (Welch PSD computed internally)
#'   bp <- eeg_band_power(eeg)
#'
#'   # From a pre-computed spectrum
#'   psd  <- eeg_psd_welch(eeg)
#'   bp   <- eeg_band_power(psd)
#'
#'   # Relative band power
#'   bp_rel <- eeg_band_power(eeg, relative = TRUE)
#'
#'   # Custom bands
#'   bp_custom <- eeg_band_power(eeg,
#'     bands = list(low_alpha = c(8, 10), high_alpha = c(10, 13)))
#' }
#'
#' @seealso \code{\link{eeg_fft}}, \code{\link{eeg_psd_welch}},
#'   \code{\link{eeg_multitaper}}
#' @export
eeg_band_power <- function(input_obj,
                           bands = list(
                             delta = c(0.5,  4),
                             theta = c(4,    8),
                             alpha = c(8,   13),
                             beta  = c(13,  30),
                             gamma = c(30, 100)
                           ),
                           method   = c("welch", "fft", "multitaper"),
                           relative = FALSE,
                           verbose  = TRUE,
                           ...) {
  
  method      <- match.arg(method)
  valid_raw   <- c("eeg", "eeg_epochs", "eeg_evoked")
  input_class <- class(input_obj)[1L]
  
  # ========== INPUT VALIDATION ==========
  if (!input_class %in% c("eeg_spectrum", valid_raw)) {
    stop("Input must be an 'eeg_spectrum' object or one of: ",
         paste(valid_raw, collapse = ", "), ".", call. = FALSE)
  }
  
  if (!is.list(bands) || length(bands) == 0) {
    stop("bands must be a non-empty named list of c(low_hz, high_hz) pairs.",
         call. = FALSE)
  }
  
  for (nm in names(bands)) {
    b <- bands[[nm]]
    if (length(b) != 2L || !is.numeric(b) || b[1L] >= b[2L] || any(b < 0)) {
      stop("Band '", nm, "' must be c(low_hz, high_hz) with 0 <= low < high.",
           call. = FALSE)
    }
  }
  
  # ========== COMPUTE SPECTRUM IF GIVEN RAW INPUT ==========
  if (input_class != "eeg_spectrum") {
    if (verbose) cat("Computing", method, "PSD for band power extraction...\n")
    spectrum_obj <- switch(method,
                           welch      = eeg_psd_welch(input_obj, verbose = verbose, ...),
                           fft        = eeg_fft(input_obj, verbose = verbose, ...),
                           multitaper = eeg_multitaper(input_obj, verbose = verbose, ...)
    )
  } else {
    spectrum_obj <- input_obj
  }
  
  freqs      <- spectrum_obj$frequencies
  band_names <- names(bands)
  
  # ========== TRAPEZOIDAL INTEGRATION HELPER ==========
  .integrate_band <- function(pwr_vec, f_low, f_high) {
    idx <- which(freqs >= f_low & freqs <= f_high)
    if (length(idx) < 2L) return(NA_real_)
    # trapz: sum of trapezoids with base = freq_step, heights = adjacent PSD values
    f_sub <- freqs[idx]
    p_sub <- pwr_vec[idx]
    sum(diff(f_sub) * (p_sub[-length(p_sub)] + p_sub[-1L]) / 2)
  }
  
  # ========== BUILD BAND POWER DATA FRAME ==========
  .make_bp_df <- function(pwr_mat) {
    # pwr_mat: [channels x frequencies]
    bp <- data.frame(channel = rownames(pwr_mat), stringsAsFactors = FALSE)
    for (nm in band_names) {
      b       <- bands[[nm]]
      bp[[nm]] <- apply(pwr_mat, 1L, .integrate_band,
                        f_low = b[1L], f_high = b[2L])
    }
    
    if (relative) {
      totals <- rowSums(bp[, band_names, drop = FALSE], na.rm = TRUE)
      for (nm in band_names) {
        bp[[nm]] <- bp[[nm]] / totals
      }
    }
    bp
  }
  
  # ========== EVOKED: one data frame per condition ==========
  if (spectrum_obj$input_class == "eeg_evoked") {
    result <- lapply(spectrum_obj$conditions, function(cond) {
      .make_bp_df(spectrum_obj$power[[cond]])
    })
    names(result) <- spectrum_obj$conditions
    return(result)
  }
  
  # ========== EEG / EPOCHS: single data frame ==========
  bp_df <- .make_bp_df(spectrum_obj$power)
  
  if (verbose) {
    units <- if (relative) "(relative)" else "(uV^2)"
    cat("\nBand Power Summary ", units, "\n", sep = "")
    cat(strrep("-", 60), "\n")
    print(bp_df, digits = 4)
    cat("\n")
  }
  
  invisible(bp_df)
}


# ============================================================================
#                     S3 PRINT METHOD  (print.eeg_spectrum)
# ============================================================================

#' Print Method for eeg_spectrum Objects
#'
#' Displays a formatted summary of an \code{eeg_spectrum} object, including
#' method parameters, frequency range, and spectral resolution.
#'
#' @param x   An object of class 'eeg_spectrum'.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns x.
#' @export
print.eeg_spectrum <- function(x, ...) {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("EEG Spectrum Object\n")
  cat(strrep("=", 60), "\n\n")
  
  cat("Method:            ", toupper(x$method), "\n")
  cat("Input class:       ", x$input_class, "\n")
  cat("Channels:          ", length(x$channels), "\n")
  
  n_f     <- length(x$frequencies)
  f_res   <- if (n_f > 1) round(x$frequencies[2L] - x$frequencies[1L], 4) else NA
  cat("Frequency bins:    ", n_f, "\n")
  cat("Frequency range:   ", round(min(x$frequencies), 2), "-",
      round(max(x$frequencies), 2), "Hz\n")
  cat("Freq. resolution:  ", f_res, "Hz\n")
  cat("Sampling rate:     ", x$sampling_rate, "Hz\n")
  
  if (x$method == "welch") {
    cat("\nWelch parameters:\n")
    cat("  Window length:   ", x$window_length, "s\n")
    cat("  Overlap:         ", x$overlap * 100, "%\n")
    cat("  Segments:        ", x$n_segments, "\n")
    cat("  Window function: ", x$window_type, "\n")
    cat("  FFT length:      ", x$n_fft, "\n")
  } else if (x$method == "multitaper") {
    h_bw <- round(x$time_bandwidth / x$n_fft * x$sampling_rate, 3)
    cat("\nMultitaper parameters:\n")
    cat("  Time-bandwidth:  ", x$time_bandwidth, "\n")
    cat("  Tapers (K):      ", x$n_tapers, "\n")
    cat("  Half-bandwidth:  +-", h_bw, "Hz\n")
    cat("  FFT length:      ", x$n_fft, "\n")
  } else if (x$method == "fft") {
    cat("\nFFT parameters:\n")
    cat("  FFT length:      ", x$n_fft, "\n")
    cat("  Window function: ", x$window_type, "\n")
    cat("  Phase stored:    ", !is.null(x$phase), "\n")
  }
  
  if (!is.null(x$conditions)) {
    cat("\nConditions:        ", paste(x$conditions, collapse = ", "), "\n")
  }
  
  if (!is.null(x$n_epochs)) {
    cat("Epochs:            ", x$n_epochs, "\n")
    cat("Per-epoch spectra: ", x$per_epoch, "\n")
  }
  
  cat("\n")
  cat(strrep("=", 60), "\n\n")
  invisible(x)
}

# ============================================================================
# END OF FILE
# ============================================================================