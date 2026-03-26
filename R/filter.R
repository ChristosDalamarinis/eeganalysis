#' ============================================================================
#'                          EEG Filtering Functions
#' ============================================================================
#'
#' This module provides functions for filtering EEG data to isolate frequency
#' bands of interest and remove noise. Filtering is a critical preprocessing
#' step that removes low-frequency drifts, high-frequency noise, and line
#' interference (50/60 Hz) before further analysis.
#'
#' The module implements Butterworth filters (highpass, lowpass, bandpass, and
#' notch) with zero-phase forward-backward filtering to preserve temporal
#' accuracy of neural events. Filter parameters follow standard EEG analysis
#' conventions and integrate seamlessly with the eeganalysis pipeline.
#'
#' Author: Christos Dalamarinis
#' Date: 2026
#' Status: In development!!!!!!! New filtering algorithm is on the making, will
#'         deemed legacy script soon.
#' ============================================================================
#'
#' Apply a Bandpass Filter to EEG Data
#'
#' Filters EEG data by attenuating frequencies outside a specified passband.
#' Supports highpass-only, lowpass-only, or true bandpass filtering using a
#' zero-phase Butterworth filter applied with forward-backward (filtfilt)
#' convolution, which preserves the temporal alignment of neural events.
#'
#' @param eeg_obj An object of class 'eeg' containing EEG data.
#'
#' @param low_freq Numeric. Lower edge of the passband in Hz (highpass cutoff).
#'                 Frequencies below this value are attenuated.
#'                 Set to \code{NULL} to skip highpass filtering and apply
#'                 lowpass only.
#'                 Common values: 0.1 (DC drift removal), 1 (ERP), 0.5 (default).
#'                 Default: 0.5
#'
#' @param high_freq Numeric. Upper edge of the passband in Hz (lowpass cutoff).
#'                  Frequencies above this value are attenuated.
#'                  Set to \code{NULL} to skip lowpass filtering and apply
#'                  highpass only.
#'                  Must be lower than the Nyquist frequency
#'                  (sampling_rate / 2).
#'                  Common values: 40 (ERP), 100 (time-frequency), 30 (default).
#'                  Default: 30
#'
#' @param filter_order Integer. Order of the Butterworth filter.
#'                     Higher values create a steeper roll-off but increase
#'                     computation time and risk of edge artefacts.
#'                     For bandpass, the effective order is doubled
#'                     (highpass + lowpass applied sequentially).
#'                     Typical range: 2–8. Default: 4.
#'
#' @param channels Character or integer vector. Channel names or indices to
#'                 filter. If \code{NULL} (default), all channels are filtered.
#'                 Useful for excluding external channels (e.g., EOG, EMG)
#'                 from filtering while keeping them in the object.
#'                 Example: \code{c("Cz", "Pz", "Oz")} or \code{1:64}.
#'
#' @param zero_phase Logical. If \code{TRUE} (default), uses forward-backward
#'                   filtering (\code{filtfilt}) which produces zero phase
#'                   distortion — essential for ERP latency measurements.
#'                   If \code{FALSE}, uses one-pass causal filtering
#'                   (\code{filter}), which is faster but introduces phase
#'                   delay. Only use \code{FALSE} for exploratory inspection.
#'
#' @param method Character string specifying the filter design method:
#'               \describe{
#'                 \item{"butter"}{(default) Butterworth filter — maximally
#'                                 flat passband, no ripple. Recommended for
#'                                 all standard EEG analyses.}
#'                 \item{"fir"}{Finite Impulse Response filter — linear phase
#'                              by design, no phase distortion even in
#'                              one-pass mode. Computationally heavier.
#'                              Requires the \code{signal} package.}
#'               }
#'
#' @param padding Integer. Number of samples to mirror-pad at each end of
#'                each channel's time series before filtering to reduce
#'                edge artefacts caused by filter ringing.
#'                Set to \code{0} to disable padding.
#'                Default: \code{3 * filter_order} (automatically scaled
#'                to filter length).
#'
#' @param notch_freq Numeric. Frequency in Hz for an additional notch (band-stop)
#'                   filter to remove line noise. Typically \code{50} (Europe)
#'                   or \code{60} (North America).
#'                   Set to \code{NULL} (default) to skip notch filtering.
#'                   Applied after bandpass filtering.
#'
#' @param notch_bandwidth Numeric. Bandwidth in Hz around \code{notch_freq}
#'                        to attenuate. Default: 2 (i.e., notch spans
#'                        \code{notch_freq ± 1 Hz}).
#'                        Only used when \code{notch_freq} is not \code{NULL}.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a summary of
#'                filter parameters and processing progress.
#'                Set to \code{FALSE} for silent batch processing.
#'
#' @return An object of class 'eeg' with the same structure as the input,
#'   with the following fields updated:
#'   \describe{
#'     \item{data}{Filtered signal matrix (channels x time points).
#'                 Unfiltered channels (if \code{channels} was specified)
#'                 retain their original values.}
#'     \item{preprocessing_history}{Appended with a string describing the
#'                                   filter parameters applied.}
#'     \item{channels}{Unchanged.}
#'     \item{sampling_rate}{Unchanged.}
#'     \item{times}{Unchanged.}
#'     \item{events}{Unchanged.}
#'     \item{metadata}{Unchanged.}
#'     \item{reference}{Unchanged.}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates all input parameters and checks Nyquist constraints.
#'   \item Resolves target channel indices (all or subset).
#'   \item Mirror-pads each channel to reduce edge artefacts.
#'   \item Designs and applies the Butterworth highpass filter (if
#'         \code{low_freq} is not NULL).
#'   \item Designs and applies the Butterworth lowpass filter (if
#'         \code{high_freq} is not NULL).
#'   \item Applies a notch filter (if \code{notch_freq} is not NULL).
#'   \item Removes padding and writes filtered data back to the matrix.
#'   \item Appends a preprocessing history entry.
#' }
#'
#' \strong{Why zero-phase filtering matters:}
#' Standard causal filters introduce phase shifts that smear the timing of
#' neural responses. Forward-backward filtering cancels phase distortion,
#' making it the standard method for ERP research. Always use the default
#' \code{zero_phase = TRUE} when peak latency is of interest.
#'
#' \strong{Recommended filter settings by paradigm:}
#' \itemize{
#'   \item ERP (e.g., P300, N170): \code{low_freq = 0.1}, \code{high_freq = 40}
#'   \item Oscillations / ERSP: \code{low_freq = 1}, \code{high_freq = 100}
#'   \item Slow cortical potentials: \code{low_freq = NULL}, \code{high_freq = 10}
#'   \item DC drift removal only: \code{low_freq = 0.5}, \code{high_freq = NULL}
#'   \item SSVEP / RIFT: \code{low_freq = 1}, \code{high_freq = 120}
#' }
#'
#' \strong{Edge artefacts:}
#' Filtering always produces distortion at the very start and end of the
#' signal. The \code{padding} argument mitigates this by reflecting the
#' signal before filtering and discarding the mirrored segments afterwards.
#' Increase padding when using high filter orders or very low cutoffs.
#'
#' @examples
#' \dontrun{
#'   # Standard ERP bandpass (0.1 – 40 Hz)
#'   eeg_filt <- eeg_bandpass(eeg, low_freq = 0.1, high_freq = 40)
#'
#'   # Highpass only — remove DC drift
#'   eeg_filt <- eeg_bandpass(eeg, low_freq = 0.5, high_freq = NULL)
#'
#'   # Lowpass only — smooth signal
#'   eeg_filt <- eeg_bandpass(eeg, low_freq = NULL, high_freq = 40)
#'
#'   # Bandpass + 50 Hz notch for European line noise
#'   eeg_filt <- eeg_bandpass(eeg,
#'                            low_freq       = 0.1,
#'                            high_freq      = 100,
#'                            notch_freq     = 50,
#'                            notch_bandwidth = 2)
#'
#'   # Filter only EEG channels, leave EOG untouched
#'   eeg_filt <- eeg_bandpass(eeg,
#'                            low_freq  = 0.1,
#'                            high_freq = 40,
#'                            channels  = 1:64)
#'
#'   # Steeper roll-off for oscillation analysis
#'   eeg_filt <- eeg_bandpass(eeg,
#'                            low_freq     = 1,
#'                            high_freq    = 100,
#'                            filter_order = 8)
#'
#'   # Silent batch processing
#'   eeg_filt <- eeg_bandpass(eeg,
#'                            low_freq  = 0.1,
#'                            high_freq = 40,
#'                            verbose   = FALSE)
#' }
#'
#' @seealso \code{\link{downsample}} for reducing sampling rate before
#'   filtering, \code{\link{eeg_rereference}} for changing the reference
#'   scheme, \code{\link{new_eeg}} for the EEG object structure.
#'
#' @importFrom signal butter filtfilt filter fir1
#' @export
eeg_bandpass <- function(eeg_obj,
                         low_freq        = 0.5,
                         high_freq       = 30,
                         filter_order    = 4,
                         channels        = NULL,
                         zero_phase      = TRUE,
                         method          = c("butter", "fir"),
                         padding         = NULL,
                         notch_freq      = NULL,
                         notch_bandwidth = 2,
                         verbose         = TRUE) {
  
  # ========== MATCH ARGUMENTS ==========
  
  method <- match.arg(method)
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("eeg_obj must be an object of class 'eeg'.", call. = FALSE)
  }
  
  if (is.null(low_freq) && is.null(high_freq)) {
    stop("At least one of 'low_freq' or 'high_freq' must be specified.",
         call. = FALSE)
  }
  
  nyquist <- eeg_obj$sampling_rate / 2
  
  if (!is.null(low_freq)) {
    if (!is.numeric(low_freq) || length(low_freq) != 1 || low_freq <= 0) {
      stop("'low_freq' must be a single positive numeric value in Hz.",
           call. = FALSE)
    }
    if (low_freq >= nyquist) {
      stop("'low_freq' (", low_freq, " Hz) must be less than the Nyquist ",
           "frequency (", nyquist, " Hz).", call. = FALSE)
    }
  }
  
  if (!is.null(high_freq)) {
    if (!is.numeric(high_freq) || length(high_freq) != 1 || high_freq <= 0) {
      stop("'high_freq' must be a single positive numeric value in Hz.",
           call. = FALSE)
    }
    if (high_freq >= nyquist) {
      stop("'high_freq' (", high_freq, " Hz) must be less than the Nyquist ",
           "frequency (", nyquist, " Hz). Consider downsampling first.",
           call. = FALSE)
    }
  }
  
  if (!is.null(low_freq) && !is.null(high_freq)) {
    if (low_freq >= high_freq) {
      stop("'low_freq' (", low_freq, " Hz) must be less than 'high_freq' (",
           high_freq, " Hz).", call. = FALSE)
    }
  }
  
  if (!is.numeric(filter_order) || length(filter_order) != 1 ||
      filter_order < 1 || filter_order %% 1 != 0) {
    stop("'filter_order' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(zero_phase) || length(zero_phase) != 1) {
    stop("'zero_phase' must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
  }
  
  if (!is.null(notch_freq)) {
    if (!is.numeric(notch_freq) || length(notch_freq) != 1 || notch_freq <= 0) {
      stop("'notch_freq' must be a single positive numeric value in Hz.",
           call. = FALSE)
    }
    if (notch_freq >= nyquist) {
      stop("'notch_freq' (", notch_freq, " Hz) must be less than the Nyquist ",
           "frequency (", nyquist, " Hz).", call. = FALSE)
    }
    if (!is.numeric(notch_bandwidth) || length(notch_bandwidth) != 1 ||
        notch_bandwidth <= 0) {
      stop("'notch_bandwidth' must be a single positive numeric value.",
           call. = FALSE)
    }
    notch_low  <- notch_freq - notch_bandwidth / 2
    notch_high <- notch_freq + notch_bandwidth / 2
    if (notch_low <= 0 || notch_high >= nyquist) {
      stop("Notch filter edges (", notch_low, " – ", notch_high, " Hz) fall ",
           "outside the valid frequency range (0 – ", nyquist, " Hz).",
           call. = FALSE)
    }
  }
  
  if (any(is.na(eeg_obj$data))) {
    stop("EEG data contains NA values. Please handle missing data before ",
         "filtering.", call. = FALSE)
  }
  
  if (any(!is.finite(eeg_obj$data))) {
    stop("EEG data contains non-finite values (Inf or -Inf). Please clean ",
         "data before filtering.", call. = FALSE)
  }
  
  # ========== RESOLVE CHANNEL INDICES ==========
  
  n_channels  <- nrow(eeg_obj$data)
  chan_names  <- eeg_obj$channels
  
  if (is.null(channels)) {
    # Filter all channels
    ch_idx <- seq_len(n_channels)
  } else if (is.character(channels)) {
    ch_idx <- match(channels, chan_names)
    if (any(is.na(ch_idx))) {
      bad <- channels[is.na(ch_idx)]
      stop("Channel(s) not found in eeg_obj: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
  } else if (is.numeric(channels)) {
    ch_idx <- as.integer(channels)
    if (any(ch_idx < 1) || any(ch_idx > n_channels)) {
      stop("Channel indices in 'channels' are out of range [1, ", n_channels, "].",
           call. = FALSE)
    }
  } else {
    stop("'channels' must be NULL, a character vector of names, or an ",
         "integer vector of indices.", call. = FALSE)
  }
  
  # ========== AUTO-CALCULATE FIR ORDER (MNE-STYLE) ==========
  # For FIR, override filter_order with MNE's automatic rule:
  #   order = round(3.3 * sampling_rate / min_cutoff_freq)
  # This ensures the filter window is long enough to resolve the lowest
  # cutoff frequency, matching MNE's firwin behaviour exactly.
  # For Butterworth, filter_order is used as-is (small values are fine).
  
  if (method == "fir") {
    sr <- eeg_obj$sampling_rate
    
    # ---- MNE transition-bandwidth rules (per-edge) ----
    #   highpass: l_trans = min(max(l_freq * 0.25, 2.0), l_freq)
    #   lowpass:  h_trans = min(max(h_freq * 0.25, 2.0), nyquist - h_freq)
    # Filter length: N = round(3.3 * sr / min_transition_bw)
    #
    # Crucially, MNE uses the SAME (longest) kernel for the entire bandpass.
    # The narrowest transition bandwidth drives the filter length for BOTH edges.
    
    hp_trans_bw <- NULL
    lp_trans_bw <- NULL
    
    if (!is.null(low_freq)) {
      hp_trans_bw <- min(max(low_freq * 0.25, 2.0), low_freq)
    }
    if (!is.null(high_freq)) {
      lp_trans_bw <- min(max(high_freq * 0.25, 2.0), nyquist - high_freq)
    }
    
    # Use the narrowest transition BW → longest kernel (matches MNE)
    min_trans_bw <- min(c(hp_trans_bw, lp_trans_bw), na.rm = TRUE)
    fir_order    <- round(3.3 * sr / min_trans_bw)
    # Ensure odd length (required for symmetric linear-phase FIR)
    if (fir_order %% 2 == 0) fir_order <- fir_order + 1L
    
    if (verbose) {
      if (!is.null(hp_trans_bw)) {
        cat("  HP transition BW:     ", round(hp_trans_bw, 4), " Hz\n", sep = "")
      }
      if (!is.null(lp_trans_bw)) {
        cat("  LP transition BW:     ", round(lp_trans_bw, 4), " Hz\n", sep = "")
      }
      cat("  FIR order (shared):   ", fir_order,
          "  [3.3 x sr / min_trans_bw = ", round(min_trans_bw, 4), " Hz]\n",
          sep = "")
    }
  }
  
  # ========== RESOLVE PADDING ==========
  
  n_timepoints <- ncol(eeg_obj$data)
  
  if (is.null(padding)) {
    if (method == "fir") {
      # For FIR: pad by half the filter length (standard overlap-add convention)
      # Cap at 10% of signal length to avoid excessive memory use
      fir_half   <- fir_order %/% 2
      max_pad    <- as.integer(floor(n_timepoints * 0.10))
      padding    <- min(fir_half, max_pad)
    } else {
      padding <- 3L * as.integer(filter_order)
    }
  } else {
    if (!is.numeric(padding) || length(padding) != 1 || padding < 0 ||
        padding %% 1 != 0) {
      stop("'padding' must be a single non-negative integer.", call. = FALSE)
    }
    padding <- as.integer(padding)
  }
  
  if (padding > 0 && padding >= n_timepoints) {
    warning("'padding' (", padding, ") is >= the number of time points (",
            n_timepoints, "). Disabling padding.", call. = FALSE, immediate. = TRUE)
    padding <- 0L
  }
  
  # ========== PRINT HEADER ==========
  
  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("EEG Bandpass Filtering\n")
    cat(strrep("=", 70), "\n")
    cat("  Sampling rate:        ", eeg_obj$sampling_rate, " Hz\n")
    cat("  Nyquist frequency:    ", nyquist, " Hz\n")
    
    if (!is.null(low_freq) && !is.null(high_freq)) {
      cat("  Filter type:          Bandpass (", low_freq, " – ", high_freq,
          " Hz)\n", sep = "")
    } else if (!is.null(low_freq)) {
      cat("  Filter type:          Highpass (", low_freq, " Hz)\n", sep = "")
    } else {
      cat("  Filter type:          Lowpass (", high_freq, " Hz)\n", sep = "")
    }
    
    cat("  Filter design:        ",
        if (method == "butter") "Butterworth" else "FIR (firwin / Hamming)", "\n")
    cat("  Filter order:         ",
        if (method == "fir") fir_order else filter_order, "\n")
    cat("  Zero-phase:           ", if (zero_phase) "Yes (filtfilt)" else "No (causal)", "\n")
    cat("  Padding:              ", padding, " samples\n")
    
    if (!is.null(notch_freq)) {
      cat("  Notch filter:         ", notch_freq, " Hz (±",
          notch_bandwidth / 2, " Hz)\n", sep = "")
    }
    
    cat("  Channels to filter:   ", length(ch_idx), " / ", n_channels, "\n")
    cat(strrep("-", 70), "\n")
    cat("  Processing...\n")
  }
  
  # ========== DESIGN FILTERS ==========
  
  # Normalise frequencies to [0, 1] where 1 = Nyquist
  W_low  <- if (!is.null(low_freq))  low_freq  / nyquist else NULL
  W_high <- if (!is.null(high_freq)) high_freq / nyquist else NULL
  
  if (method == "butter") {
    
    # Highpass filter object
    if (!is.null(W_low)) {
      hp_filt <- signal::butter(filter_order, W_low, type = "high")
    }
    
    # Lowpass filter object
    if (!is.null(W_high)) {
      lp_filt <- signal::butter(filter_order, W_high, type = "low")
    }
    
    # Notch (band-stop) filter object
    if (!is.null(notch_freq)) {
      W_notch_low  <- (notch_freq - notch_bandwidth / 2) / nyquist
      W_notch_high <- (notch_freq + notch_bandwidth / 2) / nyquist
      notch_filt   <- signal::butter(filter_order,
                                     c(W_notch_low, W_notch_high),
                                     type = "stop")
    }
    
  } else {
    # ---- FIR design: firwin-style Hamming window (matches MNE) ----
    # firwin builds a linear-phase FIR kernel by multiplying an ideal
    # sinc impulse response with a Hamming window, then normalising.
    # This is identical to scipy.signal.firwin (MNE's default).
    
    .firwin_hamming <- function(n_taps, cutoff_norm, pass_type = "low") {
      # n_taps      : filter length (must be odd for linear phase)
      # cutoff_norm : cutoff frequency normalised to [0, 1] (1 = Nyquist)
      # pass_type   : "low" or "high"
      m      <- (n_taps - 1) / 2
      k      <- seq(-m, m)
      # Ideal sinc lowpass kernel
      h <- ifelse(k == 0,
                  2 * cutoff_norm,
                  sin(2 * pi * cutoff_norm * k) / (pi * k))
      # Hamming window
      w <- 0.54 - 0.46 * cos(2 * pi * seq(0, n_taps - 1) / (n_taps - 1))
      h <- h * w
      # Normalise to unit gain at DC (lowpass) or Nyquist (highpass)
      h <- h / sum(h)
      # Spectral inversion for highpass
      if (pass_type == "high") {
        h        <- -h
        h[m + 1] <- h[m + 1] + 1
      }
      h
    }
    
    if (!is.null(W_low)) {
      hp_filt <- .firwin_hamming(fir_order, W_low,  pass_type = "high")
    }
    if (!is.null(W_high)) {
      lp_filt <- .firwin_hamming(fir_order, W_high, pass_type = "low")
    }
    if (!is.null(notch_freq)) {
      # Band-stop = lowpass + highpass combined
      W_nl   <- (notch_freq - notch_bandwidth / 2) / nyquist
      W_nh   <- (notch_freq + notch_bandwidth / 2) / nyquist
      h_lp   <- .firwin_hamming(fir_order, W_nl,  pass_type = "low")
      h_hp   <- .firwin_hamming(fir_order, W_nh,  pass_type = "high")
      notch_filt <- h_lp + h_hp
    }
  }
  
  # ========== FFT-BASED CONVOLUTION HELPER ==========
  # Single-pass FFT convolution matching MNE's phase='zero' approach.
  
  .fft_convolve_fir <- function(x, h, zero_ph) {
    # Single-pass FFT convolution with delay compensation.
    #
    # zero_ph = TRUE  (MNE phase='zero'):
    #   Apply symmetric FIR ONCE, shift left by group delay (N-1)/2.
    #   Zero-phase because kernel is linear-phase.  Does NOT square
    #   the magnitude response (unlike filtfilt / double-pass).
    #
    # zero_ph = FALSE (causal):
    #   Apply once, no delay compensation.
    
    n_sig <- length(x)
    n_h   <- length(h)
    n_fft <- nextn(n_sig + n_h - 1, factors = 2)
    
    X      <- fft(c(x, rep(0, n_fft - n_sig)))
    H      <- fft(c(h, rep(0, n_fft - n_h)))
    y_full <- Re(fft(X * H, inverse = TRUE)) / n_fft
    
    if (zero_ph) {
      # Group delay of symmetric FIR = (n_h - 1) / 2 samples
      delay <- (n_h - 1) %/% 2
      y     <- y_full[seq(delay + 1, delay + n_sig)]
    } else {
      y <- y_full[seq_len(n_sig)]
    }
    
    y
  }
  
  # ========== UNIFIED APPLY HELPER ==========
  # Routes to FFT convolution (FIR) or filtfilt/filter (Butterworth)
  
  .apply_filter <- function(signal_vec, filt_obj) {
    if (method == "fir") {
      .fft_convolve_fir(signal_vec, filt_obj, zero_phase)
    } else if (zero_phase) {
      signal::filtfilt(filt_obj, signal_vec)
    } else {
      as.numeric(signal::filter(filt_obj, signal_vec))
    }
  }
  
  # ========== APPLY FILTERING CHANNEL BY CHANNEL ==========
  
  filtered_data <- eeg_obj$data  # copy; unselected channels untouched
  
  # Initialise progress bar (only when verbose)
  if (verbose) {
    cat("\n")
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }
  
  for (idx in seq_along(ch_idx)) {
    
    i <- ch_idx[idx]
    
    # Per-channel message above progress bar
    if (verbose) {
      cat(sprintf("\r  Channel %d / %d  [ %s ]%s",
                  idx, length(ch_idx),
                  chan_names[i],
                  strrep(" ", 20)))
      cat("\n")
      setTxtProgressBar(pb, idx)
    }
    
    x <- filtered_data[i, ]
    
    # ---- Reflect padding (matches numpy pad mode='reflect' / MNE) ----
    if (padding > 0) {
      left_pad  <- rev(x[seq(2, min(padding + 1, n_timepoints))])
      right_pad <- rev(x[seq(max(n_timepoints - padding, 1),
                             n_timepoints - 1)])
      x <- c(left_pad, x, right_pad)
    }
    
    # ---- Highpass ----
    if (!is.null(W_low)) {
      x <- .apply_filter(x, hp_filt)
    }
    
    # ---- Lowpass ----
    if (!is.null(W_high)) {
      x <- .apply_filter(x, lp_filt)
    }
    
    # ---- Notch ----
    if (!is.null(notch_freq)) {
      x <- .apply_filter(x, notch_filt)
    }
    
    # ---- Remove padding ----
    if (padding > 0) {
      x <- x[seq(padding + 1, padding + n_timepoints)]
    }
    
    filtered_data[i, ] <- x
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========
  
  filter_desc <- if (!is.null(low_freq) && !is.null(high_freq)) {
    paste0("Bandpass filtered: ", low_freq, " – ", high_freq, " Hz")
  } else if (!is.null(low_freq)) {
    paste0("Highpass filtered: ", low_freq, " Hz")
  } else {
    paste0("Lowpass filtered: ", high_freq, " Hz")
  }
  
  history_entry <- paste0(
    filter_desc,
    " | method: ", method,
    ", order: ", filter_order,
    ", zero-phase: ", if (zero_phase) "yes" else "no",
    ", channels: ", length(ch_idx), "/", n_channels,
    if (!is.null(notch_freq)) paste0(", notch: ", notch_freq, " Hz") else ""
  )
  
  new_history <- c(eeg_obj$preprocessing_history, list(history_entry))
  
  # ========== CONSTRUCT OUTPUT EEG OBJECT ==========
  
  filtered_eeg <- new_eeg(
    data                 = filtered_data,
    channels             = eeg_obj$channels,
    sampling_rate        = eeg_obj$sampling_rate,
    times                = eeg_obj$times,
    events               = eeg_obj$events,
    metadata             = eeg_obj$metadata,
    reference            = eeg_obj$reference,
    preprocessing_history = new_history
  )
  
  # ========== PRINT SUMMARY ==========
  
  if (verbose) {
    cat("  Done.\n")
    cat(strrep("=", 70), "\n")
    cat("Filtering Complete\n")
    cat(strrep("=", 70), "\n")
    cat("  Filter applied:       ", filter_desc, "\n")
    cat("  Channels filtered:    ", length(ch_idx), " / ", n_channels, "\n")
    if (!is.null(notch_freq)) {
      cat("  Notch applied:        ", notch_freq, " Hz\n")
    }
    cat("  History entry added:  Yes\n")
    cat(strrep("=", 70), "\n\n")
  }
  
  return(filtered_eeg)
}

# ============================================================================
#
#' Apply a Notch Filter to EEG Data
#'
#' Removes a narrow frequency band from EEG data to eliminate electrical line
#' noise interference. The most common use is removing 50 Hz noise (Europe)
#' or 60 Hz noise (North America) introduced by the power grid into the
#' recording equipment. Uses a zero-phase Butterworth band-stop filter applied
#' with forward-backward (filtfilt) convolution to preserve event timing.
#'
#' @param eeg_obj An object of class 'eeg' containing EEG data.
#'
#' @param freq Numeric. The centre frequency of the notch in Hz. This is the
#'             frequency you want to remove.
#'             Use \code{50} for Europe, Asia, Africa, and Australia.
#'             Use \code{60} for North America.
#'             Default: \code{50}.
#'
#' @param bandwidth Numeric. The width of the notch in Hz, centred on
#'                  \code{freq}. A value of \code{2} means the filter removes
#'                  frequencies from \code{freq - 1} to \code{freq + 1} Hz.
#'                  Keeping this narrow preserves as much real brain signal as
#'                  possible. Typical range: 1–4 Hz.
#'                  Default: \code{2}.
#'
#' @param filter_order Integer. Order of the Butterworth band-stop filter.
#'                     Higher values create a steeper notch but may introduce
#'                     more edge artefacts. Typical range: 2–6.
#'                     Default: \code{4}.
#'
#' @param harmonics Logical. If \code{TRUE}, also applies the notch filter to
#'                  harmonics of \code{freq} — i.e. 100 Hz, 150 Hz, etc. for
#'                  a 50 Hz notch. Harmonics are integer multiples of the line
#'                  frequency that can appear in recordings with strong
#'                  interference. Only harmonics below the Nyquist frequency
#'                  are filtered. Default: \code{FALSE}.
#'
#' @param channels Character or integer vector. Channel names or indices to
#'                 filter. If \code{NULL} (default), all channels are filtered.
#'                 Useful for skipping external channels (e.g. EOG, EMG) that
#'                 do not need notch filtering.
#'                 Example: \code{c("Cz", "Pz")} or \code{1:64}.
#'
#' @param zero_phase Logical. If \code{TRUE} (default), uses forward-backward
#'                   filtering (\code{filtfilt}) which produces zero phase
#'                   distortion, preserving the timing of neural events.
#'                   If \code{FALSE}, uses one-pass causal filtering, which is
#'                   faster but shifts the signal slightly in time.
#'                   Default: \code{TRUE}.
#'
#' @param padding Integer. Number of samples to mirror-pad at each end of
#'                each channel before filtering to reduce edge artefacts.
#'                Set to \code{0} to disable. Default: \code{3 * filter_order}
#'                (automatically scaled to filter length).
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a summary of the
#'                filter parameters and processing progress. Set to \code{FALSE}
#'                for silent batch processing.
#'
#' @return An object of class 'eeg' with the same structure as the input,
#'   with the following fields updated:
#'   \describe{
#'     \item{data}{Notch-filtered signal matrix (channels x time points).
#'                 Unfiltered channels retain their original values.}
#'     \item{preprocessing_history}{Appended with a string describing the
#'                                   notch filter parameters applied.}
#'     \item{channels}{Unchanged.}
#'     \item{sampling_rate}{Unchanged.}
#'     \item{times}{Unchanged.}
#'     \item{events}{Unchanged.}
#'     \item{metadata}{Unchanged.}
#'     \item{reference}{Unchanged.}
#'   }
#'
#' @details
#' The notch filter is a band-stop filter — it passes all frequencies except
#' those within the notch window. It is applied after bandpass filtering in a
#' typical EEG pipeline.
#'
#' When \code{harmonics = TRUE}, the function automatically calculates all
#' integer multiples of \code{freq} that fall below the Nyquist frequency and
#' applies a separate notch at each one. For example, with \code{freq = 50}
#' and a sampling rate of 512 Hz (Nyquist = 256 Hz), notches are applied at
#' 50, 100, 150, and 200 Hz.
#'
#' \strong{When to use this function vs the notch argument in eeg_bandpass():}
#' The \code{notch_freq} argument inside \code{eeg_bandpass()} is a convenience
#' shortcut for applying a single notch in the same call as the bandpass.
#' Use \code{eeg_notch()} as a standalone function when you want more control,
#' need to filter harmonics, or want to apply notch filtering independently
#' from bandpass filtering.
#'
#' @examples
#' \dontrun{
#'   # Remove 50 Hz line noise (Europe)
#'   eeg_clean <- eeg_notch(eeg, freq = 50)
#'
#'   # Remove 60 Hz line noise (North America)
#'   eeg_clean <- eeg_notch(eeg, freq = 60)
#'
#'   # Remove 50 Hz and all harmonics (100, 150, 200 Hz...)
#'   eeg_clean <- eeg_notch(eeg, freq = 50, harmonics = TRUE)
#'
#'   # Wider notch for strong interference
#'   eeg_clean <- eeg_notch(eeg, freq = 50, bandwidth = 4)
#'
#'   # Filter only EEG channels, leave EOG untouched
#'   eeg_clean <- eeg_notch(eeg, freq = 50, channels = 1:64)
#'
#'   # Silent batch processing
#'   eeg_clean <- eeg_notch(eeg, freq = 50, verbose = FALSE)
#' }
#'
#' @seealso \code{\link{eeg_bandpass}} for bandpass filtering with an optional
#'   built-in notch, \code{\link{new_eeg}} for the EEG object structure.
#'
#' @importFrom signal butter filtfilt filter
#' @export
eeg_notch <- function(eeg_obj,
                      freq         = 50,
                      bandwidth    = 2,
                      filter_order = 4,
                      harmonics    = FALSE,
                      channels     = NULL,
                      zero_phase   = TRUE,
                      padding      = NULL,
                      verbose      = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  
  nyquist <- eeg_obj$sampling_rate / 2
  
  if (!is.numeric(freq) || length(freq) != 1 || freq <= 0) {
    stop("'freq' must be a single positive numeric value in Hz.", call. = FALSE)
  }
  if (freq >= nyquist) {
    stop("'freq' (", freq, " Hz) must be less than the Nyquist frequency (",
         nyquist, " Hz).", call. = FALSE)
  }
  
  if (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0) {
    stop("'bandwidth' must be a single positive numeric value.", call. = FALSE)
  }
  
  notch_low  <- freq - bandwidth / 2
  notch_high <- freq + bandwidth / 2
  
  if (notch_low <= 0) {
    stop("Notch lower edge (", notch_low, " Hz) must be above 0 Hz. ",
         "Reduce 'bandwidth' or increase 'freq'.", call. = FALSE)
  }
  if (notch_high >= nyquist) {
    stop("Notch upper edge (", notch_high, " Hz) must be below the Nyquist ",
         "frequency (", nyquist, " Hz). Reduce 'bandwidth'.", call. = FALSE)
  }
  
  if (!is.numeric(filter_order) || length(filter_order) != 1 ||
      filter_order < 1 || filter_order %% 1 != 0) {
    stop("'filter_order' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(harmonics) || length(harmonics) != 1) {
    stop("'harmonics' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(zero_phase) || length(zero_phase) != 1) {
    stop("'zero_phase' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (any(is.na(eeg_obj$data))) {
    stop("EEG data contains NA values. Please handle missing data before ",
         "filtering.", call. = FALSE)
  }
  
  if (any(!is.finite(eeg_obj$data))) {
    stop("EEG data contains non-finite values (Inf or -Inf). Please clean ",
         "data before filtering.", call. = FALSE)
  }
  
  # ========== RESOLVE CHANNEL INDICES ==========
  
  n_channels   <- nrow(eeg_obj$data)
  n_timepoints <- ncol(eeg_obj$data)
  chan_names   <- eeg_obj$channels
  
  if (is.null(channels)) {
    ch_idx <- seq_len(n_channels)
  } else if (is.character(channels)) {
    ch_idx <- match(channels, chan_names)
    if (any(is.na(ch_idx))) {
      bad <- channels[is.na(ch_idx)]
      stop("Channel(s) not found in eeg_obj: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
  } else if (is.numeric(channels)) {
    ch_idx <- as.integer(channels)
    if (any(ch_idx < 1) || any(ch_idx > n_channels)) {
      stop("Channel indices in 'channels' are out of range [1, ",
           n_channels, "].", call. = FALSE)
    }
  } else {
    stop("'channels' must be NULL, a character vector of names, or an ",
         "integer vector of indices.", call. = FALSE)
  }
  
  # ========== BUILD LIST OF NOTCH FREQUENCIES ==========
  
  # Start with the fundamental frequency
  notch_freqs <- freq
  
  # Add harmonics below Nyquist if requested
  if (harmonics) {
    harmonic_mult <- 2
    while ((freq * harmonic_mult) < nyquist) {
      h_freq      <- freq * harmonic_mult
      h_low       <- h_freq - bandwidth / 2
      h_high      <- h_freq + bandwidth / 2
      if (h_low > 0 && h_high < nyquist) {
        notch_freqs <- c(notch_freqs, h_freq)
      }
      harmonic_mult <- harmonic_mult + 1
    }
  }
  
  # ========== RESOLVE PADDING ==========
  
  if (is.null(padding)) {
    padding <- 3L * as.integer(filter_order)
  } else {
    if (!is.numeric(padding) || length(padding) != 1 || padding < 0 ||
        padding %% 1 != 0) {
      stop("'padding' must be a single non-negative integer.", call. = FALSE)
    }
    padding <- as.integer(padding)
  }
  
  if (padding > 0 && padding >= n_timepoints) {
    warning("'padding' (", padding, ") is >= the number of time points (",
            n_timepoints, "). Disabling padding.",
            call. = FALSE, immediate. = TRUE)
    padding <- 0L
  }
  
  # ========== PRINT HEADER ==========
  
  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("EEG Notch Filtering\n")
    cat(strrep("=", 70), "\n")
    cat("  Sampling rate:        ", eeg_obj$sampling_rate, " Hz\n")
    cat("  Nyquist frequency:    ", nyquist, " Hz\n")
    cat("  Notch frequency:      ", freq, " Hz\n")
    cat("  Bandwidth:            ", bandwidth, " Hz  (",
        notch_low, " – ", notch_high, " Hz)\n", sep = "")
    cat("  Filter order:         ", filter_order, "\n")
    cat("  Zero-phase:           ",
        if (zero_phase) "Yes (filtfilt)" else "No (causal)", "\n")
    cat("  Harmonics:            ",
        if (harmonics) paste(notch_freqs, collapse = ", ") else "No", " Hz\n",
        sep = "")
    cat("  Padding:              ", padding, " samples\n")
    cat("  Channels to filter:   ", length(ch_idx), " / ", n_channels, "\n")
    cat(strrep("-", 70), "\n")
    cat("  Processing...\n")
  }
  
  # ========== DESIGN FILTER OBJECTS FOR EACH NOTCH FREQUENCY ==========
  # Notch always uses Butterworth band-stop (best choice for narrow notches)
  
  .make_notch_filter <- function(nf) {
    W_low  <- (nf - bandwidth / 2) / nyquist
    W_high <- (nf + bandwidth / 2) / nyquist
    signal::butter(filter_order, c(W_low, W_high), type = "stop")
  }
  
  notch_filters <- lapply(notch_freqs, .make_notch_filter)
  
  # ========== FFT-BASED CONVOLUTION HELPER ==========
  # Used when the filter object is a raw numeric FIR kernel.
  # For Butterworth (IIR), routes to filtfilt/filter as before.
  
  .fft_convolve_fir <- function(x, h, zero_ph) {
    n_sig   <- length(x)
    n_h     <- length(h)
    n_fft   <- nextn(n_sig + n_h - 1, factors = 2)
    X       <- fft(c(x, rep(0, n_fft - n_sig)))
    H       <- fft(c(h, rep(0, n_fft - n_h)))
    y_full  <- Re(fft(X * H, inverse = TRUE)) / n_fft
    delay   <- (n_h - 1) %/% 2
    y       <- y_full[seq(delay + 1, delay + n_sig)]
    if (zero_ph) {
      y_rev   <- rev(y)
      Xr      <- fft(c(y_rev, rep(0, n_fft - n_sig)))
      yr_full <- Re(fft(Xr * H, inverse = TRUE)) / n_fft
      y       <- rev(yr_full[seq(delay + 1, delay + n_sig)])
    }
    y
  }
  
  .apply_notch_filter <- function(signal_vec, filt_obj) {
    if (zero_phase) {
      signal::filtfilt(filt_obj, signal_vec)
    } else {
      as.numeric(signal::filter(filt_obj, signal_vec))
    }
  }
  
  # ========== APPLY FILTERING ==========
  
  filtered_data <- eeg_obj$data
  
  # Initialise progress bar (only when verbose)
  if (verbose) {
    cat("\n")
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }
  
  for (idx in seq_along(ch_idx)) {
    
    i <- ch_idx[idx]
    
    # Per-channel message above progress bar
    if (verbose) {
      cat(sprintf("\r  Channel %d / %d  [ %s ]%s",
                  idx, length(ch_idx),
                  chan_names[i],
                  strrep(" ", 20)))
      cat("\n")
      setTxtProgressBar(pb, idx)
    }
    
    x <- filtered_data[i, ]
    
    # ---- Reflect padding (matches numpy pad mode='reflect' / MNE) ----
    if (padding > 0) {
      left_pad  <- rev(x[seq(2, min(padding + 1, n_timepoints))])
      right_pad <- rev(x[seq(max(n_timepoints - padding, 1),
                             n_timepoints - 1)])
      x <- c(left_pad, x, right_pad)
    }
    
    # Apply each notch sequentially
    for (filt_obj in notch_filters) {
      x <- .apply_notch_filter(x, filt_obj)
    }
    
    # Remove padding
    if (padding > 0) {
      x <- x[seq(padding + 1, padding + n_timepoints)]
    }
    
    filtered_data[i, ] <- x
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========
  
  history_entry <- paste0(
    "Notch filtered: ", freq, " Hz",
    " (", notch_low, " – ", notch_high, " Hz)",
    " | bandwidth: ", bandwidth, " Hz",
    ", order: ", filter_order,
    ", zero-phase: ", if (zero_phase) "yes" else "no",
    if (harmonics) paste0(", harmonics: ",
                          paste(notch_freqs[-1], collapse = ", "), " Hz") else "",
    ", channels: ", length(ch_idx), "/", n_channels
  )
  
  new_history <- c(eeg_obj$preprocessing_history, list(history_entry))
  
  # ========== CONSTRUCT OUTPUT EEG OBJECT ==========
  
  notched_eeg <- new_eeg(
    data                  = filtered_data,
    channels              = eeg_obj$channels,
    sampling_rate         = eeg_obj$sampling_rate,
    times                 = eeg_obj$times,
    events                = eeg_obj$events,
    metadata              = eeg_obj$metadata,
    reference             = eeg_obj$reference,
    preprocessing_history = new_history
  )
  
  # ========== PRINT SUMMARY ==========
  
  if (verbose) {
    cat("  Done.\n")
    cat(strrep("=", 70), "\n")
    cat("Notch Filtering Complete\n")
    cat(strrep("=", 70), "\n")
    cat("  Notch applied:        ", freq, " Hz\n")
    if (harmonics && length(notch_freqs) > 1) {
      cat("  Harmonics removed:    ",
          paste(notch_freqs[-1], collapse = ", "), " Hz\n")
    }
    cat("  Channels filtered:    ", length(ch_idx), " / ", n_channels, "\n")
    cat("  History entry added:  Yes\n")
    cat(strrep("=", 70), "\n\n")
  }
  
  return(notched_eeg)
}

# ============================================================================
#
#' Apply Bandpass and Optional Notch Filtering to EEG Data
#'
#' A unified wrapper that combines \code{eeg_bandpass()} and \code{eeg_notch()}
#' into a single call. Applies bandpass filtering first (to set the frequency
#' window of interest), then optionally applies a notch filter (to remove line
#' noise), in the correct order. All arguments are passed directly to the
#' underlying functions.
#'
#' @param eeg_obj An object of class 'eeg' containing EEG data.
#'
#' @param low_freq Numeric. Highpass cutoff in Hz. Frequencies below this are
#'                 removed. Set to \code{NULL} for lowpass only.
#'                 Default: \code{0.5}.
#'
#' @param high_freq Numeric. Lowpass cutoff in Hz. Frequencies above this are
#'                  removed. Set to \code{NULL} for highpass only.
#'                  Default: \code{40}.
#'
#' @param notch_freq Numeric. Centre frequency of the notch filter in Hz.
#'                   Use \code{50} (Europe) or \code{60} (North America).
#'                   Set to \code{NULL} (default) to skip notch filtering.
#'
#' @param notch_bandwidth Numeric. Width of the notch in Hz, centred on
#'                        \code{notch_freq}. Default: \code{2}.
#'                        Only used when \code{notch_freq} is not \code{NULL}.
#'
#' @param notch_harmonics Logical. If \code{TRUE}, also removes harmonics of
#'                        \code{notch_freq} (e.g. 100, 150, 200 Hz for a
#'                        50 Hz notch). Default: \code{FALSE}.
#'                        Only used when \code{notch_freq} is not \code{NULL}.
#'
#' @param filter_order Integer. Order of the Butterworth filter applied for
#'                     both bandpass and notch stages. Default: \code{4}.
#'
#' @param method Character string specifying the filter design method.
#'               \code{"butter"} (default) uses a Butterworth filter —
#'               smooth, flat passband, recommended for all standard EEG work.
#'               \code{"fir"} uses a Finite Impulse Response filter.
#'               Only applies to the bandpass stage — the notch always uses
#'               Butterworth.
#'
#' @param channels Character or integer vector. Channel names or indices to
#'                 filter. \code{NULL} (default) filters all channels.
#'
#' @param zero_phase Logical. If \code{TRUE} (default), uses zero-phase
#'                   forward-backward filtering (\code{filtfilt}).
#'
#' @param padding Integer. Mirror-padding samples added at each end before
#'                filtering to reduce edge artefacts. Default: \code{NULL}
#'                (auto: \code{3 * filter_order}).
#'
#' @param verbose Logical. If \code{TRUE} (default), prints progress for both
#'                filtering stages. Set to \code{FALSE} for silent processing.
#'
#' @return An object of class 'eeg' with filtered data and updated
#'   \code{preprocessing_history}. If \code{notch_freq} is not \code{NULL},
#'   the history will contain two entries — one for the bandpass and one for
#'   the notch.
#'
#' @examples
#' \dontrun{
#'   # Standard ERP pipeline in one line
#'   eeg_clean <- eeg_filter(data4,
#'                           low_freq   = 1,
#'                           high_freq  = 40,
#'                           notch_freq = 50)
#'
#'   # Bandpass only — no notch
#'   eeg_clean <- eeg_filter(data4, low_freq = 0.1, high_freq = 40)
#'
#'   # Bandpass + notch + harmonics, EEG channels only
#'   eeg_clean <- eeg_filter(data4,
#'                           low_freq         = 1,
#'                           high_freq        = 100,
#'                           notch_freq       = 50,
#'                           notch_harmonics  = TRUE,
#'                           channels         = 1:64)
#'
#'   # Silent batch mode
#'   eeg_clean <- eeg_filter(data4,
#'                           low_freq   = 1,
#'                           high_freq  = 40,
#'                           notch_freq = 50,
#'                           verbose    = FALSE)
#' }
#'
#' @seealso \code{\link{eeg_bandpass}} for bandpass-only filtering,
#'   \code{\link{eeg_notch}} for notch-only filtering,
#'   \code{\link{new_eeg}} for the EEG object structure.
#'
#' @export
eeg_filter <- function(eeg_obj,
                       low_freq         = 0.5,
                       high_freq        = 40,
                       notch_freq       = NULL,
                       notch_bandwidth  = 2,
                       notch_harmonics  = FALSE,
                       filter_order     = 4,
                       method           = c("butter", "fir"),
                       channels         = NULL,
                       zero_phase       = TRUE,
                       padding          = NULL,
                       verbose          = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  
  method <- match.arg(method)
  
  if (is.null(low_freq) && is.null(high_freq) && is.null(notch_freq)) {
    stop("At least one of 'low_freq', 'high_freq', or 'notch_freq' must be ",
         "specified.", call. = FALSE)
  }
  
  # ========== STAGE 1: BANDPASS ==========
  # Only run if at least one of low_freq / high_freq is specified
  
  if (!is.null(low_freq) || !is.null(high_freq)) {
    eeg_obj <- eeg_bandpass(
      eeg_obj      = eeg_obj,
      low_freq     = low_freq,
      high_freq    = high_freq,
      filter_order = filter_order,
      method       = method,
      channels     = channels,
      zero_phase   = zero_phase,
      padding      = padding,
      verbose      = verbose
    )
  }
  
  # ========== STAGE 2: NOTCH ==========
  # Only run if notch_freq is specified
  
  if (!is.null(notch_freq)) {
    eeg_obj <- eeg_notch(
      eeg_obj      = eeg_obj,
      freq         = notch_freq,
      bandwidth    = notch_bandwidth,
      filter_order = filter_order,
      harmonics    = notch_harmonics,
      channels     = channels,
      zero_phase   = zero_phase,
      padding      = padding,
      verbose      = verbose
    )
  }
  
  return(eeg_obj)
}

# End of filter.R