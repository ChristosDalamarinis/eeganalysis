#' ============================================================================
#'                       FIR Filtering Functions
#' ============================================================================
#'
#' This module implements FIR (Finite Impulse Response) filtering for EEG
#' data, designed to match MNE-Python's default filtering pipeline exactly.
#'
#' Both bandpass and notch filtering are implemented using:
#'   - firwin / Hamming window kernel design (matches scipy.signal.firwin)
#'   - reflect_limited padding  (matches MNE's default pad='reflect_limited')
#'   - Overlap-add FFT convolution (matches MNE's _overlap_add_filter)
#'   - Single-pass zero-phase correction (matches MNE's phase='zero')
#'
#' Depends on: eeg_class.R (for eeg object structure)
#'
#' Author: Christos Dalamarinis
#' Date: March 2026
#' Status: In development!!!!!!!!!! - Not tested
#' An iteration of this script is being used to finalise bugs
#' ============================================================================
#
#
# ============================================================================
#                        PRIVATE HELPER FUNCTIONS
# ============================================================================
#
# ----------------------------------------------------------------------------
# .next_fast_len() - preparation helper
# ----------------------------------------------------------------------------
#' Find the next fast FFT length (5-smooth number)
#'
#' Returns the smallest integer >= target whose only prime factors are
#' 2, 3, and 5 (a "5-smooth" or "regular" number). These are the sizes
#' at which FFT is most efficient. Matches MNE's next_fast_len() exactly.
#'
#' @param target Positive integer. The minimum size required.
#' @return Integer. The smallest 5-smooth number >= target.
#' @keywords internal
.next_fast_len <- function(target) {
  
  # Lookup table for targets <= 10000 (same as MNE / SciPy)
  hams <- c(
    8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48,
    50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128,
    135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250,
    256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450,
    480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729,
    750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080, 1125,
    1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536,
    1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160,
    2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916,
    3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840,
    3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800, 4860, 5000,
    5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144, 6250, 6400,
    6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776, 8000, 8100,
    8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000
  )
  
  # Trivial case: target <= 6, return as-is (MNE line 239)
  if (target <= 6L) return(target)
  
  # Already a power of 2? (MNE line 243: not (target & (target-1)))
  if (bitwAnd(target, target - 1L) == 0L) return(target)
  
  # Small target: binary search in lookup table (MNE line 247-248)
  if (target <= hams[length(hams)]) {
    idx <- which(hams >= target)[1]
    return(hams[idx])
  }
  
  # Large target: algorithmic search (MNE lines 250-276)
  # Find smallest N = 2^a * 3^b * 5^c >= target
  
  # Helper: number of bits to represent n (Python's int.bit_length())
  # bit_length(0) = 0, bit_length(n) = floor(log2(n)) + 1 for n >= 1
  .bit_length <- function(n) {
    if (n == 0L) return(0L)
    floor(log2(n)) + 1L
  }
  
  best <- Inf
  p5 <- 1L
  
  while (p5 < target) {
    p35 <- p5
    
    while (p35 < target) {
      # Ceiling integer division: ceil(target / p35)
      quotient <- ((target - 1L) %/% p35) + 1L
      
      # Next power of 2 >= quotient
      p2 <- 2L ^ .bit_length(quotient - 1L)
      
      N <- p2 * p35
      
      if (N == target) return(N)
      if (N < best)   best <- N
      
      p35 <- p35 * 3L
      if (p35 == target) return(p35)
    }
    
    if (p35 < best) best <- p35
    
    p5 <- p5 * 5L
    if (p5 == target) return(p5)
  }
  
  if (p5 < best) best <- p5
  
  return(best)
}

# ----------------------------------------------------------------------------
# .auto_trans_bandwidth() - preparation helper
# ----------------------------------------------------------------------------
#
#' Compute automatic transition bandwidths (MNE-style)
#'
#' Computes the transition bandwidth for highpass and/or lowpass edges
#' using MNE's automatic rule (from _triage_filter_params, lines 2228, 2267).
#'
#' For the highpass edge:
#'   l_trans_bw = min(max(0.25 * l_freq, 2.0), l_freq)
#'
#' For the lowpass edge:
#'   h_trans_bw = min(max(0.25 * h_freq, 2.0), nyquist - h_freq)
#'
#' @param l_freq  Numeric or NULL. Highpass cutoff frequency in Hz.
#' @param h_freq  Numeric or NULL. Lowpass cutoff frequency in Hz.
#' @param sfreq   Numeric. Sampling rate in Hz.
#'
#' @return A named list with:
#'   \item{l_trans_bw}{Highpass transition bandwidth in Hz, or NULL}
#'   \item{h_trans_bw}{Lowpass transition bandwidth in Hz, or NULL}
#' @keywords internal
.auto_trans_bandwidth <- function(l_freq = NULL, h_freq = NULL, sfreq) {
  
  nyquist <- sfreq / 2.0
  
  # Highpass transition bandwidth (MNE line 2228)
  l_trans_bw <- if (!is.null(l_freq)) {
    min(max(0.25 * l_freq, 2.0), l_freq)
  } else {
    NULL
  }
  
  # Lowpass transition bandwidth (MNE line 2267)
  h_trans_bw <- if (!is.null(h_freq)) {
    min(max(0.25 * h_freq, 2.0), nyquist - h_freq)
  } else {
    NULL
  }
  
  list(l_trans_bw = l_trans_bw, h_trans_bw = h_trans_bw)
}

# ----------------------------------------------------------------------------
# .auto_filter_length() - preparation helper
# ----------------------------------------------------------------------------
#
#' Compute automatic FIR filter length (MNE-style)
#'
#' Computes the filter kernel length from transition bandwidths using MNE's
#' automatic rule (_triage_filter_params lines 2300-2318, _to_samples lines
#' 2099-2125).
#'
#' Full chain:
#'   1. min_trans_bw = min(l_trans_bw, h_trans_bw)   [narrowest edge drives length]
#'   2. duration_s   = 3.3 / min_trans_bw             [3.3 = hamming length factor]
#'   3. n_samples    = ceil(duration_s x sfreq)
#'   4. force odd    = if even, add 1
#'
#' @param l_trans_bw  Numeric or NULL. Highpass transition bandwidth in Hz.
#' @param h_trans_bw  Numeric or NULL. Lowpass transition bandwidth in Hz.
#' @param sfreq       Numeric. Sampling rate in Hz.
#'
#' @return Integer. Filter kernel length in samples (always odd).
#' @keywords internal
.auto_filter_length <- function(l_trans_bw = NULL, h_trans_bw = NULL, sfreq) {
  
  # Collect non-NULL bandwidths and take the minimum (MNE lines 2302-2306)
  bws <- c(l_trans_bw, h_trans_bw)
  if (length(bws) == 0) {
    stop("At least one of l_trans_bw or h_trans_bw must be non-NULL.",
         call. = FALSE)
  }
  min_trans_bw <- min(bws)
  
  # Length factor for Hamming window + firwin (MNE line 40 + 2308)
  # _length_factors["hamming"] = 3.3, mult_fact = 1.0 for firwin
  length_factor <- 3.3
  
  # Convert to duration in seconds, then to samples (MNE lines 2308, 2121)
  duration_s <- length_factor / min_trans_bw
  n_samples  <- ceiling(duration_s * sfreq)
  
  # Force odd length (MNE lines 2122-2123 in _to_samples)
  # (n - 1) %% 2 == 1 only when n is even → adds 1 to make it odd
  n_samples <- n_samples + (n_samples - 1L) %% 2L
  
  as.integer(n_samples)
}

# ----------------------------------------------------------------------------
# .firwin_kernel() - preparation helper
# ----------------------------------------------------------------------------

#' Normalized sinc function
#'
#' sinc(0) = 1, sinc(x) = sin(pi*x) / (pi*x) for x != 0.
#' Matches numpy.sinc() exactly.
#'
#' @param x Numeric vector.
#' @return Numeric vector of same length as x.
#' @keywords internal
.sinc_norm <- function(x) {
  ifelse(x == 0, 1.0, sin(pi * x) / (pi * x))
}


#' Single lowpass FIR kernel via firwin (scipy.signal.firwin equivalent)
#'
#' Builds one lowpass kernel of length N with normalized cutoff cutoff_norm,
#' using a Hamming window. Matches scipy.signal.firwin(N, cutoff_norm,
#' window='hamming', pass_zero=True, fs=2.0) exactly.
#'
#' @param N          Integer. Number of taps (must be odd).
#' @param cutoff_norm Numeric. Normalized cutoff in [0, 1] where 1 = Nyquist.
#' @return Numeric vector of length N.
#' @keywords internal
.firwin_lowpass <- function(N, cutoff_norm) {
  
  M     <- (N - 1) / 2
  n     <- seq(0L, N - 1L)
  alpha <- n - M
  
  # Ideal sinc impulse response (lowpass)
  h <- cutoff_norm * .sinc_norm(cutoff_norm * alpha)
  
  # Hamming window: 0.54 - 0.46 * cos(2*pi*n / (N-1))
  w <- 0.54 - 0.46 * cos(2 * pi * n / (N - 1L))
  
  # Apply window and normalize so DC gain = 1
  h <- h * w
  h <- h / sum(h)
  
  h
}


#' Build FIR kernel from freq/gain arrays (MNE's _firwin_design equivalent)
#'
#' Constructs a FIR kernel by combining lowpass firwin kernels, one per
#' gain transition in the freq/gain specification. Matches MNE's
#' _firwin_design() exactly.
#'
#' Algorithm (right-to-left sweep):
#'   - When gain changes 0→1: ADD a lowpass kernel at the transition midpoint
#'   - When gain changes 1→0: SUBTRACT a lowpass kernel at the transition midpoint
#'   - No change: skip
#'
#' @param N    Integer. Total kernel length (must be odd).
#' @param freq Numeric vector. Normalized frequencies in [0, 1] starting at 0
#'             and ending at 1 (Nyquist). E.g. c(0, 0.1, 0.2, 0.9, 1).
#' @param gain Numeric vector. Gain at each frequency (0 or 1 only),
#'             same length as freq.
#'
#' @return Numeric vector of length N — the FIR kernel coefficients.
#' @keywords internal
.firwin_kernel <- function(N, freq, gain) {
  
  stopifnot(freq[1] == 0)
  stopifnot(length(freq) == length(gain))
  stopifnot(N %% 2L == 1L)
  
  h      <- numeric(N)
  center <- N %/% 2L + 1L  # 1-indexed center (= Python's N//2)
  
  # If gain at Nyquist == 1, start with "all pass" delta (MNE line 431)
  if (gain[length(gain)] == 1) {
    h[center] <- 1.0
  }
  
  prev_freq <- freq[length(freq)]
  prev_gain <- gain[length(gain)]
  
  # Walk right to left, skipping the last element (MNE line 433)
  for (i in seq(length(freq) - 1L, 1L)) {
    
    this_freq <- freq[i]
    this_gain <- gain[i]
    
    if (this_gain != prev_gain) {
      
      # Half-width of transition band (MNE line 437)
      transition <- (prev_freq - this_freq) / 2.0
      
      # Per-transition kernel length — may be shorter than N (MNE line 438-439)
      this_N <- as.integer(round(3.3 / transition))
      this_N <- this_N + 1L - this_N %% 2L  # force odd
      
      # Cutoff at midpoint of transition band (MNE line 449)
      cutoff <- (prev_freq + this_freq) / 2.0
      
      # Build lowpass kernel of length this_N (MNE line 447-453)
      this_h <- .firwin_lowpass(this_N, cutoff)
      
      # Embed this_h centered in h (MNE line 455-459)
      offset <- (N - this_N) %/% 2L
      idx    <- seq(offset + 1L, offset + this_N)
      
      if (this_gain == 0L) {
        h[idx] <- h[idx] - this_h  # subtract → removes that band
      } else {
        h[idx] <- h[idx] + this_h  # add → passes that band
      }
    }
    
    prev_gain <- this_gain
    prev_freq <- this_freq
  }
  
  h
}

# ----------------------------------------------------------------------------
# .reflect_limited_pad() - preparation helper
# ----------------------------------------------------------------------------

#' Pad a signal using MNE's reflect_limited mode
#'
#' Replicates MNE's _smart_pad() from mne/cuda.py with pad='reflect_limited'.
#'
#' reflect_limited pads with an odd reflection anchored on the edge values,
#' i.e.:  left pad  = 2*x[1]   - x[reversed inner samples]
#'         right pad = 2*x[end] - x[reversed inner samples]
#' If n_pad exceeds len(x)-1, zeros fill the remaining space.
#'
#' This differs from standard numpy 'reflect' (which excludes the edge value)
#' and 'symmetric' (which includes the edge value but does not subtract).
#'
#' Python equivalent (mne/cuda.py):
#'   l_z_pad = zeros(max(n_pad[0] - len(x) + 1, 0))
#'   r_z_pad = zeros(max(n_pad[1] - len(x) + 1, 0))
#'   out = [l_z_pad,
#'          2*x[0] - x[n_pad[0]:0:-1],
#'          x,
#'          2*x[-1] - x[-2:-n_pad[1]-2:-1],
#'          r_z_pad]
#'
#' @param x       Numeric vector. The signal to pad.
#' @param n_left  Non-negative integer. Number of samples to add at the left.
#' @param n_right Non-negative integer. Number of samples to add at the right.
#'
#' @return Numeric vector of length length(x) + n_left + n_right.
#' @keywords internal
.reflect_limited_pad <- function(x, n_left, n_right) {
  
  n <- length(x)
  
  # ---- Left padding ----
  # Zero fill when n_left exceeds what reflection can provide (MNE line: l_z_pad)
  l_z_pad <- rep(0.0, max(n_left - n + 1L, 0L))
  
  # Reflected samples: Python x[n_left:0:-1] → R indices min(n_left+1,n) down to 2
  if (n_left > 0L) {
    left_idx     <- seq(from = min(n_left + 1L, n), to = 2L, by = -1L)
    left_reflect <- 2.0 * x[1L] - x[left_idx]
  } else {
    left_reflect <- numeric(0)
  }
  
  # ---- Right padding ----
  # Zero fill when n_right exceeds what reflection can provide (MNE line: r_z_pad)
  r_z_pad <- rep(0.0, max(n_right - n + 1L, 0L))
  
  # Reflected samples: Python x[-2:-n_right-2:-1] → R indices (n-1) down to max(n-n_right, 1)
  if (n_right > 0L) {
    right_idx     <- seq(from = n - 1L, to = max(n - n_right, 1L), by = -1L)
    right_reflect <- 2.0 * x[n] - x[right_idx]
  } else {
    right_reflect <- numeric(0)
  }
  
  # ---- Assemble: [l_z_pad | left_reflect | x | right_reflect | r_z_pad] ----
  c(l_z_pad, left_reflect, x, right_reflect, r_z_pad)
}

# ----------------------------------------------------------------------------
# .overlap_add_filter() - takes all above helper and does the filtering
# ----------------------------------------------------------------------------

#' Apply FIR filter via overlap-add FFT convolution (MNE equivalent)
#'
#' Replicates MNE's _overlap_add_filter() + _1d_overlap_filter() from
#' mne/filter.py exactly, for a single 1-D signal vector.
#'
#' Pipeline:
#'   1. Compute padding length: n_edge = max(min(n_h, n_sig) - 1, 0)
#'   2. Select optimal FFT block size (cost-minimising power of 2)
#'   3. Pad signal with reflect_limited by n_edge on each side
#'   4. Overlap-add FFT convolution loop (block by block)
#'   5. shift = (n_h-1)//2 + n_edge (zero-phase correction + unpadding in one)
#'   6. Trim output back to original signal length
#'
#' @param x     Numeric vector. The signal to filter (one channel).
#' @param h     Numeric vector. The FIR kernel coefficients.
#' @param phase Character. "zero" (default) for zero-phase single-pass.
#' @param pad   Character. Padding mode. Default "reflect_limited".
#'
#' @return Numeric vector of same length as x — the filtered signal.
#' @keywords internal
.overlap_add_filter <- function(x, h, phase = "zero", pad = "reflect_limited") {
  
  n_h     <- length(h)
  n_sig   <- length(x)
  
  # ---- Trivial case: kernel of length 1 (MNE line 295-296) ----
  if (n_h == 1L) return(x * h)
  
  # ---- Padding length (MNE line 297) ----
  # n_edge = max(min(n_h, n_sig) - 1, 0)
  n_edge <- max(min(n_h, n_sig) - 1L, 0L)
  n_x    <- n_sig + 2L * n_edge   # length of padded signal
  
  # ---- Select FFT block size (MNE lines 305-326) ----
  # Use cost-minimising power of 2 in range [2*n_h-1, n_x]
  min_fft <- 2L * n_h - 1L
  
  if (n_x >= min_fft) {
    # Candidate powers of 2 from ceil(log2(min_fft)) to ceil(log2(n_x))
    exp_lo <- ceiling(log2(min_fft))
    exp_hi <- ceiling(log2(n_x))
    N      <- 2L ^ seq(exp_lo, exp_hi)
    
    # Cost function: multiplications + heuristic penalty for long FFTs
    cost   <- ceiling(n_x / (N - n_h + 1)) * N * (log2(N) + 1) +
      4e-5 * N * n_x
    n_fft  <- N[which.min(cost)]
  } else {
    # Signal too short for multiple blocks — use one block
    n_fft <- .next_fast_len(min_fft)
  }
  
  # ---- Precompute FFT of kernel (MNE: cuda_dict["h_fft"] = rfft(h, n_fft)) ----
  h_padded <- c(h, rep(0.0, n_fft - n_h))
  H_fft    <- fft(h_padded)
  
  # ---- Pad signal (MNE line 359) ----
  x_ext     <- .reflect_limited_pad(x, n_edge, n_edge)
  x_filtered <- rep(0.0, n_x)
  
  # ---- Overlap-add parameters (MNE lines 363-365) ----
  n_seg      <- n_fft - n_h + 1L
  n_segments <- ceiling(n_x / n_seg)
  
  # shift combines: zero-phase group delay + unpadding offset (MNE line 365)
  # phase="zero" → (n_h-1)//2, else 0
  group_delay <- if (startsWith(phase, "zero")) (n_h - 1L) %/% 2L else 0L
  shift       <- group_delay + n_edge
  
  # ---- Overlap-add loop (MNE lines 369-381) ----
  # All indexing kept in 0-based to match MNE exactly, converted to R at access
  for (seg_idx in seq(0L, n_segments - 1L)) {
    
    # 0-indexed segment boundaries in x_ext
    start_0 <- seg_idx * n_seg
    stop_0  <- (seg_idx + 1L) * n_seg   # exclusive
    
    # Extract segment and zero-pad to n_fft (MNE lines 372-373)
    seg_r   <- x_ext[(start_0 + 1L) : min(stop_0, n_x)]
    seg_pad <- c(seg_r, rep(0.0, n_fft - length(seg_r)))
    
    # FFT multiply (MNE line 375: _fft_multiply_repeated)
    prod <- Re(fft(fft(seg_pad) * H_fft, inverse = TRUE)) / n_fft
    
    # Output placement with shift (MNE lines 377-381)
    start_filt_0 <- max(0L, start_0 - shift)
    stop_filt_0  <- min(start_0 - shift + n_fft, n_x)  # exclusive
    start_prod_0 <- max(0L, shift - start_0)
    stop_prod_0  <- start_prod_0 + (stop_filt_0 - start_filt_0)  # exclusive
    
    if (stop_filt_0 > start_filt_0) {
      filt_r <- (start_filt_0 + 1L) : stop_filt_0
      prod_r <- (start_prod_0 + 1L) : stop_prod_0
      x_filtered[filt_r] <- x_filtered[filt_r] + prod[prod_r]
    }
  }
  
  # ---- Remove padding and return (MNE line 384) ----
  # x_filtered[:n_x - 2*n_edge] = x_filtered[:n_sig]
  x_filtered[seq_len(n_x - 2L * n_edge)]
}

# ============================================================================
#                        MAIN EXPORTED FUNCTIONS
# ============================================================================

# ----------------------------------------------------------------------------
# eeg_bandpass()
# ----------------------------------------------------------------------------

#' EEG Bandpass / Highpass / Lowpass FIR Filter
#'
#' Applies a zero-phase FIR filter to EEG data, designed to match MNE-Python's
#' default filtering pipeline exactly (method='fir', phase='zero',
#' fir_window='hamming', fir_design='firwin', pad='reflect_limited').
#'
#' Supports highpass-only (h_freq=NULL), lowpass-only (l_freq=NULL),
#' or bandpass (both specified).
#'
#' @param eeg_obj          An object of class 'eeg'.
#' @param l_freq           Numeric or NULL. Highpass cutoff in Hz.
#'                         NULL = no highpass (lowpass only).
#' @param h_freq           Numeric or NULL. Lowpass cutoff in Hz.
#'                         NULL = no lowpass (highpass only).
#' @param l_trans_bandwidth Numeric or "auto". Highpass transition bandwidth
#'                         in Hz. "auto" uses MNE's rule:
#'                         min(max(0.25*l_freq, 2.0), l_freq).
#'                         Default: "auto".
#' @param h_trans_bandwidth Numeric or "auto". Lowpass transition bandwidth
#'                         in Hz. "auto" uses MNE's rule:
#'                         min(max(0.25*h_freq, 2.0), nyquist-h_freq).
#'                         Default: "auto".
#' @param filter_length    Integer or "auto". FIR kernel length in samples.
#'                         "auto" uses MNE's formula:
#'                         ceil(3.3 / min_trans_bw * sfreq), forced odd.
#'                         Default: "auto".
#' @param fir_window       Character. Window function for kernel design.
#'                         Default: "hamming" (matches MNE).
#' @param phase            Character. "zero" (default) applies single-pass
#'                         with group delay correction (matches MNE phase='zero').
#' @param pad              Character. Padding mode. Default: "reflect_limited"
#'                         (matches MNE default).
#' @param channels         Character or integer vector, or NULL. Channel names
#'                         or indices to filter. NULL = all channels.
#' @param verbose          Logical. Print filter parameters and progress.
#'                         Default: TRUE.
#'
#' @return An object of class 'eeg' with filtered data and updated
#'   preprocessing_history.
#'
#' @export
eeg_bandpass <- function(eeg_obj,
                         l_freq              = NULL,
                         h_freq              = NULL,
                         l_trans_bandwidth   = "auto",
                         h_trans_bandwidth   = "auto",
                         filter_length       = "auto",
                         fir_window          = "hamming",
                         phase               = "zero",
                         pad                 = "reflect_limited",
                         channels            = NULL,
                         verbose             = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (is.null(l_freq) && is.null(h_freq)) {
    stop("At least one of 'l_freq' or 'h_freq' must be specified.", call. = FALSE)
  }
  if (!is.null(l_freq) && (!is.numeric(l_freq) || l_freq <= 0)) {
    stop("'l_freq' must be a positive numeric value.", call. = FALSE)
  }
  
  sfreq   <- eeg_obj$sampling_rate
  nyquist <- sfreq / 2.0
  
  if (!is.null(h_freq) && (!is.numeric(h_freq) || h_freq <= 0 || h_freq >= nyquist)) {
    stop("'h_freq' must be a positive numeric value less than the Nyquist frequency (",
         nyquist, " Hz).", call. = FALSE)
  }
  if (!is.null(l_freq) && !is.null(h_freq) && l_freq >= h_freq) {
    stop("'l_freq' (", l_freq, ") must be less than 'h_freq' (", h_freq, ").",
         call. = FALSE)
  }
  
  # ========== RESOLVE CHANNEL INDICES ==========
  
  n_channels <- nrow(eeg_obj$data)
  chan_names <- eeg_obj$channels
  
  if (is.null(channels)) {
    ch_idx <- seq_len(n_channels)
  } else if (is.character(channels)) {
    ch_idx <- match(channels, chan_names)
    if (any(is.na(ch_idx))) {
      stop("Channel(s) not found: ",
           paste(channels[is.na(ch_idx)], collapse = ", "), call. = FALSE)
    }
  } else if (is.numeric(channels)) {
    ch_idx <- as.integer(channels)
    if (any(ch_idx < 1L) || any(ch_idx > n_channels)) {
      stop("Channel indices out of range [1, ", n_channels, "].", call. = FALSE)
    }
  } else {
    stop("'channels' must be NULL, a character vector, or an integer vector.",
         call. = FALSE)
  }
  
  # ========== COMPUTE TRANSITION BANDWIDTHS ==========
  
  if (identical(l_trans_bandwidth, "auto")) {
    bw <- .auto_trans_bandwidth(l_freq, h_freq, sfreq)
    l_trans_bw <- bw$l_trans_bw
    h_trans_bw <- bw$h_trans_bw
  } else {
    l_trans_bw <- if (!is.null(l_freq)) as.numeric(l_trans_bandwidth) else NULL
    h_trans_bw <- if (!is.null(h_freq)) as.numeric(h_trans_bandwidth) else NULL
  }
  
  # ========== COMPUTE FILTER LENGTH ==========
  
  if (identical(filter_length, "auto")) {
    N <- .auto_filter_length(l_trans_bw, h_trans_bw, sfreq)
  } else {
    N <- as.integer(filter_length)
    if (N %% 2L == 0L) N <- N + 1L  # force odd
  }
  
  # ========== BUILD FREQ / GAIN ARRAYS ==========
  # MNE constructs these from the passband and stopband edges.
  # All frequencies normalised to [0, 1] where 1 = Nyquist.
  
  if (!is.null(l_freq) && !is.null(h_freq)) {
    # Bandpass
    f_s1 <- (l_freq - l_trans_bw) / nyquist
    f_p1 <- l_freq / nyquist
    f_p2 <- h_freq / nyquist
    f_s2 <- (h_freq + h_trans_bw) / nyquist
    freq  <- c(0, f_s1, f_p1, f_p2, f_s2, 1)
    gain  <- c(0,    0,    1,    1,    0, 0)
    filter_type <- paste0("Bandpass (", l_freq, " - ", h_freq, " Hz)")
  } else if (!is.null(l_freq)) {
    # Highpass only
    f_s1 <- (l_freq - l_trans_bw) / nyquist
    f_p1 <- l_freq / nyquist
    freq  <- c(0, f_s1, f_p1, 1)
    gain  <- c(0,    0,    1, 1)
    filter_type <- paste0("Highpass (", l_freq, " Hz)")
  } else {
    # Lowpass only
    f_p  <- h_freq / nyquist
    f_s  <- (h_freq + h_trans_bw) / nyquist
    freq <- c(0,  f_p, f_s, 1)
    gain <- c(1,    1,   0, 0)
    filter_type <- paste0("Lowpass (", h_freq, " Hz)")
  }
  
  # Clamp to [0, 1] to avoid floating point edge violations
  freq <- pmax(pmin(freq, 1.0), 0.0)
  
  # ========== DESIGN KERNEL ==========
  
  h <- .firwin_kernel(N, freq, gain)
  
  # ========== VERBOSE HEADER ==========
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("EEG FIR Filtering\n")
    cat(strrep("=", 70), "\n")
    cat("  Filter type:          ", filter_type, "\n")
    cat("  Sampling rate:        ", sfreq, " Hz\n")
    cat("  Nyquist:              ", nyquist, " Hz\n")
    if (!is.null(l_trans_bw))
      cat("  HP transition BW:     ", round(l_trans_bw, 4), " Hz\n")
    if (!is.null(h_trans_bw))
      cat("  LP transition BW:     ", round(h_trans_bw, 4), " Hz\n")
    cat("  Filter length:        ", N, " samples (", round(N/sfreq, 3), " s)\n")
    cat("  Window:               ", fir_window, "\n")
    cat("  Phase:                ", phase, "\n")
    cat("  Padding:              ", pad, "\n")
    cat("  Channels:             ", length(ch_idx), "/", n_channels, "\n")
    cat(strrep("-", 70), "\n")
  }
  
  # ========== APPLY FILTER CHANNEL BY CHANNEL ==========
  
  filtered_data <- eeg_obj$data
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }
  
  for (idx in seq_along(ch_idx)) {
    i <- ch_idx[idx]
    filtered_data[i, ] <- .overlap_add_filter(
      x     = filtered_data[i, ],
      h     = h,
      phase = phase,
      pad   = pad
    )
    if (verbose) setTxtProgressBar(pb, idx)
  }
  
  if (verbose) {
    close(pb)
    cat("\n", strrep("=", 70), "\n")
  }
  
  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========
  
  history_entry <- paste0(
    "FIR bandpass: ", filter_type,
    " | window: ", fir_window,
    " | length: ", N, " samples",
    " | phase: ", phase,
    " | pad: ", pad,
    " | channels: ", length(ch_idx), "/", n_channels
  )
  
  # ========== RETURN UPDATED EEG OBJECT ==========
  
  new_eeg(
    data                 = filtered_data,
    channels             = eeg_obj$channels,
    sampling_rate        = eeg_obj$sampling_rate,
    times                = eeg_obj$times,
    events               = eeg_obj$events,
    metadata             = eeg_obj$metadata,
    reference            = eeg_obj$reference,
    preprocessing_history = c(eeg_obj$preprocessing_history,
                              list(history_entry))
  )
}



# ----------------------------------------------------------------------------
# eeg_notch()
# ----------------------------------------------------------------------------
#
#' EEG Notch (Band-Stop) FIR Filter
#'
#' Removes narrow frequency bands from EEG data to eliminate electrical line
#' noise. Designed to match MNE-Python's notch_filter() exactly
#' (method='fir', phase='zero', fir_window='hamming', fir_design='firwin',
#' pad='reflect_limited').
#'
#' Internally implements notch as a bandstop FIR filter, exactly as MNE does:
#'   lows  = freq - notch_width/2 - trans_bandwidth/2
#'   highs = freq + notch_width/2 + trans_bandwidth/2
#' Then builds a bandstop kernel (gain=0 in the notch, gain=1 elsewhere).
#'
#' Multiple frequencies (e.g. harmonics) are each filtered sequentially.
#'
#' @param eeg_obj          An object of class 'eeg'.
#' @param freqs            Numeric vector. Frequencies to notch in Hz.
#'                         E.g. c(50, 100, 150) for line noise + harmonics.
#'                         Default: 50.
#' @param notch_widths     Numeric or NULL. Width of the stop band in Hz
#'                         centred at each freq. NULL uses MNE's default:
#'                         freq / 200. Default: NULL.
#' @param trans_bandwidth  Numeric. Width of each transition band in Hz.
#'                         Default: 1 (matches MNE default).
#' @param filter_length    Integer or "auto". Kernel length in samples.
#'                         "auto" uses MNE's formula. Default: "auto".
#' @param fir_window       Character. Window function. Default: "hamming".
#' @param phase            Character. "zero" for zero-phase single-pass.
#'                         Default: "zero".
#' @param pad              Character. Padding mode. Default: "reflect_limited".
#' @param channels         Character or integer vector, or NULL. Channels to
#'                         filter. NULL = all channels. Default: NULL.
#' @param verbose          Logical. Print parameters and progress.
#'                         Default: TRUE.
#'
#' @return An object of class 'eeg' with filtered data and updated
#'   preprocessing_history.
#'
#' @export
eeg_notch <- function(eeg_obj,
                      freqs            = 50,
                      notch_widths     = NULL,
                      trans_bandwidth  = 1,
                      filter_length    = "auto",
                      fir_window       = "hamming",
                      phase            = "zero",
                      pad              = "reflect_limited",
                      channels         = NULL,
                      verbose          = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!is.numeric(freqs) || length(freqs) == 0 || any(freqs <= 0)) {
    stop("'freqs' must be a non-empty numeric vector of positive values.",
         call. = FALSE)
  }
  if (!is.numeric(trans_bandwidth) || trans_bandwidth <= 0) {
    stop("'trans_bandwidth' must be a positive numeric value.", call. = FALSE)
  }
  
  sfreq   <- eeg_obj$sampling_rate
  nyquist <- sfreq / 2.0
  
  if (any(freqs >= nyquist)) {
    stop("All 'freqs' must be less than the Nyquist frequency (",
         nyquist, " Hz).", call. = FALSE)
  }
  
  # ========== RESOLVE NOTCH WIDTHS (MNE line 1536-1537) ==========
  # Default: freq / 200 per frequency
  
  if (is.null(notch_widths)) {
    notch_widths <- freqs / 200.0
  } else {
    notch_widths <- rep(as.numeric(notch_widths), length.out = length(freqs))
  }
  
  # ========== RESOLVE CHANNEL INDICES ==========
  
  n_channels <- nrow(eeg_obj$data)
  chan_names <- eeg_obj$channels
  
  if (is.null(channels)) {
    ch_idx <- seq_len(n_channels)
  } else if (is.character(channels)) {
    ch_idx <- match(channels, chan_names)
    if (any(is.na(ch_idx))) {
      stop("Channel(s) not found: ",
           paste(channels[is.na(ch_idx)], collapse = ", "), call. = FALSE)
    }
  } else if (is.numeric(channels)) {
    ch_idx <- as.integer(channels)
    if (any(ch_idx < 1L) || any(ch_idx > n_channels)) {
      stop("Channel indices out of range [1, ", n_channels, "].", call. = FALSE)
    }
  } else {
    stop("'channels' must be NULL, a character vector, or an integer vector.",
         call. = FALSE)
  }
  
  # ========== VERBOSE HEADER ==========
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("EEG FIR Notch Filtering\n")
    cat(strrep("=", 70), "\n")
    cat("  Notch frequencies:    ", paste(freqs, collapse=", "), " Hz\n")
    cat("  Notch widths:         ", paste(round(notch_widths,4), collapse=", "), " Hz\n")
    cat("  Transition bandwidth: ", trans_bandwidth, " Hz\n")
    cat("  Sampling rate:        ", sfreq, " Hz\n")
    cat("  Window:               ", fir_window, "\n")
    cat("  Phase:                ", phase, "\n")
    cat("  Padding:              ", pad, "\n")
    cat("  Channels:             ", length(ch_idx), "/", n_channels, "\n")
    cat(strrep("-", 70), "\n")
  }
  
  # ========== APPLY ONE NOTCH PER FREQUENCY ==========
  
  filtered_data  <- eeg_obj$data
  tb_2           <- trans_bandwidth / 2.0
  history_parts  <- character(0)
  
  for (fi in seq_along(freqs)) {
    
    freq <- freqs[fi]
    nw   <- notch_widths[fi]
    
    # MNE lines 1552-1553: bandstop edges
    low_edge  <- freq - nw / 2.0 - tb_2   # lower passband edge (below notch)
    high_edge <- freq + nw / 2.0 + tb_2   # upper passband edge (above notch)
    
    # Validate edges
    if (low_edge <= 0) {
      stop("Notch lower edge (", round(low_edge, 4), " Hz) is <= 0 for freq=",
           freq, " Hz. Reduce trans_bandwidth or notch_width.", call. = FALSE)
    }
    if (high_edge >= nyquist) {
      stop("Notch upper edge (", round(high_edge, 4), " Hz) >= Nyquist (",
           nyquist, " Hz) for freq=", freq, " Hz.", call. = FALSE)
    }
    
    # Transition bandwidths for this notch (both sides equal = tb_2)
    l_trans_bw <- tb_2
    h_trans_bw <- tb_2
    
    # Compute filter length
    if (identical(filter_length, "auto")) {
      N <- .auto_filter_length(l_trans_bw, h_trans_bw, sfreq)
    } else {
      N <- as.integer(filter_length)
      if (N %% 2L == 0L) N <- N + 1L
    }
    
    # Build freq/gain arrays for bandstop (MNE lines 1552-1558)
    # Bandstop = gain 1 outside notch, gain 0 inside notch
    #
    # Correct monotonically-increasing freq/gain:
    #   [0, f_p1, f_s1, f_s2, f_p2, nyquist]
    #   [1,    1,    0,    0,    1,       1]
    # where:
    #   f_p1 = low_edge           (lower passband→stopband edge)
    #   f_s1 = low_edge  + tb_2  (lower stopband start, transition inward)
    #   f_s2 = high_edge - tb_2  (upper stopband end,   transition inward)
    #   f_p2 = high_edge          (upper stopband→passband edge)
    #
    # This produces: h = delta - lowpass(midpoint_upper) + lowpass(midpoint_lower)
    # = allpass - lowpass(50.375Hz) + lowpass(49.625Hz) = bandstop ✓
    
    f_p1 <- low_edge / nyquist
    f_s1 <- (low_edge  + tb_2) / nyquist
    f_s2 <- (high_edge - tb_2) / nyquist
    f_p2 <- high_edge / nyquist
    
    freq_arr <- c(0, f_p1, f_s1, f_s2, f_p2, 1)
    gain_arr <- c(1,    1,    0,    0,    1, 1)
    
    # Clamp to [0, 1]
    freq_arr <- pmax(pmin(freq_arr, 1.0), 0.0)
    
    # Design kernel
    h <- .firwin_kernel(N, freq_arr, gain_arr)
    
    if (verbose) {
      cat(sprintf("  Notch %d/%d: %.1f Hz  [%.4f - %.4f Hz stop]  kernel=%d samples\n",
                  fi, length(freqs), freq, low_edge, high_edge, N))
      pb <- txtProgressBar(min = 0, max = length(ch_idx),
                           style = 3, width = 50, char = "=")
    }
    
    # Apply to each channel
    for (idx in seq_along(ch_idx)) {
      i <- ch_idx[idx]
      filtered_data[i, ] <- .overlap_add_filter(
        x     = filtered_data[i, ],
        h     = h,
        phase = phase,
        pad   = pad
      )
      if (verbose) setTxtProgressBar(pb, idx)
    }
    
    if (verbose) {
      close(pb)
      cat("\n")
    }
    
    history_parts <- c(history_parts,
                       paste0(freq, " Hz (±", round(nw/2 + tb_2, 4), " Hz)"))
  }
  
  if (verbose) cat(strrep("=", 70), "\n")
  
  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========
  
  history_entry <- paste0(
    "FIR notch: ", paste(history_parts, collapse = ", "),
    " | trans_bw: ", trans_bandwidth, " Hz",
    " | window: ", fir_window,
    " | phase: ", phase,
    " | pad: ", pad,
    " | channels: ", length(ch_idx), "/", n_channels
  )
  
  # ========== RETURN UPDATED EEG OBJECT ==========
  
  new_eeg(
    data                  = filtered_data,
    channels              = eeg_obj$channels,
    sampling_rate         = eeg_obj$sampling_rate,
    times                 = eeg_obj$times,
    events                = eeg_obj$events,
    metadata              = eeg_obj$metadata,
    reference             = eeg_obj$reference,
    preprocessing_history = c(eeg_obj$preprocessing_history,
                              list(history_entry))
  )
}

# End of filter1.R
