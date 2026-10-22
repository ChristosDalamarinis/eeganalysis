#' ============================================================================
#'                       FIR Filtering Functions
#' ============================================================================
#'
#' This module implements FIR (Finite Impulse Response) filtering for EEGndata.
#'
#' Both bandpass and notch filtering are implemented using:
#'   - firwin / Hamming window kernel design (matching scipy.signal.firwin)
#'   - reflect_limited padding
#'   - Overlap-add FFT convolution 
#'   - Single-pass zero-phase correction 
#'
#' Depends on: eeg_class.R (for eeg object structure)
#'
#' Author: Christos Dalamarinis
#' Date: March 2026
#' Status: Tested
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
#' at which FFT is most efficient.
#'
#' @param target Positive integer. The minimum size required.
#' @return Integer. The smallest 5-smooth number >= target.
#' @keywords internal
.next_fast_len <- function(target) {
  
  # Lookup table for targets <= 10000 (same as SciPy)
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
  
  # Trivial case: target <= 6, return as-is 
  if (target <= 6L) return(target)
  
  # Already a power of 2
  if (bitwAnd(target, target - 1L) == 0L) return(target)
  
  # Small target: binary search in lookup table 
  if (target <= hams[length(hams)]) {
    idx <- which(hams >= target)[1]
    return(hams[idx])
  }
  
  # Large target: algorithmic search 
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
#' Compute automatic transition bandwidths 
#'
#' Computes the transition bandwidth for highpass and/or lowpass edges
#' using an automatic rule.
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
  
  # Highpass transition bandwidth
  l_trans_bw <- if (!is.null(l_freq)) {
    min(max(0.25 * l_freq, 2.0), l_freq)
  } else {
    NULL
  }
  
  # Lowpass transition bandwidth
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
#' Compute automatic FIR filter length 
#'
#' Computes the filter kernel length from transition bandwidths using an
#' automatic rule.
#'
#' Full chain:
#'   1. min_trans_bw = min(l_trans_bw, h_trans_bw)   (narrowest edge drives length)
#'   2. duration_s   = 3.3 / min_trans_bw             (3.3 = hamming length factor)
#'   3. n_samples    = ceil(duration_s x sfreq)
#'   4. force odd    = if even, add 1
#'
#' @param l_trans_bw  Numeric or NULL. Highpass transition bandwidth in Hz.
#' @param h_trans_bw  Numeric or NULL. Lowpass transition bandwidth in Hz.
#' @param sfreq       Numeric. Sampling rate in Hz.
#' @param window      Character. Window function name. One of "hamming",
#'                    "hann", "blackman". Default: "hamming".
#'
#' @return Integer. Filter kernel length in samples (always odd).
#' @keywords internal
.auto_filter_length <- function(l_trans_bw = NULL, h_trans_bw = NULL, sfreq,
                                window = "hamming") {

  # Collect non-NULL bandwidths and take the minimum
  bws <- c(l_trans_bw, h_trans_bw)
  if (length(bws) == 0) {
    stop("At least one of l_trans_bw or h_trans_bw must be non-NULL.",
         call. = FALSE)
  }
  min_trans_bw <- min(bws)

  # Per-window length factors
  # hamming=3.3, hann=3.1, blackman=5.0 - mult_fact=1.0 for firwin
  length_factor <- c(hamming = 3.3, hann = 3.1, blackman = 5.0)[window]
  if (is.na(length_factor)) {
    stop("Unknown window '", window, "'. Use 'hamming', 'hann', or 'blackman'.",
         call. = FALSE)
  }

  # Convert to duration in seconds, then to samples
  duration_s <- length_factor / min_trans_bw
  n_samples  <- ceiling(duration_s * sfreq)

  # Force odd length 
  # (n - 1) %% 2 == 1 only when n is even -> adds 1 to make it odd
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
#' using the specified window. Matches scipy.signal.firwin(N, cutoff_norm,
#' window=window, pass_zero=True, fs=2.0) exactly.
#'
#' @param N           Integer. Number of taps (must be odd).
#' @param cutoff_norm Numeric. Normalized cutoff in \[0, 1\] where 1 = Nyquist.
#' @param window      Character. One of "hamming", "hann", or "blackman".
#' @return Numeric vector of length N.
#' @keywords internal
.firwin_lowpass <- function(N, cutoff_norm, window = "hamming") {

  M     <- (N - 1) / 2
  n     <- seq(0L, N - 1L)
  alpha <- n - M

  # Ideal sinc impulse response (lowpass)
  h <- cutoff_norm * .sinc_norm(cutoff_norm * alpha)

  # Window function (coefficients from Ifeachor & Jervis)
  w <- switch(window,
    hamming  = 0.54 - 0.46 * cos(2 * pi * n / (N - 1L)),
    hann     = 0.50 - 0.50 * cos(2 * pi * n / (N - 1L)),
    blackman = 0.42 - 0.50 * cos(2 * pi * n / (N - 1L)) +
                      0.08 * cos(4 * pi * n / (N - 1L)),
    stop("'window' must be one of \"hamming\", \"hann\", or \"blackman\".",
         call. = FALSE)
  )

  # Apply window and normalize so DC gain = 1
  h <- h * w
  h <- h / sum(h)

  h
}


#' Build FIR kernel from freq/gain arrays
#'
#' Constructs a FIR kernel by combining lowpass firwin kernels, one per
#' gain transition in the freq/gain specification.
#'
#' Algorithm (right-to-left sweep):
#'   - When gain changes 0->1: ADD a lowpass kernel at the transition midpoint
#'   - When gain changes 1->0: SUBTRACT a lowpass kernel at the transition midpoint
#'   - No change: skip
#'
#' Every transition needs a kernel of a certain length (longer for a narrower
#' transition band); if that is more than \code{N} samples the function stops
#' with a message that says how many samples are needed.
#'
#' @param N      Integer. Total kernel length (must be odd).
#' @param freq   Numeric vector. Normalized frequencies in \{0, 1\} starting at 0
#'               and ending at 1 (Nyquist). E.g. c(0, 0.1, 0.2, 0.9, 1).
#' @param gain   Numeric vector. Gain at each frequency (0 or 1 only),
#'               same length as freq.
#' @param window Character. Window function. One of "hamming", "hann",
#'               or "blackman". Default: "hamming".
#' @param sfreq  Numeric or NULL. Sampling rate in Hz. Only used to write the
#'               width of the transition band in Hz into the message when
#'               \code{N} is too short; with NULL the width is left out.
#'
#' @return Numeric vector of length N - the FIR kernel coefficients.
#' @keywords internal
.firwin_kernel <- function(N, freq, gain, window = "hamming", sfreq = NULL) {
  
  stopifnot(freq[1] == 0)
  stopifnot(length(freq) == length(gain))
  stopifnot(N %% 2L == 1L)
  
  h      <- numeric(N)
  center <- N %/% 2L + 1L  # 1-indexed center (= Python's N//2)
  
  # If gain at Nyquist == 1, start with "all pass" delta
  if (gain[length(gain)] == 1) {
    h[center] <- 1.0
  }
  
  prev_freq <- freq[length(freq)]
  prev_gain <- gain[length(gain)]
  
  # Walk right to left, skipping the last element 
  for (i in seq(length(freq) - 1L, 1L)) {
    
    this_freq <- freq[i]
    this_gain <- gain[i]
    
    if (this_gain != prev_gain) {
      
      # Half-width of transition band 
      transition <- (prev_freq - this_freq) / 2.0
      
      # Per-transition kernel length - may be shorter than N 
      length_factor <- switch(window,
                              hamming  = 3.3,
                              hann     = 3.1,
                              blackman = 5.0,
                              3.3
      )
      this_N <- as.integer(round(length_factor / transition))
      this_N <- this_N + 1L - this_N %% 2L  # force odd
      
      # The kernel of this transition has to fit into N samples. If it does
      # not, the width asked for cannot be had: say so (before, a negative
      # position further down gave a message that pointed nowhere).
      if (this_N > N) {
        width <- if (is.null(sfreq)) {
          ""
        } else {
          paste0(formatC((prev_freq - this_freq) * sfreq / 2, format = "f",
                         digits = 2), " Hz ")
        }
        stop("The filter length of ", N, " samples is too short for the ",
             width, "transition band, which needs at least ", this_N,
             " samples. Use a longer 'filter_length' or a wider transition ",
             "band.", call. = FALSE)
      }

      # Cutoff at midpoint of transition band 
      cutoff <- (prev_freq + this_freq) / 2.0
      
      # Build lowpass kernel of length this_N 
      this_h <- .firwin_lowpass(this_N, cutoff, window)
      
      # Embed this_h centered in h 
      offset <- (N - this_N) %/% 2L
      idx    <- seq(offset + 1L, offset + this_N)
      
      if (this_gain == 0L) {
        h[idx] <- h[idx] - this_h  # subtract -> removes that band
      } else {
        h[idx] <- h[idx] + this_h  # add -> passes that band
      }
    }
    
    prev_gain <- this_gain
    prev_freq <- this_freq
  }
  
  h
}

# ----------------------------------------------------------------------------
# .minimum_phase() - preparation helper
# ----------------------------------------------------------------------------

#' Turn a symmetric FIR kernel into a minimum-phase kernel
#'
#' A symmetric (linear-phase) kernel looks equally far into the past and into
#' the future of every sample, which is why a zero-phase filter can leave a
#' trace before an event. This function builds the kernel with the same length
#' and the same magnitude response that puts all of its action on the past
#' side (a causal, minimum-phase filter), so nothing appears before an event.
#' It uses the cepstrum (homomorphic) method: take the logarithm of the
#' magnitude response, fold the cepstrum onto positive time, and go back with
#' the exponential. The FFT is zero-padded to a power of two of about 200
#' times the kernel length, so that the folding does not wrap around.
#'
#' @references
#' Oppenheim, A. V., & Schafer, R. W. (2010). Discrete-Time Signal Processing
#' (3rd ed.). Pearson. (chapter 13: cepstrum analysis and homomorphic
#' deconvolution)
#'
#' @param h Numeric vector. A symmetric FIR kernel, as built by
#'   \code{.firwin_kernel()}.
#' @return Numeric vector, same length as \code{h}.
#' @keywords internal
.minimum_phase <- function(h) {

  n     <- length(h)
  n_fft <- 2^ceiling(log2(2 * (n - 1) / 0.01))

  # magnitude response on a fine frequency grid (zero-padded FFT); a tiny floor
  # keeps the logarithm finite where the response is exactly zero
  mag <- Mod(fft(c(h, numeric(n_fft - n))))
  mag <- mag + 1e-7 * min(mag[mag > 0])

  # cepstrum = inverse FFT of the log magnitude (fft(inverse = TRUE) does not
  # divide by the length, so that is done here)
  cep <- Re(fft(log(mag), inverse = TRUE)) / n_fft

  # fold onto positive time: keep time 0, double the positive times, drop the
  # negative times (and the single sample in the middle: n_fft is a power of
  # two, so it is even)
  fold                  <- numeric(n_fft)
  fold[1]               <- 1
  fold[2:(n_fft %/% 2)] <- 2

  # back to the time domain: exp of the folded spectrum, inverse FFT, real part
  h_min <- Re(fft(exp(fft(cep * fold)), inverse = TRUE)) / n_fft

  h_min[seq_len(n)]                                  # same length as h
}

# ----------------------------------------------------------------------------
# .reflect_limited_pad() - preparation helper
# ----------------------------------------------------------------------------

#' 
#' reflect_limited pads with an odd reflection anchored on the edge values,
#' i.e.:  left pad  = 2*x\[1\]   - x\[reversed inner samples\]
#'         right pad = 2*x\[end\] - x\[reversed inner samples\]
#' If n_pad exceeds len(x)-1, zeros fill the remaining space.
#'
#' This differs from standard numpy 'reflect' (which excludes the edge value)
#' and 'symmetric' (which includes the edge value but does not subtract).
#'
#' Python equivalent (mne/cuda.py):
#'   l_z_pad = zeros(max(n_pad\[0\] - len(x) + 1, 0))
#'   r_z_pad = zeros(max(n_pad\[1\] - len(x) + 1, 0))
#'   out = [l_z_pad,
#'          2*x\[0\] - x\[n_pad\[0\]:0:-1\],
#'          x,
#'          2*x\[-1\] - x\[-2:-n_pad\[1\]-2:-1\],
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
  # Zero fill when n_left exceeds what reflection can provide 
  l_z_pad <- rep(0.0, max(n_left - n + 1L, 0L))
  
  # Reflected samples: Python x[n_left:0:-1] -> R indices min(n_left+1,n) down to 2
  if (n_left > 0L) {
    left_idx     <- seq(from = min(n_left + 1L, n), to = 2L, by = -1L)
    left_reflect <- 2.0 * x[1L] - x[left_idx]
  } else {
    left_reflect <- numeric(0)
  }
  
  # ---- Right padding ----
  # Zero fill when n_right exceeds what reflection can provide 
  r_z_pad <- rep(0.0, max(n_right - n + 1L, 0L))
  
  # Reflected samples: Python x[-2:-n_right-2:-1] -> R indices (n-1) down to max(n-n_right, 1)
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

#' Apply FIR filter via overlap-add FFT convolution 
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
#' @return Numeric vector of same length as x - the filtered signal.
#' @keywords internal
.overlap_add_filter <- function(x, h, phase = "zero", pad = "reflect_limited") {
  
  n_h     <- length(h)
  n_sig   <- length(x)
  
  # ---- Trivial case: kernel of length 1 ----
  if (n_h == 1L) return(x * h)
  
  # ---- Padding length ----
  # n_edge = max(min(n_h, n_sig) - 1, 0)
  n_edge <- max(min(n_h, n_sig) - 1L, 0L)
  n_x    <- n_sig + 2L * n_edge   # length of padded signal
  
  # ---- Select FFT block size ----
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
    # Signal too short for multiple blocks - use one block
    n_fft <- .next_fast_len(min_fft)
  }
  
  # ---- Precompute FFT of kernel (MNE: cuda_dict["h_fft"] = rfft(h, n_fft)) ----
  h_padded <- c(h, rep(0.0, n_fft - n_h))
  H_fft    <- fft(h_padded)

  # ---- Pad signal ----
  x_ext      <- .reflect_limited_pad(x, n_edge, n_edge)
  x_filtered <- rep(0.0, n_x)

  # ---- Overlap-add parameters ----
  n_seg      <- n_fft - n_h + 1L
  n_segments <- ceiling(n_x / n_seg)

  # shift combines: zero-phase group delay + unpadding offset
  # phase="zero" -> (n_h-1)//2, else 0
  group_delay <- if (startsWith(phase, "zero")) (n_h - 1L) %/% 2L else 0L
  shift       <- group_delay + n_edge

  # ---- Overlap-add loop ----
  # All indexing kept in 0-based 
  for (seg_idx in seq(0L, n_segments - 1L)) {

    # 0-indexed segment boundaries in x_ext
    start_0 <- seg_idx * n_seg
    stop_0  <- (seg_idx + 1L) * n_seg   # exclusive

    # Extract segment and zero-pad to n_fft 
    seg_r   <- x_ext[(start_0 + 1L) : min(stop_0, n_x)]
    seg_pad <- c(seg_r, rep(0.0, n_fft - length(seg_r)))

    # FFT multiply 
    prod <- Re(fft(fft(seg_pad) * H_fft, inverse = TRUE)) / n_fft
    
    # Output placement with shift 
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
  
  # ---- Remove padding and return ----
  # x_filtered[:n_x - 2*n_edge] = x_filtered[:n_sig]
  x_filtered[seq_len(n_x - 2L * n_edge)]
}


# ============================================================================
#                        IIR FILTERING (BUTTERWORTH)
# ============================================================================
#
# An IIR filter feeds some of its own past output back in, so a short recipe
# gives a long, smooth response. Here the filter is built from its poles and
# zeros and then cut into small "second-order sections" (two poles and two
# zeros each). The signal goes through the sections one after the other.
# Keeping the filter in small pieces keeps every number accurate; one long
# recipe (two lists of numbers) loses digits for slow filters, which makes the
# filter inaccurate or even unstable.
#
# References:
#   Butterworth, S. (1930). On the theory of filter amplifiers.
#     Wireless Engineer, 7, 536-541.
#   Oppenheim, A. V., & Schafer, R. W. (2010). Discrete-Time Signal
#     Processing (3rd ed.). Pearson. (frequency transformations, bilinear
#     transform, cascade of second-order sections)

# ----------------------------------------------------------------------------
# .check_phase() / .check_pad() / .resolve_iir_params() - argument checks
# ----------------------------------------------------------------------------

#' Check the 'phase' argument of a filter function
#'
#' Accepts only the words that are really implemented for the chosen method and
#' stops with a clear message for any other word. (A word that is not backed
#' by the matching filter would give a result that is shifted in time, without
#' any warning.)
#'
#' @param phase  The value given by the user.
#' @param method Character. "fir" or "iir".
#' @return \code{phase}, unchanged, if it is allowed.
#' @keywords internal
.check_phase <- function(phase, method) {
  allowed <- if (identical(method, "iir")) {
    c("zero", "forward")
  } else {
    c("zero", "minimum")
  }
  if (!is.character(phase) || length(phase) != 1L || is.na(phase) ||
      !(phase %in% allowed)) {
    stop("'phase' must be one of ", paste0("\"", allowed, "\"", collapse = ", "),
         " for an ", toupper(method), " filter.", call. = FALSE)
  }
  phase
}

#' Check the 'pad' argument of a filter function
#'
#' Only "reflect_limited" is implemented, so any other word is an error
#' instead of being silently ignored.
#'
#' @param pad The value given by the user.
#' @return \code{pad}, unchanged, if it is allowed.
#' @keywords internal
.check_pad <- function(pad) {
  if (!identical(pad, "reflect_limited")) {
    stop("'pad' must be \"reflect_limited\" (the only padding implemented).",
         call. = FALSE)
  }
  pad
}

#' Fill in and check the 'iir_params' argument
#'
#' @param iir_params NULL or a named list with the elements \code{order} and
#'   \code{ftype}. Missing elements take the defaults (order 4, "butter").
#' @return A list with the elements \code{order} (integer) and \code{ftype}.
#' @keywords internal
.resolve_iir_params <- function(iir_params) {
  out <- list(order = 4L, ftype = "butter")
  if (is.null(iir_params)) return(out)

  if (!is.list(iir_params) ||
      (length(iir_params) > 0L &&
       (is.null(names(iir_params)) || any(names(iir_params) == "")))) {
    stop("'iir_params' must be NULL or a named list, ",
         "e.g. list(order = 4, ftype = \"butter\").", call. = FALSE)
  }
  unknown <- setdiff(names(iir_params), names(out))
  if (length(unknown) > 0L) {
    stop("Unknown element(s) in 'iir_params': ", paste(unknown, collapse = ", "),
         ". Allowed: order, ftype.", call. = FALSE)
  }
  for (nm in names(iir_params)) out[[nm]] <- iir_params[[nm]]

  if (!is.numeric(out$order) || length(out$order) != 1L || is.na(out$order) ||
      out$order < 1 || out$order != round(out$order)) {
    stop("'iir_params$order' must be a single whole number >= 1.", call. = FALSE)
  }
  if (!identical(out$ftype, "butter")) {
    stop("'iir_params$ftype' must be \"butter\" ",
         "(the only filter family implemented so far).", call. = FALSE)
  }
  out$order <- as.integer(out$order)
  out
}

#' Check the cutoff frequencies of a band-pass, high-pass or low-pass filter
#'
#' The checks that \code{eeg_bandpass()} and \code{filter_data()} have in
#' common. They need only the cutoffs and the sampling rate, nothing from an
#' eeg object, which is why both functions can call this one function instead
#' of each carrying its own copy. It is the first thing either of them checks
#' about the filter itself (see the overview before
#' \code{.make_bandpass_filter()}).
#'
#' @param l_freq,h_freq Cutoffs in Hz as given by the user (NULL = not used).
#' @param sfreq Numeric. Sampling rate in Hz.
#' @return NULL, invisibly. Stops with a message if a cutoff is not allowed.
#' @keywords internal
.check_band_edges <- function(l_freq, h_freq, sfreq) {

  # The Nyquist frequency is half the sampling rate: the highest frequency a
  # signal sampled at 'sfreq' can contain. A low-pass cutoff has to lie below it.
  nyquist <- sfreq / 2.0

  # Which kind of filter is wanted depends on which cutoffs are given:
  #   l_freq only -> high-pass,  h_freq only -> low-pass,  both -> band-pass.
  # So at least one of the two is needed.
  if (is.null(l_freq) && is.null(h_freq)) {
    stop("At least one of 'l_freq' or 'h_freq' must be specified.", call. = FALSE)
  }

  # High-pass cutoff: a positive number of Hz.
  if (!is.null(l_freq) && (!is.numeric(l_freq) || l_freq <= 0)) {
    stop("'l_freq' must be a positive numeric value.", call. = FALSE)
  }

  # Low-pass cutoff: a positive number of Hz, strictly below Nyquist.
  if (!is.null(h_freq) && (!is.numeric(h_freq) || h_freq <= 0 || h_freq >= nyquist)) {
    stop("'h_freq' must be a positive numeric value less than the Nyquist frequency (",
         nyquist, " Hz).", call. = FALSE)
  }

  # In a band-pass the high-pass edge has to lie below the low-pass edge,
  # otherwise nothing would be let through.
  if (!is.null(l_freq) && !is.null(h_freq) && l_freq >= h_freq) {
    stop("'l_freq' (", l_freq, ") must be less than 'h_freq' (", h_freq, ").",
         call. = FALSE)
  }

  # Nothing to hand back: the function either stops with a message or lets the
  # caller carry on.
  invisible(NULL)
}

#' Check the method-related settings of a band-pass filter call
#'
#' Checks 'phase', 'pad' and 'iir_params' for the chosen method and warns about
#' settings that the method does not use. The warning about FIR-only settings
#' has to know which of them the user typed on purpose; only the calling
#' function can see that (with \code{missing()}), so it is passed in.
#' \code{eeg_bandpass()} and \code{filter_data()} both call this function right
#' after \code{.check_band_edges()} and use the values it returns from there on.
#'
#' @param method Character. "fir" or "iir" (already matched).
#' @param phase,pad,iir_params The values given by the user.
#' @param fir_only_given Named logical vector with one element for each
#'   FIR-only setting (l_trans_bandwidth, h_trans_bandwidth, filter_length,
#'   fir_window): TRUE if the user set it.
#' @return A list with the checked \code{phase}, \code{pad} and
#'   \code{iir_params} (filled in for "iir", as given for "fir").
#' @keywords internal
.check_filter_settings <- function(method, phase, pad, iir_params, fir_only_given) {

  # 'phase' and 'pad' must be words that are really implemented for this
  # method (see .check_phase() and .check_pad()). Anything else stops here
  # instead of silently giving a result that is shifted in time.
  phase <- .check_phase(phase, method)
  pad   <- .check_pad(pad)

  if (method == "iir") {
    # Fill in the IIR settings that were not given (order 4, Butterworth) and
    # check the ones that were.
    iir_params <- .resolve_iir_params(iir_params)

    # The four FIR-only settings describe an FIR kernel and mean nothing for an
    # IIR filter. Only a setting the user typed on purpose deserves a warning (a
    # default does not). Telling the two apart needs missing(), which only the
    # function that owns the argument can ask: that function does it and passes
    # the answer in as 'fir_only_given'.
    if (any(fir_only_given)) {
      warning("Not used with method = \"iir\": ",
              paste(names(fir_only_given)[fir_only_given], collapse = ", "), ".",
              call. = FALSE)
    }
  } else if (!is.null(iir_params)) {
    # The other way round: 'iir_params' only means something for method = "iir".
    warning("'iir_params' is not used with method = \"fir\".", call. = FALSE)
  }

  # Hand the checked values back; the caller uses these from here on. For FIR,
  # 'iir_params' is returned as it was given (normally NULL).
  list(phase = phase, pad = pad, iir_params = iir_params)
}

# ----------------------------------------------------------------------------
# .butter_zpk() / .zpk_to_sos() / .design_iir_sos() - the design
# ----------------------------------------------------------------------------

#' Poles, zeros and gain of a digital Butterworth filter
#'
#' Builds the filter directly from its poles and zeros (no long polynomials) in
#' three steps: (1) the textbook analog low-pass filter has its poles on a
#' circle; (2) a frequency transformation moves them to the wanted cutoff(s)
#' and turns the low-pass into a high-pass, band-pass or band-stop; (3) the
#' bilinear transform converts the result to a digital filter.
#'
#' @param order Integer. Order of the low-pass prototype.
#' @param Wn    Numeric. Cutoff(s) as a fraction of the Nyquist frequency, in
#'   (0, 1): one value for "low" and "high", two (lower, upper) for "pass" and
#'   "stop".
#' @param type  Character. "low", "high", "pass" or "stop".
#' @return A list with \code{z} (zeros), \code{p} (poles, complex) and
#'   \code{k} (gain).
#' @keywords internal
.butter_zpk <- function(order, Wn, type) {

  # (1) analog low-pass prototype: cutoff 1 rad/s, poles on the left half of
  #     the unit circle, no finite zeros
  m <- seq(-order + 1, order - 1, by = 2)
  p <- -exp(1i * pi * m / (2 * order))
  z <- complex(0)
  k <- 1

  # pre-warp the cutoff(s): the bilinear transform squeezes the frequency axis,
  # so the analog cutoff is stretched first (sampling frequency taken as 2)
  fs     <- 2
  warped <- 2 * fs * tan(pi * Wn / fs)
  degree <- length(p) - length(z)            # zeros that sit at infinity

  # (2) move the poles to the cutoff(s) / change the filter type
  if (type == "low") {
    p <- p * warped
    k <- k * warped^degree
  } else if (type == "high") {
    k <- k * Re(prod(-z) / prod(-p))
    p <- warped / p
    z <- rep(0, degree)
  } else if (type == "pass") {
    bw <- warped[2] - warped[1]              # width of the band
    wo <- sqrt(warped[1] * warped[2])        # centre of the band
    p  <- p * bw / 2
    sp <- sqrt(p^2 - wo^2 + 0i)
    p  <- c(p + sp, p - sp)                  # every pole becomes two
    z  <- rep(0, degree)
    k  <- k * bw^degree
  } else {                                   # "stop"
    bw <- warped[2] - warped[1]
    wo <- sqrt(warped[1] * warped[2])
    k  <- k * Re(prod(-z) / prod(-p))
    p  <- (bw / 2) / p
    sp <- sqrt(p^2 - wo^2 + 0i)
    p  <- c(p + sp, p - sp)
    z  <- c(rep(1i * wo, degree), rep(-1i * wo, degree))
  }

  # (3) bilinear transform to the digital plane
  fs2    <- 2 * fs
  degree <- length(p) - length(z)            # zeros still at infinity land at -1
  list(z = c((fs2 + z) / (fs2 - z), rep(-1, degree)),
       p = (fs2 + p) / (fs2 - p),
       k = k * Re(prod(fs2 - z) / prod(fs2 - p)))
}

#' Pair poles and zeros two by two (second-order sections)
#'
#' Complex poles come in conjugate pairs and each pair becomes one section.
#' Real poles are paired two by two (a single left-over real pole becomes a
#' first-order section). The zeros are paired the same way. The total gain goes
#' into the first section.
#'
#' @param z,p Complex vectors. Zeros and poles.
#' @param k   Numeric. Gain.
#' @return A matrix with one row per section and the columns b0 b1 b2 a0 a1 a2
#'   (a0 is always 1).
#' @keywords internal
.zpk_to_sos <- function(z, p, k) {

  tol <- 1e-8

  pair_up <- function(r) {
    upper <- r[Im(r) > tol]                          # one of each conjugate pair
    real  <- sort(Re(r[abs(Im(r)) <= tol]))
    secs  <- lapply(upper, function(rt) c(-2 * Re(rt), Mod(rt)^2))
    i <- 1L
    j <- length(real)
    while (i < j) {                                  # smallest with largest
      secs[[length(secs) + 1L]] <- c(-(real[i] + real[j]), real[i] * real[j])
      i <- i + 1L
      j <- j - 1L
    }
    if (i == j) {                                    # one left over
      secs[[length(secs) + 1L]] <- c(-real[i], 0)
    }
    secs
  }

  p_secs <- pair_up(p)
  z_secs <- pair_up(z)
  stopifnot(length(p_secs) == length(z_secs))

  sos <- matrix(0, nrow = length(p_secs), ncol = 6L)
  for (s in seq_along(p_secs)) {
    sos[s, ] <- c(1, z_secs[[s]], 1, p_secs[[s]])
  }
  sos[1L, 1:3] <- sos[1L, 1:3] * k
  sos
}

#' Design the IIR filter asked for by eeg_bandpass()
#'
#' @param l_freq,h_freq Cutoffs in Hz, or NULL: both = band-pass, only
#'   \code{l_freq} = high-pass, only \code{h_freq} = low-pass.
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param order Integer. Filter order.
#' @return A matrix of second-order sections (see \code{.zpk_to_sos()}).
#' @keywords internal
.design_iir_sos <- function(l_freq, h_freq, sfreq, order = 4L) {
  nyquist <- sfreq / 2
  if (!is.null(l_freq) && !is.null(h_freq)) {
    type <- "pass"
    Wn   <- c(l_freq, h_freq) / nyquist
  } else if (!is.null(l_freq)) {
    type <- "high"
    Wn   <- l_freq / nyquist
  } else {
    type <- "low"
    Wn   <- h_freq / nyquist
  }
  zpk <- .butter_zpk(order, Wn, type)
  .zpk_to_sos(zpk$z, zpk$p, zpk$k)
}

# ----------------------------------------------------------------------------
# .sos_pass() / .iir_ringing_samples() / .iir_filter() - applying the filter
# ----------------------------------------------------------------------------

#' Run a signal once through all second-order sections
#'
#' Each section looks at the last two inputs and the last two outputs:
#' \code{y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]}.
#' Base R's \code{stats::filter()} does the two halves of that rule (the input
#' part and the feedback part).
#'
#' @param sos    Matrix of second-order sections (b0 b1 b2 a0 a1 a2 per row).
#' @param x      Numeric vector. The signal.
#' @param steady Logical. FALSE: every section starts with an empty memory (what
#'   a causal filter does at the start of a recording). TRUE: every section
#'   starts as if the signal had been flat, equal to its first sample, for ever
#'   before; this removes the switch-on jolt of the forward-backward method.
#' @return Numeric vector, same length as \code{x}. A missing value (NA or NaN)
#'   in \code{x} makes every value from there on missing, because the filter
#'   feeds its own output back.
#' @keywords internal
.sos_pass <- function(sos, x, steady = FALSE) {

  y     <- x
  level <- x[1]                        # the flat level the sections have "seen"

  for (s in seq_len(nrow(sos))) {
    b <- sos[s, 1:3]
    a <- sos[s, 4:6]

    if (steady) {
      gain0   <- sum(b) / sum(a)       # what this section does to a flat signal
      mem_in  <- level
      mem_out <- gain0 * level
      level   <- gain0 * level         # the next section sees this flat level
    } else {
      mem_in  <- 0
      mem_out <- 0
    }

    # input part: b0*x[n] + b1*x[n-1] + b2*x[n-2]  (two remembered inputs first)
    # The first two places have no full window yet: cut off those two and only
    # those, so that a missing value in the signal stays in the result (and is
    # carried on by the feedback part) instead of shortening it.
    fir <- stats::filter(c(rep(mem_in, 2L), y), b, sides = 1)
    fir <- as.numeric(fir)[-(1:2)]

    # feedback part: minus a1*y[n-1] minus a2*y[n-2]  (two remembered outputs)
    y <- as.numeric(stats::filter(fir, -a[2:3], method = "recursive",
                                  init = rep(mem_out, 2L)))
  }
  y
}

#' How many samples an IIR filter keeps ringing after one spike
#'
#' Sends one spike through the filter and looks at the response in blocks of
#' 1000 samples. The ringing ends at the last sample of the last block that
#' still holds a value larger than 0.1% of the largest value seen so far. The
#' result is used as the length of the padding for zero-phase filtering. A slow
#' high-pass whose response after the first sample stays below 0.1% of its
#' peak gets 0 (no padding).
#'
#' @param sos      Matrix of second-order sections.
#' @param max_try  Integer. Longest response to look at, in samples.
#' @return Integer-valued number of samples.
#' @keywords internal
.iir_ringing_samples <- function(sos, max_try = 100000L) {

  per_chunk  <- 1000L
  n_chunks   <- ceiling(max_try / per_chunk)

  spike      <- numeric(per_chunk * n_chunks)
  spike[1L]  <- 1
  h          <- abs(.sos_pass(sos, spike))     # the whole response in one go

  last_good  <- per_chunk
  thresh_val <- 0

  for (ii in seq_len(n_chunks) - 1L) {           # block number, counted from 0
    block      <- h[(ii * per_chunk + 1L):((ii + 1L) * per_chunk)]
    thresh_val <- max(0.001 * max(block), thresh_val)
    big        <- which(block > thresh_val)
    if (length(big) > 0L) {
      last_good <- big[length(big)] - 1L         # position in the block, from 0
    } else {
      return((ii - 1L) * per_chunk + last_good)
    }
  }

  warning("could not properly estimate the ringing of the filter",
          call. = FALSE)
  per_chunk * n_chunks
}

#' Filter one channel with second-order sections
#'
#' \code{phase = "forward"}: one pass, nothing can appear before an event
#' (causal); the filter starts with an empty memory and nothing is padded.
#' \code{phase = "zero"}: the signal is padded at both ends (odd reflection),
#' run forward, the result is run backward, and the padding is cut off.
#'
#' @param x      Numeric vector. One channel.
#' @param sos    Matrix of second-order sections.
#' @param phase  "zero" or "forward".
#' @param padlen Integer. Samples of padding at each end ("zero" only).
#' @return Numeric vector, same length as \code{x}.
#' @keywords internal
.iir_filter <- function(x, sos, phase, padlen) {

  if (phase == "forward") {
    return(.sos_pass(sos, x))
  }

  padlen <- min(padlen, length(x) - 1L)     # never pad by more than the signal
  x_ext  <- if (padlen > 0L) .reflect_limited_pad(x, padlen, padlen) else x

  y <- .sos_pass(sos, x_ext, steady = TRUE)             # forward
  y <- rev(.sos_pass(sos, rev(y), steady = TRUE))       # backward

  if (padlen > 0L) y[(padlen + 1L):(length(y) - padlen)] else y
}

# ----------------------------------------------------------------------------
# THE BAND-PASS RECIPE, SHARED BY eeg_bandpass() AND filter_data()
# ----------------------------------------------------------------------------
#
# WHAT THIS SECTION DOES
#   Everything that turns "a band-pass, high-pass or low-pass filter with these
#   settings" into filtered numbers is here, once, for FIR and for IIR filters.
#   Two public functions use it:
#
#     filter_data(x, sfreq, ...)    numbers in, numbers out. No eeg object.
#     eeg_bandpass(eeg_obj, ...)    a whole eeg object in, a whole one out.
#
#   Through filter_data(), the places in the package that have to filter a loose
#   row of numbers (one that is not a channel of a recording) use it as well:
#   annotate_muscle() in annotations.R, and find_bads_eog() / find_bads_ecg() in
#   ica_detect.R (they filter the ICA components and the reference signal).
#
# THE TWO STEPS: DESIGN, THEN RUN
#   1. DESIGN  .make_bandpass_filter()   settings -> a description of the filter
#   2. RUN     .apply_bandpass_filter()  description + ONE row of numbers
#                                          -> that row, filtered
#
#   They are two functions, not one, because eeg_bandpass() needs the design
#   BEFORE it filters anything: it prints the filter's details in its header and
#   writes them in the history entry. filter_data() just does both steps, one
#   after the other.
#
# WHO CALLS WHOM
#   filter_data()  --> .check_band_edges()          are the cutoffs allowed?
#                  --> .check_filter_settings()     are phase / pad / iir_params ok?
#                  --> .make_bandpass_filter()      design, once
#                  --> .apply_bandpass_filter()     run, once per row of x
#
#   eeg_bandpass() --> the same four, in the same order, plus the work that
#                      belongs to an eeg object: the channel choice, the header
#                      and progress bar, the history entry, the new eeg object.
#
#   .make_bandpass_filter()  --> FIR: .auto_filter_length(), .firwin_kernel(),
#                                     .minimum_phase() (only for phase = "minimum")
#                            --> IIR: .design_iir_sos(), .iir_ringing_samples()
#   .apply_bandpass_filter() --> FIR: .overlap_add_filter()
#                            --> IIR: .iir_filter()
#
#   The functions on the right are older and were not changed. Only the two
#   engines, .overlap_add_filter() and .iir_filter() (with the small helpers
#   they call), ever see the signal values; everything above them works on
#   settings and decides what to run.
#
# THE "spec" LIST
#   .make_bandpass_filter() returns a plain list that describes the designed
#   filter. .apply_bandpass_filter() reads it, and eeg_bandpass() reads it for
#   its header and history entry.
#     always:  method ("fir"/"iir"), filter_type (label for display), phase, pad
#     FIR:     kernel (the filter coefficients), N (their number), l_trans_bw and
#              h_trans_bw (transition widths in Hz, NULL where not used), window
#     IIR:     order, sos (the chain of second-order sections), padlen (samples
#              of padding at each end), pad_text (the padding as a sentence)
#
# WHAT IS NOT HERE
#   eeg_notch() designs a band-STOP filter with its own cutoffs and its own code;
#   it shares only the engines and the phase / pad / iir_params checks.
#   downsample() calls .design_iir_sos(), .iir_ringing_samples() and
#   .iir_filter() directly for its anti-aliasing filter.

#' Design a band-pass, high-pass or low-pass filter
#'
#' Works out the filter from its settings, without touching any data. For an
#' FIR filter that means the transition widths, the length, the wanted gain at
#' each frequency and the kernel (the symmetric one, or the minimum-phase one
#' for \code{phase = "minimum"}); for an IIR filter it means the chain of
#' second-order sections and the padding length. The result goes to
#' \code{.apply_bandpass_filter()}.
#'
#' All arguments must have been checked before (see
#' \code{.check_band_edges()} and \code{.check_filter_settings()}). This is the
#' first of the two steps (design, then run) explained in the overview above.
#'
#' @param sfreq Numeric. Sampling rate in Hz.
#' @param l_freq,h_freq Cutoffs in Hz (NULL = not used).
#' @param method Character. "fir" or "iir".
#' @param iir_params List with \code{order}, as returned by
#'   \code{.resolve_iir_params()}. Used only for "iir".
#' @param l_trans_bandwidth,h_trans_bandwidth,filter_length,fir_window See
#'   \code{\link{eeg_bandpass}}. Used only for "fir".
#' @param phase Character. "zero" or "minimum" (FIR), "zero" or "forward" (IIR).
#' @param pad Character. "reflect_limited".
#' @return A list that describes the filter. Always: \code{method},
#'   \code{filter_type} (text such as "Bandpass (1 - 40 Hz)"), \code{phase} and
#'   \code{pad}. For FIR also \code{kernel}, \code{N} (its length),
#'   \code{l_trans_bw} and \code{h_trans_bw} (NULL where not used) and
#'   \code{window}. For IIR also \code{order}, \code{sos}, \code{padlen} and
#'   \code{pad_text}.
#' @keywords internal
.make_bandpass_filter <- function(sfreq, l_freq, h_freq, method, iir_params,
                                  l_trans_bandwidth, h_trans_bandwidth,
                                  filter_length, fir_window, phase, pad) {

  # ========== SETUP, THE SAME FOR FIR AND IIR ==========

  # The Nyquist frequency (half the sampling rate). Cutoffs and transition
  # widths are given in Hz; the FIR design below works with frequencies as a
  # fraction of this one (1 = Nyquist).
  nyquist <- sfreq / 2.0

  # A short label of the filter, such as "Bandpass (1 - 40 Hz)". It is the same
  # for FIR and IIR and is only for display: eeg_bandpass() prints it in its
  # header and writes it in the history entry.
  filter_type <- if (!is.null(l_freq) && !is.null(h_freq)) {
    paste0("Bandpass (", l_freq, " - ", h_freq, " Hz)")
  } else if (!is.null(l_freq)) {
    paste0("Highpass (", l_freq, " Hz)")
  } else {
    paste0("Lowpass (", h_freq, " Hz)")
  }

  # ========== IIR: SECTIONS, AND THE PADDING LENGTH ==========
  # An IIR (Butterworth) filter is designed as a chain of small second-order
  # sections, which stays accurate even for slow filters. Nothing else has to be
  # prepared, so the IIR design is finished here and the function returns; the
  # FIR code further down only runs for method = "fir".

  if (method == "iir") {
    # The chain of second-order sections. The cutoffs are the points where ONE
    # pass of the filter is 3 dB down; a band-pass has twice as many poles as
    # 'order' says.
    sos    <- .design_iir_sos(l_freq, h_freq, sfreq, order = iir_params$order)

    # Padding. With phase = "zero" the signal is run forward and then backward,
    # and both ends are padded first so that the first and last moments of the
    # recording are not distorted. The padding has to be as long as the filter
    # keeps ringing, which .iir_ringing_samples() works out from the sections.
    # With phase = "forward" (causal, one pass) nothing is padded.
    padlen <- if (phase == "zero") .iir_ringing_samples(sos) else 0L

    # The padding as a sentence, for the header and the history entry.
    pad_text <- if (padlen > 0L) {
      paste0("reflect_limited, ", padlen, " samples at each end")
    } else {
      "none"
    }

    return(list(method = "iir", filter_type = filter_type, phase = phase,
                pad = pad, order = iir_params$order, sos = sos,
                padlen = padlen, pad_text = pad_text))
  }

  # ========== FIR: TRANSITION BANDWIDTHS ==========
  # An FIR filter cannot jump from "blocked" to "let through" at the cutoff; it
  # needs a slope. The transition bandwidth is the width of that slope in Hz.
  # Typed values are used as they are. With "auto" (the default):
  #   high-pass slope: a quarter of the cutoff, at least 2 Hz, but never wider
  #                    than the cutoff itself;
  #   low-pass slope:  a quarter of the cutoff, at least 2 Hz, but never wider
  #                    than the room left between the cutoff and Nyquist.
  # (.auto_trans_bandwidth() holds the same rule; this code has always written
  # it out here, and it is used unchanged.)

  l_trans_bw <- if (is.null(l_freq)) {
    NULL
  } else if (identical(l_trans_bandwidth, "auto")) {
    min(max(0.25 * l_freq, 2.0), l_freq)
  } else {
    as.numeric(l_trans_bandwidth)
  }

  h_trans_bw <- if (is.null(h_freq)) {
    NULL
  } else if (identical(h_trans_bandwidth, "auto")) {
    min(max(0.25 * h_freq, 2.0), nyquist - h_freq)
  } else {
    as.numeric(h_trans_bandwidth)
  }

  # ========== FIR: FILTER LENGTH ==========
  # The kernel is the list of N coefficients that makes the filter. The narrower
  # the slope, the longer it has to be. "auto" works N out from the narrowest
  # slope (.auto_filter_length()). A typed length is made odd as well, because a
  # symmetric kernel needs a centre sample. (.firwin_kernel() stops with a clear
  # message if the length is too short for the slope.)

  if (identical(filter_length, "auto")) {
    N <- .auto_filter_length(l_trans_bw, h_trans_bw, sfreq, window = fir_window)
  } else {
    N <- as.integer(filter_length)
    if (N %% 2L == 0L) N <- N + 1L  # force odd
  }

  # ========== FIR: FREQ / GAIN ARRAYS ==========
  # The wanted shape of the filter as a list of points: 'freq' are frequencies
  # as a fraction of Nyquist (0 to 1), 'gain' is what should get through there
  # (0 = blocked, 1 = untouched). Between the points the gain changes along a
  # straight line, which is the slope of the transition band. The list must
  # start at frequency 0 and end at 1.

  if (!is.null(l_freq) && !is.null(h_freq)) {
    # Band-pass: blocked up to (l_freq - slope), rising to 1 at l_freq, 1 up to
    # h_freq, falling to 0 at (h_freq + slope).
    f_s1 <- (l_freq - l_trans_bw) / nyquist
    f_p1 <- l_freq / nyquist
    f_p2 <- h_freq / nyquist
    f_s2 <- (h_freq + h_trans_bw) / nyquist
    freq  <- c(f_s1, f_p1, f_p2, f_s2)
    gain  <- c(   0,    1,    1,    0)
    # Add the two end points (0 Hz and Nyquist) when the slopes do not already
    # reach them; blocked at both.
    if (f_s1 != 0) { freq <- c(0, freq); gain <- c(0, gain) }
    if (f_s2 != 1) { freq <- c(freq, 1); gain <- c(gain, 0) }
  } else if (!is.null(l_freq)) {
    # High-pass only: blocked up to (l_freq - slope), rising to 1 at l_freq,
    # and 1 all the way to Nyquist.
    f_s1 <- (l_freq - l_trans_bw) / nyquist
    f_p1 <- l_freq / nyquist
    freq  <- c(f_s1, f_p1, 1)
    gain  <- c(   0,    1, 1)
    if (f_s1 != 0) { freq <- c(0, freq); gain <- c(0, gain) }
  } else {
    # Low-pass only: 1 from 0 Hz up to h_freq, falling to 0 at (h_freq + slope).
    f_p  <- h_freq / nyquist
    f_s  <- (h_freq + h_trans_bw) / nyquist
    freq <- c(0, f_p, f_s)
    gain <- c(1,   1,   0)
    if (f_s != 1) { freq <- c(freq, 1); gain <- c(gain, 0) }
  }

  # Rounding can push an edge a hair outside [0, 1] (for example 1.0000000000000002),
  # which the kernel design refuses. Clamp to the allowed range.
  freq <- pmax(pmin(freq, 1.0), 0.0)

  # ========== FIR: KERNEL ==========

  # The kernel: N symmetric coefficients with the wanted shape. 'sfreq' is only
  # used to put widths in Hz into the error message when N is too short.
  h <- .firwin_kernel(N, freq, gain, window = fir_window, sfreq = sfreq)

  # "minimum": the same filter shape, rebuilt so that all of its action is on
  # the past side of every sample (causal). No delay correction is made later,
  # because .overlap_add_filter() only does that for words starting with "zero".
  if (phase == "minimum") h <- .minimum_phase(h)

  list(method = "fir", filter_type = filter_type, phase = phase, pad = pad,
       kernel = h, N = N, l_trans_bw = l_trans_bw, h_trans_bw = h_trans_bw,
       window = fir_window)
}

#' Run a designed band-pass filter over one row of numbers
#'
#' The second of the two steps (design, then run). It does no filtering
#' itself: it reads the method from the description made by
#' \code{.make_bandpass_filter()} and hands the row to the engine that does
#' the work, \code{.overlap_add_filter()} for FIR or \code{.iir_filter()} for
#' IIR. Both engines were there before and are unchanged.
#'
#' @param x Numeric vector. One signal.
#' @param spec The list made by \code{.make_bandpass_filter()}.
#' @return Numeric vector, same length as \code{x}.
#' @keywords internal
.apply_bandpass_filter <- function(x, spec) {
  if (spec$method == "iir") {
    # IIR engine: runs the row through the chain of second-order sections.
    # phase "zero": padded at both ends by spec$padlen samples, forward then
    # backward (no delay). phase "forward": one pass, no padding (causal).
    .iir_filter(x, spec$sos, spec$phase, spec$padlen)
  } else {
    # FIR engine: convolves the row with the kernel (overlap-add FFT), pads
    # both ends, and for phase "zero" shifts the result back by the kernel's
    # delay so that nothing is moved in time.
    .overlap_add_filter(x = x, h = spec$kernel, phase = spec$phase, pad = spec$pad)
  }
}

# ----------------------------------------------------------------------------
# .notch_iir() - the IIR half of eeg_notch()
# ----------------------------------------------------------------------------

#' IIR half of eeg_notch()
#'
#' Called by \code{\link{eeg_notch}} when \code{method = "iir"}, after the
#' common arguments have been checked. Every frequency gets its own Butterworth
#' band-stop filter. The filters are applied one after the other, from the
#' lowest frequency to the highest, so the result is the same as calling
#' \code{eeg_notch()} once per frequency.
#'
#' @param eeg_obj         An object of class 'eeg'.
#' @param freqs           Numeric. Notch frequencies in Hz.
#' @param notch_widths    Numeric, same length as \code{freqs}. Width of each
#'   stop band in Hz.
#' @param trans_bandwidth Numeric. Transition width in Hz; the cutoffs lie half
#'   of it further out on each side of the stop band.
#' @param iir_params      List with \code{order} (see
#'   \code{.resolve_iir_params()}).
#' @param phase           "zero" or "forward".
#' @param ch_idx          Integer vector. Rows of the data to filter.
#' @param verbose         Logical. Print parameters and progress.
#' @return An object of class 'eeg'.
#' @keywords internal
.notch_iir <- function(eeg_obj, freqs, notch_widths, trans_bandwidth,
                       iir_params, phase, ch_idx, verbose) {

  sfreq      <- eeg_obj$sampling_rate
  nyquist    <- sfreq / 2
  n_channels <- nrow(eeg_obj$data)

  # Go through the frequencies from low to high, whatever order they were typed
  # in (each width stays with its frequency). The filters could be applied in
  # any order, but the very first and last samples would differ a little.
  ord          <- order(freqs)
  freqs        <- freqs[ord]
  notch_widths <- notch_widths[ord]

  # ========== CUTOFFS ==========
  # The two points where ONE pass of the filter is 3 dB down: the stop band,
  # widened by half the transition width on each side.

  tb_2     <- trans_bandwidth / 2
  low_cut  <- freqs - notch_widths / 2 - tb_2
  high_cut <- freqs + notch_widths / 2 + tb_2

  for (fi in seq_along(freqs)) {
    if (low_cut[fi] <= 0) {
      stop("Notch lower edge (", round(low_cut[fi], 4), " Hz) is <= 0 for freq=",
           freqs[fi], " Hz. Reduce trans_bandwidth or notch_width.", call. = FALSE)
    }
    if (high_cut[fi] >= nyquist) {
      stop("Notch upper edge (", round(high_cut[fi], 4), " Hz) >= Nyquist (",
           nyquist, " Hz) for freq=", freqs[fi], " Hz.", call. = FALSE)
    }
  }

  # ========== DESIGN: SECTIONS AND PADDING LENGTH, ONE SET PER FREQUENCY ==========

  sos_list <- vector("list", length(freqs))
  padlen   <- numeric(length(freqs))
  for (fi in seq_along(freqs)) {
    zpk <- .butter_zpk(iir_params$order,
                       c(low_cut[fi], high_cut[fi]) / nyquist, "stop")
    sos_list[[fi]] <- .zpk_to_sos(zpk$z, zpk$p, zpk$k)
    padlen[fi]     <- if (phase == "zero") .iir_ringing_samples(sos_list[[fi]]) else 0
  }
  pad_text <- if (any(padlen > 0)) {
    paste0("reflect_limited, ",
           paste(format(padlen, scientific = FALSE, trim = TRUE), collapse = " / "),
           " samples at each end")
  } else {
    "none"
  }

  # ========== VERBOSE HEADER ==========

  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("EEG IIR Notch Filtering\n")
    cat(strrep("=", 70), "\n")
    cat("  Notch frequencies:    ", paste(freqs, collapse = ", "), " Hz\n")
    cat("  Notch widths:         ", paste(round(notch_widths, 4), collapse = ", "), " Hz\n")
    cat("  Transition bandwidth: ", trans_bandwidth, " Hz\n")
    cat("  Sampling rate:        ", sfreq, " Hz\n")
    cat("  Family / order:        Butterworth /", iir_params$order, "\n")
    cat("  Phase:                ", phase, "\n")
    cat("  Padding:              ", pad_text, "\n")
    cat("  Channels:             ", length(ch_idx), "/", n_channels, "\n")
    cat(strrep("-", 70), "\n")
    for (fi in seq_along(freqs)) {
      cat(sprintf("  Notch %d/%d: %.1f Hz  [cutoffs %.4f - %.4f Hz, %d sections]\n",
                  fi, length(freqs), freqs[fi], low_cut[fi], high_cut[fi],
                  nrow(sos_list[[fi]])))
    }
  }

  # ========== APPLY THE FILTERS CHANNEL BY CHANNEL ==========

  filtered_data <- eeg_obj$data

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }

  for (idx in seq_along(ch_idx)) {
    i <- ch_idx[idx]
    x <- filtered_data[i, ]
    for (fi in seq_along(freqs)) {
      x <- .iir_filter(x, sos_list[[fi]], phase, padlen[fi])
    }
    filtered_data[i, ] <- x
    if (verbose) setTxtProgressBar(pb, idx)
  }

  if (verbose) {
    close(pb)
    cat("\n", strrep("=", 70), "\n")
  }

  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========

  history_parts <- paste0(freqs, " Hz (\u00B1",     # plus/minus symbol
                          round(notch_widths / 2 + tb_2, 4), " Hz)")
  history_entry <- paste0(
    "IIR notch: ", paste(history_parts, collapse = ", "),
    " | Butterworth order: ", iir_params$order,
    " | phase: ", phase,
    " | pad: ", pad_text,
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
                              list(history_entry)),
    montage               = eeg_obj$montage,
    bads                  = eeg_obj$bads,
    annotations           = eeg_obj$annotations,
    channel_types         = eeg_obj$channel_types
  )
}







# ============================================================================
#            SPECTRUM-FIT LINE-NOISE REMOVAL ("THE FINER RULER")
# ============================================================================
#
# This block holds everything behind eeg_notch(method = "spectrum_fit").
# A longer write-up with the pictures is in notes/finer-ruler-spectrum-fit.pdf.
#
# WHAT IT DOES
#   A notch filter cuts a slot out of the frequencies. Spectrum-fit works the
#   other way round: it LOOKS for steady sine waves (called "lines") in the
#   recording, measures the height and the timing of each one, and subtracts
#   exactly that wave. Nothing else in the recording is touched.
#
# THE STEPS, PER CHANNEL
#   1. Cut the recording into overlapping pieces (10 s each by default; every
#      piece starts half a piece after the one before it).
#   2. In every piece, look at the data through K smooth windows (Slepian
#      windows, also called tapers; K = bandwidth in Hz x piece length in s).
#      For every frequency on a RULER (see below) this gives two numbers:
#        - the height and timing of a steady wave at that frequency, and
#        - an F score: how much stronger the steady wave is than what is left
#          over after it has been taken out. A high score means "this is a
#          real steady wave", a low score means "this is just noise".
#   3. Frequencies whose score is above a line (the threshold, set by p_value)
#      are lines. Each line is subtracted from the piece as a sine wave.
#   4. Glue the cleaned pieces back together with a smooth cross-fade (a Hann
#      window). The fades of overlapping pieces add up to exactly 1, so data
#      that no piece changed come back unchanged.
#
# THE RULER (the "finer ruler", grid_factor)
#   The ruler is the list of frequencies that are tested. The standard ruler
#   has one mark every 1/T Hz, T being the piece length: 0.1 Hz for 10 s
#   pieces. A line that sits BETWEEN two marks gets a low score at both of
#   them and is not found at all:
#
#      marks of the standard ruler:  |        |        |        |        |
#                                   50.2     50.3     50.4     50.5     50.6
#      hum at 50.30 Hz:                       ^         exactly on a mark: found
#      hum at 50.37 Hz:                            ^    between two marks: MISSED
#
#   The finer ruler zero-pads every piece (it adds grid_factor - 1 times the
#   piece length of zeros before the Fourier transform), which gives
#   grid_factor times more marks: 0.0125 Hz apart for grid_factor = 8, so a
#   line is never more than 0.00625 Hz away from a mark. This is the default.
#   grid_factor = 1 is the standard ruler.
#
#   Measured on a made-up recording (5 uV of "brain" noise, a 20 uV hum, 60 s
#   at 512 Hz, 10 s pieces, automatic search). The numbers are the part of the
#   hum that is still there afterwards (rms of the error / rms of the hum):
#
#      hum at                     standard ruler   finer ruler (grid_factor 8)
#      50.00 Hz  (on a mark)           0.8 %              0.8 %
#      50.05 Hz  (between marks)     100   %              0.6 %
#      50.30 Hz  (on a mark)           0.6 %              0.6 %
#      50.37 Hz  (between marks)     100   %              4.0 %
#      50.45 Hz  (between marks)     100   %              0.4 %
#      49.8 -> 50.5 Hz, drifting      78   %              5.6 %   (5 s pieces)
#
#   Worst case of the finer ruler over 41 hum frequencies between two marks:
#   4.2 % left. Ruler factor 4 / 16 / 32 gave 12.5 / 3.3 / 1.9 % and took
#   0.08 / 1.1 / 2.5 s per minute of one channel at 512 Hz (factor 8: 0.18 s,
#   factor 1: 0.02 s; the windows are worked out once, about 0.3 s, for all
#   channels together). Eight is the knee of that curve. For a hum that
#   drifts, shorter pieces (filter_length = "5s") follow it better than long
#   ones: 22 % of the hum left with 10 s pieces, 5.6 % with 5 s, 3 % with 2 s.
#
#   Two things come with having more marks, and the code deals with both:
#
#   (a) Several neighbouring marks pass the test for the SAME line. Subtracting
#       a wave at every one of them would subtract the line again and again
#       (in a test, 1.9 to 2.8 times the hum was left, instead of 0.4 to 4 %).
#       So of every run of neighbouring marks that passed, only the strongest
#       one is used.
#
#   (b) More marks means more chances for noise to pass the test by luck. On
#       hum-free noise the standard ruler gave a false alarm in 3.9 % of the
#       pieces and the finer ruler, with nothing else changed, in 13.6 %. So
#       p_value is shared over ALL marks (n x grid_factor of them instead of
#       n), which keeps it meaning the same at every ruler: 2.4 % of the pieces
#       with a false alarm at p_value = 0.05 and grid_factor = 8. With
#       grid_factor = 1 this is exactly the standard rule. A hum of 1 uV or
#       more is still found in 98 to 100 % of the pieces, on a mark or between
#       two marks; only hums weaker than about 0.8 uV are found less often.
#
# WHEN freqs IS GIVEN (no automatic search)
#   There is no test then. In every piece the mark nearest to each frequency,
#   and every mark strictly inside the notch width, is removed. The finer ruler
#   is NOT used (grid_factor is ignored): every extra mark inside the zone would
#   subtract the same hum once more (the error was 6 to 7 times the hum).
#   p_value is not used either. If the user typed grid_factor or p_value
#   together with freqs, eeg_notch() warns that they have no effect.
#
# LIMITS
#   * One steady wave per mark and piece. A hum whose frequency moves during a
#     piece is smeared and only partly removed: use shorter pieces.
#   * Any wave that is steady for a whole piece is a "line", mains or not (a
#     very regular rhythm too). Look at the frequencies that were reported.
#   * A piece that contains missing values (NA) is left unchanged, and the
#     user is told how many.
#   * It is slower than a notch filter (about 0.2 s per minute of one channel
#     at 512 Hz, plus 0.3 s once for the windows); downsample first if the
#     sampling rate is higher.
#
# HOW THE CODE IS SPLIT (the same two steps as the band-pass recipe)
#   DESIGN  .make_spectrum_fit_plan(): piece length, cross-fade window, the
#           Slepian windows and the score threshold, worked out ONCE for all
#           channels.
#   RUN     .apply_spectrum_fit_plan(): runs the plan on one channel, using
#           .spectrum_fit_piece() for each piece.
#   TOP     .notch_spectrum_fit(): called by eeg_notch(); loops over the
#           channels, prints the header, reports the lines that were found
#           and builds the returned eeg object.
#   BELOW   .sturm_count(), .tridiag_top_eigenvalues() and .dpss_tapers(): the
#           Slepian windows, computed with a small solver in plain R (no
#           extra package is needed).
#
# NOTE ON .compute_dpss() IN fourier.R
#   The package already has a function for Slepian windows, .compute_dpss(),
#   used by eeg_multitaper(). It builds the whole n x n matrix and takes a
#   dense eigen() of it, which is fine for a short epoch (its help says up to
#   about 4096 samples) but takes minutes for a piece of 5120 samples (10 s at
#   512 Hz) and needs 200 MB. It also cuts a symmetric window of n points, while
#   the usual window for spectra is the periodic one (n + 1 symmetric points
#   without the last), which is what the numbers of the reference
#   implementation need. .dpss_tapers() below therefore does not use it: it finds
#   the eigenvalues with a Sturm search and builds the eigenvectors with a
#   recurrence, both in O(n) per step, so a 5120-sample window takes 0.5 s.


#' Count the eigenvalues of a tridiagonal matrix below some values
#'
#' Sturm sequence count for a symmetric tridiagonal matrix: how many of its
#' eigenvalues are smaller than each of the trial values in \code{lam}. All
#' trial values are handled in the same loop.
#'
#' @param d      Numeric. The diagonal of the matrix.
#' @param e2     Numeric, one shorter than \code{d}. The SQUARED off-diagonal.
#' @param lam    Numeric. Trial values.
#' @param pivmin Numeric. A pivot closer to zero than this is replaced by
#'   \code{-pivmin}, so that the loop never divides by zero.
#' @return An integer vector, one count per trial value.
#' @keywords internal
.sturm_count <- function(d, e2, lam, pivmin = 1e-150) {
  M <- length(d)
  q <- d[1L] - lam
  small <- abs(q) < pivmin
  if (any(small)) q[small] <- -pivmin
  cnt <- as.integer(q < 0)
  for (i in 2:M) {
    q <- d[i] - lam - e2[i - 1L] / q
    small <- abs(q) < pivmin
    if (any(small)) q[small] <- -pivmin
    cnt <- cnt + (q < 0)
  }
  cnt
}

#' The largest eigenvalues of a tridiagonal matrix
#'
#' Finds the \code{K} largest eigenvalues of a symmetric tridiagonal matrix by
#' repeated multi-section: every pass tries \code{grid} trial values inside the
#' search interval of each eigenvalue (see \code{.sturm_count()}) and keeps the
#' small piece in which the count changes, so each interval shrinks by a factor
#' of \code{grid + 1} per pass, until it is as narrow as the numbers allow.
#'
#' @param d    Numeric. The diagonal of the matrix.
#' @param e    Numeric, one shorter than \code{d}. The off-diagonal.
#' @param K    Integer. How many eigenvalues to find.
#' @param grid Integer. Trial values per interval and pass.
#' @return A numeric vector of \code{K} eigenvalues, the largest first.
#' @keywords internal
.tridiag_top_eigenvalues <- function(d, e, K, grid = 24L) {
  M  <- length(d)
  e2 <- e^2
  rad <- c(0, abs(e)) + c(abs(e), 0)     # Gershgorin bounds: every eigenvalue lies inside
  lo <- rep(min(d - rad), K)
  hi <- rep(max(d + rad), K)
  target <- M - seq_len(K) + 1L          # the k-th largest is the (M - k + 1)-th smallest
  frac   <- seq_len(grid) / (grid + 1)
  for (pass in 1:100) {
    lam   <- rep(lo, each = grid) + rep(frac, K) * rep(hi - lo, each = grid)
    cnt   <- matrix(.sturm_count(d, e2, lam), nrow = grid)
    below <- cnt < rep(target, each = grid)   # trial value still at or below the wanted eigenvalue
    j     <- colSums(below)
    lam_m <- matrix(lam, nrow = grid)
    new_lo <- ifelse(j == 0L, lo, lam_m[cbind(pmax(j, 1L), seq_len(K))])
    new_hi <- ifelse(j == grid, hi, lam_m[cbind(pmin(j + 1L, grid), seq_len(K))])
    if (all(new_lo == lo & new_hi == hi)) break
    lo <- new_lo
    hi <- new_hi
    if (all(hi - lo <= 2 * .Machine$double.eps * pmax(abs(lo), abs(hi)))) break
  }
  (lo + hi) / 2
}

#' Slepian (multitaper) windows
#'
#' The \code{K} smoothest windows of length \code{N} whose energy is most
#' concentrated in a frequency band of half-width \code{half_nbw / N}
#' (Slepian, 1978; Percival and Walden, 1993). They are the eigenvectors of a
#' tridiagonal matrix that commutes with the (much harder) band-concentration
#' problem. The eigenvalues come from \code{.tridiag_top_eigenvalues()}; each
#' eigenvector is then built by running the recurrence of the matrix from one
#' end to the middle, where the windows are mirror-symmetric (even orders) or
#' anti-symmetric (odd orders), and mirroring the second half. A periodic
#' window is made the usual way: the symmetric window of \code{N + 1} points
#' (unit length over those points) with its last point dropped. Signs follow
#' the usual convention: even orders have a positive sum, odd orders start
#' with a positive lobe.
#'
#' The package has an older function for Slepian windows,
#' \code{.compute_dpss()} (used by \code{eeg_multitaper()}). It takes a dense
#' eigen-decomposition of the full matrix, which is fine for short epochs
#' (up to about 4096 samples) but far too slow for the pieces used here, and
#' it cuts a symmetric window of \code{N} points instead of the periodic one.
#'
#' @param N        Integer. Window length in samples.
#' @param half_nbw Numeric. Time-bandwidth product NW: the bandwidth in Hz times
#'   \code{N} divided by twice the sampling rate.
#' @param K        Integer. Number of windows (orders 0 to \code{K - 1}).
#' @return A \code{K} by \code{N} matrix, one window per row.
#' @references
#' Slepian, D. (1978). Prolate spheroidal wave functions, Fourier analysis, and
#' uncertainty V. Bell System Technical Journal, 57, 1371-1430.
#'
#' Percival, D. B., & Walden, A. T. (1993). Spectral Analysis for Physical
#' Applications. Cambridge University Press.
#' @keywords internal
.dpss_tapers <- function(N, half_nbw, K) {
  M <- N + 1L
  W <- half_nbw / M
  i <- 0:(M - 1L)
  d <- ((M - 1 - 2 * i) / 2)^2 * cos(2 * pi * W)
  e <- i[-1L] * (M - i[-1L]) / 2
  lam  <- .tridiag_top_eigenvalues(d, e, K)
  half <- (M + 1L) %/% 2L                # centre row (M odd) or last row of the left half (M even)
  V <- matrix(0, M, K)
  V[1L, ] <- 1
  V[2L, ] <- -(d[1L] - lam) / e[1L]
  for (r in 2:(half - 1L)) {
    V[r + 1L, ] <- -((d[r] - lam) * V[r, ] + e[r - 1L] * V[r - 1L, ]) / e[r]
  }
  s <- rep(c(1, -1), length.out = K)     # even orders mirror-symmetric, odd orders anti-symmetric
  for (r in seq_len(M - half)) V[M + 1L - r, ] <- s * V[r, ]
  if (M %% 2L == 1L) V[half, s < 0] <- 0
  if (!all(is.finite(V))) {
    stop("The Slepian windows could not be computed (the bandwidth is too ",
         "large for this piece length).", call. = FALSE)
  }
  V <- sweep(V, 2L, sqrt(colSums(V^2)), "/")
  thr <- max(1e-7, 1 / M)
  for (k in seq_len(K)) {
    if (k %% 2L == 1L) {
      if (sum(V[, k]) < 0) V[, k] <- -V[, k]
    } else {
      first <- V[V[, k]^2 > thr, k][1L]
      if (first < 0) V[, k] <- -V[, k]
    }
  }
  t(V[seq_len(N), , drop = FALSE])
}

#' Piece length of the spectrum-fit method, in samples
#'
#' @param filter_length "auto" (10 seconds), a number of samples, or text such
#'   as "5s" or "500ms".
#' @param sfreq         Numeric. Sampling rate in Hz.
#' @param n_total       Integer. Length of the recording in samples; a piece is
#'   never longer than this.
#' @return An integer number of samples.
#' @keywords internal
.piece_length_samples <- function(filter_length, sfreq, n_total) {
  one <- length(filter_length) == 1L
  n <- if (one && is.character(filter_length) && identical(filter_length, "auto")) {
    round(10 * sfreq)
  } else if (one && is.character(filter_length) &&
             grepl("^[0-9]+([.][0-9]+)? *(s|ms)$", filter_length)) {
    value <- as.numeric(sub(" *(s|ms)$", "", filter_length))
    round(value * sfreq * if (grepl("ms$", filter_length)) 0.001 else 1)
  } else if (one && is.numeric(filter_length) && is.finite(filter_length) &&
             filter_length >= 1) {
    round(filter_length)
  } else {
    stop("For method = \"spectrum_fit\", 'filter_length' must be \"auto\" ",
         "(10 seconds), a number of samples, or text like \"5s\" or \"500ms\".",
         call. = FALSE)
  }
  n <- as.integer(min(n, n_total))
  if (n < 32L) {
    stop("The pieces would be only ", n, " samples long; the recording or ",
         "'filter_length' is too short for method = \"spectrum_fit\".",
         call. = FALSE)
  }
  n
}

#' Slepian windows and score threshold for one piece length
#'
#' The windows depend only on the piece length and the bandwidth, so they are
#' worked out once and reused for every piece (and every channel). The number
#' of windows is \code{K = floor(2 * NW)}, with \code{NW = 4} when no bandwidth
#' is given, or \code{NW = mt_bandwidth * len / (2 * sfreq)}. The threshold is
#' the F value that a score must exceed: the upper \code{p_value / (len * grid)}
#' point of the F distribution with 2 and \code{2K - 2} degrees of freedom,
#' i.e. \code{p_value} shared over all \code{len * grid} marks that are tested.
#'
#' @param len          Integer. Piece length in samples.
#' @param sfreq        Numeric. Sampling rate in Hz.
#' @param mt_bandwidth Numeric or NULL. Bandwidth in Hz.
#' @param p_value      Numeric. See \code{eeg_notch()}.
#' @param grid         Integer. The finer-ruler factor that is in use (1 when
#'   the frequencies are given).
#' @return A list with the \code{tapers} (K by len), \code{K}, \code{half_nbw}
#'   and the \code{threshold}.
#' @keywords internal
.spectrum_fit_design <- function(len, sfreq, mt_bandwidth, p_value, grid) {
  half_nbw <- if (is.null(mt_bandwidth)) 4 else mt_bandwidth * len / (2 * sfreq)
  K <- as.integer(floor(2 * half_nbw))
  if (K < 3L) {
    stop("'mt_bandwidth' is too small: with pieces of ", round(len / sfreq, 3),
         " s it gives only ", K, " window(s) and at least 3 are needed. Use at ",
         "least ", round(3 / (len / sfreq), 3), " Hz (or longer pieces).",
         call. = FALSE)
  }
  if (half_nbw >= len / 4) {
    stop("'mt_bandwidth' is too large for pieces of ", round(len / sfreq, 3),
         " s. Use a smaller bandwidth or longer pieces.", call. = FALSE)
  }
  list(tapers    = .dpss_tapers(len, half_nbw, K),
       K         = K,
       half_nbw  = half_nbw,
       threshold = stats::qf(p_value / (len * grid), 2, 2 * K - 2,
                             lower.tail = FALSE))
}

#' Cross-fade weights of one piece
#'
#' Every piece is multiplied by a smooth fade-in / fade-out before the pieces
#' are added up. Overlapping fades add up to 1. The first piece has no
#' neighbour on its left and the last piece (which may be longer than the
#' others, because it takes the samples that are left over) has none on its
#' right, so their weights are completed to keep the sum at 1 everywhere.
#'
#' @param win      Numeric. The fade of an ordinary piece (length \code{n}).
#' @param n        Integer. Ordinary piece length.
#' @param step     Integer. Distance between the starts of two pieces.
#' @param len      Integer. Length of THIS piece.
#' @param is_first Logical. Is this the first piece?
#' @param is_last  Logical. Is this the last piece?
#' @return A numeric vector of length \code{len}.
#' @keywords internal
.piece_weights <- function(win, n, step, len, is_first, is_last) {
  w <- win
  if (is_last) {
    w <- c(win, numeric(len - n))
    if (step < len) {
      for (off in seq(step, len - 1L, by = step)) {
        w[(off + 1L):len] <- w[(off + 1L):len] + win[seq_len(len - off)]
      }
    }
  }
  if (is_first && n - step >= 1L) {
    for (off in seq(n - step, 1L, by = -step)) {
      w[seq_len(off)] <- w[seq_len(off)] + win[(n - off + 1L):n]
    }
  }
  w
}

#' Plan of the spectrum-fit method (the design step)
#'
#' Everything that does not depend on the data of a channel: where the pieces
#' start and stop, the cross-fade, and the Slepian windows and threshold of each
#' piece length. Worked out once and used for every channel by
#' \code{.apply_spectrum_fit_plan()}.
#'
#' @param n_total       Integer. Length of the recording in samples.
#' @param sfreq         Numeric. Sampling rate in Hz.
#' @param filter_length See \code{.piece_length_samples()}.
#' @param mt_bandwidth  Numeric or NULL. Bandwidth in Hz.
#' @param p_value       Numeric. See \code{eeg_notch()}.
#' @param grid_factor   Integer. The finer-ruler factor asked for.
#' @param searching     Logical. TRUE for the automatic search (\code{freqs} =
#'   NULL). Only then is the finer ruler used.
#' @return A list: \code{n}, \code{step}, \code{win}, \code{starts},
#'   \code{stops}, \code{designs} (one per piece length), \code{grid},
#'   \code{searching}.
#' @keywords internal
.make_spectrum_fit_plan <- function(n_total, sfreq, filter_length, mt_bandwidth,
                                    p_value, grid_factor, searching) {
  n    <- .piece_length_samples(filter_length, sfreq, n_total)
  grid <- if (searching) as.integer(grid_factor) else 1L

  # smooth cross-fade (Hann), scaled so that overlapping copies add up to 1
  step <- n - (n + 1L) %/% 2L
  win <- if (n %% 2L == 0L) {
    0.5 - 0.5 * cos(2 * pi * (0:(n - 1L)) / n)
  } else {
    0.5 - 0.5 * cos(2 * pi * (0:(n - 1L)) / (n - 1L))
  }
  nseg    <- n %/% step
  binsums <- Reduce(`+`, lapply(0:(nseg - 1L),
                                function(ii) win[(ii * step + 1L):((ii + 1L) * step)]))
  rem <- n %% step
  if (rem != 0L) binsums[seq_len(rem)] <- binsums[seq_len(rem)] + win[(n - rem + 1L):n]
  win <- win / stats::median(binsums)

  starts <- seq(0L, n_total - n, by = step)
  stops  <- starts + n
  stops[length(stops)] <- n_total          # the last piece takes the leftover samples
  lens   <- stops - starts

  # a very long padded transform needs a lot of memory: say so early
  K_guess <- if (is.null(mt_bandwidth)) 8 else max(3, floor(mt_bandwidth * max(lens) / sfreq))
  if (as.numeric(max(lens)) * grid * K_guess > 6e7) {
    stop("The pieces are too long for this 'grid_factor' (it would need more ",
         "than about 1 GB of memory). Use a smaller 'grid_factor', a shorter ",
         "'filter_length', or downsample first.", call. = FALSE)
  }

  designs <- list()
  for (len in unique(lens)) {
    designs[[as.character(len)]] <- .spectrum_fit_design(len, sfreq, mt_bandwidth,
                                                         p_value, grid)
  }
  list(n = n, step = step, win = win, starts = starts, stops = stops,
       designs = designs, grid = grid, searching = searching)
}

#' Spectrum-fit on one piece of one channel
#'
#' Looks, for every mark of the ruler, for a steady wave; decides which marks
#' are lines; and subtracts those waves. See the overview above.
#'
#' @param x            Numeric. The piece.
#' @param sfreq        Numeric. Sampling rate in Hz.
#' @param freqs        Numeric or NULL. Given frequencies (zone removal) or
#'   NULL (automatic search).
#' @param notch_widths Numeric, same length as \code{freqs}, or NULL.
#' @param design       List from \code{.spectrum_fit_design()}.
#' @param grid         Integer. Ruler factor (1 = standard ruler).
#' @return A list: \code{x} (the cleaned piece), \code{found} (the frequencies,
#'   in Hz, of the lines that were subtracted), \code{freq} (the marks of the
#'   ruler, in Hz) and \code{fstat} (the F score of every mark; NULL when
#'   frequencies are given, because nothing is tested then).
#' @references
#' Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
#' Proceedings of the IEEE, 70, 1055-1096.
#'
#' Mitra, P., & Bokil, H. (2008). Observed Brain Dynamics. Oxford University
#' Press.
#' @keywords internal
.spectrum_fit_piece <- function(x, sfreq, freqs, notch_widths, design, grid) {
  n      <- length(x)
  tapers <- design$tapers
  K      <- nrow(tapers)
  odd    <- seq(1L, K, by = 2L)          # orders 0, 2, 4, ...: they look like a hill
  even   <- seq(2L, K, by = 2L)          # orders 1, 3, 5, ...: they have a node in the middle
  H0     <- rowSums(tapers[odd, , drop = FALSE])
  H0_sq  <- sum(H0^2)
  rads   <- 2 * pi * (0:(n - 1L)) / sfreq

  # tapered spectra on the ruler (zero-padded when grid > 1: finer marks)
  nfft <- n * as.integer(grid)
  nf   <- nfft %/% 2L + 1L
  xt   <- t(tapers) * (x - mean(x))
  if (nfft > n) xt <- rbind(xt, matrix(0, nfft - n, K))
  X <- stats::mvfft(xt)[seq_len(nf), , drop = FALSE]
  X[1L, ] <- X[1L, ] / sqrt(2)
  if (nfft %% 2L == 0L) X[nf, ] <- X[nf, ] / sqrt(2)
  fr <- (0:(nf - 1L)) * sfreq / nfft
  Xo <- X[, odd, drop = FALSE]
  Xe <- X[, even, drop = FALSE]

  # complex strength (height and timing) of a steady wave at every mark
  A <- as.vector(Xo %*% H0) / H0_sq

  fstat <- NULL
  if (is.null(freqs)) {
    # automatic search: the F score of every mark against the threshold
    x_hat <- outer(A, H0)
    num   <- (K - 1) * Mod(A)^2 * H0_sq
    den   <- rowSums(Mod(Xo - x_hat)^2) + rowSums(Mod(Xe)^2)
    den[den == 0] <- Inf
    fstat <- num / den
    idx   <- which(fstat > design$threshold)
    if (grid > 1L && length(idx) > 1L) {
      # of every run of neighbouring marks that passed, keep the strongest one
      run <- cumsum(c(1L, diff(idx) > 1L))
      idx <- vapply(split(idx, run), function(g) g[which.max(fstat[g])], 1L)
    }
  } else {
    # given frequencies: the mark nearest to each, plus every mark strictly
    # inside the notch width
    idx1 <- unique(vapply(freqs, function(lf) which.min(abs(fr - lf)), 1L))
    zone <- Reduce(`|`, Map(function(lf, nw) fr > lf - nw / 2 & fr < lf + nw / 2,
                            freqs, notch_widths))
    idx  <- sort(unique(c(idx1, which(zone))))
  }

  # subtract one sine wave per chosen mark (at most 200 at a time, for memory)
  fit <- numeric(n)
  if (length(idx)) {
    cc <- 2 * A[idx]
    for (g in split(seq_along(idx), ceiling(seq_along(idx) / 200))) {
      ph  <- outer(rads, fr[idx[g]]) + matrix(Arg(cc[g]), n, length(g), byrow = TRUE)
      fit <- fit + as.vector(cos(ph) %*% Mod(cc[g]))
    }
  }
  list(x = x - fit, found = fr[idx], freq = fr, fstat = fstat)
}

#' Run a spectrum-fit plan on one channel (the run step)
#'
#' Cuts the channel into the pieces of the plan, cleans each piece with
#' \code{.spectrum_fit_piece()} and adds the cleaned pieces back together with
#' the cross-fade weights. A piece that contains missing values cannot be
#' analysed and is left as it is.
#'
#' @param x            Numeric. One channel.
#' @param sfreq        Numeric. Sampling rate in Hz.
#' @param plan         List from \code{.make_spectrum_fit_plan()}.
#' @param freqs        Numeric or NULL. See \code{.spectrum_fit_piece()}.
#' @param notch_widths Numeric or NULL. See \code{.spectrum_fit_piece()}.
#' @return A list: \code{x} (the cleaned channel), \code{found} (a list with the
#'   frequencies found in each piece) and \code{skipped} (how many pieces were
#'   left unchanged because of missing values).
#' @keywords internal
.apply_spectrum_fit_plan <- function(x, sfreq, plan, freqs, notch_widths) {
  n_pieces <- length(plan$starts)
  out     <- numeric(length(x))
  found   <- vector("list", n_pieces)
  skipped <- 0L
  for (ii in seq_len(n_pieces)) {
    s   <- plan$starts[ii]
    e   <- plan$stops[ii]
    len <- e - s
    xs  <- x[(s + 1L):e]
    w   <- .piece_weights(plan$win, plan$n, plan$step, len,
                          is_first = ii == 1L, is_last = ii == n_pieces)
    if (anyNA(xs)) {
      res <- list(x = xs, found = numeric(0))
      skipped <- skipped + 1L
    } else {
      res <- .spectrum_fit_piece(xs, sfreq, freqs, notch_widths,
                                 plan$designs[[as.character(len)]], plan$grid)
    }
    out[(s + 1L):e] <- out[(s + 1L):e] + res$x * w
    found[[ii]] <- res$found
  }
  list(x = out, found = found, skipped = skipped)
}

#' Spectrum-fit half of eeg_notch()
#'
#' Called by \code{\link{eeg_notch}} when \code{method = "spectrum_fit"}, after
#' the common arguments have been checked. Makes the plan once, runs it on every
#' selected channel, tells the user which lines were found, and returns the new
#' \code{eeg} object. The method is explained in the comment block above this
#' function and in \code{notes/finer-ruler-spectrum-fit.pdf}.
#'
#' @param eeg_obj       An object of class 'eeg'.
#' @param freqs         Numeric or NULL. Frequencies to remove (zone removal), or
#'   NULL to search for the lines.
#' @param notch_widths  Numeric (same length as \code{freqs}) or NULL.
#' @param filter_length See \code{.piece_length_samples()}.
#' @param mt_bandwidth  Numeric or NULL. Bandwidth in Hz.
#' @param p_value       Numeric. False-alarm level of the automatic search.
#' @param grid_factor   Integer. Finer-ruler factor of the automatic search.
#' @param ch_idx        Integer vector. Rows of the data to clean.
#' @param verbose       Logical. Print parameters and progress.
#' @return An object of class 'eeg'.
#' @keywords internal
.notch_spectrum_fit <- function(eeg_obj, freqs, notch_widths, filter_length,
                                mt_bandwidth, p_value, grid_factor, ch_idx,
                                verbose) {

  sfreq      <- eeg_obj$sampling_rate
  n_channels <- nrow(eeg_obj$data)
  searching  <- is.null(freqs)

  # ========== DESIGN (ONCE FOR ALL CHANNELS) ==========

  plan    <- .make_spectrum_fit_plan(ncol(eeg_obj$data), sfreq, filter_length,
                                     mt_bandwidth, p_value, grid_factor, searching)
  d_main  <- plan$designs[[as.character(plan$n)]]
  piece_s <- plan$n / sfreq
  band_hz <- 2 * d_main$half_nbw / piece_s
  mark_hz <- sfreq / (plan$n * plan$grid)

  # ========== VERBOSE HEADER ==========

  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("EEG Spectrum-Fit Notch Filtering\n")
    cat(strrep("=", 70), "\n")
    if (searching) {
      cat("  Lines:                 searched automatically (F test)\n")
      cat("  p_value:              ", p_value, "\n")
    } else {
      cat("  Frequencies:          ", paste(freqs, collapse = ", "), " Hz\n")
      cat("  Notch widths:         ", paste(round(notch_widths, 4), collapse = ", "), " Hz\n")
    }
    cat("  Sampling rate:        ", sfreq, " Hz\n")
    cat("  Piece length:         ", round(piece_s, 3), " s,", length(plan$starts),
        " overlapping pieces\n")
    cat("  Bandwidth / windows:  ", round(band_hz, 3), " Hz /", d_main$K, "\n")
    cat("  Ruler:                 a mark every", signif(mark_hz, 4), "Hz",
        if (plan$grid > 1L) paste0("(", plan$grid, " times finer than the standard ruler)")
        else "(the standard ruler)", "\n")
    cat("  Channels:             ", length(ch_idx), "/", n_channels, "\n")
    cat(strrep("-", 70), "\n")
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }

  # ========== RUN, CHANNEL BY CHANNEL ==========

  filtered_data <- eeg_obj$data
  all_found     <- numeric(0)
  skipped       <- 0L
  for (idx in seq_along(ch_idx)) {
    i   <- ch_idx[idx]
    res <- .apply_spectrum_fit_plan(filtered_data[i, ], sfreq, plan, freqs, notch_widths)
    filtered_data[i, ] <- res$x
    all_found <- c(all_found, unlist(res$found))
    skipped   <- skipped + res$skipped
    if (verbose) setTxtProgressBar(pb, idx)
  }

  # ========== REPORT WHAT WAS FOUND ==========
  # Lines are counted in 0.1 Hz groups: (channel, piece) pairs in which a line
  # was subtracted.

  counts <- sort(table(sprintf("%.1f", round(all_found, 1))), decreasing = TRUE)
  found_text <- if (length(counts) == 0L) {
    "none"
  } else {
    shown <- utils::head(counts, 6L)
    paste0(paste0(names(shown), " Hz (x", shown, ")", collapse = ", "),
           if (length(counts) > 6L) paste0(", ... ", length(counts) - 6L, " more") else "")
  }

  if (verbose) {
    close(pb)
    cat("\n")
    cat(if (searching) "  Lines found" else "  Waves removed",
        " (frequency, number of channel-pieces): ", found_text, "\n", sep = "")
    cat(strrep("=", 70), "\n")
  }
  if (skipped > 0L) {
    warning(skipped, " piece(s) contained missing values and were left ",
            "unchanged (no line removal there).", call. = FALSE)
  }

  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========

  what <- if (searching) {
    paste0("lines searched (p = ", p_value, ")")
  } else {
    paste0(paste0(freqs, " Hz (\u00B1", round(notch_widths / 2, 4), " Hz)",
                  collapse = ", "))
  }
  history_entry <- paste0(
    "Spectrum-fit notch: ", what,
    " | pieces: ", round(piece_s, 3), " s",
    " | bandwidth: ", round(band_hz, 3), " Hz",
    " | ruler: ", if (plan$grid > 1L) paste0("x", plan$grid, " finer") else "standard",
    " | found: ", found_text,
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
                              list(history_entry)),
    montage               = eeg_obj$montage,
    bads                  = eeg_obj$bads,
    annotations           = eeg_obj$annotations,
    channel_types         = eeg_obj$channel_types
  )
}



# ============================================================================
#                        MAIN EXPORTED FUNCTIONS
# ============================================================================

# ----------------------------------------------------------------------------
# eeg_bandpass()
# ----------------------------------------------------------------------------

#' EEG Bandpass / Highpass / Lowpass Filter (FIR or IIR)
#'
#' Applies a zero-phase FIR filter (the default) or an IIR Butterworth filter
#' to EEG data. With the defaults, (method='fir', phase='zero',
#' fir_window='hamming', fir_design='firwin', pad='reflect_limited'), it is a
#' zero-phase FIR filter.
#'
#' Supports highpass-only (h_freq=NULL), lowpass-only (l_freq=NULL),
#' or bandpass (both specified).
#'
#' \strong{IIR filters} (\code{method = "iir"}) are Butterworth filters, built
#' from their poles and zeros and run as a chain of second-order sections, so
#' they stay accurate even for slow filters (a low cutoff compared with the
#' sampling rate). The cutoff frequencies are the points where ONE pass of the
#' filter is 3 dB down. With \code{phase = "zero"} the filter is run forward and
#' then backward: there is no delay, the response is 6 dB down at the cutoff,
#' and, because the backward pass looks at the signal before an event as well,
#' an event can leave a trace before it happened. With \code{phase = "forward"}
#' the filter is run once: nothing can appear before an event (the filter is
#' causal), but different frequencies are delayed by different amounts. Only
#' the Butterworth family is implemented so far.
#'
#' \strong{Causal FIR filter} (\code{phase = "minimum"}): the usual symmetric
#' filter is rebuilt as a minimum-phase filter of the same length and with the
#' same magnitude response (cepstrum method), so that all of its action lies on
#' the past side of every sample and no delay correction is made. Nothing can
#' appear before an event, but different frequencies are delayed by different
#' amounts. For a filter that has no gain at 0 Hz (a high-pass or band-pass)
#' the conversion is very sensitive to rounding: two correct implementations
#' agree to about 1e-4 of the signal size, mostly below the lowest cutoff. The
#' conversion uses a long FFT, so for a very long filter, such as a 0.1 Hz
#' high-pass at 2048 Hz (67585 taps), it takes about 40 seconds and 1 GB of
#' memory; downsampling first makes it much cheaper.
#'
#' @references
#' Butterworth, S. (1930). On the theory of filter amplifiers. Wireless
#' Engineer, 7, 536-541.
#'
#' Oppenheim, A. V., & Schafer, R. W. (2010). Discrete-Time Signal Processing
#' (3rd ed.). Pearson.
#'
#' @param eeg_obj          An object of class 'eeg'.
#' @param l_freq           Numeric or NULL. Highpass cutoff in Hz.
#'                         NULL = no highpass (lowpass only).
#' @param h_freq           Numeric or NULL. Lowpass cutoff in Hz.
#'                         NULL = no lowpass (highpass only).
#' @param method           Character. "fir" (default) or "iir".
#' @param iir_params       NULL or a named list, used only with
#'                         \code{method = "iir"}: \code{order} (default 4; a
#'                         band-pass has twice as many poles) and \code{ftype}
#'                         (only "butter" so far).
#' @param l_trans_bandwidth Numeric or "auto". Highpass transition bandwidth
#'                         in Hz. "auto" uses MNE's rule:
#'                         min(max(0.25*l_freq, 2.0), l_freq).
#'                         Default: "auto". Used only with method = "fir".
#' @param h_trans_bandwidth Numeric or "auto". Lowpass transition bandwidth
#'                         in Hz. "auto" uses MNE's rule:
#'                         min(max(0.25*h_freq, 2.0), nyquist-h_freq).
#'                         Default: "auto". Used only with method = "fir".
#' @param filter_length    Integer or "auto". FIR kernel length in samples.
#'                         "auto" uses MNE's formula:
#'                         ceil(3.3 / min_trans_bw * sfreq), forced odd.
#'                         Default: "auto". Used only with method = "fir".
#' @param fir_window       Character. Window function for kernel design.
#'                         Default: "hamming". Used only with method = "fir".
#' @param phase            Character. The allowed words depend on
#'                         \code{method}. "fir": "zero" (default), one pass
#'                         with the group delay corrected; or "minimum", one
#'                         pass with a filter that only looks at the past
#'                         (causal). "iir": "zero" (default), forward then
#'                         backward, no delay; or "forward", one pass, causal.
#'                         Any other word is an error.
#' @param pad              Character. Padding at both ends of the recording.
#'                         Only "reflect_limited" (default) is implemented.
#' @param channels         Character or integer vector, or NULL. Channel names
#'                         or indices to filter. NULL = all channels.
#' @param verbose          Logical. Print filter parameters and progress.
#'                         Default: TRUE.
#'
#' @return An object of class 'eeg' with filtered data and updated
#'   preprocessing_history.
#'
#' @seealso \code{\link{filter_data}} to run the same filter on numbers that are
#'   not in an \code{eeg} object.
#'
#' @export
eeg_bandpass <- function(eeg_obj,
                         l_freq              = NULL,
                         h_freq              = NULL,
                         method              = c("fir", "iir"),
                         iir_params          = NULL,
                         l_trans_bandwidth   = "auto",
                         h_trans_bandwidth   = "auto",
                         filter_length       = "auto",
                         fir_window          = "hamming",
                         phase               = "zero",
                         pad                 = "reflect_limited",
                         channels            = NULL,
                         verbose             = TRUE) {

  # HOW THIS FUNCTION IS BUILT
  # eeg_bandpass() looks after everything that belongs to an eeg object: the
  # channel choice, the printed header and progress bar, the history entry and
  # the new eeg object. The filtering itself (FIR and IIR) is done by the
  # band-pass recipe that is shared with filter_data(); see the overview above
  # .make_bandpass_filter(). In order:
  #   1. check the input          .check_band_edges(), .check_filter_settings()
  #   2. choose the channels      (here)
  #   3. design the filter        .make_bandpass_filter()
  #   4. print the header         (here)
  #   5. filter each channel      .apply_bandpass_filter()
  #   6. write the history entry  (here)
  #   7. build the new object     (here)

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }

  # The sampling rate comes from the object; filter_data() gets it as an argument.
  sfreq   <- eeg_obj$sampling_rate
  nyquist <- sfreq / 2.0

  # Are the cutoffs allowed? (shared with filter_data())
  .check_band_edges(l_freq, h_freq, sfreq)

  # ========== CHECK METHOD, PHASE, PAD AND IIR SETTINGS ==========

  method <- match.arg(method)       # "fir" unless the user asks for "iir"

  # missing() is TRUE when an argument was left at its default, so the check
  # can warn about FIR-only settings that the user set on purpose. Only this
  # function can ask that about its own arguments, so it does it here and passes
  # the answer on.
  chk <- .check_filter_settings(
    method, phase, pad, iir_params,
    fir_only_given = c(l_trans_bandwidth = !missing(l_trans_bandwidth),
                       h_trans_bandwidth = !missing(h_trans_bandwidth),
                       filter_length     = !missing(filter_length),
                       fir_window        = !missing(fir_window)))
  # From here on the checked values are used.
  phase      <- chk$phase
  pad        <- chk$pad
  iir_params <- chk$iir_params

  # ========== RESOLVE CHANNEL INDICES ==========
  # Which rows of the data table to filter: all of them (NULL), the ones with
  # these names, or the ones at these positions. This is eeg-object work, so it
  # stays here; filter_data() has no channels.

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

  # ========== DESIGN THE FILTER ==========
  # Step 1 of the shared recipe, done once for all channels. 'spec' describes
  # the finished filter (FIR kernel or IIR sections, lengths, padding, ...); it
  # is used below for the header and the history entry and, in the loop, to run
  # the filter. For method = "iir" the same call designs the IIR filter.

  spec <- .make_bandpass_filter(sfreq, l_freq, h_freq, method, iir_params,
                                l_trans_bandwidth, h_trans_bandwidth,
                                filter_length, fir_window, phase, pad)

  # ========== VERBOSE HEADER ==========
  # What the filter is, printed before the work starts. The first lines are the
  # same for FIR and IIR; the details differ (a kernel length and window for
  # FIR, the order and number of sections for IIR).

  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    if (method == "fir") {
      cat("EEG FIR Filtering\n")
    } else {
      cat("EEG IIR Filtering\n")
    }
    cat(strrep("=", 70), "\n")
    cat("  Filter type:          ", spec$filter_type, "\n")
    cat("  Sampling rate:        ", sfreq, " Hz\n")
    cat("  Nyquist:              ", nyquist, " Hz\n")
    if (method == "fir") {
      if (!is.null(spec$l_trans_bw))
        cat("  HP transition BW:     ", round(spec$l_trans_bw, 4), " Hz\n")
      if (!is.null(spec$h_trans_bw))
        cat("  LP transition BW:     ", round(spec$h_trans_bw, 4), " Hz\n")
      cat("  Filter length:        ", spec$N, " samples (", round(spec$N/sfreq, 3), " s)\n")
      cat("  Window:               ", spec$window, "\n")
      cat("  Phase:                ", phase, "\n")
      cat("  Padding:              ", pad, "\n")
    } else {
      cat("  Family / order:        Butterworth /", spec$order, "\n")
      cat("  Sections:             ", nrow(spec$sos), "\n")
      cat("  Phase:                ", phase, "\n")
      cat("  Padding:              ", spec$pad_text, "\n")
    }
    cat("  Channels:             ", length(ch_idx), "/", n_channels, "\n")
    cat(strrep("-", 70), "\n")
  }

  # ========== APPLY FILTER CHANNEL BY CHANNEL ==========
  # Step 2 of the shared recipe: every chosen row of the data table goes through
  # .apply_bandpass_filter() on its own. The rows that were not chosen are left
  # exactly as they are. filtered_data starts as a copy of the data table.

  filtered_data <- eeg_obj$data

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }

  for (idx in seq_along(ch_idx)) {
    i <- ch_idx[idx]
    filtered_data[i, ] <- .apply_bandpass_filter(filtered_data[i, ], spec)
    if (verbose) setTxtProgressBar(pb, idx)
  }

  if (verbose) {
    close(pb)
    cat("\n", strrep("=", 70), "\n")
  }

  # ========== BUILD PREPROCESSING HISTORY ENTRY ==========
  # One line that records what was done to the data, taken from the description
  # of the filter. It is added to the object's preprocessing history below.

  history_entry <- if (method == "fir") {
    paste0(
      "FIR bandpass: ", spec$filter_type,
      " | window: ", spec$window,
      " | length: ", spec$N, " samples",
      " | phase: ", phase,
      " | pad: ", pad,
      " | channels: ", length(ch_idx), "/", n_channels
    )
  } else {
    paste0(
      "IIR bandpass: ", spec$filter_type,
      " | Butterworth order: ", spec$order,
      " | phase: ", phase,
      " | pad: ", spec$pad_text,
      " | channels: ", length(ch_idx), "/", n_channels
    )
  }

  # ========== RETURN UPDATED EEG OBJECT ==========
  # A new eeg object: the filtered data table and the history entry added;
  # everything else (channels, times, events, bad channels, montage,
  # annotations, channel types, ...) is carried over as it was.

  new_eeg(
    data                  = filtered_data,
    channels              = eeg_obj$channels,
    sampling_rate         = eeg_obj$sampling_rate,
    times                 = eeg_obj$times,
    events                = eeg_obj$events,
    metadata              = eeg_obj$metadata,
    reference             = eeg_obj$reference,
    preprocessing_history = c(eeg_obj$preprocessing_history,
                              list(history_entry)),
    montage               = eeg_obj$montage,
    bads                  = eeg_obj$bads,
    annotations           = eeg_obj$annotations,
    channel_types         = eeg_obj$channel_types
  )
}


# ----------------------------------------------------------------------------
# filter_data()
# ----------------------------------------------------------------------------

#' Filter Plain Numbers (Band-pass, High-pass or Low-pass)
#'
#' Filters numbers that are not inside an \code{eeg} object and gives them back
#' filtered: one row of numbers (a vector), or several (a matrix with one row
#' per signal, like the data table of an \code{eeg} object). It runs exactly
#' the same filter as \code{\link{eeg_bandpass}} (FIR or IIR, zero-phase or
#' causal), with the same arguments and the same numbers, but it needs no
#' \code{eeg} object: there are no channels, events or history to look after,
#' and nothing is printed.
#'
#' Use \code{\link{eeg_bandpass}} for the channels of a recording. This function
#' is for other rows of numbers: a trace you computed from channels (for
#' example the strength of a signal over time), the time course of an ICA
#' component, a signal recorded by another device, or numbers made up to test a
#' method. \code{eeg_bandpass()} and \code{filter_data()} build and run the
#' filter with the same inner functions, so the two always agree.
#'
#' A missing value (NA or NaN) spoils the row: the filter cannot run past it,
#' so most or all of the row comes back missing. Find and deal with gaps first
#' (see \code{\link{annotate_nan}}).
#'
#' @param x Numeric vector (one signal) or numeric matrix with one row per
#'   signal and one column per time point.
#' @param sfreq Numeric. Sampling rate of the numbers, in Hz.
#' @inheritParams eeg_bandpass
#'
#' @return Numbers of the same shape as \code{x} (a vector, or a matrix with the
#'   same dimension names), filtered.
#'
#' @seealso \code{\link{eeg_bandpass}} for the channels of an \code{eeg} object.
#'
#' @examples
#' # 10 seconds of a 10 Hz wave with 50 Hz hum, sampled at 256 Hz
#' fs <- 256
#' t  <- (0:2559) / fs
#' x  <- sin(2 * pi * 10 * t) + sin(2 * pi * 50 * t)
#'
#' # band-pass 1-40 Hz: the wave stays, the hum goes
#' y <- filter_data(x, sfreq = fs, l_freq = 1, h_freq = 40)
#'
#' # a causal IIR low-pass instead (nothing appears before it happened)
#' y2 <- filter_data(x, sfreq = fs, h_freq = 40, method = "iir", phase = "forward")
#'
#' # several signals at once: one row each
#' m  <- rbind(x, 0.5 * x)
#' ym <- filter_data(m, sfreq = fs, l_freq = 1, h_freq = 40)
#'
#' @export
filter_data <- function(x,
                        sfreq,
                        l_freq              = NULL,
                        h_freq              = NULL,
                        method              = c("fir", "iir"),
                        iir_params          = NULL,
                        l_trans_bandwidth   = "auto",
                        h_trans_bandwidth   = "auto",
                        filter_length       = "auto",
                        fir_window          = "hamming",
                        phase               = "zero",
                        pad                 = "reflect_limited") {

  # HOW THIS FUNCTION WORKS
  # filter_data() is the front door of the shared band-pass recipe (see the
  # overview above .make_bandpass_filter()) for numbers that are not in an eeg
  # object. It does the same steps as eeg_bandpass(), without the eeg-object
  # work (no channels, no header, no history, no new object):
  #   1. check x and sfreq                   (here)
  #   2. check the cutoffs and the settings  .check_band_edges(),
  #                                          .check_filter_settings()
  #   3. design the filter, once             .make_bandpass_filter()
  #   4. run it over every row of x          .apply_bandpass_filter()
  # The other places in the package that filter a loose row of numbers
  # (annotate_muscle(), find_bads_eog(), find_bads_ecg()) call this function.

  # ========== INPUT VALIDATION ==========

  # x: a vector, or a matrix with one row per signal; numbers only, and not empty
  # (a 3-dimensional array has no meaning here).
  if (!is.numeric(x) || length(x) == 0L || length(dim(x)) > 2L) {
    stop("'x' must be a non-empty numeric vector, or a numeric matrix with one ",
         "row per signal.", call. = FALSE)
  }
  # eeg_bandpass() takes the sampling rate from its eeg object; here it has to
  # be given, and it must be a single positive number.
  if (!is.numeric(sfreq) || length(sfreq) != 1L || is.na(sfreq) || sfreq <= 0) {
    stop("'sfreq' must be a single positive number (the sampling rate in Hz).",
         call. = FALSE)
  }

  # Are the cutoffs allowed? (shared with eeg_bandpass())
  .check_band_edges(l_freq, h_freq, sfreq)

  # ========== CHECK METHOD, PHASE, PAD AND IIR SETTINGS ==========

  method <- match.arg(method)       # "fir" unless the user asks for "iir"

  # Same check as in eeg_bandpass(); the missing() answers have to come from
  # here because these are this function's own arguments.
  chk <- .check_filter_settings(
    method, phase, pad, iir_params,
    fir_only_given = c(l_trans_bandwidth = !missing(l_trans_bandwidth),
                       h_trans_bandwidth = !missing(h_trans_bandwidth),
                       filter_length     = !missing(filter_length),
                       fir_window        = !missing(fir_window)))

  # ========== DESIGN THE FILTER ONCE, RUN IT OVER EVERY ROW ==========

  # Step 3: design. One design is enough for any number of rows; it uses the
  # checked values in 'chk'. 'spec' describes the finished filter.
  spec <- .make_bandpass_filter(sfreq, l_freq, h_freq, method, chk$iir_params,
                                l_trans_bandwidth, h_trans_bandwidth,
                                filter_length, fir_window, chk$phase, chk$pad)

  # Step 4: run. 'out' starts as a copy of x, so a vector stays a vector, a
  # matrix keeps its shape, and names / dimension names are kept. Whole numbers
  # (integers) become ordinary numbers when the filtered values are written in.
  out <- x
  if (length(dim(x)) == 2L) {
    # a matrix: every row is a signal of its own and is filtered on its own
    for (i in seq_len(nrow(x))) {
      out[i, ] <- .apply_bandpass_filter(out[i, ], spec)
    }
  } else {
    # a vector: a single signal
    out[] <- .apply_bandpass_filter(as.vector(x), spec)
  }

  out
}



# ----------------------------------------------------------------------------
# eeg_notch()
# ----------------------------------------------------------------------------
#
#' EEG Notch (Band-Stop) Filter (FIR, IIR or Spectrum-Fit)
#'
#' Removes narrow frequency bands from EEG data to eliminate electrical line
#' noise, with a zero-phase FIR filter (the default), an IIR Butterworth
#' band-stop filter, or by finding the steady sine waves of the line noise and
#' subtracting exactly those (spectrum-fit). With the defaults, (method='fir',
#' phase='zero', fir_window='hamming', fir_design='firwin',
#' pad='reflect_limited'), it is a zero-phase FIR filter.
#'
#' Internally implements notch as a bandstop FIR filter, exactly as MNE does:
#'   lows  = freq - notch_width/2 - trans_bandwidth/2
#'   highs = freq + notch_width/2 + trans_bandwidth/2
#' Then builds a bandstop kernel (gain=0 in the notch, gain=1 elsewhere).
#'
#' Multiple frequencies (e.g. harmonics) are combined into a single multi-band
#' FIR kernel and applied in one pass.
#'
#' \strong{IIR filters} (\code{method = "iir"}) are Butterworth band-stop
#' filters, built from their poles and zeros and run as a chain of
#' second-order sections, so they stay accurate at high sampling rates (see
#' \code{\link{eeg_bandpass}}). The two cutoffs of a notch are the \code{lows}
#' and \code{highs} given above: the points where ONE pass of the filter is
#' 3 dB down. So \code{trans_bandwidth} only moves the cutoffs outwards (by
#' half of it on each side); unlike for the FIR filter it is not the width of
#' a transition band. With \code{phase = "zero"} the filter is run forward and
#' then backward: there is no delay, the response is 6 dB down at the cutoffs,
#' and, because the backward pass looks at the signal before an event as well,
#' an event can leave a trace before it happened. With \code{phase = "forward"}
#' the filter is run once: nothing can appear before an event (the filter is
#' causal), but different frequencies are delayed by different amounts. Every
#' frequency in \code{freqs} gets its own band-stop filter, and the filters are
#' applied one after the other, from the lowest frequency to the highest: the
#' result is the same as calling \code{eeg_notch()} once for each frequency.
#' Only the Butterworth family is implemented so far.
#'
#' \strong{Spectrum-fit} (\code{method = "spectrum_fit"}) does not cut a slot
#' out of the frequencies. It looks for steady sine waves ("lines") in the
#' recording, measures the height and the timing of each one, and subtracts
#' exactly that wave; the rest of the recording is not touched. Every channel
#' is cut into overlapping pieces (10 s by default, see \code{filter_length}).
#' In each piece the data are looked at through several smooth windows
#' (Slepian windows), and for every frequency of a ruler an F score says how
#' much stronger a steady wave is than the noise around it. Frequencies whose
#' score is above a threshold (set by \code{p_value}) are lines. The cleaned
#' pieces are joined again with a smooth cross-fade. With \code{freqs = NULL},
#' which is also what you get when you leave \code{freqs} out for this method,
#' the lines are searched automatically, harmonics included, and the
#' frequencies that were found are reported. With \code{freqs} given no test is
#' made: in every piece, the frequencies of the ruler that lie inside
#' \code{notch_widths} around each of them are removed.
#'
#' \strong{The finer ruler.} The standard ruler has one mark every 1/T Hz, T
#' being the piece length (0.1 Hz for 10 s pieces), and a line that falls
#' between two marks is not found at all. \code{grid_factor} (default 8) makes
#' the ruler that many times finer, by zero-padding, so that a line is never
#' more than 1 / (2 T \code{grid_factor}) Hz away from a mark. Example, on
#' made-up noise with a 20 uV hum: a hum at 50.37 Hz, between two marks, was
#' 100 percent still there with the standard ruler (\code{grid_factor = 1}) and
#' 4 percent with the default one; the worst case over hums between two marks is
#' 4.2 percent. Of every run of neighbouring marks that pass the test only the
#' strongest one is used, and \code{p_value} is shared over all marks, so false
#' alarms stay as rare as with the standard ruler. The finer ruler is used only
#' when the lines are searched (\code{freqs = NULL}); with \code{freqs} given
#' \code{grid_factor} has no effect (a warning says so if you set it). It costs time: about 0.2 s per minute of
#' one channel at 512 Hz (0.02 s with \code{grid_factor = 1}), plus about 0.3 s
#' once for the windows, so downsample first if the sampling rate is higher.
#'
#' Limits of spectrum-fit: a hum whose frequency moves during a piece is only
#' partly removed (use shorter pieces, for example \code{filter_length =
#' "5s"}); any wave that stays steady for a whole piece counts as a line, so
#' look at the frequencies that are reported; a piece that contains missing
#' values is left unchanged, with a warning. A longer explanation with figures
#' is in the file \code{notes/finer-ruler-spectrum-fit.pdf} of the package
#' sources.
#'
#' @references
#' Butterworth, S. (1930). On the theory of filter amplifiers. Wireless
#' Engineer, 7, 536-541.
#'
#' Oppenheim, A. V., & Schafer, R. W. (2010). Discrete-Time Signal Processing
#' (3rd ed.). Pearson.
#'
#' Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
#' Proceedings of the IEEE, 70, 1055-1096.
#'
#' Percival, D. B., & Walden, A. T. (1993). Spectral Analysis for Physical
#' Applications. Cambridge University Press.
#'
#' Mitra, P., & Bokil, H. (2008). Observed Brain Dynamics. Oxford University
#' Press.
#'
#' @param eeg_obj          An object of class 'eeg'.
#' @param freqs            Numeric vector or NULL. Frequencies to notch in Hz.
#'                         E.g. c(50, 100, 150) for line noise + harmonics.
#'                         NULL, only with \code{method = "spectrum_fit"},
#'                         searches for the lines automatically; for that
#'                         method it is also what you get when you leave
#'                         \code{freqs} out. Default: 50.
#' @param notch_widths     Numeric or NULL. Width of the stop band in Hz
#'                         centred at each freq. NULL uses:
#'                         freq / 200. Give one number for all frequencies or
#'                         one per frequency; each must be 0 or more (anything
#'                         else is an error). With \code{method =
#'                         "spectrum_fit"} it is the width of the zone that is
#'                         removed, and it is used only when \code{freqs} is
#'                         given. Default: NULL.
#' @param trans_bandwidth  Numeric. Width of each transition band in Hz.
#'                         Default: 1. With method = "iir" it only moves the
#'                         two cutoffs of each notch outwards, by half of it.
#' @param method           Character. "fir" (default), "iir" or
#'                         "spectrum_fit".
#' @param iir_params       NULL or a named list, used only with
#'                         \code{method = "iir"}: \code{order} (default 4; a
#'                         band-stop has twice as many poles) and \code{ftype}
#'                         (only "butter" so far).
#' @param filter_length    Integer or "auto". With method = "fir": kernel
#'                         length in samples; "auto" uses MNE's formula. With
#'                         method = "spectrum_fit": the length of each piece,
#'                         as a number of samples or as text such as "5s" or
#'                         "500ms"; "auto" is 10 seconds. Not used with
#'                         method = "iir". Default: "auto".
#' @param fir_window       Character. Window function. Default: "hamming".
#'                         Used only with method = "fir".
#' @param phase            Character. The allowed words depend on
#'                         \code{method}. "fir": "zero" (default), one pass
#'                         with the group delay corrected; or "minimum", one
#'                         pass with a filter that only looks at the past
#'                         (causal). "iir": "zero" (default), forward then
#'                         backward, no delay; or "forward", one pass, causal.
#'                         Any other word is an error.
#' @param pad              Character. Padding mode. Default: "reflect_limited".
#' @param mt_bandwidth     Numeric or NULL. Only for \code{method =
#'                         "spectrum_fit"}. The bandwidth in Hz of the smooth
#'                         windows: how wide a blur one steady wave is seen
#'                         with. The number of windows is this bandwidth times
#'                         the piece length in seconds (at least 3 are
#'                         needed). Wider finds weaker lines and costs more
#'                         time, but two lines closer together than the
#'                         bandwidth cannot be told apart. NULL = 8 windows,
#'                         which is 0.8 Hz for 10 s pieces. Default: NULL.
#' @param p_value          Numeric between 0 and 1. Only for \code{method =
#'                         "spectrum_fit"} with \code{freqs = NULL}. How easily
#'                         a frequency is called a line: roughly the chance that
#'                         a piece of pure noise gets a false line. Lower is
#'                         stricter (fewer false lines, weak lines are missed
#'                         more often). With \code{freqs} given it has no effect
#'                         (a warning says so if you set it). Default: 0.05.
#' @param grid_factor      Whole number of 1 or more. Only for \code{method =
#'                         "spectrum_fit"} with \code{freqs = NULL}. How many
#'                         times finer than the standard ruler the frequencies
#'                         are tested (see Details). 1 is the standard ruler,
#'                         which misses a line that lies between two marks.
#'                         With \code{freqs} given it has no effect (a warning
#'                         says so if you set it). Default: 8.
#' @param channels         Character or integer vector, or NULL. Channels to
#'                         filter. NULL = all channels. Default: NULL.
#' @param verbose          Logical. Print parameters and progress.
#'                         Default: TRUE.
#'
#' @return An object of class 'eeg' with filtered data and updated
#'   preprocessing_history. With the spectrum-fit method, the lines that were
#'   found are also printed (with \code{verbose = TRUE}) and written into the
#'   preprocessing_history entry.
#'
#' @export
eeg_notch <- function(eeg_obj,
                      freqs            = 50,
                      notch_widths     = NULL,
                      trans_bandwidth  = 1,
                      method           = c("fir", "iir", "spectrum_fit"),
                      iir_params       = NULL,
                      filter_length    = "auto",
                      fir_window       = "hamming",
                      phase            = "zero",
                      pad              = "reflect_limited",
                      mt_bandwidth     = NULL,
                      p_value          = 0.05,
                      grid_factor      = 8L,
                      channels         = NULL,
                      verbose          = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  method <- match.arg(method)       # "fir" unless the user asks for another one

  # freqs = NULL means "search for the lines yourself", and only spectrum-fit
  # can do that. For spectrum-fit, leaving freqs out gives the same (missing()
  # is TRUE when an argument was left at its default).
  if (method == "spectrum_fit") {
    if (missing(freqs)) freqs <- NULL
  } else if (is.null(freqs)) {
    stop("'freqs' = NULL (search for the lines automatically) needs ",
         "method = \"spectrum_fit\".", call. = FALSE)
  }
  if (!is.null(freqs) &&
      (!is.numeric(freqs) || length(freqs) == 0 || any(freqs <= 0))) {
    stop("'freqs' must be a non-empty numeric vector of positive values.",
         call. = FALSE)
  }
  if (!is.numeric(trans_bandwidth) || trans_bandwidth <= 0) {
    stop("'trans_bandwidth' must be a positive numeric value.", call. = FALSE)
  }
  if (method != "spectrum_fit") {
    phase <- .check_phase(phase, method)
    pad   <- .check_pad(pad)
  }

  if (method == "spectrum_fit") {
    # Arguments that only describe a filter are not used here: if one was set
    # on purpose, say so (missing() is TRUE when an argument was left at its
    # default).
    unused <- c(trans_bandwidth = !missing(trans_bandwidth),
                iir_params      = !is.null(iir_params),
                fir_window      = !missing(fir_window),
                phase           = !missing(phase),
                pad             = !missing(pad))
    if (is.null(freqs)) unused <- c(unused, notch_widths = !is.null(notch_widths))
    if (any(unused)) {
      warning("Not used with method = \"spectrum_fit\"",
              if (is.null(freqs)) " and freqs = NULL" else "", ": ",
              paste(names(unused)[unused], collapse = ", "), ".",
              call. = FALSE)
    }
    # The finer ruler and the p-value belong to the automatic search. With
    # freqs given there is no search, so if the user typed one of them, say so.
    if (!is.null(freqs)) {
      search_only <- c(grid_factor = !missing(grid_factor),
                       p_value     = !missing(p_value))
      if (any(search_only)) {
        warning("Only used when the lines are searched (freqs = NULL), not ",
                "when 'freqs' is given: ",
                paste(names(search_only)[search_only], collapse = ", "), ".",
                call. = FALSE)
      }
    }
    if (!is.null(mt_bandwidth) &&
        (!is.numeric(mt_bandwidth) || length(mt_bandwidth) != 1L ||
         !is.finite(mt_bandwidth) || mt_bandwidth <= 0)) {
      stop("'mt_bandwidth' must be NULL or one positive number (Hz).",
           call. = FALSE)
    }
    if (!is.numeric(p_value) || length(p_value) != 1L || !is.finite(p_value) ||
        p_value <= 0 || p_value >= 1) {
      stop("'p_value' must be one number between 0 and 1.", call. = FALSE)
    }
    if (!is.numeric(grid_factor) || length(grid_factor) != 1L ||
        !is.finite(grid_factor) || grid_factor < 1 ||
        grid_factor != round(grid_factor)) {
      stop("'grid_factor' must be a whole number of 1 or more.", call. = FALSE)
    }
    grid_factor <- as.integer(grid_factor)
  } else {
    sf_only <- c(mt_bandwidth = !is.null(mt_bandwidth),
                 p_value      = !missing(p_value),
                 grid_factor  = !missing(grid_factor))
    if (any(sf_only)) {
      warning("Only used with method = \"spectrum_fit\": ",
              paste(names(sf_only)[sf_only], collapse = ", "), ".",
              call. = FALSE)
    }
  }

  if (method == "iir") {
    iir_params <- .resolve_iir_params(iir_params)

    # These two only describe an FIR kernel. If the user set one of them on
    # purpose, say that it is not used (missing() is TRUE when an argument was
    # left at its default). trans_bandwidth is used by both methods.
    fir_only <- c(filter_length = !missing(filter_length),
                  fir_window    = !missing(fir_window))
    if (any(fir_only)) {
      warning("Not used with method = \"iir\": ",
              paste(names(fir_only)[fir_only], collapse = ", "), ".",
              call. = FALSE)
    }
  } else if (method == "fir" && !is.null(iir_params)) {
    warning("'iir_params' is not used with method = \"fir\".", call. = FALSE)
  }

  sfreq   <- eeg_obj$sampling_rate
  nyquist <- sfreq / 2.0

  if (!is.null(freqs) && any(freqs >= nyquist)) {
    stop("All 'freqs' must be less than the Nyquist frequency (",
         nyquist, " Hz).", call. = FALSE)
  }
  
  # ========== RESOLVE NOTCH WIDTHS ==========
  # Default: freq / 200 per frequency
  
  if (is.null(freqs)) {
    notch_widths <- NULL            # searching for the lines: there is no zone
  } else if (is.null(notch_widths)) {
    notch_widths <- freqs / 200.0
  } else {
    if (!is.numeric(notch_widths) || anyNA(notch_widths) ||
        any(notch_widths < 0)) {
      stop("'notch_widths' must be NULL or non-negative numbers.", call. = FALSE)
    }
    if (!(length(notch_widths) %in% c(1L, length(freqs)))) {
      stop("'notch_widths' must be NULL, one number, or one number per ",
           "frequency in 'freqs'.", call. = FALSE)
    }
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

  # ========== SPECTRUM-FIT: HAND OVER TO ITS OWN FUNCTION ==========

  if (method == "spectrum_fit") {
    return(.notch_spectrum_fit(eeg_obj, freqs, notch_widths, filter_length,
                               mt_bandwidth, p_value, grid_factor, ch_idx,
                               verbose))
  }

  # ========== IIR: HAND OVER TO ITS OWN FUNCTION ==========
  # Everything below this block is the FIR path, unchanged.

  if (method == "iir") {
    return(.notch_iir(eeg_obj, freqs, notch_widths, trans_bandwidth,
                      iir_params, phase, ch_idx, verbose))
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
  
  # ========== BUILD ONE MULTI-BAND FREQ/GAIN ARRAY ==========
  # MNE computes lows/highs for all notch frequencies, passes them all to
  # filter_data() at once, which builds a single multi-band FIR kernel.
  # We replicate that here: collect all passband/stopband edges, sort them,
  # and build one freq/gain spec covering every notch band simultaneously.

  filtered_data <- eeg_obj$data
  tb_2          <- trans_bandwidth / 2.0
  history_parts <- character(0)

  # --- Compute and validate all band edges ---
  low_edges  <- freqs - notch_widths / 2.0 - tb_2   # lower passband edge per notch
  high_edges <- freqs + notch_widths / 2.0 + tb_2   # upper passband edge per notch

  for (fi in seq_along(freqs)) {
    if (low_edges[fi] <= 0) {
      stop("Notch lower edge (", round(low_edges[fi], 4), " Hz) is <= 0 for freq=",
           freqs[fi], " Hz. Reduce trans_bandwidth or notch_width.", call. = FALSE)
    }
    if (high_edges[fi] >= nyquist) {
      stop("Notch upper edge (", round(high_edges[fi], 4), " Hz) >= Nyquist (",
           nyquist, " Hz) for freq=", freqs[fi], " Hz.", call. = FALSE)
    }
    history_parts <- c(history_parts,
                       paste0(freqs[fi], " Hz (\u00B1", # plus/minus symbol
                              round(notch_widths[fi] / 2 + tb_2, 4), " Hz)"))
  }

  # --- Build interleaved passband/stopband edges for all notches ---
  # Per notch:  f_p1 = low_edge,        f_s1 = low_edge  + tb_2  (stop starts)
  #             f_s2 = high_edge - tb_2  (stop ends),    f_p2 = high_edge
  # Gains:      f_p1->1, f_s1->0, f_s2->0, f_p2->1
  # Then sort all edges and prepend 0 / append Nyquist 

  f_p1_all <- low_edges / nyquist
  f_s1_all <- (low_edges  + tb_2) / nyquist
  f_s2_all <- (high_edges - tb_2) / nyquist
  f_p2_all <- high_edges / nyquist

  inner_freq <- c(f_p1_all, f_s1_all, f_s2_all, f_p2_all)
  inner_gain <- c(rep(1, length(freqs)), rep(0, length(freqs)),
                  rep(0, length(freqs)), rep(1, length(freqs)))

  # Sort by frequency 
  ord        <- order(inner_freq)
  inner_freq <- inner_freq[ord]
  inner_gain <- inner_gain[ord]

  # Prepend 0 and append 1 (Nyquist) with gain = 1 
  freq_arr <- c(0, inner_freq, 1)
  gain_arr <- c(1, inner_gain, 1)

  # Clamp to [0, 1] for floating-point safety
  freq_arr <- pmax(pmin(freq_arr, 1.0), 0.0)

  # --- Single filter length based on global minimum trans bandwidth ---
  # MNE uses min across all bands (both sides = tb_2 here, so min = tb_2)
  if (identical(filter_length, "auto")) {
    N <- .auto_filter_length(tb_2, tb_2, sfreq, window = fir_window)
  } else {
    N <- as.integer(filter_length)
    if (N %% 2L == 0L) N <- N + 1L
  }

  # --- Design one kernel for all notch bands ---
  h <- .firwin_kernel(N, freq_arr, gain_arr, window = fir_window, sfreq = sfreq)

  # "minimum": same as in eeg_bandpass(): a causal kernel, no delay correction
  if (phase == "minimum") h <- .minimum_phase(h)

  if (verbose) {
    for (fi in seq_along(freqs)) {
      cat(sprintf("  Notch %d/%d: %.1f Hz  [%.4f - %.4f Hz stop]\n",
                  fi, length(freqs), freqs[fi], low_edges[fi], high_edges[fi]))
    }
    cat(sprintf("  Combined kernel: %d samples (single pass)\n", N))
    pb <- txtProgressBar(min = 0, max = length(ch_idx),
                         style = 3, width = 50, char = "=")
  }

  # --- Apply the single multi-band kernel to each channel ---
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
  
  out <- new_eeg(
    data                  = filtered_data,
    channels              = eeg_obj$channels,
    sampling_rate         = eeg_obj$sampling_rate,
    times                 = eeg_obj$times,
    events                = eeg_obj$events,
    metadata              = eeg_obj$metadata,
    reference             = eeg_obj$reference,
    preprocessing_history = c(eeg_obj$preprocessing_history,
                              list(history_entry)),
    montage               = eeg_obj$montage,
    bads                  = eeg_obj$bads,
    annotations           = eeg_obj$annotations,
    channel_types         = eeg_obj$channel_types
  )

  out
}

# End of filter1.R
