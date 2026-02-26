#' ============================================================================
#'                        EEG Signal Visualisation
#' ============================================================================
#'
#' This module provides a function for visualising the raw continuous EEG
#' signal from a loaded eeg object. It allows the user to inspect any single
#' electrode over any time window of the recording in milliseconds, making it
#' suitable for quick visual checks before and after preprocessing steps such
#' as filtering, re-referencing, or downsampling.
#'
#' The function draws inspiration from the plotting conventions established in
#' epoch2.R and follows the same code style and structure used across the
#' eeganalysis package.
#'
#' Author: Christos Dalamarinis
#' Date: 2026
#' ============================================================================
#'
#' Plot a Continuous EEG Signal from a Single Electrode
#'
#' Visualises the raw or preprocessed EEG time series for a chosen electrode
#' over a specified time window. Time is displayed in milliseconds on the
#' x-axis and amplitude in microvolts on the y-axis. Event markers present
#' in the EEG object are optionally drawn as vertical lines.
#'
#' @param eeg_obj An object of class 'eeg' containing continuous EEG data.
#'
#' @param channel Character string or integer. The electrode to plot.
#'                 Either the channel name (e.g. \code{"Cz"}) or its index
#'                 (e.g. \code{48}). Default: \code{"Cz"}.
#'
#' @param tmin Numeric. Start of the time window to display, in milliseconds.
#'             Default: \code{NULL} — plots from the very beginning of the
#'             recording.
#'
#' @param tmax Numeric. End of the time window to display, in milliseconds.
#'             Default: \code{NULL} — plots to the very end of the recording.
#'
#' @param demean Logical. If \code{TRUE} (default), subtracts the mean of the
#'               plotted window from the signal before drawing, centering the
#'               trace around 0. This mirrors the baseline correction applied
#'               in \code{epoch2.R} and makes the y-axis scale comparable to
#'               evoked response plots. Set to \code{FALSE} to display the raw
#'               amplifier values.
#'
#' @param show_events Logical. If \code{TRUE} (default), draws a vertical
#'                    dashed line at each event onset that falls within the
#'                    time window, labelled with the trigger code.
#'
#' @param event_color Character string. Colour used for event marker lines
#'                    and labels. Default: \code{"red"}.
#'
#' @param line_color Character string. Colour of the EEG signal trace.
#'                   Default: \code{"black"}.
#'
#' @param line_width Numeric. Line width of the EEG signal trace.
#'                   Default: \code{1}.
#'
#' @param title Character string. Custom plot title. If \code{NULL} (default),
#'              a title is generated automatically from the channel name and
#'              time window.
#'
#' @param y_lim Numeric vector of length 2: \code{c(min, max)} amplitude
#'              in microvolts for the y-axis. If \code{NULL} (default),
#'              the y-axis is scaled automatically to the data range.
#'
#' @return Invisibly returns a list with:
#'   \describe{
#'     \item{channel}{Channel name plotted.}
#'     \item{tmin_ms}{Actual start of the plotted window in ms.}
#'     \item{tmax_ms}{Actual end of the plotted window in ms.}
#'     \item{n_samples}{Number of samples in the plotted window.}
#'     \item{n_events}{Number of event markers drawn.}
#'   }
#'
#' @details
#' The function converts the \code{tmin} and \code{tmax} arguments (in
#' milliseconds) to sample indices using the recording's sampling rate, then
#' slices the data matrix accordingly. The x-axis always displays time in
#' milliseconds for consistency with standard EEG analysis conventions.
#'
#' If \code{show_events = TRUE}, events are drawn as vertical dashed lines.
#' Trigger labels are printed at the top of the plot to avoid obscuring the
#' signal. If many events are present in the window, only the trigger code
#' is shown to keep the plot readable.
#'
#' @examples
#' \dontrun{
#'   # Plot Cz for the full recording
#'   plot_eeg_signal(eeg)
#'
#'   # Plot Cz between 1000 ms and 3000 ms
#'   plot_eeg_signal(eeg, channel = "Cz", tmin = 1000, tmax = 3000)
#'
#'   # Plot channel 48 without event markers
#'   plot_eeg_signal(eeg, channel = 48, show_events = FALSE)
#'
#'   # Custom colour and fixed y-axis
#'   plot_eeg_signal(eeg,
#'                   channel    = "Oz",
#'                   tmin       = 500,
#'                   tmax       = 2500,
#'                   line_color = "steelblue",
#'                   y_lim      = c(-100, 100))
#' }
#'
#' @seealso \code{\link{plot_epochs}} for visualising epoched data,
#'   \code{\link{eeg_bandpass}} for filtering before plotting,
#'   \code{\link{new_eeg}} for the EEG object structure.
#'
#' @importFrom graphics plot lines abline text title mtext
#' @importFrom grDevices adjustcolor
#' @export
plot_eeg_signal <- function(eeg_obj,
                            channel     = "Cz",
                            tmin        = NULL,
                            tmax        = NULL,
                            demean      = TRUE,
                            show_events = TRUE,
                            event_color = "red",
                            line_color  = "black",
                            line_width  = 1,
                            title       = NULL,
                            y_lim       = NULL) {
  
  # ========== INPUT VALIDATION ==========
  
  if (!inherits(eeg_obj, "eeg")) {
    stop("'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  
  if (is.null(eeg_obj$data)) {
    stop("'eeg_obj$data' is NULL. No signal data to plot.", call. = FALSE)
  }
  
  n_channels   <- nrow(eeg_obj$data)
  n_timepoints <- ncol(eeg_obj$data)
  sr           <- eeg_obj$sampling_rate
  times_sec    <- eeg_obj$times                          # time vector in seconds
  times_ms     <- times_sec * 1000                       # convert to milliseconds
  
  # ---- Resolve channel ----
  if (is.character(channel)) {
    ch_idx <- match(channel, eeg_obj$channels)
    if (is.na(ch_idx)) {
      stop("Channel '", channel, "' not found in eeg_obj. ",
           "Available channels: ", paste(head(eeg_obj$channels, 10), collapse = ", "),
           if (n_channels > 10) paste0(" ... (", n_channels, " total)") else "",
           call. = FALSE)
    }
    ch_name <- channel
  } else if (is.numeric(channel)) {
    ch_idx <- as.integer(channel)
    if (ch_idx < 1 || ch_idx > n_channels) {
      stop("Channel index ", ch_idx, " is out of range [1, ", n_channels, "].",
           call. = FALSE)
    }
    ch_name <- eeg_obj$channels[ch_idx]
  } else {
    stop("'channel' must be a channel name (character) or index (integer).",
         call. = FALSE)
  }
  
  # ---- Resolve time window (ms → sample indices) ----
  rec_start_ms <- times_ms[1]
  rec_end_ms   <- times_ms[n_timepoints]
  
  if (is.null(tmin)) {
    tmin <- rec_start_ms
  }
  if (is.null(tmax)) {
    tmax <- rec_end_ms
  }
  
  if (!is.numeric(tmin) || length(tmin) != 1) {
    stop("'tmin' must be a single numeric value in milliseconds.", call. = FALSE)
  }
  if (!is.numeric(tmax) || length(tmax) != 1) {
    stop("'tmax' must be a single numeric value in milliseconds.", call. = FALSE)
  }
  if (tmin >= tmax) {
    stop("'tmin' (", tmin, " ms) must be less than 'tmax' (", tmax, " ms).",
         call. = FALSE)
  }
  if (tmin < rec_start_ms || tmax > rec_end_ms) {
    warning("Requested time window [", tmin, ", ", tmax, "] ms extends beyond ",
            "the recording [", round(rec_start_ms, 1), ", ",
            round(rec_end_ms, 1), "] ms. Clipping to recording bounds.",
            call. = FALSE, immediate. = TRUE)
    tmin <- max(tmin, rec_start_ms)
    tmax <- min(tmax, rec_end_ms)
  }
  
  # Find the sample indices closest to tmin and tmax
  samp_start <- which.min(abs(times_ms - tmin))
  samp_end   <- which.min(abs(times_ms - tmax))
  
  if (samp_start >= samp_end) {
    stop("Time window too narrow — fewer than 2 samples selected. ",
         "Increase the range between 'tmin' and 'tmax'.", call. = FALSE)
  }
  
  # ---- Validate y_lim ----
  if (!is.null(y_lim)) {
    if (!is.numeric(y_lim) || length(y_lim) != 2 || y_lim[1] >= y_lim[2]) {
      stop("'y_lim' must be a numeric vector of length 2: c(min, max).",
           call. = FALSE)
    }
  }
  
  # ---- Validate logical arguments ----
  if (!is.logical(show_events) || length(show_events) != 1) {
    stop("'show_events' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(demean) || length(demean) != 1) {
    stop("'demean' must be TRUE or FALSE.", call. = FALSE)
  }
  
  # ========== EXTRACT SIGNAL SLICE ==========
  
  signal_slice <- eeg_obj$data[ch_idx, samp_start:samp_end]
  
  # Subtract the mean of the plotted window so the trace is centred around 0,
  # matching the appearance of baseline-corrected epoch plots from epoch2.R.
  if (demean) {
    signal_slice <- signal_slice - mean(signal_slice, na.rm = TRUE)
  }
  time_slice   <- times_ms[samp_start:samp_end]
  n_samples    <- length(signal_slice)
  
  # ========== IDENTIFY EVENTS IN WINDOW ==========
  
  events_in_window <- data.frame(
    onset_ms = numeric(0),
    type     = character(0),
    stringsAsFactors = FALSE
  )
  
  if (show_events && !is.null(eeg_obj$events) && nrow(eeg_obj$events) > 0) {
    ev            <- eeg_obj$events
    ev_ms         <- ev$onset_time * 1000
    in_window     <- ev_ms >= tmin & ev_ms <= tmax
    events_in_window <- data.frame(
      onset_ms = ev_ms[in_window],
      type     = ev$type[in_window],
      stringsAsFactors = FALSE
    )
  }
  
  n_events <- nrow(events_in_window)
  
  # ========== COMPUTE PLOT LIMITS ==========
  
  if (is.null(y_lim)) {
    data_range <- range(signal_slice, na.rm = TRUE)
    # Add 10% padding above and below
    padding <- (data_range[2] - data_range[1]) * 0.10
    y_lim   <- c(data_range[1] - padding, data_range[2] + padding)
  }
  
  # ========== BUILD TITLE ==========
  
  if (is.null(title)) {
    title <- paste0(
      "EEG Signal  |  Channel: ", ch_name,
      "  |  ", round(tmin), " – ", round(tmax), " ms"
    )
  }
  
  # ========== DRAW PLOT ==========
  
  # Save and restore par on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  # Add top margin when events will be labelled
  if (show_events && n_events > 0) {
    par(mar = c(4.5, 4.5, 4, 1))
  } else {
    par(mar = c(4.5, 4.5, 3, 1))
  }
  
  # Empty plot frame
  plot(
    x    = range(time_slice),
    y    = y_lim,
    type = "n",
    xlab = "Time (ms)",
    ylab = "Amplitude (microV)",
    main = title,
    xaxt = "s",
    las  = 1
  )
  
  # Zero-amplitude reference line
  abline(h = 0, col = "gray70", lty = 2, lwd = 0.8)
  
  # Event markers
  if (show_events && n_events > 0) {
    for (i in seq_len(n_events)) {
      abline(v   = events_in_window$onset_ms[i],
             col = adjustcolor(event_color, alpha.f = 0.7),
             lty = 2,
             lwd = 1.2)
    }
    # Trigger labels just above the top of the plot (in outer margin)
    mtext(
      text  = events_in_window$type,
      side  = 3,
      at    = events_in_window$onset_ms,
      col   = event_color,
      cex   = 0.65,
      line  = 0.3
    )
  }
  
  # EEG signal trace (drawn last so it sits on top)
  lines(
    x   = time_slice,
    y   = signal_slice,
    col = line_color,
    lwd = line_width
  )
  
  # ========== INVISIBLE RETURN ==========
  
  invisible(list(
    channel  = ch_name,
    tmin_ms  = round(time_slice[1], 2),
    tmax_ms  = round(time_slice[n_samples], 2),
    n_samples = n_samples,
    n_events  = n_events
  ))
}

# End of plot_signal.R