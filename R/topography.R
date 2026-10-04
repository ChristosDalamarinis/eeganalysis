#' ============================================================================
#'                  Scalp Topography (Topoplot) Plotting
#' ============================================================================
#'
#' Renders a 2D interpolated scalp heatmap of one value per channel (e.g.
#' band power from eeg_band_power()), using electrode positions from a
#' montage attached to an eeg object (see set_montage(), create_montage()).
#'
#' Author: Christos Dalamarinis
#' Date: July 2026
#' ============================================================================

#' Plot a Scalp Topography (Topoplot)
#'
#' Draws a 2D interpolated heatmap of per-channel values (e.g. band power)
#' over the scalp, using electrode positions from an attached
#' \code{\link{create_montage}}-built montage. Each channel's 3D scalp
#' position is projected to 2D with the standard EEG azimuthal-equidistant
#' topomap projection (distance from the plot center is proportional to the
#' angular distance from the vertex, Cz), then the scattered per-channel
#' values are interpolated onto a regular grid with \code{akima::interp()}.
#'
#' @param eeg_obj An object of class \code{'eeg'}. Used to read
#'   \code{eeg_obj$montage} when \code{montage} is not supplied directly, and
#'   to read \code{eeg_obj$bads} so channels marked bad are always excluded.
#' @param values Named numeric vector of one value per channel (names must
#'   match montage channel names, e.g. \code{"Cz"}). A natural source is
#'   \code{\link{eeg_band_power}}, e.g.
#'   \code{setNames(power_df$alpha, power_df$channel)}.
#' @param montage An object of class \code{'montage'} giving channel
#'   positions. If \code{NULL} (default), \code{eeg_obj$montage} is used.
#' @param interpolate_res Integer, resolution (grid points per side) of the
#'   interpolation grid. Default \code{100}.
#' @param contour Logical, whether to draw contour isolines over the heatmap.
#'   Default \code{TRUE}.
#' @param show_electrodes Logical, whether to mark electrode positions (and
#'   names) on the plot. Default \code{TRUE}.
#' @param col_palette Character vector of colors passed to
#'   \code{grDevices::colorRampPalette()} for the heatmap. Default a
#'   blue-white-red diverging palette.
#' @param main Character string, plot title. Default \code{NULL}, which uses
#'   \code{"Scalp Topography"}.
#' @param zero_centered Logical, whether to make the colour scale symmetric
#'   around zero. The strongest plotted value of either sign then gets the
#'   full colour at its end of \code{col_palette} (the first colour for the
#'   most negative value, the last for the most positive), and the middle
#'   colour means exactly zero. Default \code{FALSE}, which stretches the
#'   palette from the lowest to the highest plotted value, so the middle
#'   colour is the middle of that range and not zero. Turn it on for signed
#'   maps such as the output of \code{\link{compute_csd}}; leave it off for
#'   maps that are always positive, such as band power, where half of the
#'   palette would go unused.
#' @param extrapolate Logical, whether to colour the area outside the outer
#'   boundary of the electrodes. There is no data there, so the colours are
#'   extrapolated from the electrodes, and on sharp maps such as the output
#'   of \code{\link{compute_csd}} that guess can show values larger than any
#'   electrode measured. Default \code{TRUE} colours the whole head disc, as
#'   before. With \code{FALSE} that area is left empty (\code{NA} in
#'   \code{grid_z}).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with:
#'  \describe{
#'    \item{grid_x}{Numeric vector, interpolation grid x-axis}
#'    \item{grid_y}{Numeric vector, interpolation grid y-axis}
#'    \item{grid_z}{Numeric matrix, the interpolated values (\code{NA} outside
#'      the head disc, and also outside the electrodes when
#'      \code{extrapolate = FALSE})}
#'    \item{zlim}{Numeric vector of length 2, the lowest and highest value the
#'      colour scale runs between (symmetric around zero when
#'      \code{zero_centered = TRUE}); handy for drawing a colour bar}
#'    \item{channel_positions}{Data frame of the 2D-projected electrode
#'      positions and values actually plotted}
#'  }
#'
#' @details
#' This function requires an attached montage: build one with
#' \code{\link{create_montage}} and attach it with \code{\link{set_montage}},
#' or supply the \code{montage} argument directly.
#'
#' Channels present in \code{values} but absent from the montage, and montage
#' channels with no supplied value, are both dropped with a warning. Channels
#' listed in \code{eeg_obj$bads} (see \code{\link{new_eeg}}) are also dropped
#' with a warning, the same shared exclude list used elsewhere in the
#' package, so a channel marked bad once stays out of the topography too. At
#' least 3 channels with both a position and a value are required to
#' interpolate a surface.
#'
#' Maps of signed quantities, such as current source density from
#' \code{\link{compute_csd}} (positive and negative values mean different
#' things), are easier to read with \code{zero_centered = TRUE} and
#' \code{extrapolate = FALSE} together: the first makes the middle colour
#' mean zero, the second keeps the picture from showing values that were
#' only guessed outside the electrodes.
#'
#' @examples
#' \dontrun{
#'   eeg <- set_montage(eeg, create_montage())
#'   power <- eeg_band_power(eeg)
#'   plot_topography(eeg, setNames(power$alpha, power$channel))
#'
#'   # A signed map such as CSD: white means zero, and nothing is painted
#'   # beyond the electrodes
#'   csd <- compute_csd(eeg)
#'   plot_topography(csd, setNames(csd$data[, 1], csd$channels),
#'                   zero_centered = TRUE, extrapolate = FALSE)
#' }
#'
#' @seealso \code{\link{create_montage}}, \code{\link{set_montage}},
#'   \code{\link{eeg_band_power}}, \code{\link{compute_csd}},
#'   \code{\link{plot_eeg_signal}}
#'
#' @importFrom graphics image contour polygon points text lines par
#' @importFrom grDevices colorRampPalette
#' @importFrom akima interp
#'
#' @export
plot_topography <- function(eeg_obj,
                             values,
                             montage = NULL,
                             interpolate_res = 100,
                             contour = TRUE,
                             show_electrodes = TRUE,
                             col_palette = c("blue", "white", "red"),
                             main = NULL,
                             zero_centered = FALSE,
                             extrapolate = TRUE,
                             ...) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }

  if (is.null(montage)) {
    montage <- eeg_obj$montage
  }
  if (is.null(montage) || !inherits(montage, "montage")) {
    stop("ERROR: No montage available. Attach one with set_montage(), or ",
         "pass one explicitly via the 'montage' argument.", call. = FALSE)
  }

  if (is.null(values) || is.null(names(values))) {
    stop("ERROR: 'values' must be a named numeric vector ",
         "(names = channel names).", call. = FALSE)
  }

  # ========== MATCH VALUES TO MONTAGE CHANNELS ==========

  common         <- intersect(names(values), montage$channels)
  missing_values <- setdiff(montage$channels, names(values))
  extra_values   <- setdiff(names(values), montage$channels)

  if (length(extra_values) > 0) {
    warning("plot_topography(): ", length(extra_values),
            " value(s) have no matching montage channel and will be ",
            "ignored: ", paste(extra_values, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }
  if (length(missing_values) > 0) {
    warning("plot_topography(): ", length(missing_values),
            " montage channel(s) have no value supplied and will be ",
            "excluded: ", paste(missing_values, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  # Bad channels are a shared exclude list (eeg_obj$bads, see new_eeg()) -
  # read here rather than taking a separate per-call exclude argument, so a
  # channel marked bad once stays excluded from the topography too.
  bad_common <- intersect(common, eeg_obj$bads)
  if (length(bad_common) > 0) {
    warning("plot_topography(): ", length(bad_common),
            " channel(s) excluded because they are marked bad ",
            "(eeg_obj$bads): ", paste(bad_common, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
    common <- setdiff(common, bad_common)
  }

  if (length(common) < 3) {
    stop("ERROR: At least 3 channels with both a montage position and a ",
         "value are required to interpolate a topography (found ",
         length(common), ").", call. = FALSE)
  }

  pos  <- montage$positions[match(common, montage$positions$channel), , drop = FALSE]
  vals <- as.numeric(values[common])

  # ========== 2D AZIMUTHAL-EQUIDISTANT PROJECTION ==========
  # Derived from the montage's Cartesian x/y/z (not the theta/phi columns,
  # whose sign convention is ambiguous) - elevation is the angle from the
  # vertex (Cz, elevation = 0), azimuth is the angle around the vertical
  # axis. Distance from the plot center is proportional to elevation, which
  # is the standard EEG topomap projection.

  r         <- sqrt(pos$x^2 + pos$y^2 + pos$z^2)
  elevation <- acos(pos$z / r)
  azimuth   <- atan2(pos$y, pos$x)

  x2d <- elevation * cos(azimuth)
  y2d <- elevation * sin(azimuth)

  head_r <- max(abs(c(x2d, y2d))) * 1.1

  if (diff(range(x2d)) < 1e-8 || diff(range(y2d)) < 1e-8) {
    stop("ERROR: The selected channels are collinear (no spread in one ",
         "spatial dimension) and cannot be interpolated into a surface. ",
         "Include channels spread across both hemispheres and the ",
         "anterior-posterior axis.", call. = FALSE)
  }

  # ========== INTERPOLATE ONTO A REGULAR GRID ==========
  # Outside the outer boundary of the electrodes there is no data to
  # interpolate from. With extrapolate = TRUE the surface is extended over
  # that area anyway, and it can show values larger than any electrode
  # measured; with extrapolate = FALSE the area is left empty (NA) instead.

  fit <- akima::interp(x2d, y2d, vals,
                        xo = seq(-head_r, head_r, length.out = interpolate_res),
                        yo = seq(-head_r, head_r, length.out = interpolate_res),
                        linear = FALSE, extrap = isTRUE(extrapolate),
                        duplicate = "mean")

  grid_dist   <- sqrt(outer(fit$x^2, fit$y^2, "+"))
  fit$z[grid_dist > head_r] <- NA

  # ========== RENDER ==========

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(2, 2, 3, 2))

  palette_fn <- colorRampPalette(col_palette)

  # Colour limits. Left to itself, image() stretches the palette from the
  # lowest to the highest plotted value, which puts the middle colour at the
  # middle of that range, not at zero. zero_centered makes the range symmetric
  # instead, so the middle colour is exactly zero.
  zlim <- range(fit$z, na.rm = TRUE)
  if (isTRUE(zero_centered)) {
    zlim <- c(-1, 1) * max(abs(zlim))
  }

  image(fit$x, fit$y, fit$z,
        col = palette_fn(100), zlim = zlim,
        xlab = "", ylab = "", axes = FALSE, asp = 1,
        main = if (is.null(main)) "Scalp Topography" else main)

  if (isTRUE(contour)) {
    contour(fit$x, fit$y, fit$z, add = TRUE, col = "black", lwd = 0.5, labcex = 0.6)
  }

  # Head outline: circle + nose
  circle_theta <- seq(0, 2 * pi, length.out = 200)
  polygon(head_r * cos(circle_theta), head_r * sin(circle_theta), lwd = 1.5)
  lines(c(-0.08, 0, 0.08) * head_r, c(0.99, 1.12, 0.99) * head_r, lwd = 1.5)

  if (isTRUE(show_electrodes)) {
    points(x2d, y2d, pch = 21, bg = "black", col = "white", cex = 0.9)
    text(x2d, y2d, labels = common, pos = 3, cex = 0.6, offset = 0.4)
  }

  invisible(list(
    grid_x = fit$x,
    grid_y = fit$y,
    grid_z = fit$z,
    zlim = zlim,
    channel_positions = data.frame(channel = common, x2d = x2d, y2d = y2d,
                                    value = vals, stringsAsFactors = FALSE)
  ))
}
