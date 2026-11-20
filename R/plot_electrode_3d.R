# ============================================================================
# 3D plots of electrode positions
# ============================================================================
#
#   plot_electrode_3d()            3D plot from x/y/z positions
#   plot_electrode_3d_spherical()  3D plot from theta/phi angles
#
# By default both draw the 64 BioSemi electrodes (their numbers are kept in
# .biosemi64_table(), R/montage_biosemi.R). With the argument `layout` they
# draw any built-in layout, by name, or a montage object. The names of the
# layouts and their tables come from R/montage_templates.R, the front door of
# all layouts.
#
# The two small helpers under the banner below are used by both plots.
# ============================================================================

# ----------------------------------------------------------------------------
# Layouts for the 3D electrode plots
#
# plot_electrode_3d() and plot_electrode_3d_spherical() draw the 64 BioSemi
# electrodes. Given a layout (the name of a built-in layout, or a montage
# object) they draw that one instead. .plot_layout_table() turns a layout into
# one table (electrode, alias, x, y, z, theta, phi) and .match_layout_names()
# finds the electrode behind each channel name that was asked for.
# ----------------------------------------------------------------------------

.plot_layout_table <- function(layout) {
  if (inherits(layout, "montage")) {
    pos <- layout$positions
    return(data.frame(electrode = pos$channel, alias = NA_character_,
                      x = pos$x, y = pos$y, z = pos$z,
                      theta = pos$theta, phi = pos$phi,
                      stringsAsFactors = FALSE))
  }

  available <- .montage_templates()
  if (!is.character(layout) || length(layout) != 1 || is.na(layout) ||
      !(layout %in% available)) {
    stop("ERROR: 'layout' must be the name of a built-in layout (",
         paste(available, collapse = ", "), ") or a montage object.",
         call. = FALSE)
  }

  tab <- .montage_template_table(layout)
  data.frame(electrode = tab$name, alias = tab$alias,
             x = tab$x, y = tab$y, z = tab$z,
             theta = tab$theta, phi = tab$phi,
             stringsAsFactors = FALSE)
}

# For each channel name: the name of the electrode in the layout table, looked
# up by its name first and then by its alias, in any letter case. A name the
# layout does not know is returned as it was typed.
.match_layout_names <- function(channels, layout_table) {
  keys     <- tolower(channels)
  by_name  <- match(keys, tolower(layout_table$electrode))
  by_alias <- match(keys, tolower(layout_table$alias))
  by_alias[is.na(keys)] <- NA_integer_

  row <- ifelse(is.na(by_name), by_alias, by_name)
  ifelse(is.na(row), channels, layout_table$electrode[row])
}


#' ============================================================================
#' Plot 3D Electrode Layout on Scalp (BioSemi 64-Channel)
#'
#' Creates an interactive 3D visualization of electrode positions on a realistic
#' head model using OFFICIAL BioSemi 64-channel electrode coordinates.
#'
#' @param channels Character vector - channel names to plot. If NULL, plots all
#'   64 standard 10-20/10-10 electrodes. Can also accept BioSemi A/B naming.
#' @param highlight_channels Character vector - specific channels to highlight
#'   in a different color (optional)
#' @param show_labels Logical - whether to show electrode labels on the plot.
#'   Default is TRUE.
#' @param color_by Character - how to color electrodes. Options:
#'   - "region" = color by brain region (frontal, central, parietal, occipital, temporal)
#'   - "single" = all one color (default blue)
#'   Default is "region".
#' @param show_head Logical - whether to show the head outline. Default is TRUE.
#' @param layout The electrodes to draw. \code{NULL} (the default) draws the 64
#'   BioSemi electrodes. Otherwise the name of a built-in layout (for example
#'   \code{"standard_1005"} or \code{"biosemi128"}; see
#'   \code{\link{create_montage}} for the names) or a montage object, such as
#'   one made by \code{\link{read_custom_montage}}. With a layout,
#'   \code{channels} picks from that layout: a channel is found by its name or
#'   by its second name (the old names T3 to T6, or the BioSemi socket names),
#'   in any letter case.
#'
#' @return Returns a plotly object (interactive 3D plot viewable in RStudio)
#'
#' @details
#' This function visualizes electrode positions in 3D space using the official
#' BioSemi 64-channel electrode cap coordinates. All positions are exact Cartesian
#' coordinates (x, y, z) in millimeters from BioSemi's cap layout documentation.
#'
#' **Coordinate System:**
#' - Origin (0, 0, 0) at center of head
#' - X-axis: Left (-) to Right (+)
#' - Y-axis: Back (-) to Front (+)
#' - Z-axis: Bottom (-) to Top (+)
#' - Nose points in +Y direction
#'
#' **Other layouts:**
#' With \code{layout} the same plot is drawn for any built-in layout or montage
#' object. The positions are the layout's own x/y/z (in mm), the title names the
#' layout, and the electrodes are coloured by region as usual. For a layout of
#' hundreds of electrodes, such as \code{"standard_1005"}, use
#' \code{show_labels = FALSE}: the names overlap.
#'
#' @examples
#' \dontrun{
#'
#'   # Example 1: Plot all 64 electrodes
#'   plot_electrode_3d()
#'
#'   # Example 2: Plot specific channels
#'   plot_electrode_3d(channels = c("Cz", "Fz", "Pz", "Oz", "T7", "T8"))
#'
#'   # Example 3: Highlight midline electrodes
#'   plot_electrode_3d(highlight_channels = c("Fpz", "Fz", "Cz", "Pz", "Oz"))
#'
#'   # Example 4: No head outline, just electrodes
#'   plot_electrode_3d(show_head = FALSE)
#'
#'   # Example 5: Another layout, here the 335 electrodes of the 10-05 system
#'   plot_electrode_3d(layout = "standard_1005", show_labels = FALSE)
#'
#'   # Example 6: Some channels of a layout, given by their old names
#'   plot_electrode_3d(channels = c("Cz", "T3", "T4"), layout = "standard_1020")
#' }
#'
#' @keywords internal
#' @importFrom plotly plot_ly add_trace layout
#' @export
#'
plot_electrode_3d <- function(channels = NULL,
                              highlight_channels = NULL,
                              show_labels = TRUE,
                              color_by = "region",
                              show_head = TRUE,
                              layout = NULL) {
  
  # ========================================================================
  # LOAD REQUIRED PACKAGE
  # ========================================================================
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive 3D plots.\n",
         "Install it with: install.packages('plotly')")
  }
  
  # ========================================================================
  # THE ELECTRODES TO DRAW
  # ========================================================================
  # The electrodes and their positions come from a layout. With no layout they
  # are the 64 BioSemi electrodes: the "biosemi64" layout, the one montages use,
  # so these numbers are kept in one place (.biosemi64_table(), in
  # R/montage_biosemi.R).
  
  layout_table     <- .plot_layout_table(if (is.null(layout)) "biosemi64" else layout)
  layout_label     <- if (is.character(layout)) layout else "Custom montage"
  electrode_coords <- layout_table[, c("electrode", "x", "y", "z")]
  
  # ========================================================================
  # DETERMINE WHICH CHANNELS TO PLOT
  # ========================================================================
  
  if (is.null(channels)) {
    channels_to_plot <- electrode_coords$electrode
  } else {
    channels_to_plot <- channels
  }
  
  # Map BioSemi names to standard names if needed
  plot_data <- data.frame(
    original_name = channels_to_plot,
    standard_name = channels_to_plot,
    stringsAsFactors = FALSE
  )
  
  # The layout's own names and second names, found in any letter case: for the
  # 64 BioSemi electrodes the socket names A1 to B32, for the other layouts
  # also old names such as T3.
  plot_data$standard_name <- .match_layout_names(channels_to_plot, layout_table)
  
  # Merge with coordinates
  plot_data <- merge(plot_data, electrode_coords, 
                     by.x = "standard_name", by.y = "electrode",
                     all.x = TRUE)
  
  # Remove electrodes without coordinates
  missing_coords <- plot_data$original_name[is.na(plot_data$x)]
  if (length(missing_coords) > 0) {
    warning("No coordinates found for: ", paste(missing_coords, collapse = ", "))
    plot_data <- plot_data[!is.na(plot_data$x), ]
  }
  
  if (nrow(plot_data) == 0) {
    stop("No valid electrodes to plot.")
  }
  
  # ========================================================================
  # ASSIGN COLORS BY BRAIN REGION
  # ========================================================================
  
  if (color_by == "region") {
    plot_data$color_group <- "Other"
    plot_data$color_group[grepl("^Fp|^AF|^F[^CT]", plot_data$standard_name)] <- "Frontal"
    plot_data$color_group[grepl("^C|^FC", plot_data$standard_name)] <- "Central"
    plot_data$color_group[grepl("^P[^O]|^CP", plot_data$standard_name)] <- "Parietal"
    plot_data$color_group[grepl("^O|^PO|^I", plot_data$standard_name)] <- "Occipital"
    plot_data$color_group[grepl("^T[0-9]|^FT|^TP|^TT", plot_data$standard_name)] <- "Temporal"
    
    color_palette <- c(
      "Frontal" = "#e41a1c",
      "Central" = "#377eb8",
      "Parietal" = "#4daf4a",
      "Occipital" = "#984ea3",
      "Temporal" = "#ff7f00",
      "Other" = "#999999"
    )
    
  } else {
    plot_data$color_group <- "Electrode"
    color_palette <- c("Electrode" = "#377eb8")
  }
  
  # Highlight specific channels if requested
  if (!is.null(highlight_channels)) {
    plot_data$is_highlighted <- plot_data$original_name %in% highlight_channels |
      plot_data$standard_name %in% highlight_channels
  } else {
    plot_data$is_highlighted <- FALSE
  }
  
  # ========================================================================
  # CREATE REALISTIC HEAD OUTLINE (ELLIPSOID)
  # ========================================================================
  
  # Create ellipsoid head shape (more realistic than sphere)
  # Semi-axes: a (left-right), b (front-back), c (top-bottom)
  a <- 80   # Width (left-right)
  b <- 95   # Length (front-back)
  c <- 95   # Height (top-bottom, slightly elongated)
  
  # Parametric ellipsoid surface
  n_points <- 30
  u <- seq(0, 2*pi, length.out = n_points)
  v <- seq(0, pi, length.out = n_points)
  
  # Create wireframe circles (longitude and latitude lines)
  head_lines <- list()
  
  # Longitude lines (vertical circles from front to back)
  for (i in seq(1, length(u), by = 3)) {
    x_line <- a * sin(v) * cos(u[i])
    y_line <- b * sin(v) * sin(u[i])
    z_line <- c * cos(v)
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # Latitude lines (horizontal circles)
  for (j in seq(1, length(v), by = 3)) {
    x_line <- a * sin(v[j]) * cos(u)
    y_line <- b * sin(v[j]) * sin(u)
    z_line <- rep(c * cos(v[j]), length(u))
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # ========================================================================
  # CREATE NOSE INDICATOR
  # ========================================================================
  
  nose_x <- c(0, 8, -8, 0, 0)
  nose_y <- c(100, 88, 88, 100, 100)
  nose_z <- c(20, 12, 12, 20, 20)
  
  # ========================================================================
  # BUILD PLOTLY 3D FIGURE
  # ========================================================================
  
  p <- plotly::plot_ly()
  
  # Add head wireframe (if requested)
  if (show_head) {
    for (line in head_lines) {
      p <- plotly::add_trace(
        p,
        x = line$x, y = line$y, z = line$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = "lightgray", width = 1),
        hoverinfo = "none",
        showlegend = FALSE
      )
    }
  }
  
  # Add nose indicator
  p <- plotly::add_trace(
    p,
    x = nose_x, y = nose_y, z = nose_z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 5),
    hoverinfo = "none",
    showlegend = FALSE,
    name = "Nose"
  )
  
  # Add electrodes (by color group)
  for (group in unique(plot_data$color_group)) {
    group_data <- plot_data[plot_data$color_group == group, ]
    
    marker_sizes <- ifelse(group_data$is_highlighted, 14, 10)
    marker_borders <- ifelse(group_data$is_highlighted, "black", "white")
    border_widths <- ifelse(group_data$is_highlighted, 3, 2)
    
    p <- plotly::add_trace(
      p,
      data = group_data,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = ifelse(show_labels, "markers+text", "markers"),
      marker = list(
        size = marker_sizes,
        color = color_palette[group],
        line = list(color = marker_borders, width = border_widths)
      ),
      # The label settings are only given when labels are wanted: plotly puts
      # the text mode back whenever text or textfont are present.
      text = if (show_labels) ~original_name,
      textposition = if (show_labels) "top center",
      textfont = if (show_labels) list(size = 10, color = "black", family = "Arial Black"),
      hovertext = ~paste0("<b>", original_name, "</b><br>",
                          "Standard: ", standard_name, "<br>",
                          "Region: ", color_group, "<br>",
                          "Position (mm):<br>",
                          "  X: ", round(x, 1), " (L/R)<br>",
                          "  Y: ", round(y, 1), " (Back/Front)<br>",
                          "  Z: ", round(z, 1), " (Bottom/Top)"),
      hoverinfo = "text",
      name = group,
      showlegend = TRUE
    )
  }
  
  # ========================================================================
  # CONFIGURE LAYOUT
  # ========================================================================
  
  p <- plotly::layout(
    p,
    title = list(
      text = if (is.null(layout)) {
        "BioSemi 64-Channel Electrode Layout (3D)"
      } else {
        paste0(layout_label, " Electrode Layout (3D)")
      },
      font = list(size = 20, family = "Arial", color = "black")
    ),
    scene = list(
      xaxis = list(
        title = "Left | Right",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-110, 110)
      ),
      yaxis = list(
        title = "Back | Front (Nose)",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-110, 110)
      ),
      zaxis = list(
        title = "Bottom | Top",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-60, 110)
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      ),
      aspectmode = "cube",
      bgcolor = "white"
    ),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    legend = list(
      x = 0.02,
      y = 0.98,
      bgcolor = "rgba(255, 255, 255, 0.9)",
      bordercolor = "black",
      borderwidth = 1,
      font = list(size = 12)
    )
  )
  
  return(p)
}








#' ============================================================================
#' Plot 3D Electrode Layout on Scalp (BioSemi 64-Channel - Spherical)
#'
#' Creates an interactive 3D visualization of electrode positions on a spherical
#' head model using OFFICIAL BioSemi 64-channel SPHERICAL coordinates.
#'
#' @param channels Character vector - channel names to plot. If NULL, plots all
#'   64 standard 10-20/10-10 electrodes. Can also accept BioSemi A/B naming.
#' @param highlight_channels Character vector - specific channels to highlight
#'   in a different color (optional)
#' @param show_labels Logical - whether to show electrode labels on the plot.
#'   Default is TRUE.
#' @param color_by Character - how to color electrodes. Options:
#'   - "region" = color by brain region (frontal, central, parietal, occipital, temporal)
#'   - "single" = all one color (default blue)
#'   Default is "region".
#' @param show_head Logical - whether to show the head outline. Default is TRUE.
#' @param head_radius Numeric - radius of head sphere in mm. Default is 87.54 mm
#'   (calculated from 55 cm circumference).
#' @param layout The electrodes to draw. \code{NULL} (the default) draws the 64
#'   BioSemi electrodes. Otherwise the name of a built-in layout (for example
#'   \code{"standard_1005"} or \code{"biosemi128"}; see
#'   \code{\link{create_montage}} for the names) or a montage object, such as
#'   one made by \code{\link{read_custom_montage}}. The electrodes are placed
#'   from the theta and phi of the layout on the head sphere of radius
#'   \code{head_radius}, and are not rounded to whole millimetres. With a
#'   layout, \code{channels} picks from that layout: a channel is found by its
#'   name or by its second name (the old names T3 to T6, or the BioSemi socket
#'   names), in any letter case.
#'
#' @return Returns a plotly object (interactive 3D plot viewable in RStudio)
#'
#' @details
#' This function visualizes electrode positions in 3D space using the official
#' BioSemi 64-channel electrode cap SPHERICAL coordinates (theta, phi), which are then
#' converted to Cartesian (x, y, z) for plotting.
#'
#' **Spherical Coordinate System:**
#' - theta (inclination): range from +180deg to -180deg (at Cz, inclination = 0deg)
#'   - Positive inclinations = right side of head
#'   - Negative inclinations = left side of head
#' - phi (azimuth): range from +90deg to -90deg
#'   - For positive inclinations: phi = angle from T8 (right ear)
#'   - For negative inclinations: phi = angle from T7 (left ear)
#'   - Positive = anti-clockwise, negative = clockwise
#' - r (radius): 87.54 mm (from 55 cm head circumference)
#'
#' **Conversion to Cartesian:**
#' - x = r * sin(theta) * cos(phi)
#' - y = r * sin(theta) * sin(phi)
#' - z = r * cos(theta)
#'
#' **Coordinate System:**
#' - Origin (0, 0, 0) at center of head
#' - X-axis: Left (-) to Right (+)
#' - Y-axis: Back (-) to Front (+)
#' - Z-axis: Bottom (-) to Top (+)
#' - Nose points in +Y direction
#'
#' **Other layouts:**
#' With \code{layout} the same plot is drawn for any built-in layout or montage
#' object, from the theta and phi of its electrodes. The title names the layout,
#' and the electrodes are coloured by region as usual. For a layout of hundreds
#' of electrodes, such as \code{"standard_1005"}, use
#' \code{show_labels = FALSE}: the names overlap.
#'
#' @examples
#' \dontrun{
#'
#'   # Example 1: Plot all 64 electrodes
#'   plot_electrode_3d_spherical()
#'
#'   # Example 2: Plot specific channels
#'   plot_electrode_3d_spherical(channels = c("Cz", "Fz", "Pz", "Oz", "T7", "T8"))
#'
#'   # Example 3: Highlight midline electrodes
#'   plot_electrode_3d_spherical(highlight_channels = c("Fpz", "Fz", "Cz", "Pz", "Oz"))
#'
#'   # Example 4: Different head size
#'   plot_electrode_3d_spherical(head_radius = 90)
#'
#'   # Example 5: Another layout, here the 335 electrodes of the 10-05 system
#'   plot_electrode_3d_spherical(layout = "standard_1005", show_labels = FALSE)
#'
#'   # Example 6: Some channels of a layout, given by their old names
#'   plot_electrode_3d_spherical(channels = c("Cz", "T3", "T4"), layout = "standard_1020")
#' }
#'
#' @keywords internal
#' @importFrom plotly plot_ly add_trace layout
#' @export
#'
plot_electrode_3d_spherical <- function(channels = NULL,
                                        highlight_channels = NULL,
                                        show_labels = TRUE,
                                        color_by = "region",
                                        show_head = TRUE,
                                        head_radius = 87.54,
                                        layout = NULL) {
  
  # ========================================================================
  # LOAD REQUIRED PACKAGE
  # ========================================================================
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive 3D plots.\n",
         "Install it with: install.packages('plotly')")
  }
  
  # ========================================================================
  # THE ELECTRODES TO DRAW
  # ========================================================================
  # The electrodes and their angles come from a layout. With no layout they are
  # the 64 BioSemi electrodes: the "biosemi64" layout, the one montages use, so
  # these numbers are kept in one place (.biosemi64_table(), in
  # R/montage_biosemi.R).
  
  layout_table <- .plot_layout_table(if (is.null(layout)) "biosemi64" else layout)
  layout_label <- if (is.character(layout)) layout else "Custom montage"
  electrode_coords_spherical <- layout_table[, c("electrode", "theta", "phi")]
  
  # ========================================================================
  # CONVERT SPHERICAL TO CARTESIAN COORDINATES
  # ========================================================================
  
  # Convert degrees to radians
  deg_to_rad <- function(deg) {
    return(deg * pi / 180)
  }
  
  # Spherical to Cartesian conversion
  # x = r * sin(theta) * cos(phi)
  # y = r * sin(theta) * sin(phi)
  # z = r * cos(theta)
  
  electrode_coords_spherical$theta_rad <- deg_to_rad(electrode_coords_spherical$theta)
  electrode_coords_spherical$phi_rad <- deg_to_rad(electrode_coords_spherical$phi)
  
  electrode_coords_spherical$x <- head_radius * sin(electrode_coords_spherical$theta_rad) * 
    cos(electrode_coords_spherical$phi_rad)
  electrode_coords_spherical$y <- head_radius * sin(electrode_coords_spherical$theta_rad) * 
    sin(electrode_coords_spherical$phi_rad)
  electrode_coords_spherical$z <- head_radius * cos(electrode_coords_spherical$theta_rad)
  
  # Round to match BioSemi precision (a layout is left unrounded)
  if (is.null(layout)) {
    electrode_coords_spherical$x <- round(electrode_coords_spherical$x, 0)
    electrode_coords_spherical$y <- round(electrode_coords_spherical$y, 0)
    electrode_coords_spherical$z <- round(electrode_coords_spherical$z, 0)
  }
  
  # ========================================================================
  # DETERMINE WHICH CHANNELS TO PLOT
  # ========================================================================
  
  if (is.null(channels)) {
    channels_to_plot <- electrode_coords_spherical$electrode
  } else {
    channels_to_plot <- channels
  }
  
  # Map BioSemi names to standard names if needed
  plot_data <- data.frame(
    original_name = channels_to_plot,
    standard_name = channels_to_plot,
    stringsAsFactors = FALSE
  )
  
  # The layout's own names and second names, found in any letter case: for the
  # 64 BioSemi electrodes the socket names A1 to B32, for the other layouts
  # also old names such as T3.
  plot_data$standard_name <- .match_layout_names(channels_to_plot, layout_table)
  
  # Merge with coordinates
  plot_data <- merge(plot_data, electrode_coords_spherical, 
                     by.x = "standard_name", by.y = "electrode",
                     all.x = TRUE)
  
  # Remove electrodes without coordinates
  missing_coords <- plot_data$original_name[is.na(plot_data$x)]
  if (length(missing_coords) > 0) {
    warning("No coordinates found for: ", paste(missing_coords, collapse = ", "))
    plot_data <- plot_data[!is.na(plot_data$x), ]
  }
  
  if (nrow(plot_data) == 0) {
    stop("No valid electrodes to plot.")
  }
  
  # ========================================================================
  # ASSIGN COLORS BY BRAIN REGION
  # ========================================================================
  
  if (color_by == "region") {
    plot_data$color_group <- "Other"
    plot_data$color_group[grepl("^Fp|^AF|^F[^CT]", plot_data$standard_name)] <- "Frontal"
    plot_data$color_group[grepl("^C|^FC", plot_data$standard_name)] <- "Central"
    plot_data$color_group[grepl("^P[^O]|^CP", plot_data$standard_name)] <- "Parietal"
    plot_data$color_group[grepl("^O|^PO|^I", plot_data$standard_name)] <- "Occipital"
    plot_data$color_group[grepl("^T[0-9]|^FT|^TP|^TT", plot_data$standard_name)] <- "Temporal"
    
    color_palette <- c(
      "Frontal" = "#e41a1c",
      "Central" = "#377eb8",
      "Parietal" = "#4daf4a",
      "Occipital" = "#984ea3",
      "Temporal" = "#ff7f00",
      "Other" = "#999999"
    )
    
  } else {
    plot_data$color_group <- "Electrode"
    color_palette <- c("Electrode" = "#377eb8")
  }
  
  # Highlight specific channels if requested
  if (!is.null(highlight_channels)) {
    plot_data$is_highlighted <- plot_data$original_name %in% highlight_channels |
      plot_data$standard_name %in% highlight_channels
  } else {
    plot_data$is_highlighted <- FALSE
  }
  
  # ========================================================================
  # CREATE PERFECT SPHERE HEAD OUTLINE (from spherical coordinates)
  # ========================================================================
  
  # Create perfect sphere wireframe
  n_points <- 30
  u <- seq(0, 2*pi, length.out = n_points)
  v <- seq(0, pi, length.out = n_points)
  
  # Create wireframe circles (longitude and latitude lines)
  head_lines <- list()
  
  # Longitude lines (vertical circles from front to back)
  for (i in seq(1, length(u), by = 3)) {
    x_line <- head_radius * sin(v) * cos(u[i])
    y_line <- head_radius * sin(v) * sin(u[i])
    z_line <- head_radius * cos(v)
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # Latitude lines (horizontal circles)
  for (j in seq(1, length(v), by = 3)) {
    x_line <- head_radius * sin(v[j]) * cos(u)
    y_line <- head_radius * sin(v[j]) * sin(u)
    z_line <- rep(head_radius * cos(v[j]), length(u))
    head_lines[[length(head_lines) + 1]] <- list(x = x_line, y = y_line, z = z_line)
  }
  
  # ========================================================================
  # CREATE NOSE INDICATOR
  # ========================================================================
  
  nose_length <- head_radius * 0.15  # 15% of radius
  nose_x <- c(0, head_radius * 0.09, -head_radius * 0.09, 0, 0)
  nose_y <- c(head_radius + nose_length, head_radius, head_radius, head_radius + nose_length, head_radius + nose_length)
  nose_z <- c(head_radius * 0.23, head_radius * 0.14, head_radius * 0.14, head_radius * 0.23, head_radius * 0.23)
  
  # ========================================================================
  # BUILD PLOTLY 3D FIGURE
  # ========================================================================
  
  p <- plotly::plot_ly()
  
  # Add head wireframe (if requested)
  if (show_head) {
    for (line in head_lines) {
      p <- plotly::add_trace(
        p,
        x = line$x, y = line$y, z = line$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = "lightgray", width = 1),
        hoverinfo = "none",
        showlegend = FALSE
      )
    }
  }
  
  # Add nose indicator
  p <- plotly::add_trace(
    p,
    x = nose_x, y = nose_y, z = nose_z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 5),
    hoverinfo = "none",
    showlegend = FALSE,
    name = "Nose"
  )
  
  # Add electrodes (by color group)
  for (group in unique(plot_data$color_group)) {
    group_data <- plot_data[plot_data$color_group == group, ]
    
    marker_sizes <- ifelse(group_data$is_highlighted, 14, 10)
    marker_borders <- ifelse(group_data$is_highlighted, "black", "white")
    border_widths <- ifelse(group_data$is_highlighted, 3, 2)
    
    p <- plotly::add_trace(
      p,
      data = group_data,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = ifelse(show_labels, "markers+text", "markers"),
      marker = list(
        size = marker_sizes,
        color = color_palette[group],
        line = list(color = marker_borders, width = border_widths)
      ),
      # The label settings are only given when labels are wanted: plotly puts
      # the text mode back whenever text or textfont are present.
      text = if (show_labels) ~original_name,
      textposition = if (show_labels) "top center",
      textfont = if (show_labels) list(size = 10, color = "black", family = "Arial Black"),
      hovertext = ~paste0("<b>", original_name, "</b><br>",
                          "Standard: ", standard_name, "<br>",
                          "Region: ", color_group, "<br>",
                          "Spherical (degrees):<br>",
                          "  theta (inclination): ", round(theta, 1), "deg<br>",
                          "  phi (azimuth): ", round(phi, 1), "deg<br>",
                          "Cartesian (mm):<br>",
                          "  X: ", round(x, 1), " (L/R)<br>",
                          "  Y: ", round(y, 1), " (Back/Front)<br>",
                          "  Z: ", round(z, 1), " (Bottom/Top)"),
      hoverinfo = "text",
      name = group,
      showlegend = TRUE
    )
  }
  
  # ========================================================================
  # CONFIGURE LAYOUT
  # ========================================================================
  
  axis_range <- head_radius * 1.3
  
  p <- plotly::layout(
    p,
    title = list(
      text = if (is.null(layout)) {
        "BioSemi 64-Channel Electrode Layout (3D - Spherical)"
      } else {
        paste0(layout_label, " Electrode Layout (3D - Spherical)")
      },
      font = list(size = 20, family = "Arial", color = "black")
    ),
    scene = list(
      xaxis = list(
        title = "Left | Right",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range, axis_range)
      ),
      yaxis = list(
        title = "Back | Front (Nose)",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range, axis_range)
      ),
      zaxis = list(
        title = "Bottom | Top",
        showgrid = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        range = c(-axis_range * 0.7, axis_range)
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      ),
      aspectmode = "cube",
      bgcolor = "white"
    ),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    legend = list(
      x = 0.02,
      y = 0.98,
      bgcolor = "rgba(255, 255, 255, 0.9)",
      bordercolor = "black",
      borderwidth = 1,
      font = list(size = 12)
    )
  )
  
  return(p)
}
