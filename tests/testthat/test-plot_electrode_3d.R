# ============================================================================
#                 Test File for plot_electrode_3d() and
#                         plot_electrode_3d_spherical()
# ============================================================================
#
# This test file covers the two interactive 3D electrode plots in
# channel_info2.R. They draw the 64 BioSemi electrodes by default, and any
# built-in layout or montage object when the argument `layout` is given.
#
# Functions tested:
#   1. plot_electrode_3d()           - 3D plot from x/y/z positions
#   2. plot_electrode_3d_spherical() - 3D plot from theta/phi angles
#
# The plots are plotly objects, so the tests build them with
# plotly::plotly_build() and read the electrodes back out of the traces.
#
# Author: Christos Dalamarinis
# Date: October 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# Helper: the electrodes of a plot as a table (name as typed, name in the
# layout, region group, x, y, z, marker size, sorted by name). The names are
# read from the hover text, which every plot has whether or not it has labels.
plotted_electrodes <- function(p) {
  built  <- plotly::plotly_build(p)
  traces <- Filter(function(t) grepl("markers", t$mode), built$x$data)
  out <- do.call(rbind, lapply(traces, function(t) {
    data.frame(
      name     = sub("^<b>(.*?)</b>.*$", "\\1", t$hovertext, perl = TRUE),
      standard = sub("^.*Standard: ([^<]*)<br>.*$", "\\1", t$hovertext),
      group    = t$name,
      x = t$x, y = t$y, z = t$z,
      size = rep_len(t$marker$size, length(t$x)),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$name), ]
  rownames(out) <- NULL
  out
}

plot_title <- function(p) plotly::plotly_build(p)$x$layout$title$text

# Helper: the position of one electrode in a plotted table.
xyz_of <- function(plotted, name) {
  unlist(plotted[plotted$name == name, c("x", "y", "z")], use.names = FALSE)
}

skip_if_not_installed("plotly")


# ============================================================================
# 1. WITHOUT A LAYOUT: THE 64 BIOSEMI ELECTRODES, AS BEFORE
# ============================================================================

test_that("1.1 plot_electrode_3d() with no layout draws the 64 BioSemi electrodes", {
  p <- plot_electrode_3d()
  d <- plotted_electrodes(p)

  expect_equal(plot_title(p), "BioSemi 64-Channel Electrode Layout (3D)")
  expect_equal(nrow(d), 64)
  expect_equal(c(sum(d$x), sum(d$y), sum(d$z)), c(0, -173, 1825))
  expect_equal(xyz_of(d, "Cz"), c(0, 0, 88))
  expect_equal(xyz_of(d, "Iz"), c(0, -79, -37))
  expect_equal(as.vector(table(d$group)[c("Frontal", "Central", "Parietal",
                                          "Occipital", "Temporal")]),
               c(17, 14, 18, 9, 6))
})

test_that("1.2 plot_electrode_3d_spherical() with no layout draws the same 64", {
  p <- plot_electrode_3d_spherical()
  d <- plotted_electrodes(p)

  expect_equal(plot_title(p), "BioSemi 64-Channel Electrode Layout (3D - Spherical)")
  expect_equal(nrow(d), 64)
  expect_equal(c(sum(d$x), sum(d$y), sum(d$z)), c(0, -173, 1825))
  expect_equal(xyz_of(d, "Cz"), c(0, 0, 88))
  # Rounded to whole millimetres, like the BioSemi sheet
  expect_equal(c(d$x, d$y, d$z), round(c(d$x, d$y, d$z)))
})

test_that("1.3 BioSemi socket names still work with no layout", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    d <- plotted_electrodes(f(channels = c("A1", "B32", "Cz")))
    expect_equal(d$name, c("A1", "B32", "Cz"))
    expect_equal(d$standard, c("Fp1", "O2", "Cz"))
    expect_equal(xyz_of(d, "A1"), c(-27, 83, -3))
  }
})

test_that("1.4 an unknown channel with no layout gives the same single warning", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    expect_warning(p <- f(channels = c("Cz", "Banana")),
                   "No coordinates found for: Banana")
    expect_equal(plotted_electrodes(p)$name, "Cz")
  }
})


# ============================================================================
# 2. A LAYOUT GIVEN BY NAME
# ============================================================================

test_that("2.1 plot_electrode_3d() draws every electrode of standard_1005 at its position", {
  p <- plot_electrode_3d(layout = "standard_1005", show_labels = FALSE)
  d <- plotted_electrodes(p)
  tab <- eeganalysis:::.montage_template_table("standard_1005")

  expect_equal(plot_title(p), "standard_1005 Electrode Layout (3D)")
  expect_equal(nrow(d), 335)
  expect_setequal(d$name, tab$name)

  at <- match(d$name, tab$name)
  expect_identical(d$x, tab$x[at])
  expect_identical(d$y, tab$y[at])
  expect_identical(d$z, tab$z[at])
})

test_that("2.2 plot_electrode_3d_spherical() draws standard_1005 from its angles, unrounded", {
  p <- plot_electrode_3d_spherical(layout = "standard_1005", show_labels = FALSE)
  d <- plotted_electrodes(p)
  tab <- eeganalysis:::.montage_template_table("standard_1005")

  expect_equal(plot_title(p), "standard_1005 Electrode Layout (3D - Spherical)")
  expect_equal(nrow(d), 335)

  at <- match(d$name, tab$name)
  expect_equal(d$x, tab$x[at], tolerance = 1e-5)
  expect_equal(d$y, tab$y[at], tolerance = 1e-5)
  expect_equal(d$z, tab$z[at], tolerance = 1e-5)
  # Not rounded to whole millimetres as the BioSemi 64 default is
  expect_true(any(abs(d$x - round(d$x)) > 0.1))
})

test_that("2.3 every built-in layout is drawn with all its electrodes", {
  for (nm in eeganalysis:::.montage_templates()) {
    n <- nrow(eeganalysis:::.montage_template_table(nm))
    expect_equal(nrow(plotted_electrodes(plot_electrode_3d(layout = nm, show_labels = FALSE))), n,
                 info = nm)
    expect_equal(nrow(plotted_electrodes(plot_electrode_3d_spherical(layout = nm, show_labels = FALSE))), n,
                 info = nm)
  }
})

test_that("2.4 the layout 'biosemi64' gives the same plot as no layout", {
  a <- plotted_electrodes(plot_electrode_3d())
  b <- plotted_electrodes(plot_electrode_3d(layout = "biosemi64"))
  expect_identical(a[, c("name", "x", "y", "z")], b[, c("name", "x", "y", "z")])

  # The spherical plot rounds the default to whole mm, a layout is not rounded:
  # the two agree to within BioSemi's own rounding (about half a millimetre)
  a <- plotted_electrodes(plot_electrode_3d_spherical())
  b <- plotted_electrodes(plot_electrode_3d_spherical(layout = "biosemi64"))
  expect_identical(a$name, b$name)
  expect_lt(max(abs(as.matrix(a[, c("x", "y", "z")]) - as.matrix(b[, c("x", "y", "z")]))), 0.5)
})

test_that("2.5 head_radius sets the sphere the layout is drawn on", {
  d <- plotted_electrodes(plot_electrode_3d_spherical(layout = "standard_1020",
                                                      head_radius = 100))
  expect_equal(sqrt(d$x^2 + d$y^2 + d$z^2), rep(100, 19), tolerance = 1e-9)

  d <- plotted_electrodes(plot_electrode_3d_spherical(layout = "standard_1020"))
  expect_equal(sqrt(d$x^2 + d$y^2 + d$z^2), rep(87.54, 19), tolerance = 1e-9)
})


# ============================================================================
# 3. CHANNEL NAMES WITH A LAYOUT
# ============================================================================

test_that("3.1 old names are found by their second name and keep the name typed", {
  old <- c("T3", "T4", "T5", "T6")
  new <- c("T7", "T8", "P7", "P8")
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    d   <- plotted_electrodes(f(channels = old, layout = "standard_1020"))
    ref <- plotted_electrodes(f(channels = new, layout = "standard_1020"))
    expect_setequal(d$name, old)
    expect_equal(d$standard[match(old, d$name)], new)
    for (i in 1:4) {
      expect_equal(xyz_of(d, old[i]), xyz_of(ref, new[i]))
    }
  }
})

test_that("3.2 names are found in any letter case", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    typed <- c("cz", "FZ", "pZ")
    d <- plotted_electrodes(f(channels = typed, layout = "standard_1010"))
    expect_setequal(d$name, typed)
    expect_equal(d$standard[match(typed, d$name)], c("Cz", "Fz", "Pz"))
  }
})

test_that("3.3 BioSemi sockets are found with a BioSemi layout", {
  d <- plotted_electrodes(plot_electrode_3d(channels = c("A1", "B32"), layout = "biosemi64"))
  expect_equal(d$standard, c("Fp1", "O2"))
  expect_equal(xyz_of(d, "A1"), c(-27, 83, -3))

  expect_no_warning(
    d <- plotted_electrodes(plot_electrode_3d(channels = c("A1", "A32"), layout = "biosemi32"))
  )
  expect_equal(d$name, c("A1", "A32"))
  expect_false(anyNA(d$x))
})

test_that("3.4 the BioSemi 64 socket table is not used with another layout", {
  # A1 is Fp1 of the default plot; the 10-20 layout has no electrode A1
  expect_warning(p <- plot_electrode_3d(channels = c("A1", "Cz"), layout = "standard_1020"),
                 "No coordinates found for: A1")
  expect_equal(plotted_electrodes(p)$name, "Cz")
})

test_that("3.5 a channel the layout does not know gives one warning, the rest is drawn", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    expect_warning(p <- f(channels = c("Cz", "Banana", "T3"), layout = "standard_1020"),
                   "No coordinates found for: Banana")
    expect_equal(plotted_electrodes(p)$name, c("Cz", "T3"))
  }
})

test_that("3.6 no channel known is an error", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    expect_error(suppressWarnings(f(channels = c("Banana", "Apple"), layout = "standard_1020")),
                 "No valid electrodes to plot")
  }
})

test_that("3.7 a layout and highlight_channels highlight by the layout's name too", {
  d <- plotted_electrodes(plot_electrode_3d(channels = c("T3", "Cz", "Fz"),
                                            layout = "standard_1020",
                                            highlight_channels = "T7"))
  expect_equal(d$size[d$name == "T3"], 14)
  expect_equal(d$size[d$name != "T3"], c(10, 10))
})


# ============================================================================
# 4. A MONTAGE OBJECT AS THE LAYOUT
# ============================================================================

test_that("4.1 a montage object is drawn at its own positions", {
  mon <- create_montage(c("Cz", "Fz", "Pz"), template = "standard_1020")
  titles <- c("Custom montage Electrode Layout (3D)",
              "Custom montage Electrode Layout (3D - Spherical)")
  plots  <- list(plot_electrode_3d, plot_electrode_3d_spherical)
  for (i in 1:2) {
    p <- plots[[i]](layout = mon)
    d <- plotted_electrodes(p)
    expect_equal(plot_title(p), titles[i])
    expect_equal(d$name, c("Cz", "Fz", "Pz"))
    at <- match(d$name, mon$positions$channel)
    expect_equal(d$x, mon$positions$x[at], tolerance = 1e-5)
    expect_equal(d$y, mon$positions$y[at], tolerance = 1e-5)
    expect_equal(d$z, mon$positions$z[at], tolerance = 1e-5)
  }
})

test_that("4.2 channels pick from a montage object, in any letter case", {
  mon <- create_montage(c("Cz", "Fz", "Pz", "Fp1"), template = "standard_1020")
  typed <- c("cz", "FP1")
  d <- plotted_electrodes(plot_electrode_3d(channels = typed, layout = mon))
  expect_setequal(d$name, typed)
  expect_equal(d$standard[match(typed, d$name)], c("Cz", "Fp1"))
})

test_that("4.3 a montage read from a file is drawn", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  writeLines(c("name,x,y,z", "Cz,0,0,9", "Fz,0,6,7", "Pz,0,-6,7", "T7,-8,0,0"), f)
  mon <- read_custom_montage(f, unit = "cm", recenter = FALSE)

  d <- plotted_electrodes(plot_electrode_3d(layout = mon))
  expect_equal(d$name, c("Cz", "Fz", "Pz", "T7"))
  expect_equal(xyz_of(d, "Cz"), c(0, 0, 90))
  expect_equal(xyz_of(d, "T7"), c(-80, 0, 0))

  d <- plotted_electrodes(plot_electrode_3d_spherical(layout = mon, head_radius = 90))
  expect_equal(xyz_of(d, "Cz"), c(0, 0, 90), tolerance = 1e-9)
})


# ============================================================================
# 5. BAD LAYOUTS
# ============================================================================

test_that("5.1 a layout that is not a layout name or a montage is an error", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    expect_error(f(layout = "nope"), "'layout' must be the name of a built-in layout")
    expect_error(f(layout = c("standard_1020", "standard_1010")), "'layout' must be")
    expect_error(f(layout = NA_character_), "'layout' must be")
    expect_error(f(layout = 42), "'layout' must be")
    expect_error(f(layout = list()), "'layout' must be")
  }
  # The message lists the layouts to choose from
  expect_error(plot_electrode_3d(layout = "nope"), "standard_1005")
  expect_error(plot_electrode_3d(layout = "nope"), "biosemi128")
})


# ============================================================================
# 6. LABELS
# ============================================================================

test_that("6.1 show_labels = FALSE draws markers only, and no warning from plotly", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    p <- f(show_labels = FALSE)
    expect_no_warning(built <- plotly::plotly_build(p))
    modes <- vapply(built$x$data, function(t) t$mode, "")
    expect_false(any(grepl("text", modes)))
    # The hover text is still there
    expect_equal(nrow(plotted_electrodes(p)), 64)
  }
})

test_that("6.2 show_labels = TRUE draws the names next to the markers", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    built  <- plotly::plotly_build(f(channels = c("Cz", "Fz"), layout = "standard_1020"))
    traces <- Filter(function(t) grepl("markers", t$mode), built$x$data)
    expect_true(all(vapply(traces, function(t) grepl("text", t$mode), TRUE)))
    expect_setequal(unlist(lapply(traces, function(t) t$text)), c("Cz", "Fz"))
  }
})

test_that("6.3 labels can be switched off for a layout too", {
  built <- plotly::plotly_build(plot_electrode_3d(layout = "standard_1005", show_labels = FALSE))
  traces <- Filter(function(t) grepl("markers", t$mode), built$x$data)
  expect_true(all(vapply(traces, function(t) identical(t$mode, "markers"), TRUE)))
  expect_true(all(vapply(traces, function(t) is.null(t$text), TRUE)))
})


# ============================================================================
# 7. REGION COLOURS FOR THE LARGER LAYOUTS
# ============================================================================

test_that("7.1 the TTP and inion-side electrodes of 10-05 get a region", {
  for (f in list(plot_electrode_3d, plot_electrode_3d_spherical)) {
    d <- plotted_electrodes(f(layout = "standard_1005", show_labels = FALSE))
    expect_true(all(d$group[grepl("^TTP", d$name)] == "Temporal"))
    expect_true(all(d$group[d$name %in% c("I1", "I2", "I1h", "I2h", "Iz")] == "Occipital"))
    expect_false("Other" %in% d$group)
  }
})

test_that("7.2 color_by = 'single' draws one group for a layout", {
  d <- plotted_electrodes(plot_electrode_3d(layout = "standard_1010", color_by = "single",
                                            show_labels = FALSE))
  expect_equal(unique(d$group), "Electrode")
  expect_equal(nrow(d), 72)
})


# ============================================================================
# SUMMARY OF WHAT THIS FILE COVERS
# ============================================================================
#
# - No layout: the 64 BioSemi electrodes, titles, positions, region counts,
#   rounding of the spherical plot, BioSemi socket names and the warning for
#   an unknown channel are all as before
# - A layout by name: all 335 electrodes of standard_1005 at the positions of
#   the layout (exactly for the Cartesian plot, from the angles for the
#   spherical one, unrounded), every built-in layout drawn with all its
#   electrodes, 'biosemi64' the same as no layout, head_radius sets the
#   sphere
# - Channel names with a layout: old names T3 to T6 and BioSemi sockets found
#   by their second name, any letter case, the BioSemi 64 socket table not
#   used with another layout, one warning for an unknown channel, an error
#   when none is known, highlight_channels by the layout's name
# - A montage object as the layout: its own positions, picking channels, a
#   montage read from a file in centimetres
# - Bad layouts: five kinds of wrong input give the same clear error that
#   lists the layouts
# - Labels: show_labels = FALSE really draws markers only (plotly used to put
#   the text back), show_labels = TRUE draws the names
# - Region colours: TTP and I1/I2 electrodes of 10-05 get a region
# ============================================================================
