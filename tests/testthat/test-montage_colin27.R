# ============================================================================
#                  Test File for montage_colin27.R
# ============================================================================
#
# This test file covers the two real-head layouts "colin27_1020" (94
# electrodes) and "colin27_1005" (343 electrodes): positions measured on a real
# head, taken unchanged from the electrode tables of the FieldTrip toolbox
# (inst/extdata/fieldtrip/standard_1020.elc and standard_1005.elc).
#
# Parts tested:
#   1. The two files in the package    - shipped, unchanged, well formed
#   2. The layout names                - listed, older layouts untouched
#   3. The tables                      - names, positions, centre, rows kept
#   4. The functions that use layouts  - create_montage(), set_montage(),
#                                        plot_electrode_3d(), interpolate_bads(),
#                                        compute_csd()
#
# Author: Christos Dalamarinis
# Date: November 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# The two layouts and their files
colin27_files <- c(colin27_1020 = "standard_1020.elc",
                   colin27_1005 = "standard_1005.elc")
colin27_rows  <- c(colin27_1020 = 94, colin27_1005 = 343)

colin27_file <- function(name) {
  system.file("extdata", "fieldtrip", name, package = "eeganalysis")
}

# The checksum of a file with every carriage return taken out. On Windows Git
# may turn the line endings of a text file into CR LF when it is checked out;
# the numbers in the file do not change.
md5_without_cr <- function(path) {
  bytes <- readBin(path, "raw", n = file.size(path))
  tmp <- tempfile(fileext = ".elc")
  on.exit(unlink(tmp))
  writeBin(bytes[bytes != as.raw(13)], tmp)
  unname(tools::md5sum(tmp))
}

# Reads an .elc file without the package's own reader: the lines, the numbers
# of the Positions block (landmarks included) and the names of the Labels block.
read_elc_plain <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  i_pos <- which(lines == "Positions")
  i_lab <- which(lines == "Labels")
  rows  <- strsplit(lines[(i_pos + 1):(i_lab - 1)], "[[:space:]]+")
  names <- lines[(i_lab + 1):length(lines)]
  list(lines = lines,
       xyz   = do.call(rbind, lapply(rows, as.numeric)),
       name  = names[nzchar(names)])
}

# The table of a layout without its landmark rows, as the file gives it
raw_electrodes <- function(layout) {
  e <- read_elc_plain(colin27_file(colin27_files[[layout]]))
  xyz <- e$xyz[-(1:3), , drop = FALSE]
  rownames(xyz) <- e$name[-(1:3)]
  xyz
}

colin27_table <- function(layout) eeganalysis:::.montage_template_table(layout)

# A small recording whose channels are in the 94-electrode layout
recording_in_colin27 <- function(chans, n = 60) {
  set.seed(7)
  new_eeg(data = matrix(rnorm(length(chans) * n), nrow = length(chans)),
          channels = chans, sampling_rate = 100)
}

# ============================================================================
#      TEST SUITE 1: The Two Files in the Package
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Shipped
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies both electrode files are found in the package.
test_that("1.1 the two electrode files are shipped with the package", {
  for (f in unname(colin27_files)) {
    expect_true(nzchar(colin27_file(f)), info = f)
    expect_true(file.exists(colin27_file(f)), info = f)
  }
})

# ----------------------------------------------------------------------------
# Test 1.2: Unchanged copies
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the files are exact copies of the originals
# (checksum of the original files, line endings aside). The files are included
# unchanged, so a change in them is a mistake.
test_that("1.2 the files are unchanged copies of the originals", {
  expect_equal(md5_without_cr(colin27_file("standard_1020.elc")),
               "cc6837af23f426d7d0c1d250a7f9fe52")
  expect_equal(md5_without_cr(colin27_file("standard_1005.elc")),
               "c76c30b155e952029f64c83f8836bf05")
})

# ----------------------------------------------------------------------------
# Test 1.3: Well formed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies each file says it is in millimetres, has as many
# position rows and names as its header says, and starts with the three
# landmarks.
test_that("1.3 the files hold what their header says", {
  for (f in unname(colin27_files)) {
    e <- read_elc_plain(colin27_file(f))
    n <- as.integer(sub("^NumberPositions=[[:space:]]*", "",
                        e$lines[grepl("^NumberPositions=", e$lines)]))
    expect_true(any(grepl("^UnitPosition[[:space:]]+mm$", e$lines)), info = f)
    expect_equal(nrow(e$xyz), n, info = f)
    expect_equal(length(e$name), n, info = f)
    expect_equal(e$name[1:3], c("LPA", "RPA", "Nz"), info = f)
    expect_false(anyNA(e$xyz), info = f)
  }
})

# ----------------------------------------------------------------------------
# Test 1.4: The two files agree
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an electrode that is in both files has the same
# numbers in both. This is what lets both layouts share one shift (test 3.4).
test_that("1.4 an electrode in both files has the same numbers in both", {
  a <- raw_electrodes("colin27_1020")
  b <- raw_electrodes("colin27_1005")
  common <- intersect(rownames(a), rownames(b))
  expect_length(common, 92)
  expect_identical(a[common, ], b[common, ])
  expect_setequal(setdiff(rownames(a), rownames(b)), c("O9", "O10"))
})

# ============================================================================
#      TEST SUITE 2: The Layout Names
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Listed after the older layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the two names come last in the list of layouts.
test_that("2.1 the layout names are listed after the older layouts", {
  expect_equal(eeganalysis:::.colin27_templates(), c("colin27_1020", "colin27_1005"))
  expect_equal(eeganalysis:::.montage_templates(),
               c("biosemi16", "biosemi32", "biosemi64", "biosemi128",
                 "biosemi160", "biosemi256",
                 "standard_1020", "standard_1010", "standard_1005",
                 "colin27_1020", "colin27_1005"))
})

# ----------------------------------------------------------------------------
# Test 2.2: Older layouts untouched
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the layouts that were there before are as they
# were: their sizes, and one head radius for all electrodes of the standard
# ones.
test_that("2.2 the older layouts are as they were", {
  expect_equal(nrow(colin27_table("biosemi64")), 64)
  expect_equal(nrow(colin27_table("biosemi256")), 256)
  for (tp in c("standard_1020", "standard_1010", "standard_1005")) {
    tab <- colin27_table(tp)
    expect_equal(nrow(tab), c(standard_1020 = 19, standard_1010 = 72,
                              standard_1005 = 335)[[tp]], info = tp)
    expect_true(all(tab$radius == 87.54), info = tp)
  }
})

# ----------------------------------------------------------------------------
# Test 2.3: Unknown names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a name that is not a real-head layout is an error,
# also through create_montage().
test_that("2.3 a name that is not a real-head layout is an error", {
  expect_error(eeganalysis:::.colin27_template_table("colin27_1010"),
               "Unknown real-head layout 'colin27_1010'")
  expect_error(eeganalysis:::.colin27_template_table(c("colin27_1020", "colin27_1005")),
               "Unknown real-head layout")
  expect_error(eeganalysis:::.colin27_template_table(NA_character_),
               "Unknown real-head layout")
  expect_error(create_montage(template = "colin27_1010"),
               "Unknown montage template.*colin27_1020, colin27_1005")
})

# ============================================================================
#      TEST SUITE 3: The Tables
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Every electrode of the file, no landmarks
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies each table has the columns of the other layouts,
# one row for every electrode of its file in the same order (so 94 and 343
# rows, nothing folded into an alias), no landmark rows and no repeated name.
test_that("3.1 each table has every electrode of its file and no landmarks", {
  for (nm in names(colin27_files)) {
    tab <- colin27_table(nm)
    e   <- read_elc_plain(colin27_file(colin27_files[[nm]]))
    expect_equal(names(tab),
                 c("name", "alias", "x", "y", "z", "theta", "phi", "radius"),
                 info = nm)
    expect_equal(nrow(tab), colin27_rows[[nm]], info = nm)
    expect_equal(tab$name, e$name[-(1:3)], info = nm)
    expect_false(any(c("LPA", "RPA", "Nz") %in% tab$name), info = nm)
    expect_false(anyDuplicated(tolower(tab$name)) > 0, info = nm)
    expect_true(all(is.na(tab$alias)), info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 3.2: The numbers of the file, shifted
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the positions are the numbers of the file with one
# and the same shift taken off (only the rounding to six decimals is left), and
# that the shift is the centre of the sphere that fits the 94 electrodes of the
# first file, for both layouts.
test_that("3.2 the positions are the numbers of the file, shifted by one amount", {
  fit <- eeganalysis:::.montage_sphere_fit(
    raw_electrodes("colin27_1020")[, 1], raw_electrodes("colin27_1020")[, 2],
    raw_electrodes("colin27_1020")[, 3])

  for (nm in names(colin27_files)) {
    mine  <- as.matrix(colin27_table(nm)[, c("x", "y", "z")])
    raw   <- unname(raw_electrodes(nm))
    shift <- colMeans(mine - raw)
    expect_lt(max(abs(sweep(mine - raw, 2, shift))), 1e-5)   # the shape is the file's
    expect_equal(unname(shift), -fit$centre, tolerance = 1e-5)
  }
})

# ----------------------------------------------------------------------------
# Test 3.3: Centred, and the size of the head
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the sphere that fits the 94 electrodes best is
# centred on the origin (the other layouts are, and the functions that use
# positions expect it), and has the radius of the original head.
test_that("3.3 the best-fit sphere of colin27_1020 is centred on the origin", {
  tab <- colin27_table("colin27_1020")
  fit <- eeganalysis:::.montage_sphere_fit(tab$x, tab$y, tab$z)
  expect_lt(max(abs(fit$centre)), 1e-4)
  expect_equal(fit$radius, 99.36, tolerance = 1e-3)
})

# ----------------------------------------------------------------------------
# Test 3.4: One shift for both layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an electrode that is in both layouts has exactly
# the same position in both.
test_that("3.4 an electrode in both layouts has the same position in both", {
  a <- colin27_table("colin27_1020")
  b <- colin27_table("colin27_1005")
  common <- intersect(a$name, b$name)
  cols <- c("x", "y", "z", "theta", "phi", "radius")
  expect_length(common, 92)
  expect_identical(unname(as.matrix(a[match(common, a$name), cols])),
                   unname(as.matrix(b[match(common, b$name), cols])))
})

# ----------------------------------------------------------------------------
# Test 3.5: Angles and distance agree with the positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies theta, phi and radius are those of the x, y and z
# of the same row, and that the radius is the distance from the origin.
test_that("3.5 theta, phi and radius belong to the position of the row", {
  for (nm in names(colin27_files)) {
    tab <- colin27_table(nm)
    ang <- eeganalysis:::.montage_angles(tab$x, tab$y, tab$z)
    expect_equal(tab$theta, ang$theta, tolerance = 1e-6, info = nm)
    expect_equal(tab$phi, ang$phi, tolerance = 1e-6, info = nm)
    expect_equal(tab$radius, sqrt(tab$x^2 + tab$y^2 + tab$z^2), tolerance = 1e-6, info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 3.6: Known numbers of the file
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the reading with numbers copied by hand from the
# files: the distance between T7 and T8 does not depend on the shift, so it
# must be the distance between the two rows of the file.
test_that("3.6 the distance between T7 and T8 is the one of the file", {
  t7 <- c(-84.1611, -16.0187, -9.3460)
  t8 <- c(85.0799, -15.0203, -9.4900)
  expected <- sqrt(sum((t8 - t7)^2))
  for (nm in names(colin27_files)) {
    tab <- colin27_table(nm)
    p7 <- unlist(tab[tab$name == "T7", c("x", "y", "z")])
    p8 <- unlist(tab[tab$name == "T8", c("x", "y", "z")])
    expect_equal(sqrt(sum((p8 - p7)^2)), expected, tolerance = 1e-6, info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 3.7: The old names are rows of their own
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies T3, T4, T5 and T6 are rows of their own, at exactly
# the position of T7, T8, P7 and P8, as in the original files, and that the
# mastoid and ear electrodes are rows too.
test_that("3.7 T3 to T6 are rows of their own at the position of T7, T8, P7, P8", {
  old <- c("T3", "T4", "T5", "T6")
  new <- c("T7", "T8", "P7", "P8")
  for (nm in names(colin27_files)) {
    tab <- colin27_table(nm)
    expect_true(all(c(old, new, "M1", "M2", "A1", "A2") %in% tab$name), info = nm)
    cols <- c("x", "y", "z")
    expect_identical(unname(as.matrix(tab[match(old, tab$name), cols])),
                     unname(as.matrix(tab[match(new, tab$name), cols])),
                     info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 3.8: A real head is not round
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the radius column holds the distance of each
# electrode: it differs from electrode to electrode, between about 84 and
# 112 mm, unlike the standard layouts.
test_that("3.8 the distance from the origin differs from electrode to electrode", {
  for (nm in names(colin27_files)) {
    r <- colin27_table(nm)$radius
    expect_gt(min(r), 84)
    expect_lt(max(r), 112)
    expect_gt(diff(range(r)), 20)
  }
})

# ============================================================================
#      TEST SUITE 4: The Functions That Use Layouts
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: create_montage()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() gives all electrodes of the layout
# by default, with the positions of the table and the layout as coordinate
# frame.
test_that("4.1 create_montage() builds a montage of a real-head layout", {
  for (nm in names(colin27_files)) {
    m   <- create_montage(template = nm)
    tab <- colin27_table(nm)
    expect_s3_class(m, "montage")
    expect_equal(m$coord_frame, nm)
    expect_equal(m$channels, tab$name)
    expect_equal(m$positions$x, tab$x)
    expect_equal(m$positions$y, tab$y)
    expect_equal(m$positions$z, tab$z)
    expect_equal(m$positions$radius, tab$radius)
  }
})

# ----------------------------------------------------------------------------
# Test 4.2: Channels by name
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies channels are found by name in any letter case and
# keep the name typed, and that a name the layout does not have is dropped with
# a warning.
test_that("4.2 channels are found by name in any letter case", {
  tab <- colin27_table("colin27_1020")
  m <- create_montage(c("cz", "FZ", "pz"), template = "colin27_1020")
  expect_equal(m$channels, c("cz", "FZ", "pz"))
  expect_equal(m$positions$z, tab$z[match(c("Cz", "Fz", "Pz"), tab$name)])

  expect_warning(m <- create_montage(c("Cz", "Banana"), template = "colin27_1005"),
                 "1 channel\\(s\\).*Banana")
  expect_equal(m$channels, "Cz")
})

# ----------------------------------------------------------------------------
# Test 4.3: Old names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a recording that names its channels T3 and T4 finds
# them, with the same positions as T7 and T8.
test_that("4.3 the old names are found as electrodes of their own", {
  m <- create_montage(c("T3", "T7", "T4", "T8"), template = "colin27_1020")
  expect_equal(m$channels, c("T3", "T7", "T4", "T8"))
  p <- m$positions
  expect_identical(unname(as.matrix(p[1, c("x", "y", "z")])), unname(as.matrix(p[2, c("x", "y", "z")])))
  expect_identical(unname(as.matrix(p[3, c("x", "y", "z")])), unname(as.matrix(p[4, c("x", "y", "z")])))
})

# ----------------------------------------------------------------------------
# Test 4.4: set_montage() by name
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a layout can be attached to a recording by its
# name, and the channels keep the names they have in the recording.
test_that("4.4 set_montage() attaches a real-head layout by name", {
  eeg <- suppressMessages(set_montage(recording_in_colin27(c("Fz", "Cz", "Pz", "T3")),
                                      "colin27_1020"))
  tab <- colin27_table("colin27_1020")
  expect_equal(eeg$montage$coord_frame, "colin27_1020")
  expect_equal(eeg$montage$channels, c("Fz", "Cz", "Pz", "T3"))
  expect_equal(eeg$montage$positions$x, tab$x[match(c("Fz", "Cz", "Pz", "T3"), tab$name)])
})

# ----------------------------------------------------------------------------
# Test 4.5: The 3D plots
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies both 3D electrode plots draw a real-head layout
# (test-plot_electrode_3d.R checks that every layout is drawn with all its
# electrodes).
test_that("4.5 the 3D electrode plots draw a real-head layout", {
  expect_s3_class(plot_electrode_3d(layout = "colin27_1005", show_labels = FALSE), "plotly")
  expect_s3_class(plot_electrode_3d_spherical(layout = "colin27_1020", show_labels = FALSE), "plotly")
})

# ----------------------------------------------------------------------------
# Test 4.6: interpolate_bads()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies interpolate_bads() works with a real-head layout
# without a warning, and that it uses only the directions of the electrodes:
# the same montage with every electrode moved onto a sphere gives the same
# result.
test_that("4.6 interpolate_bads() works with a real-head layout", {
  chans <- c("Fz", "FCz", "Cz", "CPz", "Pz", "F3", "F4", "C3", "C4", "P3", "P4",
             "Fp1", "Fp2", "O1", "O2", "T7", "T8", "F7", "F8", "P7", "P8")
  eeg <- suppressMessages(set_montage(recording_in_colin27(chans), "colin27_1020"))
  eeg$bads <- "C3"

  expect_warning(out <- suppressMessages(interpolate_bads(eeg)), NA)
  expect_true(all(is.finite(out$data)))
  c3 <- match("C3", chans)
  expect_false(isTRUE(all.equal(out$data[c3, ], eeg$data[c3, ])))

  round_eeg <- eeg
  p <- round_eeg$montage$positions
  r <- sqrt(p$x^2 + p$y^2 + p$z^2)
  for (k in c("x", "y", "z")) p[[k]] <- p[[k]] / r * 87.54
  round_eeg$montage$positions <- p
  out_round <- suppressMessages(interpolate_bads(round_eeg))
  expect_equal(out$data, out_round$data, tolerance = 1e-8)
})

# ----------------------------------------------------------------------------
# Test 4.7: compute_csd()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies compute_csd() gives finite values with a real-head
# layout, uses only the directions of the electrodes (the result is the same on
# a sphere when the head radius is given), and says that the positions are not
# close to spherical, because a real head is not round.
test_that("4.7 compute_csd() works with a real-head layout and warns that it is not round", {
  chans <- c("Fz", "FCz", "Cz", "CPz", "Pz", "F3", "F4", "C3", "C4", "P3", "P4",
             "Fp1", "Fp2", "O1", "O2", "T7", "T8", "F7", "F8", "P7", "P8")
  eeg <- suppressMessages(set_montage(recording_in_colin27(chans), "colin27_1020"))

  expect_warning(out <- suppressMessages(compute_csd(eeg, head_radius = 90)),
                 "not close to spherical")
  expect_true(all(is.finite(out$data)))

  round_eeg <- eeg
  p <- round_eeg$montage$positions
  r <- sqrt(p$x^2 + p$y^2 + p$z^2)
  for (k in c("x", "y", "z")) p[[k]] <- p[[k]] / r * 87.54
  round_eeg$montage$positions <- p
  expect_warning(out_round <- suppressMessages(compute_csd(round_eeg, head_radius = 90)), NA)
  expect_equal(out$data, out_round$data, tolerance = 1e-8)
})
