# ============================================================================
#                       Test File for montage.R
# ============================================================================
#
# This test file provides comprehensive testing for the montage.R script,
# which defines the 'montage' S3 class (channel scalp positions) and the
# functions to build one from the built-in BioSemi electrode template and
# attach it to an 'eeg' object.
#
# Functions tested:
#   1. new_montage()     - Low-level montage constructor
#   2. create_montage()  - Builds a montage from the BioSemi 64 template
#   3. set_montage()     - Attaches a montage to an eeg object
#   4. print.montage()   - Print method for montage objects
#
# Author: Christos Dalamarinis
# Date: July 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#          TEST SUITE 1: new_montage() - Low-Level Constructor
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Basic object creation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that new_montage() builds a valid 'montage'
# object from arbitrary channel names and a positions data frame.
test_that("new_montage creates a valid montage object", {
  positions <- data.frame(
    channel = c("A", "B"),
    x = c(1, 2), y = c(3, 4), z = c(5, 6),
    theta = c(10, 20), phi = c(30, 40), radius = c(87.54, 87.54),
    stringsAsFactors = FALSE
  )

  m <- new_montage(channels = c("A", "B"), positions = positions,
                    coord_frame = "custom")

  expect_s3_class(m, "montage")
  expect_equal(m$channels, c("A", "B"))
  expect_equal(m$coord_frame, "custom")
  expect_equal(nrow(m$positions), 2)
})

# ----------------------------------------------------------------------------
# Test 1.2: Default coord_frame
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the default coord_frame is "biosemi64".
test_that("new_montage defaults coord_frame to biosemi64", {
  m <- new_montage(channels = "Cz",
                    positions = data.frame(channel = "Cz", x = 0, y = 0, z = 88,
                                            theta = 0, phi = 0, radius = 87.54))
  expect_equal(m$coord_frame, "biosemi64")
})

# ============================================================================
#     TEST SUITE 2: create_montage() - Build From BioSemi Template
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Default call returns all 64 standard electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that create_montage() with no arguments returns
# a montage with all 64 standard 10-20/10-10 BioSemi electrodes, and that
# external/status channels are excluded from the default set.
test_that("create_montage() with no arguments returns all 64 standard electrodes", {
  m <- create_montage()

  expect_s3_class(m, "montage")
  expect_equal(length(m$channels), 64)
  expect_equal(nrow(m$positions), 64)
  expect_true(all(c("Cz", "Fz", "Pz", "Oz", "Fp1", "Fp2") %in% m$channels))
  expect_false(any(grepl("^EXG|^GSR|^Status$", m$channels, ignore.case = TRUE)))
})

# ----------------------------------------------------------------------------
# Test 2.2: Subset selection, mixed standard/BioSemi naming, case-insensitivity
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that create_montage() accepts a subset of channels
# mixing standard 10-20 names, BioSemi A/B names, and non-standard casing -
# all should resolve to the same underlying electrode positions.
test_that("create_montage() resolves standard names, BioSemi names, and mixed case", {
  m <- create_montage(c("Cz", "cz", "CZ", "B16"))

  expect_equal(nrow(m$positions), 4)

  # Cz and B16 are the same physical electrode (see get_electrode_database())
  cz_rows <- m$positions[m$positions$channel %in% c("Cz", "cz", "CZ", "B16"), ]
  expect_equal(length(unique(cz_rows$x)), 1)
  expect_equal(length(unique(cz_rows$y)), 1)
  expect_equal(length(unique(cz_rows$z)), 1)
  expect_equal(unique(cz_rows$x), 0)
  expect_equal(unique(cz_rows$z), 88)
})

# ----------------------------------------------------------------------------
# Test 2.3: Unmatched channels are dropped with a warning
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that channels not present in the BioSemi template
# (typos, or external channels with no coordinates) are dropped with a
# warning rather than causing an error, as long as at least one channel
# resolves.
test_that("create_montage() warns and drops unmatched/external channels", {
  expect_warning(
    m <- create_montage(c("Cz", "NOT_A_CHANNEL", "EXG1")),
    "NOT_A_CHANNEL"
  )
  expect_equal(m$channels, "Cz")
})

# ----------------------------------------------------------------------------
# Test 2.4: All channels unmatched raises an error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() errors informatively when none
# of the requested channels can be resolved.
test_that("create_montage() errors when no channels resolve", {
  expect_error(
    suppressWarnings(create_montage(c("NOT_REAL_1", "NOT_REAL_2"))),
    "None of the requested channels"
  )
})

# ----------------------------------------------------------------------------
# Test 2.5: Unknown template errors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() rejects unsupported templates.
test_that("create_montage() errors on an unknown template", {
  expect_error(create_montage(template = "not_a_real_template"), "Unknown montage template")
})

# ============================================================================
#          TEST SUITE 3: set_montage() - Attach Montage to eeg Object
# ============================================================================

# Helper: build a small synthetic eeg object with real electrode names plus
# one external channel, so channel_types has both "eeg" and "external".
make_test_eeg <- function() {
  chans <- c("Cz", "Fz", "Pz", "Oz", "T7", "T8", "EXG1")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  new_eeg(data = mat, channels = chans, sampling_rate = 100)
}

# ----------------------------------------------------------------------------
# Test 3.1: Successful attachment
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() attaches a montage restricted to
# the eeg object's actual "eeg"-type channels (excluding the external
# channel), appends to preprocessing_history, and returns an 'eeg' object.
test_that("set_montage() attaches montage scoped to eeg-type channels only", {
  eeg <- make_test_eeg()
  eeg2 <- suppressWarnings(set_montage(eeg, create_montage()))

  expect_s3_class(eeg2, "eeg")
  expect_s3_class(eeg2$montage, "montage")

  # EXG1 is classified "external" and must not receive a scalp position
  expect_false("EXG1" %in% eeg2$montage$channels)
  expect_true(all(c("Cz", "Fz", "Pz", "Oz", "T7", "T8") %in% eeg2$montage$channels))

  # preprocessing_history grew by one entry
  expect_equal(length(eeg2$preprocessing_history),
               length(eeg$preprocessing_history) + 1)
  expect_true(grepl("Montage attached", eeg2$preprocessing_history[[length(eeg2$preprocessing_history)]]))
})

# ----------------------------------------------------------------------------
# Test 3.2: Input validation
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() requires an 'eeg' object and a
# 'montage' object, in the correct argument positions.
test_that("set_montage() validates input classes", {
  eeg <- make_test_eeg()
  m <- create_montage()

  expect_error(set_montage(list(), m), "class 'eeg'")
  expect_error(set_montage(eeg, list()), "class 'montage'")
})

# ----------------------------------------------------------------------------
# Test 3.3: Warnings for mismatched channel sets
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() warns (not errors) when the eeg
# object has EEG channels missing from the montage, or the montage has
# channels not present in the data.
test_that("set_montage() warns on partial channel overlap", {
  chans <- c("Cz", "Fz", "NotInMontage")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100)

  small_montage <- create_montage(c("Cz", "Fz", "Pz"))

  expect_warning(
    expect_warning(
      eeg2 <- set_montage(eeg, small_montage),
      "NotInMontage"
    ),
    "Pz"
  )
  expect_setequal(eeg2$montage$channels, c("Cz", "Fz"))
})

# ----------------------------------------------------------------------------
# Test 3.4: No overlap at all errors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() errors when there is no overlap
# at all between the montage and the eeg object's channels.
test_that("set_montage() errors when there is no channel overlap", {
  chans <- c("NotAnElectrode1", "NotAnElectrode2")
  mat <- matrix(rnorm(length(chans) * 50), nrow = length(chans))
  eeg <- new_eeg(data = mat, channels = chans, sampling_rate = 100)

  m <- create_montage(c("Cz", "Fz"))

  expect_error(suppressWarnings(set_montage(eeg, m)), "No overlap")
})

# ============================================================================
#              TEST SUITE 4: print.montage() - Print Method
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: Print output contains key summary fields
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies print.montage() runs without error and reports
# the coordinate template and channel count.
test_that("print.montage displays a formatted summary", {
  m <- create_montage(c("Cz", "Fz", "Pz"))

  output <- capture.output(print(m))
  expect_true(any(grepl("Montage Object Summary", output)))
  expect_true(any(grepl("biosemi64", output)))
  expect_true(any(grepl("3", output)))
})

# ----------------------------------------------------------------------------
# Test 4.2: Print returns its argument invisibly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies print.montage() follows the standard R print
# method convention of invisibly returning x.
test_that("print.montage returns x invisibly", {
  m <- create_montage(c("Cz", "Fz"))
  ret <- withVisible(print(m))
  expect_false(ret$visible)
  expect_identical(ret$value, m)
})

# ============================================================================
#         TEST SUITE 5: BioSemi Cap Layouts (create_montage(template = ))
# ============================================================================

# Helper: checks one electrode against the Cartesian numbers BioSemi prints in
# its own coordinate sheet (Cap_coords_all.xls). The sheet prints x and y to the
# whole millimetre and z in full; the montage uses a head radius of 87.54 mm
# where the sheet uses 87.5352 mm, so z differs by about 0.005 mm.
expect_biosemi_position <- function(m, channel, x, y, z) {
  p <- m$positions[m$positions$channel == channel, ]
  expect_equal(nrow(p), 1, info = channel)
  expect_lte(abs(p$x - x), 0.51)
  expect_lte(abs(p$y - y), 0.51)
  expect_lte(abs(p$z - z), 0.01)
}

# ----------------------------------------------------------------------------
# Test 5.1: Every layout has the right number of distinct electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies each of the six BioSemi layouts returns one row
# per electrode (16, 32, 64, 128, 160, 256), with no repeated names, finite
# positions, and the layout name stored as the coordinate frame.
test_that("every BioSemi layout has the right number of distinct electrodes", {
  sizes <- c(biosemi16 = 16, biosemi32 = 32, biosemi64 = 64,
             biosemi128 = 128, biosemi160 = 160, biosemi256 = 256)

  for (tmpl in names(sizes)) {
    m <- create_montage(template = tmpl)
    expect_s3_class(m, "montage")
    expect_equal(length(m$channels), sizes[[tmpl]], info = tmpl)
    expect_equal(nrow(m$positions), sizes[[tmpl]], info = tmpl)
    expect_false(anyDuplicated(m$channels) > 0, info = tmpl)
    expect_true(all(is.finite(as.matrix(m$positions[, c("x", "y", "z")]))), info = tmpl)
    expect_equal(m$coord_frame, tmpl)
  }
})

# ----------------------------------------------------------------------------
# Test 5.2: Every electrode lies on the head sphere
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every position is at the head radius (87.54 mm)
# from the centre of the head. The 64-channel table stores whole millimetres
# (as BioSemi prints them), so it is allowed up to 1 mm of rounding.
test_that("every electrode of every layout lies on the 87.54 mm head sphere", {
  for (tmpl in c("biosemi16", "biosemi32", "biosemi64",
                 "biosemi128", "biosemi160", "biosemi256")) {
    p <- create_montage(template = tmpl)$positions
    dist <- sqrt(p$x^2 + p$y^2 + p$z^2)
    expect_true(all(abs(dist - 87.54) < if (tmpl == "biosemi64") 1 else 1e-4),
                info = tmpl)
    expect_true(all(p$radius == 87.54), info = tmpl)
  }
})

# ----------------------------------------------------------------------------
# Test 5.3: Positions match the numbers BioSemi prints
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a few electrodes per layout against the Cartesian
# values printed in BioSemi's own sheet (the expected numbers below were read
# from that sheet, not computed here).
test_that("layout positions match BioSemi's printed coordinates", {
  m16 <- create_montage(template = "biosemi16")
  expect_biosemi_position(m16, "F4", 48, 59, 43.7676)
  expect_biosemi_position(m16, "Pz", 0, -63, 60.8071)
  expect_biosemi_position(m16, "Cz", 0, 0, 87.5352)

  m32 <- create_montage(template = "biosemi32")
  expect_biosemi_position(m32, "Fp1", -27, 83, -3.0549)
  expect_biosemi_position(m32, "Fz", 0, 63, 60.8071)
  expect_biosemi_position(m32, "Cz", 0, 0, 87.5352)

  m128 <- create_montage(template = "biosemi128")
  expect_biosemi_position(m128, "A1", 0, 0, 87.5352)
  expect_biosemi_position(m128, "A2", 0, -17, 85.7779)
  expect_biosemi_position(m128, "D32", -69, -50, -20.4347)

  m256 <- create_montage(template = "biosemi256")
  expect_biosemi_position(m256, "A2", 0, -14, 86.4092)
  expect_biosemi_position(m256, "H32", -31, -70, -43.2373)
})

# ----------------------------------------------------------------------------
# Test 5.4: The 160 and 256 layouts keep the angles BioSemi's cap uses
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the last electrode of the 160 layout and the
# first two of the 256 layout carry the angles (in degrees) from BioSemi's
# sheet, including the fractional ones of the 256 cap (9.2, not 9).
test_that("160 and 256 layouts keep the angles of BioSemi's caps", {
  m160 <- create_montage(c("A1", "A2", "E32"), template = "biosemi160")
  expect_equal(m160$positions$theta, c(0, 11.5, -115))
  expect_equal(m160$positions$phi, c(0, -90, 45))

  m256 <- create_montage(c("A1", "A2", "H32"), template = "biosemi256")
  expect_equal(m256$positions$theta, c(0, 9.2, -119.6))
  expect_equal(m256$positions$phi, c(0, -90, 66))
})

# ----------------------------------------------------------------------------
# Test 5.5: Names of the large layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 128, 160 and 256 layouts are named by BioSemi
# socket (32 per letter: A1-D32, A1-E32, A1-H32), in that order.
test_that("128, 160 and 256 layouts use BioSemi socket names in order", {
  expect_equal(create_montage(template = "biosemi128")$channels,
               paste0(rep(LETTERS[1:4], each = 32), 1:32))
  expect_equal(create_montage(template = "biosemi160")$channels,
               paste0(rep(LETTERS[1:5], each = 32), 1:32))
  expect_equal(create_montage(template = "biosemi256")$channels,
               paste0(rep(LETTERS[1:8], each = 32), 1:32))
})

# ----------------------------------------------------------------------------
# Test 5.6: 16- and 32-channel layouts accept the socket names too
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 16 and 32 channel caps answer to both the 10-20
# name and BioSemi's socket name (A<row>), in any letter case, and that both
# names give the same position. The channel column keeps the name as typed.
test_that("16 and 32 channel layouts accept 10-20 names and socket names", {
  m32 <- create_montage(c("Cz", "A32", "a1", "Fp1", "A13", "Pz"), template = "biosemi32")
  p <- m32$positions
  expect_equal(p$channel, c("Cz", "A32", "a1", "Fp1", "A13", "Pz"))
  expect_equal(p[1, c("x", "y", "z")], p[2, c("x", "y", "z")], ignore_attr = TRUE)   # Cz = A32
  expect_equal(p[3, c("x", "y", "z")], p[4, c("x", "y", "z")], ignore_attr = TRUE)   # Fp1 = A1
  expect_equal(p[5, c("x", "y", "z")], p[6, c("x", "y", "z")], ignore_attr = TRUE)   # Pz = A13

  m16 <- create_montage(c("A8", "Cz", "A15", "Oz"), template = "biosemi16")
  p <- m16$positions
  expect_equal(p[1, c("x", "y", "z")], p[2, c("x", "y", "z")], ignore_attr = TRUE)   # Cz = A8
  expect_equal(p[3, c("x", "y", "z")], p[4, c("x", "y", "z")], ignore_attr = TRUE)   # Oz = A15
})

# ----------------------------------------------------------------------------
# Test 5.7: 16 and 32 channel layouts agree with the 64 layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode of the 16 and 32 channel caps sits
# at the same angles in the 64 channel layout (they are subsets of the same
# 10-20/10-10 positions), a cross-check between two independent BioSemi tables.
test_that("16 and 32 channel layouts agree with the 64 channel layout", {
  for (tmpl in c("biosemi16", "biosemi32")) {
    small <- create_montage(template = tmpl)
    big   <- create_montage(small$channels, template = "biosemi64")
    expect_equal(big$positions$theta, small$positions$theta, info = tmpl)
    expect_equal(big$positions$phi, small$positions$phi, info = tmpl)
  }
})

# ----------------------------------------------------------------------------
# Test 5.8: The large layouts have no 10-20 names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a 10-20 name is not found in the 128 layout (which
# only knows BioSemi socket names): it is dropped with a warning naming it, and
# a missing (NA) channel name never matches anything.
test_that("128 layout drops 10-20 names and NA names with a warning", {
  expect_warning(m <- create_montage(c("A1", "Cz", "D32"), template = "biosemi128"), "Cz")
  expect_equal(m$channels, c("A1", "D32"))

  expect_warning(m <- create_montage(c("A1", NA), template = "biosemi128"), "NA")
  expect_equal(m$channels, "A1")
})

# ----------------------------------------------------------------------------
# Test 5.9: Unknown or malformed template
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an unknown template stops with a message that
# lists the layouts available, and that a template must be a single string.
test_that("create_montage() lists the available layouts for a bad template", {
  expect_error(create_montage(template = "biosemi96"),
               "Unknown montage template 'biosemi96'.*biosemi16.*biosemi256")
  expect_error(create_montage(template = c("biosemi32", "biosemi64")),
               "Unknown montage template")
  expect_error(create_montage(template = 32), "Unknown montage template")
  expect_error(create_montage(template = NULL), "Unknown montage template")
})

# ----------------------------------------------------------------------------
# Test 5.10: The 64 channel layout is unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the default layout still returns BioSemi's whole
# millimetre values for the 64 channel cap, still accepts both naming systems,
# and still labels the rows with the lower-case names.
test_that("biosemi64 (the default) is unchanged", {
  expect_identical(create_montage(), create_montage(template = "biosemi64"))

  m <- create_montage(c("Fp1", "Cz", "O2"))
  expect_equal(m$positions$x, c(-27, 0, 27))
  expect_equal(m$positions$y, c(83, 0, -83))
  expect_equal(m$positions$z, c(-3, 88, -3))
  expect_equal(rownames(m$positions), c("fp1", "cz", "o2"))

  expect_equal(create_montage("B16")$positions$z, 88)    # B16 is Cz
})

# ----------------------------------------------------------------------------
# Test 5.11: A layout attaches to a matching recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() accepts a 128 channel layout for a
# 128 channel recording named A1-D32, and a 32 channel layout for a recording
# named by BioSemi socket (A1-A32), both without any warning.
test_that("set_montage() attaches BioSemi layouts to matching recordings", {
  chans128 <- paste0(rep(LETTERS[1:4], each = 32), 1:32)
  eeg128 <- new_eeg(data = matrix(rnorm(128 * 20), nrow = 128),
                    channels = chans128, sampling_rate = 100)
  expect_no_warning(eeg128 <- set_montage(eeg128, create_montage(template = "biosemi128")))
  expect_equal(eeg128$montage$channels, chans128)
  expect_equal(eeg128$montage$coord_frame, "biosemi128")

  chans32 <- paste0("A", 1:32)
  eeg32 <- new_eeg(data = matrix(rnorm(32 * 20), nrow = 32),
                   channels = chans32, sampling_rate = 100)
  expect_no_warning(
    eeg32 <- set_montage(eeg32, create_montage(chans32, template = "biosemi32"))
  )
  expect_equal(eeg32$montage$channels, chans32)
  expect_equal(eeg32$montage$positions$theta[32], 0)     # A32 is Cz
})

# ----------------------------------------------------------------------------
# Test 5.12: No angle of any layout has been changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Guards the whole angle tables. Each layout's angles are
# summed with the electrode's row number as weight (so a typo, a swap of two
# electrodes or a shifted row all change the sum). The expected sums were
# computed from the tables when they were checked, value by value, against
# BioSemi's coordinate sheet.
test_that("the angle tables of the BioSemi layouts are unchanged", {
  fingerprint <- function(tmpl) {
    p <- create_montage(template = tmpl)$positions
    i <- seq_len(nrow(p))
    c(theta = sum(p$theta * i), phi = sum(p$phi * i))
  }

  expect_equal(fingerprint("biosemi16"),  c(theta = 2612,       phi = -2142))
  expect_equal(fingerprint("biosemi32"),  c(theta = 19072,      phi = 5539))
  expect_equal(fingerprint("biosemi128"), c(theta = -52313.5,   phi = 22266))
  expect_equal(fingerprint("biosemi160"), c(theta = -157757,    phi = 59992.5))
  expect_equal(fingerprint("biosemi256"), c(theta = -769570.8,  phi = 346766.25))
})

# ============================================================================
#        TEST SUITE 6: read_custom_montage() - Montage From a File
# ============================================================================

# Helper: writes `lines` to a temporary file with the given extension.
write_montage_file <- function(lines, ext) {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}

# Helper: comma separated lines (name, x, y, z, with a header) for a positions
# table, with every position moved by `shift` (mm) and then multiplied by
# `scale` (1 for millimetres, 0.1 for centimetres, 0.001 for metres).
positions_as_csv <- function(positions, shift = c(0, 0, 0), scale = 1) {
  c("name,x,y,z",
    sprintf("%s,%.10f,%.10f,%.10f", positions$channel,
            (positions$x + shift[1]) * scale, (positions$y + shift[2]) * scale,
            (positions$z + shift[3]) * scale))
}

# Five electrodes on a small head, plus the three head landmarks, in an old
# style .elc file (bare "x y z" lines) that states its unit. `scale` multiplies
# every number, to write the same head in another size.
elc_old_style <- function(unit = "mm", scale = 1) {
  coords <- rbind(c(-86.0, -20.0, -48.0), c(85.8, -20.0, -48.0), c(0.0, 86.8, -40.0),
                  c(0, 0, 100), c(50, 0, 80), c(-50, 0, 80),
                  c(0, 60, 80), c(0, -60, 80)) * scale
  write_montage_file(c(
    "# ASA electrode file", "ReferenceLabel\tavg", paste0("UnitPosition\t", unit),
    "NumberPositions=\t8", "Positions",
    sprintf("%.6f %.6f %.6f", coords[, 1], coords[, 2], coords[, 3]),
    "Labels", "LPA", "RPA", "Nz", "Cz", "C4", "C3", "Fz", "Pz", ""
  ), ".elc")
}

# ----------------------------------------------------------------------------
# Test 6.1: An .elc file in the old style
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an .elc file with bare "x y z" lines is read, the
# three head landmarks (LPA, RPA, Nz) are left out, positions are kept as
# written when recenter = FALSE, and the coordinate frame names the file.
test_that("read_custom_montage reads an old style .elc file", {
  f <- elc_old_style()
  m <- read_custom_montage(f, recenter = FALSE)

  expect_s3_class(m, "montage")
  expect_equal(m$channels, c("Cz", "C4", "C3", "Fz", "Pz"))
  expect_equal(m$positions$x, c(0, 50, -50, 0, 0))
  expect_equal(m$positions$y, c(0, 0, 0, 60, -60))
  expect_equal(m$positions$z, c(100, 80, 80, 80, 80))
  expect_equal(m$coord_frame, paste0("file: ", basename(f)))
})

# ----------------------------------------------------------------------------
# Test 6.2: An .elc file in the new style, in metres
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the "E01 : x y z" style, with several labels per
# line, and that positions given in metres come out in millimetres.
test_that("read_custom_montage reads a new style .elc file given in metres", {
  f <- write_montage_file(c(
    "UnitPosition m", "Positions",
    "E01 : 0.0 0.0 0.1", "E02 : 0.05 0.0 0.08", "E03 : -0.05 0.0 0.08",
    "E04 : 0.0 0.06 0.08", "E05 : 0.0 -0.06 0.08",
    "Labels", "Cz C4 C3", "Fz Pz", ""
  ), ".elc")
  m <- read_custom_montage(f, recenter = FALSE)

  expect_equal(m$channels, c("Cz", "C4", "C3", "Fz", "Pz"))
  expect_equal(m$positions$x, c(0, 50, -50, 0, 0))
  expect_equal(m$positions$z, c(100, 80, 80, 80, 80))
})

# ----------------------------------------------------------------------------
# Test 6.3: .sfp, .tsv and .txt files
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the space and tab separated formats, with landmark
# rows (any letter case) left out, comment lines and blank lines ignored, and
# T9 and T10 kept because they are electrodes in the 10-10 system.
test_that("read_custom_montage reads .sfp, .tsv and .txt files", {
  sfp <- write_montage_file(c(
    "FidNz 0.0 86.8 -40.0", "FidT9 -86.0 -20.0 -48.0", "FidT10 85.8 -20.0 -48.0",
    "Cz 0.0 0.0 100.0", "C4\t50.0\t0.0\t80.0", "C3 -50.0 0.0 80.0",
    "Fz 0.0 60.0 80.0", "Pz 0.0 -60.0 80.0"), ".sfp")
  tsv <- write_montage_file(c(
    "Cz\t0\t0\t100", "C4\t50\t0\t80", "C3\t-50\t0\t80",
    "Fz\t0\t60\t80", "Pz\t0\t-60\t80"), ".tsv")
  txt <- write_montage_file(c(
    "# a comment", "", "Cz 0 0 100", "C4 50 0 80", "C3 -50 0 80",
    "Fz 0 60 80", "Pz 0 -60 80"), ".txt")

  for (f in c(sfp, tsv, txt)) {
    m <- read_custom_montage(f, unit = "mm", recenter = FALSE)
    expect_equal(m$channels, c("Cz", "C4", "C3", "Fz", "Pz"), info = f)
    expect_equal(m$positions$x, c(0, 50, -50, 0, 0), info = f)
    expect_equal(m$positions$y, c(0, 0, 0, 60, -60), info = f)
  }

  landmarks <- write_montage_file(c(
    "nasion 0 86 -40", "lpa -86 -20 -48", "RPA 86 -20 -48",
    "T9 -80 -10 -45", "T10 80 -10 -45", "Cz 0 0 100"), ".sfp")
  m <- read_custom_montage(landmarks, unit = "mm", recenter = FALSE)
  expect_equal(m$channels, c("T9", "T10", "Cz"))
})

# ----------------------------------------------------------------------------
# Test 6.4: A .csv file with a header, quotes and a byte-order mark
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies columns are found by their header names even when
# they come in another order, quotes around names are removed, and the
# byte-order mark some programs put at the start of a text file is ignored.
test_that("read_custom_montage reads a .csv file with a header, quotes and a BOM", {
  f <- tempfile(fileext = ".csv")
  rows <- c("X,Y,Z,Label", "0,0,100,\"Cz\"", "50,0,80,\"C4\"", "-50,0,80,\"C3\"",
            "0,60,80,\"Fz\"", "0,-60,80,\"Pz\"")
  writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)),
             charToRaw(paste0(paste(rows, collapse = "\n"), "\n"))), f)

  m <- read_custom_montage(f, unit = "mm", recenter = FALSE)

  expect_equal(m$channels, c("Cz", "C4", "C3", "Fz", "Pz"))
  expect_equal(m$positions$x, c(0, 50, -50, 0, 0))
  expect_equal(m$positions$y, c(0, 0, 0, 60, -60))
  expect_equal(m$positions$z, c(100, 80, 80, 80, 80))
})

# ----------------------------------------------------------------------------
# Test 6.5: A head that is not centred on the origin is recentred
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 32 channel BioSemi head, moved by (5, -12, 4)
# mm and written to a file in millimetres, centimetres or metres, comes back
# exactly as the template (positions, and the angles and radius calculated
# from them). The template is a perfect sphere, so the answer is known in
# advance. Writing it in all three units also checks the shift that puts the
# centre back is converted to millimetres like the positions.
test_that("read_custom_montage puts the centre of a shifted head back at the origin", {
  p <- create_montage(template = "biosemi32")$positions

  for (unit in c("mm", "cm", "m")) {
    scale <- c(mm = 1, cm = 0.1, m = 0.001)[[unit]]
    f <- write_montage_file(
      positions_as_csv(p, shift = c(5, -12, 4), scale = scale), ".csv")

    m <- read_custom_montage(f, unit = unit)

    expect_equal(m$positions$x, p$x, tolerance = 1e-6, info = unit)
    expect_equal(m$positions$y, p$y, tolerance = 1e-6, info = unit)
    expect_equal(m$positions$z, p$z, tolerance = 1e-6, info = unit)
    expect_equal(m$positions$theta, p$theta, tolerance = 1e-6, info = unit)
    expect_equal(m$positions$phi, p$phi, tolerance = 1e-6, info = unit)
    expect_equal(m$positions$radius, rep(87.54, 32), tolerance = 1e-6, info = unit)
  }
})

# ----------------------------------------------------------------------------
# Test 6.6: recenter = FALSE keeps the file's positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the same shifted head is left exactly as written
# when recenter = FALSE.
test_that("read_custom_montage with recenter = FALSE keeps the positions as written", {
  p <- create_montage(template = "biosemi32")$positions
  f <- write_montage_file(positions_as_csv(p, shift = c(5, -12, 4)), ".csv")

  m <- read_custom_montage(f, unit = "mm", recenter = FALSE)

  expect_equal(m$positions$x, p$x + 5, tolerance = 1e-6)
  expect_equal(m$positions$y, p$y - 12, tolerance = 1e-6)
  expect_equal(m$positions$z, p$z + 4, tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 6.7: Recentring works for a head that is not a sphere
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that, for a head whose electrodes sit at different
# distances (alternately 10 percent nearer and farther) and are moved off the
# origin, the best fitting sphere of the montage that comes out is centred on
# the origin.
test_that("read_custom_montage centres a head that is not a sphere", {
  p <- create_montage(template = "biosemi64")$positions
  stretch <- rep(c(0.9, 1.1), length.out = nrow(p))
  odd <- p
  odd$x <- p$x * stretch; odd$y <- p$y * stretch; odd$z <- p$z * stretch
  f <- write_montage_file(positions_as_csv(odd, shift = c(6, -9, 2)), ".csv")

  m <- read_custom_montage(f, unit = "mm")
  fit <- eeganalysis:::.montage_sphere_fit(m$positions$x, m$positions$y,
                                          m$positions$z)

  expect_false(is.null(fit))
  expect_lt(max(abs(fit$centre)), 1e-6)
})

# ----------------------------------------------------------------------------
# Test 6.8: The angles follow BioSemi's convention
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies theta and phi, worked out from x, y and z, equal the
# angles of the 128 channel BioSemi layout (whose positions are not rounded),
# for all 128 electrodes, including the negative inclinations on the left side
# and the ones below the ear line.
test_that("read_custom_montage works out theta and phi like the built-in layouts", {
  p <- create_montage(template = "biosemi128")$positions
  f <- write_montage_file(positions_as_csv(p), ".csv")

  m <- read_custom_montage(f, unit = "mm", recenter = FALSE)

  expect_equal(m$positions$theta, p$theta, tolerance = 1e-6)
  expect_equal(m$positions$phi, p$phi, tolerance = 1e-6)
  expect_equal(m$positions$radius, rep(87.54, 128), tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 6.9: The unit
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the unit is worked out from the head size when the
# file does not state one (with a message), the function stops when the size fits
# no unit, a unit that is given wins over the file and over the head size, and
# a bad unit is rejected.
test_that("read_custom_montage works out, takes or rejects the unit", {
  p <- create_montage(template = "biosemi32")$positions
  scaled_file <- function(scale) {
    write_montage_file(sprintf("%s %.10f %.10f %.10f", p$channel,
                               p$x * scale, p$y * scale, p$z * scale), ".sfp")
  }

  expect_message(m_mm <- read_custom_montage(scaled_file(1)), "taken as mm")
  expect_message(m_cm <- read_custom_montage(scaled_file(0.1)), "taken as cm")
  expect_message(m_m  <- read_custom_montage(scaled_file(0.001)), "taken as m\\.")
  for (m in list(m_mm, m_cm, m_m)) {
    expect_equal(m$positions$x, p$x, tolerance = 1e-6)
    expect_equal(m$positions$z, p$z, tolerance = 1e-6)
  }

  expect_error(read_custom_montage(scaled_file(100)), "Cannot tell the unit")

  expect_no_message(m <- read_custom_montage(scaled_file(0.1), unit = "cm"))
  expect_equal(m$positions$x, p$x, tolerance = 1e-6)

  # a unit that is given wins over the one the file states
  m <- read_custom_montage(elc_old_style("mm"), unit = "cm", recenter = FALSE)
  expect_equal(m$positions$x, c(0, 50, -50, 0, 0) * 10)

  # a unit the file states is used as it is, even if the head size would
  # suggest another one: no message, no guessing
  expect_no_message(
    m <- read_custom_montage(elc_old_style("mm", scale = 0.1), recenter = FALSE)
  )
  expect_equal(m$positions$x, c(0, 5, -5, 0, 0))

  expect_error(read_custom_montage(elc_old_style(), unit = "km"), "'unit' must be")
  expect_error(read_custom_montage(elc_old_style("inches")), "unknown unit")
})

# ----------------------------------------------------------------------------
# Test 6.10: Files that cannot be read
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies each kind of unusable file stops with a message that
# says what is wrong.
test_that("read_custom_montage stops with a clear message for a bad file", {
  expect_error(read_custom_montage(file.path(tempdir(), "no_such_file.elc")),
               "existing file")
  expect_error(read_custom_montage(write_montage_file("Cz 0 0 100", ".xyz")),
               "Unsupported file type '.xyz'")

  expect_error(
    read_custom_montage(write_montage_file(c("Cz 0 0 100", "Cz 1 2 3", "Fz 0 60 80"), ".sfp"),
                        unit = "mm"),
    "more than once in the file: Cz")
  expect_error(
    read_custom_montage(write_montage_file(c("Cz 0 0 100", "Fz abc 60 80"), ".sfp"),
                        unit = "mm"),
    "position of Fz is not made of three numbers")
  expect_error(
    read_custom_montage(write_montage_file(c("Cz 0 0 100", "Fz 0 60"), ".sfp"),
                        unit = "mm"),
    "fewer than four values")
  expect_error(
    read_custom_montage(write_montage_file(c("Nz 0 86 -40", "LPA -86 -20 -48"), ".sfp"),
                        unit = "mm"),
    "only head landmarks")
  expect_error(
    read_custom_montage(write_montage_file(character(0), ".csv"), unit = "mm"),
    "no electrode positions")

  expect_error(
    read_custom_montage(write_montage_file(c(
      "UnitPosition mm", "Positions", "0 0 100", "50 0 80", "-50 0 80",
      "Labels", "Cz", "C4", ""), ".elc")),
    "3 positions but 2 labels")
  expect_error(
    read_custom_montage(write_montage_file(c("UnitPosition mm", "Positions", "0 0 100"), ".elc")),
    "Not a valid .elc file")

  expect_error(read_custom_montage(elc_old_style(), recenter = NA),
               "'recenter' must be")
})

# ----------------------------------------------------------------------------
# Test 6.11: Too few electrodes to find the centre
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when no head centre can be found (fewer than 5
# electrodes, a flat arrangement, or electrodes nowhere near a sphere) the
# positions are kept as written and a warning says the centre could not be found.
test_that("read_custom_montage warns when the head centre cannot be found", {
  f <- write_montage_file(c("Cz 0 0 100", "C4 50 0 80", "C3 -50 0 80"), ".sfp")

  expect_warning(m <- read_custom_montage(f, unit = "mm"),
                 "centre of the head could not be found")
  expect_equal(m$positions$x, c(0, 50, -50))

  # four electrodes that are not on one plane fit a sphere exactly, with nothing
  # left over to tell whether it is a head, so five are needed
  four <- write_montage_file(c("Cz 0 0 100", "C4 50 0 80", "Fz 0 60 80",
                               "Pz 0 -60 80"), ".sfp")
  expect_warning(m <- read_custom_montage(four, unit = "mm"),
                 "centre of the head could not be found")
  expect_equal(m$positions$y, c(0, 0, 60, -60))

  # six electrodes on a flat grid do not determine a head centre either
  flat <- write_montage_file(c("A 0 0 0", "B 10 0 0", "C 20 0 0",
                               "D 0 10 0", "E 10 10 0", "F 20 10 0"), ".sfp")
  expect_warning(m2 <- read_custom_montage(flat, unit = "mm"),
                 "centre of the head could not be found")
  expect_equal(m2$positions$x, c(0, 10, 20, 0, 10, 20))

  # electrodes alternately 40 mm and 100 mm from the origin do give a sphere,
  # but one that is nowhere near them, so it is not used either
  p <- create_montage(template = "biosemi32")$positions
  odd <- p
  near_far <- rep(c(40, 100) / 87.54, length.out = nrow(p))
  odd$x <- p$x * near_far; odd$y <- p$y * near_far; odd$z <- p$z * near_far
  not_a_head <- write_montage_file(positions_as_csv(odd), ".csv")
  expect_warning(m3 <- read_custom_montage(not_a_head, unit = "mm"),
                 "centre of the head could not be found")
  expect_equal(m3$positions$x, odd$x, tolerance = 1e-6)

  # a nearly flat bowl (z = (x^2 + y^2) / 2000) fits a sphere of about a metre
  # radius with almost no error; that would move every electrode by about a
  # metre, so it is not used either
  grid <- expand.grid(x = seq(-40, 40, by = 20), y = seq(-30, 30, by = 20))
  bowl <- sprintf("E%02d %.6f %.6f %.6f", seq_len(nrow(grid)), grid$x, grid$y,
                  (grid$x^2 + grid$y^2) / 2000)
  expect_warning(m4 <- read_custom_montage(write_montage_file(bowl, ".sfp"),
                                           unit = "mm"),
                 "centre of the head could not be found")
  expect_equal(m4$positions$x, grid$x)
})

# ----------------------------------------------------------------------------
# Test 6.12: The montage is built like the one from create_montage()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the object has the same parts, the same position
# columns (all numbers), the lower-case names as row names when they are all
# different, and prints like any montage.
test_that("read_custom_montage returns a montage built like create_montage()", {
  ref <- create_montage(template = "biosemi32")
  f <- write_montage_file(positions_as_csv(ref$positions), ".csv")

  m <- read_custom_montage(f, unit = "mm", recenter = FALSE)

  expect_s3_class(m, "montage")
  expect_identical(names(m), names(ref))
  expect_identical(names(m$positions), names(ref$positions))
  expect_identical(m$channels, m$positions$channel)
  expect_identical(rownames(m$positions), tolower(m$channels))
  expect_true(all(vapply(m$positions[-1], is.numeric, logical(1))))
  expect_true(any(grepl("file: ", capture.output(print(m)), fixed = TRUE)))

  clash <- write_montage_file(c("Cz 0 0 100", "CZ 5 0 99", "Fz 0 60 80"), ".sfp")
  m2 <- read_custom_montage(clash, unit = "mm", recenter = FALSE)
  expect_equal(rownames(m2$positions), c("1", "2", "3"))
})

# ----------------------------------------------------------------------------
# Test 6.13: A montage from a file attaches to a recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() accepts a montage read from a file for
# a recording with the same channel names, with no warning, and keeps its
# coordinate frame.
test_that("a montage read from a file attaches with set_montage()", {
  ref <- create_montage(template = "biosemi32")
  f <- write_montage_file(positions_as_csv(ref$positions, shift = c(3, -5, 2)), ".csv")
  m <- read_custom_montage(f, unit = "mm")

  eeg <- new_eeg(data = matrix(rnorm(32 * 20), nrow = 32),
                 channels = m$channels, sampling_rate = 100)
  expect_no_warning(eeg <- set_montage(eeg, m))

  expect_equal(eeg$montage$channels, m$channels)
  expect_equal(eeg$montage$coord_frame, m$coord_frame)
})

# ============================================================================
#     TEST SUITE 7: standard_1010 - the standard 10-10 electrode system
# ============================================================================

# Helper: the positions of the standard_1010 layout.
standard_1010_positions <- function() {
  create_montage(template = "standard_1010")$positions
}

# Helper: signed angles of the named electrodes, in the order of `names`.
angles_of <- function(positions, names, column) {
  positions[[column]][match(names, positions$channel)]
}

# ----------------------------------------------------------------------------
# Test 7.1: Which electrodes the layout has
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies standard_1010 has the 72 electrodes of the 10-10
# standard (its 74 positions without the ear electrodes A1 and A2), each name
# once, in front-to-back order, and leaves out the ear electrodes, the nasion
# landmark and the eight AF/PO positions the rules give but the standard does
# not include.
test_that("standard_1010 has the 72 electrodes of the 10-10 standard without the ears", {
  m <- create_montage(template = "standard_1010")

  expect_s3_class(m, "montage")
  expect_equal(length(m$channels), 72)
  expect_false(anyDuplicated(m$channels) > 0)
  expect_equal(m$coord_frame, "standard_1010")
  expect_equal(head(m$channels, 5), c("Fp1", "Fpz", "Fp2", "AF7", "AF3"))
  expect_equal(tail(m$channels, 4), c("O1", "Oz", "O2", "Iz"))
  expect_true(all(c("Cz", "Pz", "AFz", "FT9", "TP10", "P9", "F10", "T9", "T10",
                    "PO8") %in% m$channels))
  expect_false(any(c("A1", "A2", "Nz", "AF1", "AF2", "AF5", "AF6", "PO1", "PO2",
                     "PO5", "PO6") %in% m$channels))
  expect_false(anyNA(m$positions))
})

# ----------------------------------------------------------------------------
# Test 7.2: Same shape and scale as the other layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode lies on the 87.54 mm head sphere,
# the positions table has the same columns as the BioSemi layouts, and the
# lower-case names are the row names.
test_that("standard_1010 is built like the other layouts and sits on the head sphere", {
  m <- create_montage(template = "standard_1010")
  p <- m$positions

  expect_equal(sqrt(p$x^2 + p$y^2 + p$z^2), rep(87.54, 72), tolerance = 1e-6)
  expect_true(all(p$radius == 87.54))
  expect_identical(names(p), names(create_montage(template = "biosemi32")$positions))
  expect_identical(rownames(p), tolower(m$channels))

  # no value is a negative zero (printed and compared like 0, but not identical)
  numbers <- unlist(p[, c("x", "y", "z", "theta", "phi")])
  expect_false(any(numbers == 0 & 1 / numbers < 0))
})

# ----------------------------------------------------------------------------
# Test 7.3: The midline, the ear line, the outer ring and the inferior row
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the electrodes the 10-10 rules place at 10 percent
# steps. The expected angles were worked out by hand from the rules, not by the
# code: the nasion-inion arc is 230 degrees (landmarks 115 degrees from the
# vertex), so each 10 percent step is 23 degrees; the outer ring electrodes
# are 18 degrees apart around a circle 92 degrees from the vertex; the
# inferior row is on the circle 115 degrees from the vertex.
test_that("standard_1010 places the 10 percent electrodes where the rules put them", {
  p <- standard_1010_positions()

  midline <- c("Fpz", "AFz", "Fz", "FCz", "Cz", "CPz", "Pz", "POz", "Oz", "Iz")
  expect_equal(angles_of(p, midline, "theta"), c(92, 69, 46, 23, 0, 23, 46, 69, 92, 115))
  expect_equal(angles_of(p, midline, "phi"), c(90, 90, 90, 90, 0, -90, -90, -90, -90, -90))

  ear_line <- c("T9", "T7", "C5", "C3", "C1", "C2", "C4", "C6", "T8", "T10")
  expect_equal(angles_of(p, ear_line, "theta"),
               c(-115, -92, -69, -46, -23, 23, 46, 69, 92, 115))
  expect_equal(angles_of(p, ear_line, "phi"), rep(0, 10), tolerance = 1e-6)

  left_ring <- c("Fp1", "AF7", "F7", "FT7", "TP7", "P7", "PO7", "O1")
  expect_equal(angles_of(p, left_ring, "theta"), rep(-92, 8))
  expect_equal(angles_of(p, left_ring, "phi"), c(-72, -54, -36, -18, 18, 36, 54, 72))
  right_ring <- c("Fp2", "AF8", "F8", "FT8", "TP8", "P8", "PO8", "O2")
  expect_equal(angles_of(p, right_ring, "theta"), rep(92, 8))
  expect_equal(angles_of(p, right_ring, "phi"), c(72, 54, 36, 18, -18, -36, -54, -72))

  inferior_left <- c("F9", "FT9", "TP9", "P9")
  expect_equal(angles_of(p, inferior_left, "theta"), rep(-115, 4))
  expect_equal(angles_of(p, inferior_left, "phi"), c(-36, -18, 18, 36))
  inferior_right <- c("F10", "FT10", "TP10", "P10")
  expect_equal(angles_of(p, inferior_right, "theta"), rep(115, 4))
  expect_equal(angles_of(p, inferior_right, "phi"), c(36, 18, -18, -36))
})

# ----------------------------------------------------------------------------
# Test 7.4: Rows are cut into equal parts along one circle
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the property that defines the rows, independent of
# how the code gets there: the electrodes of a row lie in one plane (one
# circle on the sphere) and consecutive electrodes are the same distance
# apart. The same holds for the outer ring and for the midline.
test_that("standard_1010 rows are equally spaced along one circle", {
  p <- standard_1010_positions()
  xyz <- function(names) {
    as.matrix(p[match(names, p$channel), c("x", "y", "z")])
  }

  rows <- list(
    AF   = c("AF7", "AF3", "AFz", "AF4", "AF8"),
    F    = c("F7", "F5", "F3", "F1", "Fz", "F2", "F4", "F6", "F8"),
    FC   = c("FT7", "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6", "FT8"),
    C    = c("T7", "C5", "C3", "C1", "Cz", "C2", "C4", "C6", "T8"),
    CP   = c("TP7", "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6", "TP8"),
    P    = c("P7", "P5", "P3", "P1", "Pz", "P2", "P4", "P6", "P8"),
    PO   = c("PO7", "PO3", "POz", "PO4", "PO8"),
    ring = c("Fp1", "AF7", "F7", "FT7", "T7", "TP7", "P7", "PO7", "O1"),
    mid  = c("Fpz", "AFz", "Fz", "FCz", "Cz", "CPz", "Pz", "POz", "Oz")
  )
  for (row in names(rows)) {
    pts <- xyz(rows[[row]])
    steps <- sqrt(rowSums((pts[-1, ] - pts[-nrow(pts), ])^2))
    expect_lt(max(steps) / min(steps), 1 + 1e-5, label = paste(row, "spacing"))

    centred <- sweep(pts, 2, colMeans(pts))
    flatness <- svd(centred)$d[3] / svd(centred)$d[1]
    expect_lt(flatness, 1e-6, label = paste(row, "flatness"))
  }
})

# ----------------------------------------------------------------------------
# Test 7.5: Left and right mirror each other
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode with an odd number (1, 3, 5, 7, 9)
# has a partner on the other side (2, 4, 6, 8, 10) at the mirrored position.
test_that("standard_1010 is mirror symmetric", {
  p <- standard_1010_positions()
  left <- p$channel[grepl("[13579]$", p$channel)]
  partner <- function(name) {
    if (grepl("9$", name)) return(sub("9$", "10", name))
    chartr("1357", "2468", name)
  }
  right <- vapply(left, partner, character(1))

  expect_length(left, 31)
  expect_true(all(right %in% p$channel))

  a <- p[match(left, p$channel), c("x", "y", "z")]
  b <- p[match(right, p$channel), c("x", "y", "z")]
  expect_equal(a$x, -b$x, tolerance = 1e-6)
  expect_equal(a$y, b$y, tolerance = 1e-6)
  expect_equal(a$z, b$z, tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 7.6: Agreement with the BioSemi 64 channel layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 64 electrodes standard_1010 shares with the
# BioSemi 64 channel layout are at the same places up to BioSemi's rounding
# (it prints whole millimetres). Sixteen are known to differ by more than
# 1.3 mm: AF3, AF4, PO3 and PO4 by about 5.7 mm, FC5, FC6, CP5 and CP6 by
# about 2.1 mm, F3, F4, P3 and P4 by about 1.7 mm, F5, F6, P5 and P6 by about
# 1.5 mm. The other 48 are within 1.3 mm. Nothing else may differ more.
test_that("standard_1010 agrees with the BioSemi 64 layout where the rules are unambiguous", {
  p <- standard_1010_positions()
  b <- create_montage(template = "biosemi64")$positions
  shared <- intersect(b$channel, p$channel)
  expect_length(shared, 64)

  a <- p[match(shared, p$channel), c("x", "y", "z")]
  z <- b[match(shared, b$channel), c("x", "y", "z")]
  distance <- setNames(sqrt(rowSums((as.matrix(a) - as.matrix(z))^2)), shared)

  expect_lt(median(distance), 0.7)
  expect_lt(max(distance), 6.5)
  expect_setequal(names(distance)[distance > 1.3],
                  c("AF3", "AF4", "PO3", "PO4", "FC5", "FC6", "CP5", "CP6",
                    "F3", "F4", "P3", "P4", "F5", "F6", "P5", "P6"))
})

# ----------------------------------------------------------------------------
# Test 7.7: No number of the table has been changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Guards the whole table. Each column is summed with the row
# number as weight, so a changed value, a swapped pair or a shifted row alters
# the sum. The expected sums were computed after the table had been checked
# against the rules and the BioSemi layout (tests 7.3 to 7.6).
test_that("the standard_1010 table is unchanged", {
  p <- standard_1010_positions()
  i <- seq_len(nrow(p))
  sums <- vapply(c("x", "y", "z", "theta", "phi"),
                 function(col) sum(p[[col]] * i), numeric(1))

  expect_equal(sums, c(x = 10983.471656, y = -80079.227446, z = 53235.56376,
                       theta = 39294.388892, phi = -23940), tolerance = 1e-8)
})

# ----------------------------------------------------------------------------
# Test 7.8: Names that are not in the layout, and the template list
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies names the 10-10 standard does not have (the ear
# electrodes, an AF position it leaves out) are dropped with a warning naming
# them, names match in any letter case, and an unknown template lists
# standard_1010 among the ones available.
test_that("create_montage() handles names outside standard_1010 and lists it", {
  expect_warning(m <- create_montage(c("Cz", "A1", "AF1", "pz"), template = "standard_1010"),
                 "A1, AF1")
  expect_equal(m$channels, c("Cz", "pz"))

  expect_error(create_montage(template = "standard_1030"),
               "Available templates: .*biosemi256, standard_1020, standard_1010")
})

# ----------------------------------------------------------------------------
# Test 7.9: The layout attaches to a matching recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() accepts the layout for a recording
# whose 72 channels carry its names, without any warning.
test_that("standard_1010 attaches to a recording with its channel names", {
  m <- create_montage(template = "standard_1010")
  eeg <- new_eeg(data = matrix(rnorm(72 * 20), nrow = 72),
                 channels = m$channels, sampling_rate = 100)

  expect_no_warning(eeg <- set_montage(eeg, m))
  expect_equal(eeg$montage$channels, m$channels)
  expect_equal(eeg$montage$coord_frame, "standard_1010")
})

# ============================================================================
#     TEST SUITE 8: standard_1005 - the standard 10-05 electrode system
# ============================================================================

# Helper: the positions of the standard_1005 layout.
standard_1005_positions <- function() {
  create_montage(template = "standard_1005")$positions
}

# Helpers: how the 15 rows between Fp and O are named. `medial` names a row
# from column 5 inwards and on the midline, `lateral` from column 7 outwards
# (the rows between F and P run along the temporal electrodes there).
row_medial  <- c("AFp", "AF", "AFF", "F", "FFC", "FC", "FCC", "C", "CCP", "CP",
                 "CPP", "P", "PPO", "PO", "POO")
row_lateral <- c("AFp", "AF", "AFF", "F", "FFT", "FT", "FTT", "T", "TTP", "TP",
                 "TPP", "P", "PPO", "PO", "POO")

# The 21 electrodes of row i, from the far left to the far right.
row_of_21 <- function(i) {
  c(paste0(row_lateral[i], c("9", "9h", "7", "7h")),
    paste0(row_medial[i], c("5", "5h", "3", "3h", "1", "1h", "z", "2h", "2",
                            "4h", "4", "6h", "6")),
    paste0(row_lateral[i], c("8h", "8", "10h", "10")))
}

# ----------------------------------------------------------------------------
# Test 8.1: Which electrodes the layout has
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies standard_1005 has exactly the 335 electrodes of the
# 10-05 system, written out here row by row: 15 rows of 21 electrodes between
# the five Fp and the five O, OI and I electrodes, front to back and left to
# right. The list was also compared with an independent list of the system's
# names. Every electrode of standard_1010 is among them; the landmarks, the
# mastoid and ear electrodes and the old names T3 to T6 are not (the old names
# are aliases of T7, T8, P7 and P8, see test 9.5).
test_that("standard_1005 has the 335 electrodes of the 10-05 system", {
  rows <- c(
    "Fp1 Fp1h Fpz Fp2h Fp2",
    "AFp9 AFp9h AFp7 AFp7h AFp5 AFp5h AFp3 AFp3h AFp1 AFp1h AFpz AFp2h AFp2 AFp4h AFp4 AFp6h AFp6 AFp8h AFp8 AFp10h AFp10",
    "AF9 AF9h AF7 AF7h AF5 AF5h AF3 AF3h AF1 AF1h AFz AF2h AF2 AF4h AF4 AF6h AF6 AF8h AF8 AF10h AF10",
    "AFF9 AFF9h AFF7 AFF7h AFF5 AFF5h AFF3 AFF3h AFF1 AFF1h AFFz AFF2h AFF2 AFF4h AFF4 AFF6h AFF6 AFF8h AFF8 AFF10h AFF10",
    "F9 F9h F7 F7h F5 F5h F3 F3h F1 F1h Fz F2h F2 F4h F4 F6h F6 F8h F8 F10h F10",
    "FFT9 FFT9h FFT7 FFT7h FFC5 FFC5h FFC3 FFC3h FFC1 FFC1h FFCz FFC2h FFC2 FFC4h FFC4 FFC6h FFC6 FFT8h FFT8 FFT10h FFT10",
    "FT9 FT9h FT7 FT7h FC5 FC5h FC3 FC3h FC1 FC1h FCz FC2h FC2 FC4h FC4 FC6h FC6 FT8h FT8 FT10h FT10",
    "FTT9 FTT9h FTT7 FTT7h FCC5 FCC5h FCC3 FCC3h FCC1 FCC1h FCCz FCC2h FCC2 FCC4h FCC4 FCC6h FCC6 FTT8h FTT8 FTT10h FTT10",
    "T9 T9h T7 T7h C5 C5h C3 C3h C1 C1h Cz C2h C2 C4h C4 C6h C6 T8h T8 T10h T10",
    "TTP9 TTP9h TTP7 TTP7h CCP5 CCP5h CCP3 CCP3h CCP1 CCP1h CCPz CCP2h CCP2 CCP4h CCP4 CCP6h CCP6 TTP8h TTP8 TTP10h TTP10",
    "TP9 TP9h TP7 TP7h CP5 CP5h CP3 CP3h CP1 CP1h CPz CP2h CP2 CP4h CP4 CP6h CP6 TP8h TP8 TP10h TP10",
    "TPP9 TPP9h TPP7 TPP7h CPP5 CPP5h CPP3 CPP3h CPP1 CPP1h CPPz CPP2h CPP2 CPP4h CPP4 CPP6h CPP6 TPP8h TPP8 TPP10h TPP10",
    "P9 P9h P7 P7h P5 P5h P3 P3h P1 P1h Pz P2h P2 P4h P4 P6h P6 P8h P8 P10h P10",
    "PPO9 PPO9h PPO7 PPO7h PPO5 PPO5h PPO3 PPO3h PPO1 PPO1h PPOz PPO2h PPO2 PPO4h PPO4 PPO6h PPO6 PPO8h PPO8 PPO10h PPO10",
    "PO9 PO9h PO7 PO7h PO5 PO5h PO3 PO3h PO1 PO1h POz PO2h PO2 PO4h PO4 PO6h PO6 PO8h PO8 PO10h PO10",
    "POO9 POO9h POO7 POO7h POO5 POO5h POO3 POO3h POO1 POO1h POOz POO2h POO2 POO4h POO4 POO6h POO6 POO8h POO8 POO10h POO10",
    "O1 O1h Oz O2h O2",
    "OI1 OI1h OIz OI2h OI2",
    "I1 I1h Iz I2h I2"
  )
  expected <- unlist(strsplit(rows, " ", fixed = TRUE))
  m <- create_montage(template = "standard_1005")

  expect_s3_class(m, "montage")
  expect_equal(m$coord_frame, "standard_1005")
  expect_length(expected, 335)
  expect_identical(m$channels, expected)
  expect_false(anyDuplicated(m$channels) > 0)
  expect_true(all(create_montage(template = "standard_1010")$channels %in% m$channels))
  expect_false(any(c("M1", "M2", "A1", "A2", "Nz", "LPA", "RPA", "T3", "T4",
                     "T5", "T6") %in% m$channels))
  expect_false(anyNA(m$positions))
})

# ----------------------------------------------------------------------------
# Test 8.2: Same shape and scale as the other layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode lies on the 87.54 mm head sphere,
# the positions table has the same columns as the other layouts, the
# lower-case names are the row names, and no value is a negative zero.
test_that("standard_1005 is built like the other layouts and sits on the head sphere", {
  m <- create_montage(template = "standard_1005")
  p <- m$positions

  expect_equal(sqrt(p$x^2 + p$y^2 + p$z^2), rep(87.54, 335), tolerance = 1e-6)
  expect_true(all(p$radius == 87.54))
  expect_identical(names(p), names(create_montage(template = "biosemi32")$positions))
  expect_identical(rownames(p), tolower(m$channels))

  numbers <- unlist(p[, c("x", "y", "z", "theta", "phi")])
  expect_false(any(numbers == 0 & 1 / numbers < 0))
})

# ----------------------------------------------------------------------------
# Test 8.3: The 10-10 electrodes did not move
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 72 electrodes of standard_1010 are in
# standard_1005 at exactly the same place: not close, identical, in every
# column of the positions table.
test_that("standard_1005 holds every standard_1010 electrode at exactly the same place", {
  a <- create_montage(template = "standard_1010")$positions
  b <- create_montage(a$channel, template = "standard_1005")$positions

  expect_identical(b, a)
})

# ----------------------------------------------------------------------------
# Test 8.4: The 5 percent electrodes are where the rules put them
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the angles of the electrodes cut at 5 percent
# steps. The expected values were worked out by hand from the rules, not by
# the code: the nasion-inion arc is 230 degrees, so 5 percent of it is 11.5
# degrees; the outer ring is 180 degrees from Fpz to Oz at 92 degrees from the
# vertex, so a step is 9 degrees; the lower row is the circle 115 degrees
# from the vertex, also 9 degrees a step. The C row is a great circle with
# the ear line, so its electrodes run in 11.5 degree steps from T9 to T10.
test_that("standard_1005 places the 5 percent electrodes where the rules put them", {
  p <- standard_1005_positions()

  midline <- c("Fpz", "AFpz", "AFz", "AFFz", "Fz", "FFCz", "FCz", "FCCz", "Cz",
               "CCPz", "CPz", "CPPz", "Pz", "PPOz", "POz", "POOz", "Oz", "OIz", "Iz")
  expect_equal(angles_of(p, midline, "theta"),
               c(92, 80.5, 69, 57.5, 46, 34.5, 23, 11.5, 0, 11.5, 23, 34.5, 46,
                 57.5, 69, 80.5, 92, 103.5, 115))
  expect_equal(angles_of(p, midline, "phi"), c(rep(90, 8), 0, rep(-90, 10)))

  left_ring <- c("Fp1h", "Fp1", "AFp7", "AF7", "AFF7", "F7", "FFT7", "FT7", "FTT7",
                 "T7", "TTP7", "TP7", "TPP7", "P7", "PPO7", "PO7", "POO7", "O1", "O1h")
  right_ring <- c("Fp2h", "Fp2", "AFp8", "AF8", "AFF8", "F8", "FFT8", "FT8", "FTT8",
                  "T8", "TTP8", "TP8", "TPP8", "P8", "PPO8", "PO8", "POO8", "O2", "O2h")
  expect_equal(angles_of(p, left_ring, "theta"), rep(-92, 19))
  expect_equal(angles_of(p, left_ring, "phi"), seq(-81, 81, by = 9), tolerance = 1e-6)
  expect_equal(angles_of(p, right_ring, "theta"), rep(92, 19))
  expect_equal(angles_of(p, right_ring, "phi"), seq(81, -81, by = -9), tolerance = 1e-6)

  inferior_left <- c("AFp9", "AF9", "AFF9", "F9", "FFT9", "FT9", "FTT9", "T9", "TTP9",
                     "TP9", "TPP9", "P9", "PPO9", "PO9", "POO9", "I1", "I1h")
  inferior_right <- c("AFp10", "AF10", "AFF10", "F10", "FFT10", "FT10", "FTT10", "T10",
                      "TTP10", "TP10", "TPP10", "P10", "PPO10", "PO10", "POO10", "I2", "I2h")
  expect_equal(angles_of(p, inferior_left, "theta"), rep(-115, 17))
  expect_equal(angles_of(p, inferior_left, "phi"), seq(-63, 81, by = 9), tolerance = 1e-6)
  expect_equal(angles_of(p, inferior_right, "theta"), rep(115, 17))
  expect_equal(angles_of(p, inferior_right, "phi"), seq(63, -81, by = -9), tolerance = 1e-6)

  c_row <- row_of_21(8)
  expect_equal(c_row[c(1, 11, 21)], c("T9", "Cz", "T10"))
  expect_equal(angles_of(p, c_row, "theta"), -115 + 11.5 * (0:20))
  expect_equal(angles_of(p, c_row, "phi"), rep(0, 21), tolerance = 1e-6)

  below_o <- c("OI1", "OI1h", "OI2", "OI2h")
  expect_equal(angles_of(p, below_o, "theta"), c(-103.5, -103.5, 103.5, 103.5))
  expect_equal(angles_of(p, below_o, "phi"), c(72, 81, -72, -81), tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 8.5: Rows are cut into equal parts along one circle
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the property that defines the rows, independent of
# how the code gets there: from column 7 to column 8, the 17 electrodes of
# each of the 15 rows lie in one plane (one circle on the sphere) and
# consecutive electrodes are the same distance apart. The same holds for the
# outer ring and for the midline.
test_that("standard_1005 rows are cut into equal parts along one circle", {
  p <- standard_1005_positions()
  xyz <- function(names) {
    as.matrix(p[match(names, p$channel), c("x", "y", "z")])
  }
  circles <- lapply(seq_along(row_medial), function(i) row_of_21(i)[3:19])
  names(circles) <- row_medial
  circles$ring <- c("Fp1h", "Fp1", "AFp7", "AF7", "AFF7", "F7", "FFT7", "FT7", "FTT7",
                    "T7", "TTP7", "TP7", "TPP7", "P7", "PPO7", "PO7", "POO7", "O1", "O1h")
  circles$mid <- c("Fpz", "AFpz", "AFz", "AFFz", "Fz", "FFCz", "FCz", "FCCz", "Cz",
                   "CCPz", "CPz", "CPPz", "Pz", "PPOz", "POz", "POOz", "Oz", "OIz")

  for (circle in names(circles)) {
    pts <- xyz(circles[[circle]])
    steps <- sqrt(rowSums((pts[-1, ] - pts[-nrow(pts), ])^2))
    expect_lt(max(steps) / min(steps), 1 + 1e-5, label = paste(circle, "spacing"))

    centred <- sweep(pts, 2, colMeans(pts))
    flatness <- svd(centred)$d[3] / svd(centred)$d[1]
    expect_lt(flatness, 1e-6, label = paste(circle, "flatness"))
  }
})

# ----------------------------------------------------------------------------
# Test 8.6: Half electrodes outside the ear line, and OI
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the places the paper leaves open and the code
# fixes: the electrode 9h (10h) is halfway along the shortest path between
# its neighbours 7 and 9 (8 and 10), and OI1 (OI2) is halfway between O1 (O2)
# and I1 (I2). "Halfway along the shortest path" means the direction of the
# sum of the two neighbours' directions.
test_that("standard_1005 puts the half electrodes outside the ear line halfway", {
  p <- standard_1005_positions()
  direction <- function(names) {
    v <- unname(as.matrix(p[match(names, p$channel), c("x", "y", "z")]))
    v / sqrt(rowSums(v^2))
  }
  halfway <- function(a, b) {
    v <- direction(a) + direction(b)
    v / sqrt(rowSums(v^2))
  }

  expect_equal(direction(paste0(row_lateral, "9h")),
               halfway(paste0(row_lateral, "9"), paste0(row_lateral, "7")),
               tolerance = 1e-6)
  expect_equal(direction(paste0(row_lateral, "10h")),
               halfway(paste0(row_lateral, "10"), paste0(row_lateral, "8")),
               tolerance = 1e-6)
  expect_equal(direction(c("OI1", "OI2")),
               halfway(c("O1", "O2"), c("I1", "I2")), tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 8.7: New rows lie halfway between the rows around them
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the paper's own example for the new rows: the
# electrode halfway between C3h and CP3h is CCP3h. It holds for the whole
# grid: at every column, each new row between two rows (AFF, FFC, FCC, CCP,
# CPP and PPO) is within a fraction of a millimetre of the halfway point of
# the two rows around it.
test_that("standard_1005 new rows lie halfway between the rows around them", {
  p <- standard_1005_positions()
  direction <- function(names) {
    v <- unname(as.matrix(p[match(names, p$channel), c("x", "y", "z")]))
    v / sqrt(rowSums(v^2))
  }
  distance_mm <- function(a, b) 87.54 * sqrt(rowSums((a - b)^2))
  halfway <- function(a, b) {
    v <- a + b
    v / sqrt(rowSums(v^2))
  }

  centre <- halfway(direction(c("C3h", "C4h")), direction(c("CP3h", "CP4h")))
  expect_lt(max(distance_mm(direction(c("CCP3h", "CCP4h")), centre)), 0.3)

  for (i in c(3, 5, 7, 9, 11, 13)) {          # AFF, FFC, FCC, CCP, CPP, PPO
    between <- halfway(direction(row_of_21(i - 1)), direction(row_of_21(i + 1)))
    gap <- distance_mm(direction(row_of_21(i)), between)
    expect_lt(max(gap), 0.6, label = paste(row_medial[i], "row"))
  }
})

# ----------------------------------------------------------------------------
# Test 8.8: Left and right mirror each other
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode on the left (columns 1, 3, 5, 7, 9,
# and the same with an h) has a partner on the right (columns 2, 4, 6, 8, 10)
# at the mirrored position. Fp1/Fp2, O1/O2, OI1/OI2 and I1/I2 count as
# columns 1 and 2.
test_that("standard_1005 is mirror symmetric", {
  p <- standard_1005_positions()
  sides <- c("1" = "2", "2" = "1", "3" = "4", "4" = "3", "5" = "6", "6" = "5",
             "7" = "8", "8" = "7", "9" = "10", "10" = "9")
  partner <- function(name) {
    part <- regmatches(name, regexec("^(.*?)(10|[1-9])(h?)$", name))[[1]]
    paste0(part[2], sides[[part[3]]], part[4])
  }
  left <- p$channel[p$x < -1e-6]
  right <- vapply(left, partner, character(1))

  expect_length(left, 158)
  expect_true(all(right %in% p$channel))
  expect_equal(sum(abs(p$x) < 1e-6), 19)      # the midline electrodes

  a <- p[match(left, p$channel), c("x", "y", "z")]
  b <- p[match(right, p$channel), c("x", "y", "z")]
  expect_equal(a$x, -b$x, tolerance = 1e-6)
  expect_equal(a$y, b$y, tolerance = 1e-6)
  expect_equal(a$z, b$z, tolerance = 1e-6)
})

# ----------------------------------------------------------------------------
# Test 8.9: No number of the table has been changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Guards the whole table. Each column is summed with the row
# number as weight, so a changed value, a swapped pair or a shifted row alters
# the sum. The expected sums were computed after the table had been checked
# against the rules and an independent list of the names (tests 8.1 to 8.8).
test_that("the standard_1005 table is unchanged", {
  p <- standard_1005_positions()
  i <- seq_len(nrow(p))
  sums <- vapply(c("x", "y", "z", "theta", "phi"),
                 function(col) sum(p[[col]] * i), numeric(1))

  expect_equal(sums, c(x = 97768.371484, y = -1979969.041776, z = 1017234.760881,
                       theta = 352566.900276, phi = -195030), tolerance = 1e-8)
})

# ----------------------------------------------------------------------------
# Test 8.10: Names that are not in the layout, and the template list
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a few electrodes can be picked by name in any
# letter case without a warning, names the 10-05 system does not have (the
# mastoid and ear electrodes) are dropped with a warning naming them, and an
# unknown template lists the standard layouts among the available ones.
test_that("create_montage() picks electrodes from standard_1005 and lists it", {
  expect_no_warning(m <- create_montage(c("Cz", "c1h", "CCP3h", "FFT7"),
                                        template = "standard_1005"))
  expect_equal(m$channels, c("Cz", "c1h", "CCP3h", "FFT7"))

  expect_warning(m <- create_montage(c("Cz", "M1", "A2"), template = "standard_1005"),
                 "M1, A2")
  expect_equal(m$channels, "Cz")

  expect_error(create_montage(template = "standard_1030"),
               "Available templates: .*standard_1020, standard_1010, standard_1005")
})

# ----------------------------------------------------------------------------
# Test 8.11: The layout attaches to a matching recording
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies set_montage() accepts the layout for a recording
# whose 335 channels carry its names, without any warning.
test_that("standard_1005 attaches to a recording with its channel names", {
  m <- create_montage(template = "standard_1005")
  eeg <- new_eeg(data = matrix(rnorm(335 * 20), nrow = 335),
                 channels = m$channels, sampling_rate = 100)

  expect_no_warning(eeg <- set_montage(eeg, m))
  expect_equal(eeg$montage$channels, m$channels)
  expect_equal(eeg$montage$coord_frame, "standard_1005")
})

# ============================================================================
#     TEST SUITE 9: standard_1020 - the original 10-20 electrode system
# ============================================================================

# Helper: the positions of the standard_1020 layout.
standard_1020_positions <- function() {
  create_montage(template = "standard_1020")$positions
}

# ----------------------------------------------------------------------------
# Test 9.1: Which electrodes the layout has
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies standard_1020 has exactly the 19 scalp electrodes
# of the original 10-20 system, front to back and left to right. Fpz and Oz
# (not in the original system), the ear and mastoid electrodes, the nasion and
# the old names T3 to T6 are not among the channel names.
test_that("standard_1020 has the 19 electrodes of the original 10-20 system", {
  m <- create_montage(template = "standard_1020")

  expect_s3_class(m, "montage")
  expect_equal(m$coord_frame, "standard_1020")
  expect_identical(m$channels,
                   c("Fp1", "Fp2",
                     "F7", "F3", "Fz", "F4", "F8",
                     "T7", "C3", "Cz", "C4", "T8",
                     "P7", "P3", "Pz", "P4", "P8",
                     "O1", "O2"))
  expect_false(any(c("Fpz", "Oz", "A1", "A2", "M1", "M2", "Nz", "T3", "T4", "T5",
                     "T6") %in% m$channels))
  expect_false(anyNA(m$positions))
})

# ----------------------------------------------------------------------------
# Test 9.2: Same shape and scale as the other layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies every electrode lies on the 87.54 mm head sphere,
# the positions table has the same columns as the other layouts, the
# lower-case names are the row names, and no value is a negative zero.
test_that("standard_1020 is built like the other layouts and sits on the head sphere", {
  m <- create_montage(template = "standard_1020")
  p <- m$positions

  expect_equal(sqrt(p$x^2 + p$y^2 + p$z^2), rep(87.54, 19), tolerance = 1e-6)
  expect_true(all(p$radius == 87.54))
  expect_identical(names(p), names(create_montage(template = "biosemi32")$positions))
  expect_identical(rownames(p), tolower(m$channels))

  numbers <- unlist(p[, c("x", "y", "z", "theta", "phi")])
  expect_false(any(numbers == 0 & 1 / numbers < 0))
})

# ----------------------------------------------------------------------------
# Test 9.3: The 10-20 electrodes are at the same places in the other layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 19 electrodes are in standard_1010 and in
# standard_1005 at exactly the same place: not close, identical, in every
# column of the positions table.
test_that("standard_1020 electrodes are at exactly the same places in the other standard layouts", {
  a <- standard_1020_positions()

  for (template in c("standard_1010", "standard_1005")) {
    b <- create_montage(a$channel, template = template)$positions
    expect_identical(b, a, label = template)
  }
})

# ----------------------------------------------------------------------------
# Test 9.4: Where the rules of the 10-20 system put the electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the electrodes against the definition of the
# system, with angles worked out by hand. The midline is cut at 30, 50 and 70
# percent of the nasion-inion arc (230 degrees: Fz 46, Cz 0, Pz 46 degrees
# from the vertex). The outer ring at 92 degrees from the vertex holds Fp, F7,
# T7, P7 and O1 at 36 degrees apart (10 percent of the circumference), and C3
# and C4 are halfway between Cz and T7 or T8. F3, F4, P3 and P4 are halfway
# between the midline electrode and the lateral electrode of their row, which
# is the same distance from both.
test_that("standard_1020 places the electrodes where the rules of the system put them", {
  p <- standard_1020_positions()

  midline <- c("Fz", "Cz", "Pz")
  expect_equal(angles_of(p, midline, "theta"), c(46, 0, 46))
  expect_equal(angles_of(p, midline, "phi"), c(90, 0, -90))

  left <- c("Fp1", "F7", "T7", "P7", "O1")
  right <- c("Fp2", "F8", "T8", "P8", "O2")
  expect_equal(angles_of(p, left, "theta"), rep(-92, 5))
  expect_equal(angles_of(p, left, "phi"), c(-72, -36, 0, 36, 72), tolerance = 1e-6)
  expect_equal(angles_of(p, right, "theta"), rep(92, 5))
  expect_equal(angles_of(p, right, "phi"), c(72, 36, 0, -36, -72), tolerance = 1e-6)

  expect_equal(angles_of(p, c("C3", "C4"), "theta"), c(-46, 46))
  expect_equal(angles_of(p, c("C3", "C4"), "phi"), c(0, 0), tolerance = 1e-6)

  xyz <- function(name) unlist(p[match(name, p$channel), c("x", "y", "z")])
  distance <- function(a, b) sqrt(sum((xyz(a) - xyz(b))^2))
  halfway <- list(c("F3", "F7", "Fz"), c("F4", "F8", "Fz"),
                  c("C3", "T7", "Cz"), c("C4", "T8", "Cz"),
                  c("P3", "P7", "Pz"), c("P4", "P8", "Pz"))
  for (h in halfway) {
    expect_equal(distance(h[1], h[2]), distance(h[1], h[3]), tolerance = 1e-5,
                 label = paste(h[1], "halfway between", h[2], "and", h[3]))
  }
})

# ----------------------------------------------------------------------------
# Test 9.5: The old names T3, T4, T5 and T6
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the three standard layouts accept the old names of
# T7, T8, P7 and P8 without a warning, in any letter case: the channel keeps
# the name it was asked for and gets the position of the present name. A
# layout still never lists the old names as its own channels, and the BioSemi
# layouts do not know them.
test_that("the standard layouts accept the old names T3 to T6", {
  columns <- c("x", "y", "z", "theta", "phi", "radius")

  for (template in c("standard_1020", "standard_1010", "standard_1005")) {
    expect_no_warning(old <- create_montage(c("T3", "t4", "T5", "T6", "Cz"),
                                            template = template))
    new <- create_montage(c("T7", "T8", "P7", "P8", "Cz"), template = template)

    expect_equal(old$channels, c("T3", "t4", "T5", "T6", "Cz"), label = template)
    expect_equal(unname(as.matrix(old$positions[, columns])),
                 unname(as.matrix(new$positions[, columns])), label = template)
  }

  expect_no_warning(both <- create_montage(c("T3", "T7"), template = "standard_1020"))
  expect_equal(both$positions$x[1], both$positions$x[2])

  expect_warning(create_montage(c("Cz", "T3"), template = "biosemi32"), "T3")
})

# ----------------------------------------------------------------------------
# Test 9.6: Left and right mirror each other
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the eight electrode pairs (Fp1/Fp2, F7/F8, F3/F4,
# T7/T8, C3/C4, P7/P8, P3/P4, O1/O2) are mirror images and Fz, Cz and Pz are
# on the midline.
test_that("standard_1020 is mirror symmetric", {
  p <- standard_1020_positions()
  left <- c("Fp1", "F7", "F3", "T7", "C3", "P7", "P3", "O1")
  right <- c("Fp2", "F8", "F4", "T8", "C4", "P8", "P4", "O2")

  a <- p[match(left, p$channel), c("x", "y", "z")]
  b <- p[match(right, p$channel), c("x", "y", "z")]
  expect_equal(a$x, -b$x, tolerance = 1e-6)
  expect_equal(a$y, b$y, tolerance = 1e-6)
  expect_equal(a$z, b$z, tolerance = 1e-6)
  expect_true(all(a$x < 0))
  expect_equal(p$x[match(c("Fz", "Cz", "Pz"), p$channel)], c(0, 0, 0))
})

# ----------------------------------------------------------------------------
# Test 9.7: No number of the table has been changed
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Guards the whole table. Each column is summed with the row
# number as weight, so a changed value, a swapped pair or a shifted row alters
# the sum. The expected sums were computed after the table had been checked
# against the rules (tests 9.1 to 9.6).
test_that("the standard_1020 table is unchanged", {
  p <- standard_1020_positions()
  i <- seq_len(nrow(p))
  sums <- vapply(c("x", "y", "z", "theta", "phi"),
                 function(col) sum(p[[col]] * i), numeric(1))

  expect_equal(sums, c(x = 1288.243486, y = -5680.0813, z = 4699.38632,
                       theta = 2544.039556, phi = -900), tolerance = 1e-8)
})

# ----------------------------------------------------------------------------
# Test 9.8: A recording that uses the old names
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the situation the old names exist for. A 19
# channel recording named the old way (T3, T4, T5, T6) gets its montage from
# create_montage() and set_montage() without any warning, the channels keep
# their names, and T3 sits where T7 sits.
test_that("a recording with the old channel names gets a standard_1020 montage", {
  chans <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "T3", "C3", "Cz", "C4",
             "T4", "T5", "P3", "Pz", "P4", "T6", "O1", "O2")
  eeg <- new_eeg(data = matrix(rnorm(19 * 20), nrow = 19), channels = chans,
                 sampling_rate = 100)

  expect_no_warning(m <- create_montage(chans, template = "standard_1020"))
  expect_no_warning(eeg <- set_montage(eeg, m))
  expect_equal(eeg$montage$channels, chans)

  new <- standard_1020_positions()
  at <- function(positions, name) unlist(positions[match(name, positions$channel), c("x", "y", "z")])
  expect_equal(at(eeg$montage$positions, "T3"), at(new, "T7"))
  expect_equal(at(eeg$montage$positions, "T6"), at(new, "P8"))
})

# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - new_montage(): basic construction, default coord_frame
# - create_montage(): default 64-electrode set, mixed/case-insensitive naming,
#   unmatched-channel warnings, all-unmatched error, unknown template error
# - set_montage(): scoping to eeg-type channels, input validation, partial
#   overlap warnings, no-overlap error, preprocessing_history update
# - print.montage(): formatted output, invisible return
# - BioSemi layouts (16, 32, 64, 128, 160, 256): electrode counts, head
#   sphere, positions against BioSemi's printed values, 160/256 angles,
#   socket names, 10-20 + socket aliases (16/32), agreement of the 16/32
#   layouts with the 64 layout, NA/10-20 names in the 128 layout, bad template
#   messages, the 64 layout unchanged, attaching to matching recordings,
#   whole-table fingerprints of the angles
# - read_custom_montage(): .elc (old and new style, any unit), .sfp, .tsv,
#   .txt and .csv files (header, quotes, byte-order mark), landmarks left out,
#   recentring (sphere and non-sphere), recenter = FALSE, theta/phi against
#   the BioSemi angles, the unit (stated, worked out, given, rejected), every
#   kind of unusable file, too few electrodes, the same object shape as
#   create_montage(), attaching with set_montage()
# - standard_1010: the 72 electrodes, the head sphere and table shape, the
#   10 percent electrodes (midline, ear line, outer ring, inferior row) against
#   angles worked out by hand, equal spacing along one circle for every row,
#   mirror symmetry, agreement with the BioSemi 64 layout (and the 16 known
#   exceptions), a fingerprint of the whole table, names outside the layout,
#   attaching with set_montage()
# - standard_1005: the 335 electrodes (written out row by row), the head
#   sphere and table shape, every standard_1010 electrode at exactly the same
#   place, the 5 percent electrodes (midline, outer ring, lower row, C row,
#   OI) against angles worked out by hand, equal spacing along one circle for
#   every row, the half electrodes outside the ear line and OI halfway between
#   their neighbours, new rows halfway between the rows around them (the
#   paper's CCP3h example and the whole grid), mirror symmetry, a fingerprint
#   of the whole table, picking electrodes by name and names outside the
#   layout, attaching with set_montage()
# - standard_1020: the 19 electrodes, the head sphere and table shape, the
#   same places as in the other two standard layouts (identical), the
#   defining angles worked out by hand and F3/C3/P3 (and the right side)
#   halfway on their row, the old names T3 to T6 in all three standard
#   layouts, mirror symmetry, a fingerprint of the whole table, a recording
#   that uses the old names
# ============================================================================
