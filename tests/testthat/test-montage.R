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
# ============================================================================
