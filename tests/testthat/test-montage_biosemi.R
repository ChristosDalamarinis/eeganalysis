# ============================================================================
#                  Test File for montage_biosemi.R
# ============================================================================
#
# This test file covers everything that lives in montage_biosemi.R: the layouts
# of the BioSemi caps (16, 32, 64, 128, 160 and 256 channels), the 64-channel
# layout in particular (the electrode database reads its names and positions
# from it), and the BioSemi channels that are not on the scalp (EXG1-EXG8,
# GSR1-GSR2, Erg1-Erg2, Resp, Plet, Temp and Status). The numbers are pinned,
# so that a value that is moved or edited by mistake is noticed.
#
# Parts tested:
#   1. The layout names and the front door - montage_templates.R puts the
#                                            families of layouts together
#   2. The 64-channel layout               - order, names, positions
#   3. The electrode database              - positions from the layout, and the
#                                            description of each electrode
#   4. The non-scalp channels              - what they are, where they are listed
#   5. The functions that use them         - create_montage(), the 3D plots,
#                                            get_electrode_position(),
#                                            detect_external_channels() and
#                                            detect_electrode_naming_system()
#
# Author: Christos Dalamarinis
# Date: November 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# The number of electrodes of every layout create_montage() can build, in the
# order of eeganalysis:::.montage_templates()
all_layouts <- c(biosemi16 = 16, biosemi32 = 32, biosemi64 = 64,
                 biosemi128 = 128, biosemi160 = 160, biosemi256 = 256,
                 standard_1020 = 19, standard_1010 = 72, standard_1005 = 335,
                 colin27_1020 = 94, colin27_1005 = 343)

layout_table   <- function(layout) eeganalysis:::.montage_template_table(layout)
layout_columns <- c("name", "alias", "x", "y", "z", "theta", "phi", "radius")

# The 64 electrodes of the BioSemi 64-channel cap, in BioSemi's order
biosemi64_names <- c(
  "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3", "FC1",
  "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7",
  "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz", "Fpz", "Fp2",
  "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4",
  "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2",
  "P2", "P4", "P6", "P8", "P10", "PO8", "PO4", "O2"
)

# What the database says about the position of each of them, in the same order
biosemi64_positions <- c(
  "Prefrontal left", "Anterior frontal left", "Anterior frontal left-center",
  "Frontal left-midline", "Frontal left", "Frontal left-center",
  "Frontal left temporal", "Fronto-temporal left", "Fronto-central left",
  "Fronto-central left-center", "Fronto-central left-midline",
  "Central left-midline", "Central left", "Central left", "Temporal left",
  "Centro-temporal left", "Centro-parietal left",
  "Centro-parietal left-center", "Centro-parietal left-midline",
  "Parietal left-midline", "Parietal left", "Parietal left-center",
  "Parietal left temporal", "Parietal left ear", "Parieto-occipital left",
  "Parieto-occipital left-center", "Occipital left", "Inion",
  "Occipital midline", "Parieto-occipital midline", "Parietal midline",
  "Centro-parietal midline", "Prefrontal midline", "Prefrontal right",
  "Anterior frontal right", "Anterior frontal right-center",
  "Anterior frontal midline", "Frontal midline", "Frontal right-midline",
  "Frontal right", "Frontal right-center", "Frontal right temporal",
  "Fronto-temporal right", "Fronto-central right",
  "Fronto-central right-center", "Fronto-central right-midline",
  "Fronto-central midline", "Central midline", "Central right-midline",
  "Central right", "Central right", "Temporal right", "Centro-temporal right",
  "Centro-parietal right", "Centro-parietal right-center",
  "Centro-parietal right-midline", "Parietal right-midline", "Parietal right",
  "Parietal right-center", "Parietal right temporal", "Parietal right ear",
  "Parieto-occipital right", "Parieto-occipital right-center",
  "Occipital right"
)

# The BioSemi channels that are not on the scalp
external_names <- c(paste0("EXG", 1:8), "GSR1", "GSR2", "Erg1", "Erg2", "Resp",
                    "Plet", "Temp", "Status")

# The 21 electrodes of the original 10-20 system among the 64
system_1020 <- c("Fp1", "Fp2", "Fpz", "F3", "F4", "F7", "F8", "Fz", "C3", "C4",
                 "Cz", "T7", "T8", "P3", "P4", "P7", "P8", "Pz", "O1", "O2",
                 "Oz")

# The fields of every entry of the electrode database
entry_fields <- c("position_name", "position_type", "region", "standard_system",
                  "biosemi_name", "standard_name", "cartesian_coords",
                  "spherical_coords")

# ============================================================================
#      TEST SUITE 1: The Layout Names and the Front Door
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: The names of the layouts
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the BioSemi caps are listed in order, and that the
# front door lists the BioSemi caps first, then the standard systems and the
# real-head layouts, as before.
test_that("1.1 the layouts are listed family by family", {
  expect_equal(eeganalysis:::.biosemi_templates(),
               c("biosemi16", "biosemi32", "biosemi64", "biosemi128",
                 "biosemi160", "biosemi256"))
  expect_equal(eeganalysis:::.montage_templates(),
               c(eeganalysis:::.biosemi_templates(),
                 eeganalysis:::.standard_templates(),
                 eeganalysis:::.colin27_templates()))
  expect_equal(eeganalysis:::.montage_templates(), names(all_layouts))
})

# ----------------------------------------------------------------------------
# Test 1.2: A table for every layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the front door finds the table of every layout,
# whichever file builds it: the same columns, the expected number of
# electrodes, a position for each, and no name twice.
test_that("1.2 the front door gives a table for every layout", {
  for (nm in names(all_layouts)) {
    tab <- layout_table(nm)
    expect_equal(names(tab), layout_columns, info = nm)
    expect_equal(nrow(tab), unname(all_layouts[nm]), info = nm)
    expect_false(any(is.na(as.matrix(tab[c("x", "y", "z", "theta", "phi",
                                           "radius")]))), info = nm)
    expect_equal(anyDuplicated(tab$name), 0L, info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 1.3: Each family builds its own tables
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the table the front door gives is the table of the
# file the layout belongs to, for the BioSemi caps, the standard systems and
# the real-head layouts.
test_that("1.3 the front door gives the table of the file of each layout", {
  for (nm in eeganalysis:::.biosemi_templates()) {
    expect_identical(layout_table(nm),
                     eeganalysis:::.biosemi_template_table(nm), info = nm)
  }
  for (nm in eeganalysis:::.standard_templates()) {
    expect_identical(layout_table(nm),
                     eeganalysis:::.standard_template_table(nm), info = nm)
  }
  for (nm in eeganalysis:::.colin27_templates()) {
    expect_identical(layout_table(nm),
                     eeganalysis:::.colin27_template_table(nm), info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 1.4: A name that is not a layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the front door stops with a message that lists the
# layouts there are, instead of an error from deep inside one of the families.
test_that("1.4 a name that is not a layout stops with the list of layouts", {
  expect_error(layout_table("nonsense"), "Unknown layout 'nonsense'")
  expect_error(layout_table("nonsense"), "biosemi64.*standard_1020.*colin27_1005")
})

# ============================================================================
#      TEST SUITE 2: The 64-Channel Layout
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Order, names and sockets
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the 64-channel layout has the 64 electrodes in
# BioSemi's order, each with its socket name (A1-A32, then B1-B32) as alias, and
# one head radius for all.
test_that("2.1 the 64-channel layout has BioSemi's order, names and sockets", {
  tab <- layout_table("biosemi64")

  expect_equal(nrow(tab), 64)
  expect_equal(tab$name, biosemi64_names)
  expect_equal(tab$alias, paste0(rep(c("A", "B"), each = 32), 1:32))
  expect_equal(tab$name[tab$alias == "B16"], "Cz")
  expect_equal(tab$radius, rep(87.54, 64))
  expect_identical(eeganalysis:::.biosemi64_table(), tab)
})

# ----------------------------------------------------------------------------
# Test 2.2: Fourteen electrodes, value by value
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies x, y, z (whole millimetres), theta and phi
# (degrees) of fourteen electrodes spread over the head against BioSemi's sheet.
test_that("2.2 the positions of fourteen electrodes are BioSemi's", {
  tab <- layout_table("biosemi64")
  pinned <- list(
    Fp1 = c(-27, 83, -3, -92, -72),
    AF7 = c(-51, 71, -3, -92, -54),
    F3  = c(-48, 59, 44, -60, -51),
    FC1 = c(-33, 33, 74, -32, -45),
    Cz  = c(0, 0, 88, 0, 0),
    T7  = c(-87, 0, -3, -92, 0),
    P9  = c(-64, -47, -37, -115, 36),
    Iz  = c(0, -79, -37, 115, -90),
    Oz  = c(0, -87, -3, 92, -90),
    O2  = c(27, -83, -3, 92, -72),
    Fpz = c(0, 87, -3, 92, 90),
    T8  = c(87, 0, -3, 92, 0),
    P10 = c(64, -47, -37, 115, -36),
    PO4 = c(36, -76, 24, 74, -65)
  )
  for (nm in names(pinned)) {
    row <- tab[tab$name == nm, c("x", "y", "z", "theta", "phi")]
    expect_equal(unlist(row, use.names = FALSE), pinned[[nm]], info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 2.3: All 64 electrodes agree with themselves
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies for every electrode, not only the pinned ones, that
# x, y and z are whole millimetres that agree with the angles to within the
# rounding (half a millimetre), and that the left, right and midline electrodes
# are on the side their names say (odd numbers left, even numbers right, names
# that end in z on the midline).
test_that("2.3 the 64 positions agree with their angles and their sides", {
  tab <- layout_table("biosemi64")
  th <- tab$theta * pi / 180
  ph <- tab$phi * pi / 180

  expect_equal(tab$x, round(tab$x))
  expect_equal(tab$y, round(tab$y))
  expect_equal(tab$z, round(tab$z))
  expect_lt(max(abs(tab$x - 87.54 * sin(th) * cos(ph))), 0.5)
  expect_lt(max(abs(tab$y - 87.54 * sin(th) * sin(ph))), 0.5)
  expect_lt(max(abs(tab$z - 87.54 * cos(th))), 0.5)

  midline <- grepl("z$", tab$name)
  odd     <- as.integer(sub("^[A-Za-z]+", "", tab$name)) %% 2 == 1
  side    <- ifelse(midline, 0, ifelse(odd, -1, 1))
  expect_equal(sign(tab$x), side)
  expect_false(any(tab$theta < 0 & tab$x > 0))
  expect_false(any(tab$theta > 0 & tab$x < 0))
})

# ============================================================================
#      TEST SUITE 3: The Electrode Database
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Positions come from the layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the database has 288 keys (each of the 64 scalp
# electrodes under four keys, each of the 16 other channels under two), and
# that every scalp electrode is the same entry under its standard name, its
# BioSemi name and both in lower case, with the position of the layout row.
test_that("3.1 the database holds the 64-channel layout under four keys each", {
  db  <- get_electrode_database()
  tab <- layout_table("biosemi64")

  expect_equal(length(db), 64 * 4 + 16 * 2)
  for (i in seq_len(nrow(tab))) {
    entry <- db[[tab$name[i]]]
    expect_equal(c(entry$cartesian_coords$x, entry$cartesian_coords$y,
                   entry$cartesian_coords$z, entry$spherical_coords$theta,
                   entry$spherical_coords$phi, entry$spherical_coords$radius),
                 c(tab$x[i], tab$y[i], tab$z[i], tab$theta[i], tab$phi[i],
                   tab$radius[i]), info = tab$name[i])
    expect_equal(entry$standard_name, tab$name[i], info = tab$name[i])
    expect_equal(entry$biosemi_name, tab$alias[i], info = tab$name[i])
    expect_identical(db[[tab$alias[i]]], entry)
    expect_identical(db[[tolower(tab$name[i])]], entry)
    expect_identical(db[[tolower(tab$alias[i])]], entry)
  }
})

# ----------------------------------------------------------------------------
# Test 3.2: The words for each position
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the name of the position of all 64 electrodes, one
# by one, in the order of the layout. The descriptions are listed in the
# database in the order of the layout, so a description that moved to another
# electrode would show here.
test_that("3.2 every electrode has the right name for its position", {
  db  <- get_electrode_database()
  tab <- layout_table("biosemi64")

  said <- vapply(tab$name, function(n) db[[n]]$position_name, character(1))
  expect_equal(unname(said), biosemi64_positions)
})

# ----------------------------------------------------------------------------
# Test 3.3: Type, region and system
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the type, the region and the system (10-20 or
# 10-10) of all 64 electrodes. The type and the region follow from the letters
# of the name (Fp, AF, F, FT, FC, C, T, TP, CP, P, PO, O, I), and exactly 21
# electrodes belong to the 10-20 system.
test_that("3.3 type, region and system follow from the name of the electrode", {
  db  <- get_electrode_database()
  tab <- layout_table("biosemi64")

  type_of   <- c(Fp = "Prefrontal", AF = "Anterior Frontal", F = "Frontal",
                 FT = "Fronto-Temporal", FC = "Fronto-Central",
                 C = "Central", T = "Temporal", TP = "Centro-Temporal",
                 CP = "Centro-Parietal", P = "Parietal",
                 PO = "Parieto-Occipital", O = "Occipital", I = "Occipital")
  region_of <- c(Fp = "Frontal", AF = "Frontal", F = "Frontal",
                 FT = "Temporal", FC = "Central", C = "Central",
                 T = "Temporal", TP = "Temporal", CP = "Parietal",
                 P = "Parietal", PO = "Occipital", O = "Occipital",
                 I = "Occipital")
  prefix <- sub("[0-9z]+$", "", tab$name)

  got <- function(field) {
    vapply(tab$name, function(n) db[[n]][[field]], character(1), USE.NAMES = FALSE)
  }
  expect_equal(got("position_type"), unname(type_of[prefix]))
  expect_equal(got("region"), unname(region_of[prefix]))
  expect_equal(got("standard_system"),
               ifelse(tab$name %in% system_1020, "10-20", "10-10"))
  expect_equal(sum(got("standard_system") == "10-20"), 21)
})

# ----------------------------------------------------------------------------
# Test 3.4: Left, right and midline
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the description of each electrode says "left" when
# the electrode is on the left of the head, "right" when on the right, and
# neither on the midline, by the sign of its x in the layout.
test_that("3.4 the descriptions agree with the side of the layout", {
  db  <- get_electrode_database()
  tab <- layout_table("biosemi64")

  said  <- vapply(tab$name, function(n) db[[n]]$position_name, character(1))
  words <- ifelse(grepl("left", said), -1, ifelse(grepl("right", said), 1, 0))
  expect_equal(unname(words), sign(tab$x))
})

# ============================================================================
#      TEST SUITE 4: The Non-Scalp Channels
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: The list and its entries
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the list has the 16 channels in order, each entry
# has the fields of a database entry, belongs to BioSemi, has no position (all
# of it NA) and carries its own name as standard name.
test_that("4.1 the 16 non-scalp channels are listed with no position", {
  ext <- eeganalysis:::.biosemi_external_channels()

  expect_equal(names(ext), external_names)
  for (nm in external_names) {
    entry <- ext[[nm]]
    expect_equal(names(entry), entry_fields, info = nm)
    expect_equal(entry$standard_name, nm, info = nm)
    expect_equal(entry$standard_system, "BioSemi", info = nm)
    expect_true(is.na(entry$biosemi_name), info = nm)
    expect_equal(names(entry$cartesian_coords), c("x", "y", "z"), info = nm)
    expect_equal(names(entry$spherical_coords), c("theta", "phi", "radius"),
                 info = nm)
    expect_true(all(is.na(unlist(entry$cartesian_coords))), info = nm)
    expect_true(all(is.na(unlist(entry$spherical_coords))), info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 4.2: What each channel is
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the type, the region and the description of each of
# the 16 channels. detect_external_channels() decides by the type.
test_that("4.2 each non-scalp channel has its type, region and description", {
  ext <- eeganalysis:::.biosemi_external_channels()
  field <- function(f) vapply(ext, function(e) e[[f]], character(1), USE.NAMES = FALSE)

  expect_equal(field("position_type"),
               c(rep("External", 8), "GSR", "GSR", "Ergo/AUX", "Ergo/AUX",
                 "Respiration", "Plethysmograph", "Temperature", "Status"))
  expect_equal(field("region"), c(rep("External", 15), "Metadata"))
  expect_equal(field("position_name"),
               c(paste0("External electrode ", 1:8, " (EOG/ECG/EMG/GSR)"),
                 "Galvanic skin response channel 1 (electrodermal activity/arousal)",
                 "Galvanic skin response channel 2 (electrodermal activity/arousal)",
                 "Ergo/AUX input 1 (Response Switch/Photocell/Microphone)",
                 "Ergo/AUX input 2 (Response Switch/Photocell/Microphone)",
                 "Respiration belt (breathing rate/pattern measurement)",
                 "Plethysmograph (blood volume/pulse measurement)",
                 "Temperature sensor (body temperature measurement)",
                 "Status channel (event markers/triggers)"))
})

# ----------------------------------------------------------------------------
# Test 4.3: In the database
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the database holds each of the 16 channels, as
# listed, under its name and under its name in lower case.
test_that("4.3 the database holds the non-scalp channels under two keys each", {
  db  <- get_electrode_database()
  ext <- eeganalysis:::.biosemi_external_channels()

  for (nm in external_names) {
    expect_identical(db[[nm]], ext[[nm]], info = nm)
    expect_identical(db[[tolower(nm)]], ext[[nm]], info = nm)
  }
})

# ----------------------------------------------------------------------------
# Test 4.4: Found as external channels
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies detect_external_channels() finds the 15 external
# channels of the list through the database, in any letter case, and leaves
# the Status channel (event markers, not a signal) out.
test_that("4.4 detect_external_channels() finds them through the database", {
  expect_equal(detect_external_channels(external_names),
               setdiff(external_names, "Status"))
  expect_equal(detect_external_channels(tolower(external_names)),
               setdiff(tolower(external_names), "status"))
  expect_equal(detect_external_channels(c("Cz", "EXG1", "Fp1", "Temp")),
               c("EXG1", "Temp"))
})

# ----------------------------------------------------------------------------
# Test 4.5: In no layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies none of the 16 channels is a row of any layout
# table, because they have no position.
test_that("4.5 no layout has a row for a non-scalp channel", {
  for (nm in names(all_layouts)) {
    expect_false(any(external_names %in% layout_table(nm)$name), info = nm)
  }
})

# ============================================================================
#      TEST SUITE 5: The Functions That Use Them
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: create_montage()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies create_montage() gives the positions of the
# 64-channel layout, in its order, by the standard name and by the BioSemi
# name, and drops a non-scalp channel with a warning.
test_that("5.1 create_montage() gives the positions of the 64-channel layout", {
  tab <- layout_table("biosemi64")
  m <- create_montage()

  expect_equal(m$channels, tab$name)
  expect_equal(m$positions$x, tab$x)
  expect_equal(m$positions$y, tab$y)
  expect_equal(m$positions$z, tab$z)
  expect_equal(m$positions$theta, tab$theta)
  expect_equal(m$positions$phi, tab$phi)
  expect_equal(m$positions$radius, tab$radius)

  by_socket <- create_montage(c("A1", "B16", "B32"))
  expect_equal(by_socket$positions$x, tab$x[match(c("A1", "B16", "B32"), tab$alias)])
  expect_warning(create_montage(c("Cz", "EXG1")), "EXG1")
})

# ----------------------------------------------------------------------------
# Test 5.2: The 3D plots
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the table the two 3D electrode plots draw from is
# the 64-channel layout, with the BioSemi names as aliases.
test_that("5.2 the 3D plots draw the 64-channel layout", {
  tab  <- layout_table("biosemi64")
  plot <- eeganalysis:::.plot_layout_table("biosemi64")

  expect_equal(plot$electrode, tab$name)
  expect_equal(plot$alias, tab$alias)
  expect_equal(plot$x, tab$x)
  expect_equal(plot$y, tab$y)
  expect_equal(plot$z, tab$z)
  expect_equal(plot$theta, tab$theta)
  expect_equal(plot$phi, tab$phi)
})

# ----------------------------------------------------------------------------
# Test 5.3: get_electrode_position()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies get_electrode_position() finds a scalp electrode by
# either of its names and in lower case, finds a non-scalp channel, and says
# "Unknown" for a name that is in neither list.
test_that("5.3 get_electrode_position() finds scalp and non-scalp channels", {
  expect_equal(get_electrode_position("B16")$standard_name, "Cz")
  expect_equal(get_electrode_position("cz")$position_name, "Central midline")
  expect_equal(get_electrode_position("A24")$standard_name, "P9")
  expect_equal(get_electrode_position("exg4")$position_type, "External")
  expect_equal(get_electrode_position("Status")$region, "Metadata")
  expect_equal(get_electrode_position("Ch99")$position_type, "Unknown")
})

# ----------------------------------------------------------------------------
# Test 5.4: detect_electrode_naming_system()
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the naming detector recognises the 64 BioSemi
# sockets, and gives the standard name behind each of them in the order of the
# layout (it gets them from the database), and still recognises the standard
# names when the non-scalp channels are in the list too.
test_that("5.4 detect_electrode_naming_system() reads the 64-channel layout", {
  sockets <- paste0(rep(c("A", "B"), each = 32), 1:32)

  by_socket <- detect_electrode_naming_system(sockets)
  expect_equal(by_socket$naming_system, "BioSemi_AB")
  expect_equal(by_socket$channels_detected$standard_equivalent, biosemi64_names)

  by_name <- detect_electrode_naming_system(c(biosemi64_names, external_names))
  expect_equal(by_name$naming_system, "10-20/10-10")
})
