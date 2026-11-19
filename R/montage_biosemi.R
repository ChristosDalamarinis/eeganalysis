# ============================================================================
#                BioSemi layouts and channels used by create_montage()
# ============================================================================
#
# Everything that belongs to BioSemi is kept in this file:
#
#   * The electrode layouts of the BioSemi caps with 16, 32, 64, 128, 160 and
#     256 channels, which create_montage() can build (.biosemi_templates()).
#   * The BioSemi channels that are not on the scalp (EXG1-EXG8, GSR1-GSR2,
#     Erg1-Erg2, Resp, Plet, Temp and Status). They have no position, so they
#     are in no layout; they are listed in .biosemi_external_channels().
#
# montage_templates.R is the front door to all layouts (this one, the standard
# layouts of montage_standard.R and the real-head layouts of montage_colin27.R).
# get_electrode_database() (R/channel_info2.R) reads the 64 channel layout and
# the list of non-scalp channels from here, and adds its own descriptions.
#
# BioSemi publishes the electrode positions of its caps with 16, 32, 64, 128,
# 160 and 256 channels in the workbook "Cap_coords_all.xls" (one sheet per
# cap, https://www.biosemi.com/download.htm). Each electrode is given by two
# angles in degrees, with the head treated as a sphere:
#
#   theta  inclination: 0 at the vertex (Cz), then growing towards the ears
#          and below; negative on the left side of the head, positive on the
#          right (range -180 to 180).
#   phi    azimuth, measured from the ear line (T8 for positive inclinations,
#          T7 for negative ones); positive is anti-clockwise seen from above
#          (range -90 to 90).
#
# From these two angles BioSemi's own sheet derives the Cartesian position,
# and so do we (x towards the right ear, y towards the nose, z towards the top
# of the head):
#
#   x = r * sin(theta) * cos(phi)
#   y = r * sin(theta) * sin(phi)
#   z = r * cos(theta)
#
# r is the head radius, 87.54 mm (a 55 cm head circumference divided by 2*pi).
# It is kept in .BIOSEMI_RADIUS_MM (montage_templates.R), because the standard
# layouts of montage_standard.R use the same radius, so that all layouts built
# on a sphere share one scale.
#
# The angles of the 16, 32, 128 and 160 channel caps are copied as printed in
# BioSemi's sheets; they reproduce every Cartesian value BioSemi prints next
# to them (z to six decimals, x and y to the printed whole millimetre).
#
# The 256-channel sheet is the exception: it prints its angles rounded to whole
# degrees, which do not reproduce its own Cartesian columns (z is off by up to
# 0.6 mm). The 256 angles used here keep the extra digits (for example 9.2
# instead of 9 degrees), and with them every printed z is reproduced to six
# decimals and every printed x and y to the whole millimetre.
#
# The 64-channel cap is written out in .biosemi64_table() exactly as BioSemi's
# sheet prints it: the angles, and also x, y and z, which are whole
# millimetres (Cz has z = 88 although the radius is 87.54), so they are not
# worked out from the angles like those of the other caps.
#
# Names:
#   16 and 32 channels   standard 10-20 names (Fp1, Cz, ...). Row i of the sheet
#                        is BioSemi socket A<i>, so A1, A2, ... are accepted
#                        as well (the 64-channel table numbers its rows the
#                        same way).
#   64 channels          standard 10-10 names (Fp1, AF7, ..., O2) in BioSemi's
#                        order; the sockets A1-A32 and B1-B32 are the alias
#                        (Cz is B16).
#   128, 160, 256        BioSemi's own names: A1-D32, A1-E32 and A1-H32 (32
#                        sockets per letter). The sheets also print the
#                        nearest 10-20 name next to some of them; those are
#                        only approximate, so they are not used.

# Names of the layouts built here.
.biosemi_templates <- function() {
  c("biosemi16", "biosemi32", "biosemi64", "biosemi128", "biosemi160",
    "biosemi256")
}

# One row per electrode for one BioSemi cap: name, alias (a second name that is
# accepted, NA if none), x, y, z (mm), theta, phi (degrees) and radius (mm).
.biosemi_template_table <- function(template) {
  if (identical(template, "biosemi64")) {
    return(.biosemi64_table())
  }

  cap <- .biosemi_cap_angles(template)
  th  <- cap$theta * pi / 180
  ph  <- cap$phi * pi / 180

  data.frame(
    name   = cap$name,
    alias  = cap$alias,
    x      = round(.BIOSEMI_RADIUS_MM * sin(th) * cos(ph), 6),
    y      = round(.BIOSEMI_RADIUS_MM * sin(th) * sin(ph), 6),
    z      = round(.BIOSEMI_RADIUS_MM * cos(th), 6),
    theta  = cap$theta,
    phi    = cap$phi,
    radius = .BIOSEMI_RADIUS_MM,
    stringsAsFactors = FALSE
  )
}

# The 64-channel cap, in BioSemi's order: the standard 10-10 name and the
# BioSemi socket name (A1-A32, B1-B32) of each electrode, its Cartesian
# position in whole millimetres, and its angles in degrees (BioSemi's
# Cap_coords_all.xls; head circumference 55 centimeters is radius 87.54
# milimeters).
.biosemi64_table <- function() {
  data.frame(
    name = c(
      "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3",
      "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1",
      "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz",
      "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4",
      "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4",
      "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8",
      "P10", "PO8", "PO4", "O2"
    ),
    alias = c(
      "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
      "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20",
      "A21", "A22", "A23", "A24", "A25", "A26", "A27", "A28", "A29", "A30",
      "A31", "A32", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
      "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18",
      "B19", "B20", "B21", "B22", "B23", "B24", "B25", "B26", "B27", "B28",
      "B29", "B30", "B31", "B32"
    ),
    # Cartesian coordinates (x, y, z) in mm
    x = c(
      -27, -51, -36, -25, -48, -64, -71, -83, -78, -59,
      -33, -34, -63, -82, -87, -83, -78, -59, -33, -25,
      -48, -64, -71, -64, -51, -36, -27, 0, 0, 0,
      0, 0, 0, 27, 51, 36, 0, 0, 25, 48,
      64, 71, 83, 78, 59, 33, 0, 0, 34, 63,
      82, 87, 83, 78, 59, 33, 25, 48, 64, 71,
      64, 51, 36, 27
    ),
    y = c(
      83, 71, 76, 62, 59, 55, 51, 27, 30, 31,
      33, 0, 0, 0, 0, -27, -30, -31, -33, -62,
      -59, -55, -51, -47, -71, -76, -83, -79, -87, -82,
      -63, -34, 87, 83, 71, 76, 82, 63, 62, 59,
      55, 51, 27, 30, 31, 33, 34, 0, 0, 0,
      0, 0, -27, -30, -31, -33, -62, -59, -55, -51,
      -47, -71, -76, -83
    ),
    z = c(
      -3, -3, 24, 56, 44, 23, -3, -3, 27, 56,
      74, 81, 61, 31, -3, -3, 27, 56, 74, 56,
      44, 23, -3, -37, -3, 24, -3, -37, -3, 31,
      61, 81, -3, -3, -3, 24, 31, 61, 56, 44,
      23, -3, -3, 27, 56, 74, 81, 88, 81, 61,
      31, -3, -3, 27, 56, 74, 56, 44, 23, -3,
      -37, -3, 24, -3
    ),
    # Spherical coordinates (theta, phi) in degrees
    theta = c(
      -92, -92, -74, -50, -60, -75, -92, -92, -72, -50,
      -32, -23, -46, -69, -92, -92, -72, -50, -32, -50,
      -60, -75, -92, -115, -92, -74, -92, 115, 92, 69,
      46, 23, 92, 92, 92, 74, 69, 46, 50, 60,
      75, 92, 92, 72, 50, 32, 23, 0, 23, 46,
      69, 92, 92, 72, 50, 32, 50, 60, 75, 92,
      115, 92, 74, 92
    ),
    phi = c(
      -72, -54, -65, -68, -51, -41, -36, -18, -21, -28,
      -45, 0, 0, 0, 0, 18, 21, 28, 45, 68,
      51, 41, 36, 36, 54, 65, 72, -90, -90, -90,
      -90, -90, 90, 72, 54, 65, 90, 90, 68, 51,
      41, 36, 18, 21, 28, 45, 90, 0, 0, 0,
      0, 0, -18, -21, -28, -45, -68, -51, -41, -36,
      -36, -54, -65, -72
    ),
    radius = .BIOSEMI_RADIUS_MM,
    stringsAsFactors = FALSE
  )
}

# Names and angles of the 16, 32, 128, 160 and 256 channel caps, in the row
# order of BioSemi's sheets.
.biosemi_cap_angles <- function(template) {
  switch(
    template,
    biosemi16 = list(
      name = c(
        "Fp1", "Fp2", "F4", "Fz", "F3", "T7", "C3", "Cz", "C4", "T8", "P4",
        "Pz", "P3", "O1", "Oz", "O2"
      ),
      alias = paste0("A", 1:16),
      theta = c(
        -92, 92, 60, 46, -60, -92, -46, 0, 46, 92, 60, 46, -60, -92, 92, 92
      ),
      phi = c(
        -72, 72, 51, 90, -51, 0, 0, 0, 0, 0, -51, -90, 51, 72, -90, -72
      )
    ),
    biosemi32 = list(
      name = c(
        "Fp1", "AF3", "F7", "F3", "FC1", "FC5", "T7", "C3", "CP1", "CP5",
        "P7", "P3", "Pz", "PO3", "O1", "Oz", "O2", "PO4", "P4", "P8", "CP6",
        "CP2", "C4", "T8", "FC6", "FC2", "F4", "F8", "AF4", "Fp2", "Fz", "Cz"
      ),
      alias = paste0("A", 1:32),
      theta = c(
        -92, -74, -92, -60, -32, -72, -92, -46, -32, -72, -92, -60, 46, -74,
        -92, 92, 92, 74, 60, 92, 72, 32, 46, 92, 72, 32, 60, 92, 74, 92, 46,
        0
      ),
      phi = c(
        -72, -65, -36, -51, -45, -21, 0, 0, 45, 21, 36, 51, -90, 65, 72, -90,
        -72, -65, -51, -36, -21, -45, 0, 0, 21, 45, 51, 36, 65, 72, 90, 0
      )
    ),
    biosemi128 = list(
      name  = paste0(rep(LETTERS[1:4], each = 32), 1:32),
      alias = rep(NA_character_, 128),
      theta = c(
        0, 11.5, 23, 34.5, -46, -46, -57.5, -69, -80.5, -92, -103.5, -115,
        -115, -103.5, -92, -80.5, -69, -57.5, 46, 57.5, 69, 80.5, 92, 103.5,
        115, 115, 103.5, 92, 80.5, 69, 57.5, 46, 11.5, 23, 46, 57.5, 69,
        80.5, 92, 103.5, 115, 103.5, 92, 80.5, 69, 92, 80.5, 69, 57.5, 46,
        34.5, 23, 34.5, 46, 57.5, 69, 80.5, 92, 92, 80.5, 69, 57.5, 46, 34.5,
        11.5, 23, 46, 57.5, 69, 80.5, 92, 92, 80.5, 69, 34.5, 46, 57.5, 69,
        80.5, 92, 92, 80.5, 69, 57.5, 46, 34.5, 23, -34.5, -46, -57.5, -69,
        -80.5, -92, -92, -80.5, -69, -11.5, -23, -46, -57.5, -69, -80.5, -92,
        -92, -80.5, -69, -57.5, -46, -34.5, -23, -11.5, -23, -34.5, -34.5,
        -46, -57.5, -69, -80.5, -92, -92, -80.5, -69, -57.5, -46, -69, -80.5,
        -92, -103.5
      ),
      phi = c(
        0, -90, -90, -90, 67.5, 45, 45, 54, 54, 54, 54, 54, 72, 72, 72, 72,
        72, 67.5, -90, -90, -90, -90, -90, -90, -90, -72, -72, -72, -72, -72,
        -67.5, -67.5, -18, -45, -45, -45, -54, -54, -54, -54, -54, -36, -36,
        -36, -36, -18, -18, -18, -22.5, -22.5, -30, 0, 0, 0, 0, 0, 0, 0, 18,
        18, 18, 22.5, 22.5, 30, 54, 45, 45, 45, 36, 36, 36, 54, 54, 54, 60,
        67.5, 67.5, 72, 72, 72, 90, 90, 90, 90, 90, 90, 90, -60, -67.5,
        -67.5, -72, -72, -72, -54, -54, -54, -54, -45, -45, -45, -36, -36,
        -36, -18, -18, -18, -22.5, -22.5, -30, 0, 18, 45, 30, 0, 0, 0, 0, 0,
        0, 18, 18, 18, 22.5, 22.5, 36, 36, 36, 36
      )
    ),
    biosemi160 = list(
      name  = paste0(rep(LETTERS[1:5], each = 32), 1:32),
      alias = rep(NA_character_, 160),
      theta = c(
        0, 11.5, 23, 34.5, -46, -46, -57.5, -69, -80.5, -92, -103.5, -115,
        -115, -103.5, -92, -80.5, -69, -57.5, 46, 57.5, 69, 80.5, 92, 103.5,
        115, 115, 103.5, 92, 80.5, 69, 57.5, 46, 11.5, 23, 46, 57.5, 69,
        80.5, 92, 103.5, 115, 115, 103.5, 92, 80.5, 69, 69, 80.5, 92, 103.5,
        92, 80.5, 69, 57.5, 46, 34.5, 23, 34.5, 46, 57.5, 57.5, 69, 80.5, 92,
        11.5, 23, 23, 34.5, 34.5, 46, 46, 57.5, 69, 80.5, 92, 92, 80.5, 69,
        69, 80.5, 92, 92, 80.5, 69, 57.5, 57.5, 46, 34.5, 23, 34.5, 46, 46,
        57.5, 69, 80.5, 92, -11.5, 23, 34.5, 46, 57.5, 69, 80.5, 92, -92,
        -80.5, -69, -57.5, -46, -34.5, -23, -23, -34.5, -46, -57.5, -69,
        -80.5, -92, -92, -80.5, -69, -57.5, -46, -34.5, -46, -69, -80.5, -92,
        -11.5, -23, -34.5, -46, -57.5, -69, -80.5, -92, -92, -80.5, -69,
        -57.5, -57.5, -69, -80.5, -92, -103.5, -92, -80.5, -69, -46, -34.5,
        -23, -23, -34.5, -46, -57.5, -69, -80.5, -92, -103.5, -115
      ),
      phi = c(
        0, -90, -90, -90, 72, 54, 54, 60, 60, 60, 60, 60, 75, 75, 75, 75, 75,
        72, -90, -90, -90, -90, -90, -90, -90, -75, -75, -75, -75, -75, -72,
        -72, -18, -60, -54, -54, -60, -60, -60, -60, -60, -45, -45, -45, -45,
        -45, -30, -30, -30, -30, -15, -15, -15, -36, -36, -45, -30, -22.5,
        -18, -18, 0, 0, 0, 0, 54, 30, 0, 0, 22.5, 18, 0, 18, 15, 15, 15, 30,
        30, 30, 45, 45, 45, 60, 60, 60, 54, 36, 36, 45, 60, 67.5, 54, 72, 72,
        75, 75, 75, -54, 90, 90, 90, 90, 90, 90, 90, -75, -75, -75, -72, -72,
        -67.5, -60, -30, -45, -54, -54, -60, -60, -60, -45, -45, -45, -36,
        -36, -22.5, -18, -30, -30, -30, 18, 0, 0, 0, -18, -15, -15, -15, 0,
        0, 0, 0, 18, 15, 15, 15, 30, 30, 30, 30, 18, 22.5, 30, 60, 45, 36,
        36, 45, 45, 45, 45, 45
      )
    ),
    biosemi256 = list(
      name  = paste0(rep(LETTERS[1:8], each = 32), 1:32),
      alias = rep(NA_character_, 256),
      theta = c(
        0, 9.2, 18.4, 27.6, 36.8, 46, -46, -55.2, -64.4, -73.6, -82.8, -92,
        -101.2, -110.4, -119.6, 119.6, 110.4, 101.2, 92, 82.8, 73.6, 64.4,
        55.2, 46, 55.2, 64.4, 73.6, 82.8, 92, 101.2, 110.4, 119.6, 18.4,
        27.6, 36.8, 46, 55.2, 64.4, 73.6, 82.8, 92, 101.2, 110.4, 119.6,
        110.4, 101.2, 92, 82.8, 73.6, 64.4, 55.2, 46, 27.6, 36.8, 46, 55.2,
        64.4, 73.6, 82.8, 92, 101.2, 110.4, 110.4, 101.2, 9.2, 18.4, 27.6,
        36.8, 46, 55.2, 64.4, 64.4, 73.6, 82.8, 92, 101.2, 92, 82.8, 73.6,
        82.8, 92, 92, 82.8, 73.6, 64.4, 55.2, 46, 36.8, 27.6, 36.8, 46, 55.2,
        64.4, 73.6, 82.8, 92, 9.2, 18.4, 27.6, 36.8, 46, 64.4, 73.6, 82.8,
        92, 101.2, 101.2, 92, 82.8, 73.6, 64.4, 55.2, 46, 55.2, 64.4, 73.6,
        82.8, 92, 101.2, 92, 82.8, 73.6, 64.4, 55.2, 64.4, 73.6, 82.8, 92,
        18.4, 18.4, 27.6, 36.8, 46, 46, 55.2, 64.4, 73.6, 82.8, 92, 92, 82.8,
        73.6, 64.4, 55.2, 46, 36.8, 36.8, 27.6, -27.6, -36.8, -46, -55.2,
        -64.4, -73.6, -82.8, -92, -92, -82.8, -73.6, -64.4, -9.2, -18.4,
        -27.6, -36.8, -46, -55.2, -64.4, -73.6, -82.8, -92, -101.2, -92,
        -82.8, -73.6, -64.4, -55.2, -46, -36.8, -27.6, -18.4, -27.6, -36.8,
        -46, -55.2, -64.4, -73.6, -82.8, -92, -101.2, -101.2, -92, -82.8,
        -9.2, -18.4, -27.6, -36.8, -46, -55.2, -64.4, -73.6, -82.8, -92, -92,
        -82.8, -73.6, -64.4, -55.2, -46, -55.2, -64.4, -73.6, -82.8, -92,
        -101.2, -92, -82.8, -73.6, -64.4, -64.4, -73.6, -82.8, -92, -101.2,
        -110.4, -18.4, -27.6, -36.8, -36.8, -46, -46, -55.2, -64.4, -73.6,
        -82.8, -92, -101.2, -110.4, -110.4, -101.2, -92, -82.8, -73.6, -64.4,
        -55.2, -46, -36.8, -27.6, -46, -55.2, -64.4, -73.6, -82.8, -92,
        -101.2, -110.4, -119.6
      ),
      phi = c(
        0, -90, -90, -90, -90, -90, 75, 75, 78, 78, 78.75, 78.75, 78.75, 78,
        78, -90, -90, -90, -90, -90, -90, -90, -90, -75, -75, -78, -78,
        -78.75, -78.75, -78.75, -78, -78, -54, -66, -54, -60, -60, -66, -66,
        -67.5, -67.5, -67.5, -66, -66, -54, -56.25, -56.25, -56.25, -54, -54,
        -45, -45, -42, -36, -30, -30, -42, -42, -45, -45, -45, -42, -30,
        -33.75, -18, -18, -18, -18, -15, -15, -18, -30, -30, -33.75, -33.75,
        -22.5, -22.5, -22.5, -18, -11.25, -11.25, 0, 0, -6, -6, 0, 0, 0, 6,
        18, 15, 15, 6, 6, 11.25, 11.25, 54, 18, 30, 36, 30, 18, 18, 22.5,
        22.5, 22.5, 33.75, 33.75, 33.75, 30, 30, 30, 45, 45, 42, 42, 45, 45,
        45, 56.25, 56.25, 54, 54, 60, 66, 66, 67.5, 67.5, 90, 54, 54, 54, 60,
        75, 75, 78, 78, 78.75, 78.75, 90, 90, 90, 90, 90, 90, 90, 72, 78,
        -78, -72, -75, -75, -78, -78, -78.75, -78.75, -67.5, -67.5, -66, -66,
        -54, -54, -54, -54, -60, -60, -54, -54, -56.25, -56.25, -45, -45,
        -45, -42, -42, -45, -45, -36, -30, -18, -6, -18, -30, -30, -30, -30,
        -33.75, -33.75, -33.75, -22.5, -22.5, -22.5, 18, 18, 18, 0, -15, -15,
        -18, -18, -11.25, -11.25, 0, 0, -6, -6, 0, 0, 15, 6, 6, 11.25, 11.25,
        22.5, 22.5, 22.5, 18, 18, 30, 30, 33.75, 33.75, 33.75, 30, 54, 42,
        36, 18, 15, 30, 30, 42, 42, 45, 45, 45, 42, 54, 56.25, 56.25, 56.25,
        54, 54, 45, 45, 54, 66, 60, 60, 66, 66, 67.5, 67.5, 67.5, 66, 66
      )
    )
  )
}

# The BioSemi channels that are not on the scalp: the external electrodes
# EXG1-EXG8 (used for EOG, ECG, EMG ...), the accessories GSR1-GSR2 (skin
# response), Erg1-Erg2 (ergo/AUX inputs), Resp (respiration belt), Plet
# (plethysmograph) and Temp (temperature sensor), and the Status channel (event
# markers). They have no position (x, y, z, theta, phi and radius are NA), so
# they are not rows of any layout table above; get_electrode_database() adds
# them to the database and detect_external_channels() finds them there.
.biosemi_external_channels <- function() {
  list(
    # ===== EXG1-EXG8 - External Electrodes =====
    "EXG1" = list(
      position_name = "External electrode 1 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG1",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG2" = list(
      position_name = "External electrode 2 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG2",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG3" = list(
      position_name = "External electrode 3 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG3",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG4" = list(
      position_name = "External electrode 4 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG4",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG5" = list(
      position_name = "External electrode 5 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG5",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG6" = list(
      position_name = "External electrode 6 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG6",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG7" = list(
      position_name = "External electrode 7 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG7",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "EXG8" = list(
      position_name = "External electrode 8 (EOG/ECG/EMG/GSR)",
      position_type = "External",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "EXG8",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== BIOSEMI ACCESSORIES - GSR1/GSR2 =====
    "GSR1" = list(
      position_name = "Galvanic skin response channel 1 (electrodermal activity/arousal)",
      position_type = "GSR",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "GSR1",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "GSR2" = list(
      position_name = "Galvanic skin response channel 2 (electrodermal activity/arousal)",
      position_type = "GSR",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "GSR2",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== BIOSEMI ACCESSORIES - Erg1/Erg2 (Ergo/AUX Input) =====
    "Erg1" = list(
      position_name = "Ergo/AUX input 1 (Response Switch/Photocell/Microphone)",
      position_type = "Ergo/AUX",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Erg1",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    "Erg2" = list(
      position_name = "Ergo/AUX input 2 (Response Switch/Photocell/Microphone)",
      position_type = "Ergo/AUX",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Erg2",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== BIOSEMI ACCESSORIES - Resp (Respiration Belt) =====
    "Resp" = list(
      position_name = "Respiration belt (breathing rate/pattern measurement)",
      position_type = "Respiration",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Resp",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== BIOSEMI ACCESSORIES - Plet (Plethysmograph) =====
    "Plet" = list(
      position_name = "Plethysmograph (blood volume/pulse measurement)",
      position_type = "Plethysmograph",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Plet",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== BIOSEMI ACCESSORIES - Temp (Temperature) =====
    "Temp" = list(
      position_name = "Temperature sensor (body temperature measurement)",
      position_type = "Temperature",
      region = "External",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Temp",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    ),
    # ===== Status =====
    "Status" = list(
      position_name = "Status channel (event markers/triggers)",
      position_type = "Status",
      region = "Metadata",
      standard_system = "BioSemi",
      biosemi_name = NA,
      standard_name = "Status",
      cartesian_coords = list(x = NA, y = NA, z = NA),
      spherical_coords = list(theta = NA, phi = NA, radius = NA)
    )
  )
}
