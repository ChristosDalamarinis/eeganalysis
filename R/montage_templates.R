# ============================================================================
#             BioSemi cap layouts used by create_montage()
# ============================================================================
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
# This is the same radius the 64-channel table in get_electrode_database()
# uses, so all templates share one scale.
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
# The 64-channel cap is not repeated here: it is read from
# get_electrode_database(), which holds BioSemi's numbers for that cap (and
# its A1-B32 names), so create_montage() gives the same result as before for
# "biosemi64".
#
# Names:
#   16 and 32 channels   standard 10-20 names (Fp1, Cz, ...). Row i of the sheet
#                        is BioSemi socket A<i>, so A1, A2, ... are accepted
#                        as well (the 64-channel database numbers its rows the
#                        same way).
#   128, 160, 256        BioSemi's own names: A1-D32, A1-E32 and A1-H32 (32
#                        sockets per letter). The sheets also print the
#                        nearest 10-20 name next to some of them; those are
#                        only approximate, so they are not used.

.BIOSEMI_RADIUS_MM <- 87.54

# Names of the layouts create_montage() can build: the BioSemi caps here and
# the standard systems of montage_standard.R.
.montage_templates <- function() {
  c("biosemi16", "biosemi32", "biosemi64", "biosemi128", "biosemi160",
    "biosemi256", .standard_templates())
}

# One row per electrode for one template: name, alias (a second name that is
# accepted, NA if none), x, y, z (mm), theta, phi (degrees) and radius (mm).
.montage_template_table <- function(template) {
  if (identical(template, "biosemi64")) {
    return(.biosemi64_table())
  }
  if (template %in% .standard_templates()) {
    return(.standard_template_table(template))
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

# The 64-channel cap, taken from the electrode database. Each scalp electrode
# sits in the database under four keys (both names, both lower-cased), so the
# repeats are dropped and the database's own order is kept.
.biosemi64_table <- function() {
  db <- get_electrode_database()

  is_scalp <- vapply(db, function(e) {
    isTRUE(e$standard_system %in% c("10-20", "10-10"))
  }, logical(1))
  scalp <- unname(db[is_scalp])
  scalp <- scalp[!duplicated(vapply(scalp, function(e) e$standard_name,
                                    character(1)))]

  data.frame(
    name   = vapply(scalp, function(e) e$standard_name, character(1)),
    alias  = vapply(scalp, function(e) e$biosemi_name, character(1)),
    x      = vapply(scalp, function(e) e$cartesian_coords$x, numeric(1)),
    y      = vapply(scalp, function(e) e$cartesian_coords$y, numeric(1)),
    z      = vapply(scalp, function(e) e$cartesian_coords$z, numeric(1)),
    theta  = vapply(scalp, function(e) e$spherical_coords$theta, numeric(1)),
    phi    = vapply(scalp, function(e) e$spherical_coords$phi, numeric(1)),
    radius = vapply(scalp, function(e) e$spherical_coords$radius, numeric(1)),
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
