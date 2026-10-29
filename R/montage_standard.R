# ============================================================================
#          Standard electrode systems, built from the rules that define them
# ============================================================================
#
# The 10-20, 10-10 and 10-05 systems are defined by rules, not by a table of
# coordinates. The rules are described in Oostenveld R, Praamstra P (2001),
# "The five percent electrode system for high-resolution EEG and ERP
# measurements", Clinical Neurophysiology 112:713-719 (Sections 2.1 and 2.2).
# The positions below are worked out from those rules; no coordinates are
# copied from anywhere.
#
# They are worked out on the same sphere as the BioSemi layouts: radius
# 87.54 mm, with the nasion, the inion and the two preauricular points 115
# degrees from the vertex, as in BioSemi's coordinate sheet. That makes the
# arc from the nasion over the vertex to the inion 230 degrees long, so 10
# percent of it is 23 degrees. The layouts then sit on one scale, need no
# recentring, and agree with the BioSemi layouts for the electrodes they
# share (see below).
#
# The rules of the 10-10 system, in the order they are applied (every arc is
# cut into equal parts by angle, and an "arc through three points" is the
# part of the circle in which the plane through those points cuts the
# sphere):
#   * Cz, the vertex, is the middle of the nasion-inion arc and of the arc
#     between the two preauricular points.
#   * Midline: Fpz, AFz, Fz, FCz, Cz, CPz, Pz, POz and Oz are at 10, 20, ...,
#     90 percent of the nasion-inion arc.
#   * T7 and T8 are at 10 and 90 percent of the arc between the preauricular
#     points, over the vertex.
#   * Outer ring: the arc from Fpz through T7 to Oz holds Fp1, AF7, F7, FT7,
#     T7, TP7, P7, PO7 and O1 at 10 percent steps; the arc from Fpz through T8
#     to Oz holds the right-hand ones (Fp2 ... O2).
#   * Rows: the arc from a left ring electrode through the midline electrode
#     of its row to the right ring electrode (AF7-AFz-AF8, F7-Fz-F8,
#     FT7-FCz-FT8, T7-Cz-T8, TP7-CPz-TP8, P7-Pz-P8, PO7-POz-PO8) is cut into 8
#     equal parts, which gives the columns 7 5 3 1 z 2 4 6 8.
#   * Inferior row: the circle through the nasion, a preauricular point and
#     the inion (115 degrees from the vertex) holds F9, FT9, T9, TP9 and P9
#     (and F10 ... P10 on the right) at 10 percent steps; T9 and T10 are the
#     preauricular points themselves.
#
# standard_1010 has 72 electrodes: the 74 positions of the 10-10 standard
# without the two ear electrodes A1 and A2, which are not scalp positions.
# The landmark Nz is not an electrode. AF1, AF2, AF5, AF6, PO1, PO2, PO5 and
# PO6 follow from the same rules but are not part of the 10-10 standard.
#
# For the 64 electrodes that standard_1010 shares with the BioSemi 64 channel
# layout, the positions are the same up to BioSemi's rounding (within about
# 1 mm) for 48 of them. The other 16 differ more: AF3, AF4, PO3 and PO4
# (about 5.7 mm), FC5, FC6, CP5 and CP6 (about 2.1 mm), F3, F4, P3 and P4
# (about 1.7 mm) and F5, F6, P5 and P6 (about 1.5 mm).
#
# The 10-05 system (standard_1005) halves every step of the 10-10 system, in
# both directions. The same rules apply, with 5 percent instead of 10:
#   * Midline: AFpz, AFFz, FFCz, FCCz, CCPz, CPPz, PPOz, POOz and OIz are
#     added, so the midline is cut at 5, 10, 15, ... 95 percent.
#   * Outer ring: the ring is cut at 5 percent steps as well. Next to the
#     10-10 electrodes it holds Fp1h, AFp7, AFF7, FFT7, FTT7, TTP7, TPP7,
#     PPO7, POO7 and O1h (and the right-hand ones).
#   * Rows: the 7 rows of the 10-10 system and the 8 new ones between them
#     (AFp, AFF, FFC, FCC, CCP, CPP, PPO, POO) are each the arc through a
#     left ring electrode, the midline electrode and the right ring electrode,
#     cut into 16 equal parts. Every second part is a column of the 10-10
#     system, the others are the new "half" columns.
#   * Inferior row: cut at 5 percent steps as well, so each of the 15 rows
#     has a column 9 and a column 10. I1 and I2 are the 90 percent
#     electrodes of this row, at the place where the 10-20 system would put
#     O9 and O10.
#
# Names. A "half" electrode lies halfway between two neighbours of its
# row and carries the name of the more lateral one with an "h" added: between
# Cz and C1 is C1h, between C3 and C5 is C5h. A new row is named after the
# two rows around it, like the points of a compass: CCP lies between C and
# CP. From column 7 outwards the rows between F and P run along the temporal
# electrodes, so there the names change from FFC, FC, FCC, C, CCP, CP, CPP
# to FFT, FT, FTT, T, TTP, TP, TPP (FFT7, FT7, FTT7, T7, TTP7, TP7, TPP7).
# The electrodes at the very front and back, which are not on a row (Fp, O,
# OI, I), keep five places each: Fp1, Fp1h, Fpz, Fp2h, Fp2 and so on.
#
# Where the paper leaves the position open the choice made here is: a half
# electrode between column 7 and column 9 (or 8 and 10) is halfway along the
# shortest path between its two neighbours, and OI1 is halfway between O1 and
# I1. Every other half electrode is a part of the arc it is on.
#
# standard_1005 has 335 electrodes: 5 + 15 * 21 + 15, front to back, with
# every position of standard_1010 at the same place. The landmark Nz, the
# mastoid electrodes (M1, M2) and the ear electrodes (A1, A2) have no rule and
# are not part of it.
#
# The 10-20 system (standard_1020) is the original international system and
# the smallest of the three. The 10-10 system was built to keep its places, so
# its 19 scalp electrodes (Fp1, Fp2, F7, F3, Fz, F4, F8, T7, C3, Cz, C4, T8,
# P7, P3, Pz, P4, P8, O1 and O2) are electrodes of the grid too, at exactly
# the same places as in the other two layouts. Fpz and Oz are not part of the
# original system, and neither are its two ear electrodes (A1, A2), which are
# not scalp positions.
#
# Old names. T7, T8, P7 and P8 were called T3, T4, T5 and T6 in the original
# 10-20 system. All three layouts accept the old names as well: the table
# gives them as the alias of the electrode, and create_montage() looks a name
# up by its alias when it has no electrode of that name.

# Names of the layouts built here.
.standard_templates <- function() {
  c("standard_1020", "standard_1010", "standard_1005")
}

# The 19 electrodes of standard_1020, front to back and left to right.
.standard_1020_names <- function() {
  c("Fp1", "Fp2",
    "F7", "F3", "Fz", "F4", "F8",
    "T7", "C3", "Cz", "C4", "T8",
    "P7", "P3", "Pz", "P4", "P8",
    "O1", "O2")
}

# The old names of four electrodes, by their present names.
.standard_old_names <- function() {
  c(T7 = "T3", T8 = "T4", P7 = "T5", P8 = "T6")
}

# The 72 electrodes of standard_1010, front to back and left to right.
.standard_1010_names <- function() {
  c("Fp1", "Fpz", "Fp2",
    "AF7", "AF3", "AFz", "AF4", "AF8",
    "F9", "F7", "F5", "F3", "F1", "Fz", "F2", "F4", "F6", "F8", "F10",
    "FT9", "FT7", "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6", "FT8",
    "FT10",
    "T9", "T7", "C5", "C3", "C1", "Cz", "C2", "C4", "C6", "T8", "T10",
    "TP9", "TP7", "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6", "TP8",
    "TP10",
    "P9", "P7", "P5", "P3", "P1", "Pz", "P2", "P4", "P6", "P8", "P10",
    "PO7", "PO3", "POz", "PO4", "PO8",
    "O1", "Oz", "O2",
    "Iz")
}

# The 15 rows of the 10-05 system between Fp and O, front to back. `medial`
# is how the row is named from column 5 inwards and on the midline, `lateral`
# from column 7 outwards.
.standard_row_medial <- function() {
  c("AFp", "AF", "AFF", "F", "FFC", "FC", "FCC", "C", "CCP", "CP", "CPP", "P",
    "PPO", "PO", "POO")
}
.standard_row_lateral <- function() {
  c("AFp", "AF", "AFF", "F", "FFT", "FT", "FTT", "T", "TTP", "TP", "TPP", "P",
    "PPO", "PO", "POO")
}

# The 21 electrodes of one row, from the far left to the far right: columns 9,
# 9h, 7, 7h, 5, 5h, 3, 3h, 1, 1h, z, 2h, 2, 4h, 4, 6h, 6, 8h, 8, 10h, 10.
.standard_row_names <- function(medial, lateral) {
  c(paste0(lateral, c("9", "9h", "7", "7h")),
    paste0(medial, c("5", "5h", "3", "3h", "1", "1h", "z", "2h", "2", "4h",
                     "4", "6h", "6")),
    paste0(lateral, c("8h", "8", "10h", "10")))
}

# The 335 electrodes of standard_1005, front to back and left to right.
.standard_1005_names <- function() {
  rows <- unlist(lapply(seq_along(.standard_row_medial()), function(i) {
    .standard_row_names(.standard_row_medial()[i], .standard_row_lateral()[i])
  }))
  c("Fp1", "Fp1h", "Fpz", "Fp2h", "Fp2",
    rows,
    "O1", "O1h", "Oz", "O2h", "O2",
    "OI1", "OI1h", "OIz", "OI2h", "OI2",
    "I1", "I1h", "Iz", "I2h", "I2")
}

.standard_cross <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

# The point at fraction f (0 to 1) of the way along the arc that starts at
# `first`, passes through `through` and ends at `last` (points on the unit
# sphere): the circle in which the plane through the three points cuts the
# sphere, cut into equal angles about its centre.
.standard_arc_point <- function(first, through, last, f) {
  normal <- .standard_cross(last - first, through - first)
  normal <- normal / sqrt(sum(normal^2))
  offset <- sum(normal * first)
  centre <- offset * normal

  e1 <- first - centre
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- .standard_cross(normal, e1)

  angle_of <- function(point) {
    angle <- atan2(sum((point - centre) * e2), sum((point - centre) * e1))
    if (angle < 0) angle + 2 * pi else angle
  }
  through_angle <- angle_of(through)
  last_angle <- angle_of(last)
  total <- if (through_angle < last_angle) last_angle else last_angle - 2 * pi

  angle <- total * f
  centre + sqrt(max(0, 1 - offset^2)) * (cos(angle) * e1 + sin(angle) * e2)
}

# Every position the rules above give, as a named list of points on the unit
# sphere (x towards the right ear, y towards the nose, z towards the top):
# the 335 electrodes of standard_1005, which include those of standard_1010.
.standard_points <- function() {
  edge <- 115 * pi / 180
  nasion    <- c(0, sin(edge), cos(edge))
  inion     <- c(0, -sin(edge), cos(edge))
  left_ear  <- c(-sin(edge), 0, cos(edge))
  right_ear <- c(sin(edge), 0, cos(edge))
  vertex    <- c(0, 0, 1)

  along <- function(first, through, last, fractions) {
    lapply(fractions, function(f) .standard_arc_point(first, through, last, f))
  }
  # halfway along the shortest path between two points
  halfway <- function(a, b) {
    middle <- a + b
    middle / sqrt(sum(middle^2))
  }
  medial  <- .standard_row_medial()
  lateral <- .standard_row_lateral()
  points <- list()

  # midline, at 5 percent steps from Fpz (10) to OIz (95); Iz is the inion
  midline <- c("Fpz", paste0(medial, "z"), "Oz", "OIz")
  points[midline] <- along(nasion, vertex, inion, (2:19) / 20)
  points$Iz <- inion
  points$T7 <- .standard_arc_point(left_ear, vertex, right_ear, 0.1)
  points$T8 <- .standard_arc_point(left_ear, vertex, right_ear, 0.9)

  # outer ring, at 5 percent steps from Fp1h (5) to O1h (95)
  left_ring  <- c("Fp1h", "Fp1", paste0(lateral, "7"), "O1", "O1h")
  right_ring <- c("Fp2h", "Fp2", paste0(lateral, "8"), "O2", "O2h")
  points[left_ring]  <- along(points$Fpz, points$T7, points$Oz, (1:19) / 20)
  points[right_ring] <- along(points$Fpz, points$T8, points$Oz, (1:19) / 20)

  # inferior row, at 5 percent steps from column 9 of AFp (15) to I1h (95)
  steps <- (3:19) / 20
  points[c(paste0(lateral, "9"), "I1", "I1h")]   <- along(nasion, left_ear, inion, steps)
  points[c(paste0(lateral, "10"), "I2", "I2h")] <- along(nasion, right_ear, inion, steps)

  # rows: columns 7 to 8 are sixteenths of the arc through the ring electrodes
  # and the midline electrode (the ends and the middle are already there);
  # the half columns 9h and 10h are halfway between 7 and 9, and 8 and 10
  inner <- c(1:7, 9:15)
  for (i in seq_along(medial)) {
    row <- .standard_row_names(medial[i], lateral[i])
    points[row[3 + inner]] <- along(points[[row[3]]], points[[row[11]]],
                                    points[[row[19]]], inner / 16)
    points[[row[2]]]  <- halfway(points[[row[1]]], points[[row[3]]])
    points[[row[20]]] <- halfway(points[[row[19]]], points[[row[21]]])
  }

  # below O1 and O2: OI1 is halfway between O1 and I1, OI1h halfway to OIz
  points$OI1 <- halfway(points$O1, points$I1)
  points$OI2 <- halfway(points$O2, points$I2)
  points[c("OI1h", "OI2h")] <- along(points$OI1, points$OIz, points$OI2,
                                     c(0.25, 0.75))

  points
}

# One row per electrode for a layout built here: the same columns as the
# tables of the BioSemi layouts (see .montage_template_table()). The alias is
# the old name of T7, T8, P7 and P8 (T3, T4, T5, T6) and NA for every other
# electrode, which has one name only.
.standard_template_table <- function(template) {
  names_in_layout <- switch(template,
                            standard_1020 = .standard_1020_names(),
                            standard_1010 = .standard_1010_names(),
                            standard_1005 = .standard_1005_names())
  points <- .standard_points()[names_in_layout]

  xyz <- .BIOSEMI_RADIUS_MM * t(vapply(points, function(p) p, numeric(3)))
  angles <- .montage_angles(xyz[, 1], xyz[, 2], xyz[, 3])

  # six decimals, and no negative zeros (adding 0 turns -0 into 0)
  clean <- function(v) unname(round(v, 6) + 0)

  data.frame(
    name   = names_in_layout,
    alias  = unname(.standard_old_names()[names_in_layout]),
    x      = clean(xyz[, 1]),
    y      = clean(xyz[, 2]),
    z      = clean(xyz[, 3]),
    theta  = clean(angles$theta),
    phi    = clean(angles$phi),
    radius = .BIOSEMI_RADIUS_MM,
    stringsAsFactors = FALSE
  )
}
