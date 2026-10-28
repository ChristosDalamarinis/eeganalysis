# ============================================================================
#          Standard electrode systems, built from the rules that define them
# ============================================================================
#
# The 10-10 system is defined by rules, not by a table of coordinates. The
# rules are described in Oostenveld R, Praamstra P (2001), "The five percent
# electrode system for high-resolution EEG and ERP measurements", Clinical
# Neurophysiology 112:713-719 (Section 2.1). The positions below are worked
# out from those rules; no coordinates are copied from anywhere.
#
# They are worked out on the same sphere as the BioSemi layouts: radius
# 87.54 mm, with the nasion, the inion and the two preauricular points 115
# degrees from the vertex, as in BioSemi's coordinate sheet. That makes the
# arc from the nasion over the vertex to the inion 230 degrees long, so 10
# percent of it is 23 degrees. The layouts then sit on one scale, need no
# recentring, and agree with the BioSemi layouts for the electrodes they
# share (see below).
#
# The rules, in the order they are applied (every arc is cut into equal parts
# by angle, and an "arc through three points" is the part of the circle in
# which the plane through those points cuts the sphere):
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

# Names of the layouts built here.
.standard_templates <- function() {
  "standard_1010"
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
# sphere (x towards the right ear, y towards the nose, z towards the top).
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
  points <- list()

  midline <- c("Fpz", "AFz", "Fz", "FCz", "Cz", "CPz", "Pz", "POz", "Oz")
  points[midline] <- along(nasion, vertex, inion, (1:9) / 10)
  points$T7 <- .standard_arc_point(left_ear, vertex, right_ear, 0.1)
  points$T8 <- .standard_arc_point(left_ear, vertex, right_ear, 0.9)

  left_ring  <- c("Fp1", "AF7", "F7", "FT7", "T7", "TP7", "P7", "PO7", "O1")
  right_ring <- c("Fp2", "AF8", "F8", "FT8", "T8", "TP8", "P8", "PO8", "O2")
  points[left_ring]  <- along(points$Fpz, points$T7, points$Oz, (1:9) / 10)
  points[right_ring] <- along(points$Fpz, points$T8, points$Oz, (1:9) / 10)

  rows <- list(AF = c("AF7", "AFz", "AF8"), F  = c("F7", "Fz", "F8"),
               FC = c("FT7", "FCz", "FT8"), C  = c("T7", "Cz", "T8"),
               CP = c("TP7", "CPz", "TP8"), P  = c("P7", "Pz", "P8"),
               PO = c("PO7", "POz", "PO8"))
  columns <- c("5", "3", "1", "2", "4", "6")
  parts <- c(1, 2, 3, 5, 6, 7) / 8
  for (row in names(rows)) {
    ends <- rows[[row]]
    points[paste0(row, columns)] <- along(points[[ends[1]]], points[[ends[2]]],
                                          points[[ends[3]]], parts)
  }

  inferior <- c("F", "FT", "T", "TP", "P")
  tenths <- (3:7) / 10
  points[paste0(inferior, "9")]  <- along(nasion, left_ear, inion, tenths)
  points[paste0(inferior, "10")] <- along(nasion, right_ear, inion, tenths)
  points$Iz <- inion

  points
}

# One row per electrode for a layout built here: the same columns as the
# tables of the BioSemi layouts (see .montage_template_table()). The alias is
# always NA: these layouts have one name per electrode.
.standard_template_table <- function(template) {
  names_in_layout <- switch(template, standard_1010 = .standard_1010_names())
  points <- .standard_points()[names_in_layout]

  xyz <- .BIOSEMI_RADIUS_MM * t(vapply(points, function(p) p, numeric(3)))
  angles <- .montage_angles(xyz[, 1], xyz[, 2], xyz[, 3])

  # six decimals, and no negative zeros (adding 0 turns -0 into 0)
  clean <- function(v) unname(round(v, 6) + 0)

  data.frame(
    name   = names_in_layout,
    alias  = NA_character_,
    x      = clean(xyz[, 1]),
    y      = clean(xyz[, 2]),
    z      = clean(xyz[, 3]),
    theta  = clean(angles$theta),
    phi    = clean(angles$phi),
    radius = .BIOSEMI_RADIUS_MM,
    stringsAsFactors = FALSE
  )
}
