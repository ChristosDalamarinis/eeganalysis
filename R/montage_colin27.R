# ============================================================================
#            Real-head layouts used by create_montage(): colin27
# ============================================================================
#
# Two more layouts for create_montage(), set_montage() and the 3D electrode
# plots. Their positions were measured on a real head, they are not worked out
# on a sphere (compare montage_standard.R):
#
#   colin27_1020   94 electrodes    file standard_1020.elc
#   colin27_1005   343 electrodes   file standard_1005.elc
#
# They are the electrode tables of the FieldTrip toolbox, made by Robert
# Oostenveld on the Colin27 template head (one person's head, scanned many
# times and averaged): millimetres, in the MNI coordinate system. The two files
# are kept exactly as FieldTrip publishes them, in inst/extdata/fieldtrip/
# (GNU General Public License, version 3 or later; the credit is in the
# DESCRIPTION file and in the help of create_montage()).
#
# What is done with a file:
#
#   - It is read with read_custom_montage(), which leaves out the three
#     landmark rows (LPA, RPA and Nz).
#   - Every position is shifted by the same amount, so that the centre of the
#     sphere that fits the electrodes best is at the origin. The other layouts
#     are centred there, and so are the functions that use positions
#     (interpolate_bads(), compute_csd() and plot_topography()). One shift for
#     all electrodes leaves the shape exactly as it is in the file.
#   - The shift is the same for both layouts: the centre of the sphere that
#     fits the 94 electrodes of colin27_1020 (the set that most caps are a part
#     of). The two files give the electrodes they share the same numbers, so an
#     electrode then has the same position in both layouts, as in the standard
#     layouts of montage_standard.R.
#   - Every other row is kept, as in the original. The old names T3, T4, T5 and
#     T6 are rows of their own, with exactly the position of T7, T8, P7 and P8,
#     and M1, M2, A1 and A2 are rows too. So the alias column is empty.
#
# A real head is not round, so the radius column holds the distance of each
# electrode from the origin (the other layouts have one head radius for all of
# them).
#
# Author: Christos Dalamarinis
# Date: Nov - 2026

# The file of each layout, in inst/extdata/fieldtrip/ of the package.
.colin27_files <- function() {
  c(colin27_1020 = "standard_1020.elc",
    colin27_1005 = "standard_1005.elc")
}

# Names of the layouts of this file.
.colin27_templates <- function() {
  names(.colin27_files())
}

# The path of the file of one layout. Stops for a name that is not a layout of
# this file, and when the file is not in the package.
.colin27_path <- function(template) {
  files <- .colin27_files()
  if (!is.character(template) || length(template) != 1 || is.na(template) ||
      !(template %in% names(files))) {
    stop("ERROR: Unknown real-head layout '", paste(template, collapse = ", "),
         "'. Available layouts: ", paste(names(files), collapse = ", "), ".",
         call. = FALSE)
  }

  path <- system.file("extdata", "fieldtrip", files[[template]],
                      package = "eeganalysis")
  if (!nzchar(path)) {
    stop("ERROR: The electrode file '", files[[template]], "' of the layout '",
         template, "' was not found in the package.", call. = FALSE)
  }
  path
}

# The centre of the sphere that fits the 94 electrodes of colin27_1020 best,
# in the coordinates of the files (mm). Both layouts are shifted by this.
.colin27_centre <- function() {
  raw <- read_custom_montage(.colin27_path("colin27_1020"),
                             recenter = FALSE)$positions
  .montage_sphere_fit(raw$x, raw$y, raw$z)$centre
}

# One row per electrode for one layout, in the form of .montage_template_table():
# name, alias (NA here), x, y, z (mm), theta, phi (degrees) and radius (mm).
.colin27_template_table <- function(template) {
  raw <- read_custom_montage(.colin27_path(template), recenter = FALSE)$positions

  centre <- .colin27_centre()
  x <- raw$x - centre[1]
  y <- raw$y - centre[2]
  z <- raw$z - centre[3]
  angles <- .montage_angles(x, y, z)

  # six decimals, and no negative zeros (adding 0 turns -0 into 0), as in the
  # other layouts
  clean <- function(v) unname(round(v, 6) + 0)

  data.frame(
    name   = raw$channel,
    alias  = NA_character_,
    x      = clean(x),
    y      = clean(y),
    z      = clean(z),
    theta  = clean(angles$theta),
    phi    = clean(angles$phi),
    radius = clean(angles$radius),
    stringsAsFactors = FALSE
  )
}
