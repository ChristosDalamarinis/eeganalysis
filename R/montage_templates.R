# ============================================================================
#                 The layouts create_montage() can build
# ============================================================================
#
# This file is the front door to the layouts. Each family of layouts has a file
# of its own:
#
#   montage_biosemi.R    the BioSemi caps with 16, 32, 64, 128, 160 and 256
#                        channels (and BioSemi's channels that are not on the
#                        scalp)
#   montage_standard.R   the standard 10-20, 10-10 and 10-05 systems, computed
#                        from the rules that define them
#   montage_colin27.R    the same systems measured on a real head
#
# Every family file has the same two functions: .<family>_templates() gives the
# names of its layouts, and .<family>_template_table(template) gives the table
# of one of them. The two functions below put the families together, and they
# are the only ones the rest of the package calls: create_montage(),
# set_montage() and the 3D electrode plots (through .plot_layout_table()).
#
# A layout table has one row per electrode and the columns name, alias (a
# second name that is accepted, NA if none), x, y, z (mm), theta, phi (degrees)
# and radius (mm). Positions are x towards the right ear, y towards the nose
# and z towards the top of the head, with the origin at the middle of the head.
#
# Author: Christos Dalamarinis
# Date: Oct - 2026

# The head radius of the layouts that are built on a sphere, 87.54 mm (a 55 cm
# head circumference divided by 2*pi). The BioSemi caps and the standard
# layouts use it, so they all share one scale; the real-head layouts have a
# distance of their own for every electrode.
.BIOSEMI_RADIUS_MM <- 87.54

# Names of the layouts create_montage() can build: the BioSemi caps of
# montage_biosemi.R, the standard systems of montage_standard.R and the
# real-head layouts of montage_colin27.R.
.montage_templates <- function() {
  c(.biosemi_templates(), .standard_templates(), .colin27_templates())
}

# One row per electrode for one layout, in the form described above. Stops for
# a name that is not a layout.
.montage_template_table <- function(template) {
  if (template %in% .biosemi_templates()) {
    return(.biosemi_template_table(template))
  }
  if (template %in% .standard_templates()) {
    return(.standard_template_table(template))
  }
  if (template %in% .colin27_templates()) {
    return(.colin27_template_table(template))
  }

  stop("ERROR: Unknown layout '", paste(template, collapse = ", "),
       "'. Available layouts: ", paste(.montage_templates(), collapse = ", "),
       ".", call. = FALSE)
}
