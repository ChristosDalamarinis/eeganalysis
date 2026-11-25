# ============================================================================
#                Test File for topomap_interpolation.R
# ============================================================================
#
# This test file provides comprehensive testing for the smooth surface behind
# plot_topography(): a mesh of triangles between the electrodes, a ring of
# extra points around the head, and a cubic (Clough-Tocher) surface over the
# mesh.
#
# Functions tested:
#   1. .topomap_grid()                - the surface read on a regular grid
#   2. .topo_surface(), .topo_read()  - the surface and a read-out at any point
#   3. .topo_ring()                   - the ring of extra points
#   4. .topo_qhull_spacing(), .topo_edge_median()
#                                     - the electrode spacing the ring follows
#   5. .topo_delaunay()               - the mesh of triangles
#   6. .topo_triangle_neighbours(), .topo_slopes(), .topo_locate()
#
# Reference numbers: the numbers written out in Suite 1 were made once by the
# interpolator that MNE-Python uses for its topographic maps (cubic, with the
# ring of extra points around the head and the mean of the neighbours as the
# value of a ring point), run on the same electrodes and the same values. They
# are plain numbers here, so the tests need no Python. If a layout of the
# package ever changes, the numbers of that layout have to be made again.
#
# Author: Christos Dalamarinis
# Date: November 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# A small, irregular set of 14 electrodes with a value each: no two triangles
# of its mesh are tied, and no three electrodes lie on a line. `mne_spacing` is
# the distance between ring points that the reference program works out for it.
small_set <- function() {
  list(
    x = c(-1.30, -0.95, -0.40, 0.15, 0.70, 1.20, -1.05, -0.35, 0.35, 1.00,
          -0.60, 0.10, 0.85, -0.05),
    y = c(0.20, 0.95, 1.25, 1.35, 1.00, 0.35, -0.60, -0.85, -0.95, -0.55,
          0.15, 0.30, 0.15, -0.20),
    v = c(0.8, -0.3, 1.1, 0.4, -0.7, 0.2, 0.9, -1.0, 0.5, 0.6,
          -0.2, 1.3, -0.4, 0.0),
    head_r = 1.8,
    mne_spacing = 0.85656691254080175
  )
}

# The 2D position of every electrode of a built-in layout, projected the way
# plot_topography() projects them, and the radius of the head disc.
layout_2d <- function(template) {
  montage <- create_montage(template = template)
  pos <- montage$positions[match(montage$channels, montage$positions$channel), ]
  r <- sqrt(pos$x^2 + pos$y^2 + pos$z^2)
  elevation <- acos(pos$z / r)
  azimuth   <- atan2(pos$y, pos$x)
  x <- elevation * cos(azimuth)
  y <- elevation * sin(azimuth)
  list(x = x, y = y, head_r = max(abs(c(x, y))) * 1.1)
}

# The values given to every layout in Suite 1.
smooth_pattern <- function(x, y) {
  sin(2.5 * x) + 0.6 * cos(3 * y) + 0.4 * cos(sqrt(x^2 + y^2))
}

# Reads cells (i, j) = (position along x, position along y) of a grid matrix.
cells_of <- function(z, cells) {
  vapply(seq_len(nrow(cells)), function(k) z[cells[k, 1], cells[k, 2]], numeric(1))
}

# Is each point (px, py) inside the convex hull of the points (x, y), and how
# far is it from the edge of the hull? Written out here on purpose, so that the
# tests do not lean on the code they test.
hull_check <- function(px, py, x, y) {
  h  <- grDevices::chull(x, y)
  hx <- x[h]
  hy <- y[h]
  n  <- length(hx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    crosses <- ((hy[i] > py) != (hy[j] > py)) &
      (px < (hx[j] - hx[i]) * (py - hy[i]) / (hy[j] - hy[i]) + hx[i])
    inside <- xor(inside, crosses)
    j <- i
  }
  edge_distance <- vapply(seq_along(px), function(q) {
    ends <- c(seq_len(n)[-1], 1)
    min(vapply(seq_len(n), function(e) {
      x1 <- hx[e]; y1 <- hy[e]; x2 <- hx[ends[e]]; y2 <- hy[ends[e]]
      t <- max(0, min(1, ((px[q] - x1) * (x2 - x1) + (py[q] - y1) * (y2 - y1)) /
                        ((x2 - x1)^2 + (y2 - y1)^2)))
      sqrt((px[q] - (x1 + t * (x2 - x1)))^2 + (py[q] - (y1 + t * (y2 - y1)))^2)
    }, numeric(1)))
  }, numeric(1))
  list(inside = inside, edge_distance = edge_distance)
}

# Checks a mesh against the definition of a Delaunay triangulation: it covers
# exactly the convex hull of the points, every point is a corner of some
# triangle, and no point lies inside the circle through the corners of a
# triangle (the number returned is the worst case of that, relative; 0 or
# below is perfect).
mesh_problems <- function(points, mesh) {
  worst <- -Inf
  for (t in seq_len(nrow(mesh))) {
    a <- points[mesh[t, 1], ]; b <- points[mesh[t, 2], ]; c <- points[mesh[t, 3], ]
    orient <- sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
    bx <- b[1] - a[1]; by <- b[2] - a[2]; cx <- c[1] - a[1]; cy <- c[2] - a[2]
    dx <- points[, 1] - a[1]; dy <- points[, 2] - a[2]
    det3 <- (bx^2 + by^2) * (cx * dy - dx * cy) -
      (cx^2 + cy^2) * (bx * dy - dx * by) + (dx^2 + dy^2) * (bx * cy - cx * by)
    size <- (bx^2 + by^2) * (abs(cx * dy) + abs(dx * cy)) +
      (cx^2 + cy^2) * (abs(bx * dy) + abs(dx * by)) +
      (dx^2 + dy^2) * (abs(bx * cy) + abs(cx * by))
    violation <- -orient * det3 / pmax(size, 1e-300)
    violation[mesh[t, ]] <- -Inf
    worst <- max(worst, violation)
  }
  area <- sum(abs((points[mesh[, 2], 1] - points[mesh[, 1], 1]) *
                    (points[mesh[, 3], 2] - points[mesh[, 1], 2]) -
                    (points[mesh[, 2], 2] - points[mesh[, 1], 2]) *
                    (points[mesh[, 3], 1] - points[mesh[, 1], 1]))) / 2
  h  <- grDevices::chull(points)
  hx <- points[h, 1] - mean(points[, 1])   # centred, so far-away points do not lose digits
  hy <- points[h, 2] - mean(points[, 2])
  hull_area <- abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy)) / 2
  list(violation = worst, area_ratio = area / hull_area,
       unused = nrow(points) - length(unique(as.vector(mesh))))
}

# ============================================================================
#     TEST SUITE 1: The Same Numbers as the Reference Program
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: A small irregular set, ring spaced as the reference program does
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the whole surface (mesh, ring, slopes, cubic
# pieces, reading) against numbers made by the reference program for 14
# hand-written electrodes, at 12 cells of an 11 x 11 grid: in the middle, between
# electrodes, and in the band between the outer electrodes and the edge of the
# head. The distance between ring points is the one the reference program
# works out, so this isolates the surface itself from the choice of that
# distance (tested in Suite 4).
test_that("a small irregular set gives the numbers of the reference program", {
  s <- small_set()
  cells <- cbind(c(6, 4, 8, 6, 6, 3, 9, 2, 6, 10, 5, 7),
                 c(6, 6, 6, 3, 9, 4, 8, 6, 10, 6, 8, 4))
  reference <- c(0.551740127728, -0.167714287829, -0.256967888975,
                 -0.516461889134, 0.960745188899, 0.96313460626,
                 -0.372522724096, 1.20970749828, 0.653934968014,
                 0.512855615446, 0.733613200218, 0.382723564708)

  fit <- eeganalysis:::.topomap_grid(s$x, s$y, s$v, s$head_r, n_grid = 11,
                                     ring_spacing = s$mne_spacing)

  expect_lt(max(abs(cells_of(fit$z, cells) - reference)), 1e-6)
})

# ----------------------------------------------------------------------------
# Test 1.2: Three BioSemi caps
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the maps of the 16, 32 and 64 electrode BioSemi
# layouts at 8 cells each of the 100 x 100 grid (middle, between electrodes,
# towards the edge of the head), with the distance between ring points worked
# out by the package itself. On these three caps the package works out the
# same distance as the reference program, so the numbers agree to 1e-6.
test_that("the 16, 32 and 64 electrode BioSemi layouts give the reference numbers", {
  cells <- cbind(c(50, 30, 50, 70, 25, 80, 12, 50),
                 c(50, 50, 30, 60, 25, 75, 50, 90))
  reference <- list(
    biosemi16 = c(0.9634790032, -0.05130858294, -0.1390814893, 1.47422389,
                  -1.399465747, 0.7185409151, 0.8761135444, -0.2128461975),
    biosemi32 = c(0.9542585381, -0.08511450186, -0.1174426786, 1.553219191,
                  -1.135819655, -0.03902377937, 0.9381698401, -0.1249909904),
    biosemi64 = c(0.9431011361, 0.0924593925, -0.3591011167, 1.233131179,
                  -0.8878931912, -0.7482324387, 1.465682408, 0.1244004951)
  )

  for (template in names(reference)) {
    p   <- layout_2d(template)
    fit <- eeganalysis:::.topomap_grid(p$x, p$y, smooth_pattern(p$x, p$y),
                                       p$head_r, n_grid = 100)
    expect_lt(max(abs(cells_of(fit$z, cells) - reference[[template]])), 1e-6)
  }
})

# ----------------------------------------------------------------------------
# Test 1.3: The 10-10 layout and a real-head layout with shared positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the maps of the standard 10-10 layout (72
# electrodes) and the real-head 10-10 layout colin27_1020 (94 electrodes, of
# which 4 share a position with another one, so the surface is made from 90).
# On these layouts the distance between ring points that the package works out
# is within 1 percent of the reference program's, so the ring has a few points
# more or less than the reference one and the band outside the electrodes
# differs a little; the numbers agree to 1e-3 on a map that spans 3.7.
test_that("the 10-10 and colin27_1020 layouts give the reference numbers closely", {
  cells <- cbind(c(50, 30, 50, 70, 25, 80, 12, 50),
                 c(50, 50, 30, 60, 25, 75, 50, 90))
  reference <- list(
    standard_1010 = c(0.9424791351, 0.08673090655, -0.3562822802, 1.229391181,
                      -0.8935123139, -0.8253877265, 1.406050948, 0.125386925),
    colin27_1020  = c(0.9327523122, 0.2244420622, -0.4417380556, 0.9771758033,
                      -0.5474906242, -1.117166072, 1.398248224, 0.2003911749)
  )

  for (template in names(reference)) {
    p   <- layout_2d(template)
    fit <- eeganalysis:::.topomap_grid(p$x, p$y, smooth_pattern(p$x, p$y),
                                       p$head_r, n_grid = 100)
    expect_lt(max(abs(cells_of(fit$z, cells) - reference[[template]])), 1e-3)
  }
})

# ============================================================================
#            TEST SUITE 2: Properties of the Surface
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: The surface passes through every electrode value
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that reading the surface at the position of an
# electrode gives that electrode's value (to 1e-9), which is what makes it an
# interpolation and not a smoothing. Checked on the small set, on the 64
# electrode layout, and on colin27_1020, where the 4 electrodes that share a
# position must give the mean of their values.
test_that("the surface passes through the value of every electrode", {
  s  <- small_set()
  p  <- layout_2d("biosemi64")
  pc <- layout_2d("colin27_1020")

  surface <- eeganalysis:::.topo_surface(s$x, s$y, s$v, s$head_r)
  expect_lt(max(abs(eeganalysis:::.topo_read(surface, surface$electrodes[, 1],
                                             surface$electrodes[, 2]) -
                      surface$values)), 1e-9)
  expect_equal(sort(surface$values), sort(s$v))

  v64 <- smooth_pattern(p$x, p$y)
  surface <- eeganalysis:::.topo_surface(p$x, p$y, v64, p$head_r)
  expect_lt(max(abs(eeganalysis:::.topo_read(surface, surface$electrodes[, 1],
                                             surface$electrodes[, 2]) -
                      surface$values)), 1e-9)
  expect_equal(nrow(surface$electrodes), 64)

  vc <- smooth_pattern(pc$x, pc$y)
  surface <- eeganalysis:::.topo_surface(pc$x, pc$y, vc, pc$head_r)
  expect_lt(max(abs(eeganalysis:::.topo_read(surface, surface$electrodes[, 1],
                                             surface$electrodes[, 2]) -
                      surface$values)), 1e-9)
  expect_equal(nrow(surface$electrodes), 90)
})

# ----------------------------------------------------------------------------
# Test 2.2: The map does not depend on the order of the electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the same electrodes in a different order give the
# same map (to 1e-10), because they are put in the order of their positions
# first. Checked with a shuffle of the small set and with the 64 electrode
# layout in reverse. The old interpolation package could hang R on some
# orders of the 128 electrode layout; here the order cannot matter.
test_that("the map is the same whatever the order of the electrodes", {
  s <- small_set()
  base <- eeganalysis:::.topomap_grid(s$x, s$y, s$v, s$head_r, n_grid = 31)
  shuffled <- c(9, 3, 14, 1, 7, 12, 5, 10, 2, 13, 6, 11, 4, 8)
  again <- eeganalysis:::.topomap_grid(s$x[shuffled], s$y[shuffled],
                                       s$v[shuffled], s$head_r, n_grid = 31)
  expect_equal(again$z, base$z, tolerance = 1e-10)

  p <- layout_2d("biosemi64")
  v <- smooth_pattern(p$x, p$y)
  a <- eeganalysis:::.topomap_grid(p$x, p$y, v, p$head_r, n_grid = 41)
  b <- eeganalysis:::.topomap_grid(rev(p$x), rev(p$y), rev(v), p$head_r, n_grid = 41)
  expect_equal(b$z, a$z, tolerance = 1e-10)
})

# ----------------------------------------------------------------------------
# Test 2.3: Electrodes at the same position count as one, with the mean value
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that two electrodes at one position give exactly
# the map of one electrode there with the mean of the two values (checked with
# values 0.2 and 0.6 against 0.4), and that this does not depend on which of
# the two comes first.
test_that("electrodes at the same position count as one with the mean value", {
  s <- small_set()
  once <- s$v
  once[3] <- 0.4
  twice_x <- c(s$x, s$x[3])
  twice_y <- c(s$y, s$y[3])
  twice_v <- c(replace(s$v, 3, 0.2), 0.6)

  reference <- eeganalysis:::.topomap_grid(s$x, s$y, once, s$head_r, n_grid = 31)
  merged <- eeganalysis:::.topomap_grid(twice_x, twice_y, twice_v, s$head_r,
                                        n_grid = 31)
  swapped <- eeganalysis:::.topomap_grid(twice_x, twice_y,
                                         c(replace(s$v, 3, 0.6), 0.2),
                                         s$head_r, n_grid = 31)

  expect_equal(merged$z, reference$z, tolerance = 1e-10)
  expect_equal(swapped$z, reference$z, tolerance = 1e-10)
})

# ----------------------------------------------------------------------------
# Test 2.4: extrapolate = FALSE only empties the cells outside the electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies with extrapolate = FALSE that every cell that is
# still painted holds exactly the value it has with extrapolate = TRUE, and
# that the empty cells inside the head disc are exactly the ones outside the
# convex hull of the electrodes (found here with an independent test; cells
# less than 1e-6 from the edge of the hull are left out of the comparison).
test_that("extrapolate = FALSE only empties the cells outside the electrodes", {
  s <- small_set()
  with_ring <- eeganalysis:::.topomap_grid(s$x, s$y, s$v, s$head_r, n_grid = 41,
                                           extrapolate = TRUE)
  without   <- eeganalysis:::.topomap_grid(s$x, s$y, s$v, s$head_r, n_grid = 41,
                                           extrapolate = FALSE)

  kept <- is.finite(without$z)
  expect_identical(without$z[kept], with_ring$z[kept])

  cx <- rep(without$x, times = 41)
  cy <- rep(without$y, each = 41)
  in_disc <- sqrt(cx^2 + cy^2) <= s$head_r
  hull <- hull_check(cx, cy, s$x, s$y)
  clear <- in_disc & hull$edge_distance > 1e-6

  expect_true(all(is.finite(with_ring$z[in_disc])))
  expect_identical(as.vector(kept)[clear], hull$inside[clear])
  expect_true(any(in_disc & !as.vector(kept)))
})

# ----------------------------------------------------------------------------
# Test 2.5: Constant values give a constant map
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that when every electrode has the same value, every
# cell of the head disc has that value (to 1e-12): the ring takes the same
# value, all slopes are zero, and nothing is added anywhere.
test_that("constant values give a constant map", {
  s <- small_set()
  fit <- eeganalysis:::.topomap_grid(s$x, s$y, rep(2.5, length(s$x)), s$head_r,
                                     n_grid = 31)
  inside <- is.finite(fit$z)

  expect_gt(sum(inside), 500)
  expect_lt(max(abs(fit$z[inside] - 2.5)), 1e-12)
})

# ----------------------------------------------------------------------------
# Test 2.6: The unit of the positions and of the values does not matter
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the map follows the numbers it is given (to 1e-8,
# compared after bringing the numbers back to the size of the plain map, so the
# comparison is on numbers of size 1): the same map when all positions and the
# head radius are multiplied by 1000 or by 0.001; the same map after
# multiplying the values by 1e-9 or 1e9 (values in volts instead of
# microvolts, where the slope estimate would otherwise stop after its first
# sweep and give a badly fitted map); and the same map plus 1e6 when 1e6 is
# added to every value.
test_that("the map does not change with the unit of the positions or of the values", {
  s <- small_set()
  base <- eeganalysis:::.topomap_grid(s$x, s$y, s$v, s$head_r, n_grid = 29)

  for (k in c(1000, 0.001)) {
    scaled <- eeganalysis:::.topomap_grid(s$x * k, s$y * k, s$v, s$head_r * k,
                                          n_grid = 29)
    expect_equal(scaled$z, base$z, tolerance = 1e-8)
  }

  for (k in c(1e-9, 1e9)) {
    scaled <- eeganalysis:::.topomap_grid(s$x, s$y, s$v * k, s$head_r, n_grid = 29)
    expect_equal(scaled$z / k, base$z, tolerance = 1e-8)
  }

  shifted <- eeganalysis:::.topomap_grid(s$x, s$y, s$v + 1e6, s$head_r, n_grid = 29)
  expect_equal(shifted$z - 1e6, base$z, tolerance = 1e-8)
})

# ----------------------------------------------------------------------------
# Test 2.7: A perfectly regular grid of electrodes
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a 5 x 5 lattice of electrodes, where many
# triangles are tied (four electrodes on one circle) and whole rows lie on a
# line, gives a map with a value in every cell of the head disc, no warning,
# and the electrode values at the electrodes.
test_that("a regular lattice of electrodes gives a complete, exact map", {
  lattice <- expand.grid(x = -2:2, y = -2:2)
  v <- sin(lattice$x) + cos(lattice$y) + 0.3 * lattice$x * lattice$y

  expect_no_warning(
    fit <- eeganalysis:::.topomap_grid(lattice$x, lattice$y, v, head_r = 3,
                                       n_grid = 41)
  )
  expect_true(all(is.finite(fit$z[sqrt(outer(fit$x^2, fit$y^2, "+")) <= 3])))

  surface <- eeganalysis:::.topo_surface(lattice$x, lattice$y, v, 3)
  expect_lt(max(abs(eeganalysis:::.topo_read(surface, surface$electrodes[, 1],
                                             surface$electrodes[, 2]) -
                      surface$values)), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.8: Many electrodes on one line, and one more beside it
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies a row of 7 electrodes on a straight line plus 2
# electrodes off it (like a midline strip of a cap) gives a map with a value
# in every cell of the head disc and passes through the electrode values (a
# row on a line must not give flat triangles).
test_that("a straight row of electrodes plus two beside it still gives a map", {
  x <- c(seq(-1.5, 1.5, length.out = 7), -0.4, 0.9)
  y <- c(rep(0, 7), 1.1, -0.8)
  v <- c(0.1, 0.5, -0.2, 0.9, 0.3, -0.6, 0.2, 1.0, -0.7)

  expect_no_warning(fit <- eeganalysis:::.topomap_grid(x, y, v, 2, n_grid = 31))
  expect_true(all(is.finite(fit$z[sqrt(outer(fit$x^2, fit$y^2, "+")) <= 2])))

  surface <- eeganalysis:::.topo_surface(x, y, v, 2)
  expect_lt(max(abs(eeganalysis:::.topo_read(surface, surface$electrodes[, 1],
                                             surface$electrodes[, 2]) -
                      surface$values)), 1e-9)
})

# ----------------------------------------------------------------------------
# Test 2.9: The grid
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the grid that is returned: n_grid values from
# -head_r to head_r in both directions, a matrix of that size in which z[i, j]
# belongs to the point (x[i], y[j]) (a value that is high on the right of the
# head sits in a high column index), NA in every corner outside the disc, and a
# value in the four cells in the middle of each side, which lie exactly on the
# edge of the disc and count as inside it.
test_that("the grid runs from -head_r to head_r and z[i, j] is the point (x[i], y[j])", {
  s <- small_set()
  tilt <- 2 * s$x + 0 * s$y

  fit <- eeganalysis:::.topomap_grid(s$x, s$y, tilt, s$head_r, n_grid = 25)

  expect_equal(fit$x, seq(-1.8, 1.8, length.out = 25))
  expect_equal(fit$y, fit$x)
  expect_equal(dim(fit$z), c(25, 25))
  expect_true(all(is.na(fit$z[c(1, 25), c(1, 25)])))
  expect_true(all(is.finite(fit$z[cbind(c(1, 25, 13, 13), c(13, 13, 1, 25))])))
  # the value follows x: it grows along the first index and not along the second
  middle <- fit$z[, 13]
  expect_gt(middle[19] - middle[7], 0)
  expect_lt(abs(fit$z[13, 19] - fit$z[13, 7]), abs(middle[19] - middle[7]))
})

# ============================================================================
#                  TEST SUITE 3: Errors
# ============================================================================

# ----------------------------------------------------------------------------
# Test 3.1: Fewer than 3 different positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an error when the electrodes sit at fewer than 3
# different positions, also when there are more than 3 electrodes (4
# electrodes in 2 places).
test_that("fewer than 3 different positions is an error", {
  expect_error(eeganalysis:::.topomap_grid(c(0, 1, 0, 1), c(0, 1, 0, 1),
                                           c(1, 2, 3, 4), 2),
               "At least 3 different electrode positions")
  expect_error(eeganalysis:::.topomap_grid(c(0, 1), c(0, 1), c(1, 2), 2),
               "At least 3 different electrode positions")
})

# ----------------------------------------------------------------------------
# Test 3.2: Electrodes on one line
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an error that says "collinear" when all electrodes
# are on one line, also a slanted line (which has a spread in both x and y).
test_that("electrodes on one line are an error", {
  expect_error(eeganalysis:::.topomap_grid(c(0, 1, 2, 3), c(0, 1, 2, 3),
                                           c(1, 2, 3, 4), 3),
               "collinear")
  expect_error(eeganalysis:::.topomap_grid(c(0, 0, 0, 0), c(-1, 0, 1, 2),
                                           c(1, 2, 3, 4), 3),
               "collinear")
})

# ----------------------------------------------------------------------------
# Test 3.3: Values or positions that are not numbers
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies an error that says "finite" when a value or a
# position is NA, NaN or infinite, instead of a map with holes in it.
test_that("a value or a position that is not a finite number is an error", {
  s <- small_set()
  bad_value <- replace(s$v, 4, NA)
  bad_inf   <- replace(s$v, 4, Inf)
  bad_pos   <- replace(s$x, 2, NaN)

  expect_error(eeganalysis:::.topomap_grid(s$x, s$y, bad_value, s$head_r), "finite")
  expect_error(eeganalysis:::.topomap_grid(s$x, s$y, bad_inf, s$head_r), "finite")
  expect_error(eeganalysis:::.topomap_grid(bad_pos, s$y, s$v, s$head_r), "finite")
})

# ============================================================================
#               TEST SUITE 4: The Ring and Its Spacing
# ============================================================================

# ----------------------------------------------------------------------------
# Test 4.1: The ring of extra points
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the ring: its first point lies on the x axis at
# 1.1 times the head radius plus the spacing, all points are on the same
# circle and equally far apart, their number is 2 pi / asin(spacing / head
# radius) rounded (here 31 for a spacing of 0.2 on a head of radius 1), never
# fewer than 12, and never more than 720 (also for a spacing so small that the
# plain number would not fit in an integer).
test_that("the ring has the right size, number of points and place", {
  ring <- eeganalysis:::.topo_ring(head_r = 1, spacing = 0.2)

  expect_equal(nrow(ring), 31)
  expect_equal(ring[1, ], c(1.1 + 0.2, 0))
  expect_equal(sqrt(rowSums(ring^2)), rep(1.3, 31))
  angles <- atan2(ring[, 2], ring[, 1]) %% (2 * pi)
  expect_equal(diff(angles), rep(2 * pi / 31, 30))

  expect_equal(nrow(eeganalysis:::.topo_ring(1, 0.8)), 12)
  expect_equal(nrow(eeganalysis:::.topo_ring(1, 5)), 12)
  expect_equal(nrow(eeganalysis:::.topo_ring(1, 1e-9)), 720)
})

# ----------------------------------------------------------------------------
# Test 4.2: The fallback spacing: every edge counted once
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the spacing worked out from our own mesh, which is
# used when Qhull cannot join the electrodes (Test 4.7), is the median of the
# lengths of the edges of the mesh with every edge counted ONCE. For a unit square the mesh
# has 4 sides of length 1 and one diagonal of length sqrt(2): counted once the
# median is 1, whereas counting the shared diagonal twice would give 1.207.
# Also checks a 4 x 4 lattice (spacing 1), where the sides outnumber the
# diagonals, and that the number scales with the size of the layout.
test_that("the electrode spacing is the median edge length, each edge once", {
  square <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  mesh <- eeganalysis:::.topo_delaunay(square)
  expect_equal(eeganalysis:::.topo_edge_median(square, mesh), 1)

  lattice <- as.matrix(expand.grid(0:3, 0:3))
  expect_equal(eeganalysis:::.topo_edge_median(lattice,
                                               eeganalysis:::.topo_delaunay(lattice)), 1)
  expect_equal(eeganalysis:::.topo_edge_median(lattice * 7,
                                               eeganalysis:::.topo_delaunay(lattice * 7)), 7)
})

# ----------------------------------------------------------------------------
# Test 4.3: The spacing on the BioSemi 64 electrode layout
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that for the 64 electrode BioSemi layout the
# package works out the spacing 0.453052 that the reference program works
# out, which gives a ring of 30 points, so that the maps of this layout agree
# with the reference numbers up to the edge of the head (Test 1.2).
test_that("the 64 electrode layout gets the spacing and ring of the reference program", {
  p <- layout_2d("biosemi64")
  surface <- eeganalysis:::.topo_surface(p$x, p$y, smooth_pattern(p$x, p$y), p$head_r)
  spacing <- eeganalysis:::.topo_edge_median(surface$electrodes,
                                             surface$electrode_mesh)

  expect_equal(spacing, 0.45305192797089022, tolerance = 1e-9)
  expect_equal(nrow(surface$points) - nrow(surface$electrodes), 30)
})

# ----------------------------------------------------------------------------
# Test 4.4: The spacing measured by Qhull is the spacing of the reference program
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that .topo_qhull_spacing() works out the same
# spacing as the reference program. The numbers were made once by the
# reference program on the electrodes of each layout. They are the layouts where
# the number does not move when the positions are nudged in their last digit
# (checked by nudging them by 1e-16 and 1e-14, 40 times each), so the test
# cannot fail on another computer by the way a tie is settled: the 16 and 64
# electrode BioSemi layouts, the 10-20 layout (where counting every edge once,
# as in Test 4.2, gives a ring with one point too few) and colin27_1020. The
# electrodes are given in the order of their positions, as .topo_surface() gives
# them.
test_that("the spacing measured by Qhull is the spacing of the reference program", {
  reference <- c(biosemi16     = 0.82643716454132354,
                 biosemi64     = 0.45305192797089022,
                 standard_1020 = 0.84040499302832883,
                 colin27_1020  = 0.42598863618569416)

  for (layout in names(reference)) {
    p <- layout_2d(layout)
    surface <- eeganalysis:::.topo_surface(p$x, p$y, smooth_pattern(p$x, p$y),
                                           p$head_r)
    expect_equal(eeganalysis:::.topo_qhull_spacing(surface$electrodes),
                 reference[[layout]], tolerance = 1e-9, label = layout)
  }
})

# ----------------------------------------------------------------------------
# Test 4.5: Every layout gets the ring of the reference program
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that on 8 layouts of the package the ring of extra
# points has as many points as the reference program puts around the same
# electrodes (the numbers were made once by the reference program). The 10-20
# layout (13 points) and the 160 electrode BioSemi layout (50 points) are the
# two where the simple rule of Test 4.2 gives one or two points too few. Left
# out are the three layouts whose ring moves by a point when the positions are
# nudged in their last digit (the 256 electrode BioSemi layout, 10-10 and
# 10-05: 61, 30 and 66 points here), because a test on them could fail on
# another computer without anything being wrong.
test_that("the layouts with a steady ring get the ring of the reference program", {
  ring_points <- c(biosemi16 = 13, biosemi32 = 17, biosemi64 = 30,
                   biosemi128 = 42, biosemi160 = 50,
                   standard_1020 = 13,
                   colin27_1020 = 36, colin27_1005 = 73)

  for (layout in names(ring_points)) {
    p <- layout_2d(layout)
    surface <- eeganalysis:::.topo_surface(p$x, p$y, smooth_pattern(p$x, p$y),
                                           p$head_r)
    expect_equal(nrow(surface$points) - nrow(surface$electrodes),
                 ring_points[[layout]], label = layout)
  }
})

# ----------------------------------------------------------------------------
# Test 4.6: The spacing follows the unit of the positions
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the spacing measured by Qhull grows and
# shrinks with the positions (the 64 electrode layout in thousand times larger
# and smaller units, such as millimetres and kilometres instead of metres), so
# the ring does not depend on the unit the positions are given in.
test_that("the spacing measured by Qhull follows the unit of the positions", {
  p <- layout_2d("biosemi64")
  electrodes <- cbind(p$x, p$y)
  spacing <- eeganalysis:::.topo_qhull_spacing(electrodes)

  expect_true(is.finite(spacing) && spacing > 0)
  expect_equal(eeganalysis:::.topo_qhull_spacing(electrodes * 1000), spacing * 1000)
  expect_equal(eeganalysis:::.topo_qhull_spacing(electrodes / 1000), spacing / 1000)
})

# ----------------------------------------------------------------------------
# Test 4.7: Electrodes that Qhull cannot join fall back on our own spacing
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Qhull refuses to join electrodes that lie exactly on one
# circle (the 4 corners of a square). Verifies that .topo_qhull_spacing() then
# gives NA instead of an error, and that the map is still drawn, with the
# spacing of Test 4.2 (1 for a unit square), so the result is the same as when
# that spacing is given.
test_that("electrodes that Qhull cannot join fall back on the edge median", {
  square <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  values <- c(1, -1, 2, 0.5)

  expect_true(is.na(eeganalysis:::.topo_qhull_spacing(square)))

  by_default <- eeganalysis:::.topomap_grid(square[, 1], square[, 2], values,
                                            head_r = 2)$z
  pinned <- eeganalysis:::.topomap_grid(square[, 1], square[, 2], values,
                                        head_r = 2, ring_spacing = 1)$z
  expect_equal(by_default, pinned)
  expect_true(any(is.finite(by_default)))
})

# ============================================================================
#                 TEST SUITE 5: The Mesh of Triangles
# ============================================================================

# ----------------------------------------------------------------------------
# Test 5.1: A unit square becomes two triangles
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the smallest mesh: 4 corners of a square give 2
# triangles that use all 4 corners and together cover the square (area 1).
test_that("the corners of a square give two triangles", {
  square <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  mesh <- eeganalysis:::.topo_delaunay(square)

  expect_equal(dim(mesh), c(2, 3))
  expect_equal(sort(unique(as.vector(mesh))), 1:4)
  p <- mesh_problems(square, mesh)
  expect_equal(p$area_ratio, 1)
})

# ----------------------------------------------------------------------------
# Test 5.2: A valid mesh on awkward sets of points
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies on 5 kinds of point sets (random points, a regular
# lattice, a row plus one point, a small cluster inside a wider scatter, and a
# small set far from the origin) that the mesh covers exactly the convex hull
# of the points (area ratio 1), uses every point, and has no point inside the
# circle through the corners of any triangle (to 1e-6 relative).
test_that("the mesh is a valid Delaunay mesh on awkward sets of points", {
  set.seed(42)
  sets <- list(
    random  = cbind(runif(40, -2, 2), runif(40, -2, 2)),
    lattice = as.matrix(expand.grid(1:6, 1:5)),
    row     = cbind(c(seq(-2, 2, length.out = 9), 0.3), c(rep(0, 9), 1)),
    cluster = rbind(cbind(rnorm(15, 0, 0.05), rnorm(15, 0, 0.05)),
                    cbind(runif(25, -2, 2), runif(25, -2, 2))),
    far     = cbind(runif(25, 1e6, 1e6 + 1), runif(25, -2e6, -2e6 + 1))
  )

  for (nm in names(sets)) {
    pts <- sets[[nm]]
    p <- mesh_problems(pts, eeganalysis:::.topo_delaunay(pts))
    expect_equal(p$unused, 0, info = nm)
    expect_equal(p$area_ratio, 1, tolerance = 1e-6, info = nm)
    expect_lt(p$violation, 1e-6)
  }
})

# ----------------------------------------------------------------------------
# Test 5.3: Neighbouring triangles
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies for the two triangles of a square that each has the
# other as its one neighbour (opposite the corner they do not share) and NA
# for the other two corners, and that on a bigger mesh the relation is
# symmetric: if triangle b is across an edge of triangle a, then a is across
# an edge of b. Finally, that none is missed: every edge that is not on the
# outer boundary of the mesh is shared by two triangles, so the number of
# neighbours found is 3 times the number of triangles minus the number of edges
# on the boundary (for random points, the number of corners of their hull).
test_that("neighbouring triangles are found across shared edges", {
  square <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  mesh <- eeganalysis:::.topo_delaunay(square)
  nb <- eeganalysis:::.topo_triangle_neighbours(mesh)

  expect_equal(dim(nb), c(2, 3))
  expect_equal(sum(!is.na(nb)), 2)
  expect_equal(nb[1, which(!is.na(nb[1, ]))], 2)
  expect_equal(nb[2, which(!is.na(nb[2, ]))], 1)

  pts <- as.matrix(expand.grid(1:5, 1:4))
  mesh <- eeganalysis:::.topo_delaunay(pts)
  nb <- eeganalysis:::.topo_triangle_neighbours(mesh)
  for (a in seq_len(nrow(mesh))) {
    for (b in stats::na.omit(nb[a, ])) {
      expect_true(a %in% nb[b, ])
      expect_length(intersect(mesh[a, ], mesh[b, ]), 2)
    }
  }

  set.seed(3)
  scatter <- cbind(runif(40), runif(40))
  mesh <- eeganalysis:::.topo_delaunay(scatter)
  nb <- eeganalysis:::.topo_triangle_neighbours(mesh)
  expect_equal(sum(!is.na(nb)),
               3 * nrow(mesh) - length(grDevices::chull(scatter)))

  # A hand-made mesh in which the shared edge (points 2 and 3) is listed in
  # opposite directions by the two triangles: it must still be found.
  by_hand <- rbind(c(1, 2, 3), c(4, 3, 2))
  nb <- eeganalysis:::.topo_triangle_neighbours(by_hand)
  expect_equal(nb[1, 1], 2)
  expect_equal(nb[2, 1], 1)
  expect_equal(sum(!is.na(nb)), 2)
})

# ----------------------------------------------------------------------------
# Test 5.4: Finding the triangle a point lies in
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies, for the two triangles of a unit square, that a
# point inside gets the triangle it lies in with weights that add up to 1 and
# put it back where it was, that a corner gets weight 1 on that corner, that a
# point on the outer edge is found, that a point 1e-10 outside the edge is still
# placed (a point that lies on an edge can end up a rounding error outside it,
# and must not be lost), and that a point 1e-6 outside, and any point outside
# the mesh, gets NA.
test_that("points are placed in the triangle they lie in", {
  square <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  mesh <- eeganalysis:::.topo_delaunay(square)

  qx <- c(0.3, 0.8, 0, 1, 0.5, 1.5, -0.2)
  qy <- c(0.2, 0.9, 0, 1, 0, 0.5, 0.5)
  found <- eeganalysis:::.topo_locate(square, mesh, qx, qy)

  expect_false(anyNA(found$triangle[1:5]))
  expect_true(all(is.na(found$triangle[6:7])))
  weights <- found$weights[1:5, ]
  expect_equal(rowSums(weights), rep(1, 5))
  for (k in 1:5) {
    corners <- square[mesh[found$triangle[k], ], ]
    expect_equal(as.vector(weights[k, ] %*% corners), c(qx[k], qy[k]))
  }
  expect_equal(max(found$weights[3, ]), 1)
  expect_equal(max(found$weights[4, ]), 1)

  near <- eeganalysis:::.topo_locate(square, mesh, c(0.5, 0.5), c(-1e-10, -1e-6))
  expect_false(is.na(near$triangle[1]))
  expect_true(is.na(near$triangle[2]))
})

# ----------------------------------------------------------------------------
# Test 5.5: The bend across an edge is never undefined
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies that the number the Clough-Tocher scheme divides by
# when it chooses how the surface bends across an edge (den) is -1 or less for
# every edge of the mesh, on a lattice (many tied triangles), on a ring of
# points on one circle, and on random points. It is 3 * c - 1 with c zero or
# negative, so it can never be zero.
test_that("the number divided by when choosing the bend across an edge is at most -1", {
  set.seed(11)
  sets <- list(
    lattice = as.matrix(expand.grid(1:7, 1:6)),
    circle  = rbind(cbind(cos(seq(0, 2 * pi, length.out = 21)[-21]),
                          sin(seq(0, 2 * pi, length.out = 21)[-21])), c(0, 0)),
    random  = cbind(runif(50), runif(50))
  )

  for (nm in names(sets)) {
    pts <- sets[[nm]]
    mesh <- eeganalysis:::.topo_delaunay(pts)
    terms <- eeganalysis:::.topo_edge_terms(
      pts, mesh, eeganalysis:::.topo_triangle_neighbours(mesh))
    expect_lte(max(terms$den, na.rm = TRUE), -1 + 1e-9)
  }
})

# ----------------------------------------------------------------------------
# Test 5.6: The slopes settle, and say so when they do not
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the slope estimate gives the slope of a plane
# exactly (a plane has no bending, so slope (a, b) everywhere, to 1e-5), and
# warns, instead of staying silent, when it is stopped after one sweep before
# it has settled.
test_that("the slopes of a plane are the plane's slope, and a warning says when they do not settle", {
  set.seed(7)
  pts <- cbind(runif(25, -2, 2), runif(25, -2, 2))
  f <- 1.5 * pts[, 1] - 0.7 * pts[, 2] + 2
  mesh <- eeganalysis:::.topo_delaunay(pts)
  joined <- eeganalysis:::.topo_vertex_neighbours(mesh, nrow(pts))

  slopes <- eeganalysis:::.topo_slopes(pts, joined, f)
  expect_lt(max(abs(slopes[, 1] - 1.5)), 1e-5)
  expect_lt(max(abs(slopes[, 2] + 0.7)), 1e-5)

  g <- sin(pts[, 1]) * cos(pts[, 2])
  expect_warning(eeganalysis:::.topo_slopes(pts, joined, g, max_sweeps = 1),
                 "did not settle")
})

# ============================================================================
#     TEST SUITE 6: What the Package Does and Does Not Depend On
# ============================================================================

# ----------------------------------------------------------------------------
# Test 6.1: akima is not used any more
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the package does not import akima (the
# interpolation package whose licence limits commercial use) in its namespace
# and does not list it as a dependency in DESCRIPTION.
test_that("the package does not depend on akima", {
  expect_false("akima" %in% names(getNamespaceImports("eeganalysis")))

  description <- system.file("DESCRIPTION", package = "eeganalysis")
  skip_if(!nzchar(description), "DESCRIPTION is not available here")
  fields <- read.dcf(description, fields = c("Imports", "Depends", "LinkingTo"))
  expect_false(any(grepl("akima", fields, fixed = TRUE)))
})

# ----------------------------------------------------------------------------
# Test 6.2: geometry is declared, because the ring spacing needs it
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Verifies the package lists geometry (the package that
# carries Qhull) as an import in DESCRIPTION, and that it can be used, since
# .topo_qhull_spacing() calls it.
test_that("the package declares geometry, which the ring spacing uses", {
  description <- system.file("DESCRIPTION", package = "eeganalysis")
  skip_if(!nzchar(description), "DESCRIPTION is not available here")
  fields <- read.dcf(description, fields = "Imports")
  expect_true(grepl("geometry", fields, fixed = TRUE))
  expect_true(requireNamespace("geometry", quietly = TRUE))
})

# ============================================================================
#                     SUMMARY OF TEST COVERAGE
# ============================================================================
# - Reference numbers from the interpolator MNE-Python uses: a small irregular
#   set (ring spacing pinned), the 16, 32 and 64 electrode BioSemi layouts, and
#   the 10-10 and colin27_1020 layouts (the latter has electrodes that share a
#   position)
# - Properties: passes through every electrode value, same map in any order,
#   shared positions averaged, extrapolate = FALSE only empties the cells
#   outside the electrodes, constant values give a constant map, no dependence
#   on the unit of the positions, regular lattices, a straight row of
#   electrodes, the grid layout
# - Errors: fewer than 3 positions, electrodes on one line, values or positions
#   that are not finite
# - The ring: size, number of points, limits; the spacing measured by Qhull
#   (the reference numbers of 4 layouts, the number of ring points of the 8
#   layouts whose ring is steady, the unit of the positions, the way out when
#   Qhull cannot join the electrodes); the fallback spacing as the median edge
#   length with every edge counted once; the spacing and ring of the 64
#   electrode layout
# - The mesh: a square, validity on awkward point sets, neighbouring triangles,
#   finding the triangle of a point, the slope estimate and its warning
# - akima is gone from the namespace and from DESCRIPTION
# ============================================================================
