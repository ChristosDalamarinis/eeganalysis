# ============================================================================
#               The Smooth Surface Behind plot_topography()
# ============================================================================
#
# plot_topography() knows one value per electrode, but it has to colour every
# pixel of the head. The functions in this file fill in the pixels: they turn
# the values at the electrodes into a smooth surface and read that surface at
# every point of a regular grid. Nothing here is exported; plot_topography()
# calls .topomap_grid(), and the other functions are its parts.
#
# The method is the one used for topographic maps in MNE-Python (mne/viz/
# topomap.py), written here in plain R. One number is not worked out in R: the
# typical distance between electrodes, which sets the number of points of the
# ring (step 2). It is measured on the triangles of the Qhull library, through
# the R package geometry, because the way MNE-Python measures it depends on the
# inner order of Qhull (see .topo_qhull_spacing()):
#
#   1. The electrodes are joined into a mesh of triangles (a Delaunay
#      triangulation: the circle through the corners of any triangle holds no
#      other electrode).
#   2. A ring of extra points is placed around the head, a little outside it.
#      Each ring point gets the mean value of the electrodes it is joined to,
#      so the colours fade out smoothly towards the edge of the head instead
#      of running away beyond the outer electrodes, and the surface bends
#      sensibly at the outer electrodes themselves.
#   3. A cubic surface is fitted over every triangle (the Clough-Tocher
#      scheme). It passes through every value, and the pieces join without a
#      kink, because the slope at every point is chosen so that the surface
#      bends as little as possible overall.
#   4. The surface is read at the points of the grid that lie inside the head.
#      With extrapolate = FALSE the cells outside the outer boundary of the
#      electrodes are then left empty; the surface itself does not change.
#
# Credits and licence: the cubic fitting (step 3, with its estimate of the
# slopes and the way a point is found in the mesh and read off) is a port of
# the Clough-Tocher interpolator of scipy (scipy/interpolate/_interpnd.pyx),
# and the ring of extra points (step 2) follows MNE-Python's mne/viz/
# topomap.py. Both are released under the BSD 3-Clause licence, which asks
# that the notice is kept. The notice for scipy is:
#
#   Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions
#   are met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
#   3. Neither the name of the copyright holder nor the names of its
#      contributors may be used to endorse or promote products derived
#      from this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The mesh itself is built here with a small triangulation routine of our own
# (.topo_delaunay()) instead of the library routine of the models. It uses a
# test for "inside the circle" that stays accurate for thin triangles, so that
# electrodes that lie exactly on a line (such as the midline of a cap) or
# exactly on a circle (such as a regular lattice) give a proper mesh without
# any special treatment.


# ============================================================================
#                  THE FUNCTION plot_topography() CALLS
# ============================================================================

# Turns the values at the electrodes into a smooth surface and reads it on a
# regular grid.
#
# x, y         the 2D positions of the electrodes (numeric vectors)
# values       one value per electrode (numeric vector, all finite)
# head_r       the radius of the head disc; the grid runs from -head_r to
#              head_r in both directions
# n_grid       the number of grid points per side
# extrapolate  TRUE gives every cell of the head disc a value; FALSE leaves the
#              cells outside the outer boundary of the electrodes empty (NA),
#              and every other cell keeps exactly the value it has with TRUE
# ring_spacing the distance between neighbouring ring points; NULL (the
#              default) uses the typical distance between electrodes, measured
#              by .topo_qhull_spacing(). Only given to fix the ring at a chosen
#              value (the tests do this).
#
# Returns a list with x and y (the grid positions) and z (a matrix with
# z[i, j] for the point (x[i], y[j])); z is NA outside the head disc and, when
# extrapolate = FALSE, outside the outer boundary of the electrodes.
#
# The result does not depend on the order of the electrodes: they are put in
# the order of their positions first, and electrodes that share a position
# are replaced by one electrode with the mean of their values.
.topomap_grid <- function(x, y, values, head_r, n_grid = 100,
                          extrapolate = TRUE, ring_spacing = NULL) {
  surface <- .topo_surface(x, y, values, head_r, ring_spacing)

  # The grid points inside the head disc. The disc is worked out on a grid
  # from -1 to 1, so a cell on the edge of the disc is in or out whatever the
  # unit of the positions is.
  unit_grid <- seq(-1, 1, length.out = n_grid)
  inside    <- which(sqrt(outer(unit_grid^2, unit_grid^2, "+")) <= 1)
  grid_x    <- unit_grid[(inside - 1) %% n_grid + 1] * head_r
  grid_y    <- unit_grid[(inside - 1) %/% n_grid + 1] * head_r

  z <- matrix(NA_real_, n_grid, n_grid)
  z[inside] <- .topo_read(surface, grid_x, grid_y)

  # Without extrapolation the cells outside the outer boundary of the
  # electrodes are left empty; the other cells keep the value they have.
  if (!isTRUE(extrapolate)) {
    z[inside[!.topo_on_electrodes(surface, grid_x, grid_y)]] <- NA
  }

  list(x = unit_grid * head_r, y = unit_grid * head_r, z = z)
}

# Builds the smooth surface from the values at the electrodes (see the steps at
# the top of the file). Arguments as for .topomap_grid().
#
# Everything inside is worked out with the positions divided by the head
# radius and the values shifted by their mean and divided by their largest
# distance from it. The slope estimate stops when its change falls below a
# fixed small number, which only means something for numbers of about size 1;
# without this, values in volts (around 1e-6) would stop it after the first
# sweep and give a badly fitted map. The surface is linear in the values, so
# the shift and the scale are undone when the surface is read.
#
# Returns a list with
#   electrodes       the electrode positions that were used (one row per
#                    position, in the order of the positions), in the
#                    original unit
#   values           their values, in the original unit
#   unit, centre, spread   the head radius, the mean value and the largest
#                    distance from it, which the surface is scaled by
#   electrode_mesh   the mesh of the electrodes alone
#   points, mesh     the electrodes and the ring (scaled), and the mesh over them
#   pieces           the cubic pieces of the surface (see .topo_pieces())
.topo_surface <- function(x, y, values, head_r, ring_spacing = NULL) {

  if (any(!is.finite(c(x, y, values)))) {
    stop("ERROR: The electrode positions and the values must be finite ",
         "numbers (no NA, NaN or Inf).", call. = FALSE)
  }

  # One point per position, in the order of the positions (left to right, then
  # front to back). Positions that agree to 8 decimals are one position.
  round_x     <- round(x, 8)
  round_y     <- round(y, 8)
  by_position <- order(round_x, round_y)
  x           <- x[by_position]
  y           <- y[by_position]
  round_x     <- round_x[by_position]
  round_y     <- round_y[by_position]
  values      <- values[by_position]

  first_at_position <- c(TRUE, diff(round_x) != 0 | diff(round_y) != 0)
  group  <- cumsum(first_at_position)
  values <- as.vector(rowsum(values, group)) / tabulate(group)
  x      <- x[first_at_position]
  y      <- y[first_at_position]
  electrodes <- cbind(x, y)

  if (nrow(electrodes) < 3) {
    stop("ERROR: At least 3 different electrode positions are required to ",
         "interpolate a topography (found ", nrow(electrodes), ").",
         call. = FALSE)
  }
  singular <- svd(sweep(electrodes, 2, colMeans(electrodes)),
                  nu = 0, nv = 0)$d
  if (singular[2] < 1e-8 * singular[1]) {
    stop("ERROR: The selected channels are collinear (all on one line) and ",
         "cannot be interpolated into a surface. Include channels spread ",
         "across both hemispheres and the anterior-posterior axis.",
         call. = FALSE)
  }

  # The scaling (see above).
  unit   <- head_r
  centre <- mean(values)
  spread <- max(abs(values - centre))
  if (spread == 0) spread <- 1
  scaled_electrodes <- electrodes / unit
  scaled_values     <- (values - centre) / spread

  # The mesh of the electrodes alone. It is used, with extrapolate = FALSE, to
  # find the cells that lie outside the electrodes, and to find the spacing for
  # the ring when Qhull cannot measure it (see .topo_qhull_spacing()).
  electrode_mesh <- .topo_delaunay(scaled_electrodes)

  # The ring is spaced like the electrodes, and the mesh is built again with
  # the ring points in it. The ring is there whatever `extrapolate` is: it
  # decides the slope of the surface at the outer electrodes, and the surface
  # is closer to the truth inside the electrodes with it than without it.
  # The spacing is measured on the electrodes in their own unit, which is what
  # the reference program does, and then scaled like everything else.
  spacing <- if (!is.null(ring_spacing)) {
    ring_spacing / unit
  } else {
    typical <- .topo_qhull_spacing(electrodes)
    if (is.na(typical)) {
      .topo_edge_median(scaled_electrodes, electrode_mesh)
    } else {
      typical / unit
    }
  }
  ring   <- .topo_ring(1, spacing)
  points <- rbind(scaled_electrodes, ring)
  mesh   <- .topo_delaunay(points)

  n_electrodes <- nrow(electrodes)
  joined_to    <- .topo_vertex_neighbours(mesh, nrow(points))

  # Each ring point gets the mean of the electrodes it is joined to.
  ring_values <- rep(NA_real_, nrow(ring))
  for (k in seq_len(nrow(ring))) {
    near <- joined_to[[n_electrodes + k]]
    near <- near[near <= n_electrodes]
    if (length(near) > 0) ring_values[k] <- mean(scaled_values[near])
  }
  no_electrode <- is.na(ring_values)
  ring_values[no_electrode] <- if (all(no_electrode)) {
    mean(scaled_values)
  } else {
    mean(ring_values[!no_electrode])
  }
  all_values <- c(scaled_values, ring_values)

  slopes <- .topo_slopes(points, joined_to, all_values)
  pieces <- .topo_pieces(points, mesh, .topo_triangle_neighbours(mesh),
                         all_values, slopes)

  list(electrodes = electrodes, values = values,
       unit = unit, centre = centre, spread = spread,
       electrode_mesh = electrode_mesh,
       points = points, mesh = mesh, pieces = pieces)
}

# Reads the surface at the points (qx[k], qy[k]), given in the original unit;
# NA for a point outside the mesh (which, with the ring, is any point beyond
# the ring).
.topo_read <- function(surface, qx, qy) {
  located <- .topo_locate(surface$points, surface$mesh,
                          qx / surface$unit, qy / surface$unit)
  .topo_evaluate(surface$pieces, located) * surface$spread + surface$centre
}

# Is each point (qx[k], qy[k]) inside the outer boundary of the electrodes?
.topo_on_electrodes <- function(surface, qx, qy) {
  electrodes <- surface$points[seq_len(nrow(surface$electrodes)), , drop = FALSE]
  located <- .topo_locate(electrodes, surface$electrode_mesh,
                          qx / surface$unit, qy / surface$unit)
  !is.na(located$triangle)
}


# ============================================================================
#                       THE RING OF EXTRA POINTS
# ============================================================================

# The typical distance between neighbouring electrodes, measured the way the
# program this follows (MNE-Python) measures it: the electrodes are joined into
# triangles by Qhull (the library inside scipy; the R package geometry carries
# the same library), with the same options, and the spacing is the median of
# the lengths of two of the three edges of every triangle, the edge from the
# first to the second corner and the edge from the second to the third corner.
#
# Which edge is left out depends on the order in which Qhull lists the corners
# of a triangle. That order comes from the inner working of Qhull (the order
# in which it created its points) and cannot be worked out from outside, which
# is why the number is measured by Qhull itself and not in R. scipy also
# swaps the first two corners of every triangle that Qhull lists clockwise, so
# that the triangle runs anticlockwise; the same is done here, because it
# decides which edge is left out.
#
# electrodes  a matrix with the x and y position of every electrode, one row
#             per position, in the order of the positions (the same order the
#             electrodes are in inside .topo_surface())
#
# Returns the spacing in the unit of `electrodes`; NA when Qhull cannot join
# the electrodes (for example 4 electrodes that lie exactly on one circle),
# and then .topo_edge_median() is used instead.
#
# A limit that comes with copying the method: where electrodes lie on a very
# regular pattern, neighbouring triangles are tied (four electrodes on one
# circle), and Qhull settles a tie by the last digit of the positions. On the
# BioSemi layout with 256 electrodes and on the 10-10 and 10-05 layouts, a
# change of one unit in the last digit of the positions can move the spacing by
# one or two percent, and so the ring by one or two points. The program this
# follows has the same limit, so on those layouts the ring can differ from its
# ring by a point on another computer. On the other layouts the number does not
# move.
.topo_qhull_spacing <- function(electrodes) {
  electrodes <- unname(electrodes)
  mesh <- tryCatch(geometry::delaunayn(electrodes, options = "Qc Qt Q11"),
                   error = function(e) NULL)
  if (!is.matrix(mesh) || ncol(mesh) != 3 || nrow(mesh) == 0) return(NA_real_)

  x <- electrodes[, 1]
  y <- electrodes[, 2]
  turn <- (x[mesh[, 2]] - x[mesh[, 1]]) * (y[mesh[, 3]] - y[mesh[, 1]]) -
    (y[mesh[, 2]] - y[mesh[, 1]]) * (x[mesh[, 3]] - x[mesh[, 1]])
  clockwise <- turn < 0
  mesh[clockwise, c(1, 2)] <- mesh[clockwise, c(2, 1)]

  edge_length <- function(a, b) sqrt((x[a] - x[b])^2 + (y[a] - y[b])^2)
  spacing <- median(c(edge_length(mesh[, 1], mesh[, 2]),
                      edge_length(mesh[, 2], mesh[, 3])))
  if (!is.finite(spacing) || spacing <= 0) return(NA_real_)
  spacing
}

# The typical distance between neighbouring electrodes, worked out from our own
# mesh: the median length of the edges of the mesh, every edge counted once.
# It is the way out for .topo_qhull_spacing() when Qhull cannot join the
# electrodes. It gives the same ring as Qhull's number on the BioSemi layouts
# with 16, 32, 64, 128 and 256 electrodes; on a few other layouts the ring has
# one or two points more or fewer, which changes the colours only in the band
# between the outer electrodes and the edge of the head.
.topo_edge_median <- function(points, mesh) {
  edges <- rbind(mesh[, c(1, 2)], mesh[, c(2, 3)], mesh[, c(1, 3)])
  edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  edges <- edges[!duplicated(edges[, 1] * (max(edges) + 1) + edges[, 2]), ,
                 drop = FALSE]
  median(sqrt(rowSums((points[edges[, 1], , drop = FALSE] -
                         points[edges[, 2], , drop = FALSE])^2)))
}

# Points on a circle a little outside the head: the head radius plus 10 percent
# plus one electrode spacing. They are placed so that neighbouring ring points
# are about one electrode spacing apart, but never fewer than 12 (and never
# more than 720, which no electrode layout comes near).
.topo_ring <- function(head_r, spacing) {
  angle  <- asin(min(spacing / head_r, 1))
  n_ring <- as.integer(min(720, max(12, round(2 * pi / angle))))
  theta  <- seq(0, 2 * pi, length.out = n_ring + 1)[seq_len(n_ring)]
  radius <- 1.1 * head_r + spacing
  cbind(radius * cos(theta), radius * sin(theta))
}


# ============================================================================
#                        THE MESH OF TRIANGLES
# ============================================================================

# Delaunay triangulation of a set of 2D points (Bowyer-Watson: the points are
# added one at a time; every triangle whose circle holds the new point is cut
# away, and the hole is filled with triangles that join the new point to the
# edge of the hole).
#
# Returns a matrix with one row per triangle and the numbers of its three
# corners (rows of `points`).
.topo_delaunay <- function(points) {
  n <- nrow(points)

  # The points are centred and scaled to about -1 to 1.
  centred <- sweep(points, 2, colMeans(points))
  scaled  <- centred / max(abs(centred))

  # A triangle far larger than all the points to start from; its three corners
  # are points n + 1 to n + 3, and every triangle that uses one is dropped at
  # the end.
  far  <- 1e5
  pts  <- rbind(scaled, c(-far, -far), c(far, -far), c(0, far))
  turn <- function(tri) {
    sign((pts[tri[, 2], 1] - pts[tri[, 1], 1]) *
           (pts[tri[, 3], 2] - pts[tri[, 1], 2]) -
         (pts[tri[, 2], 2] - pts[tri[, 1], 2]) *
           (pts[tri[, 3], 1] - pts[tri[, 1], 1]))
  }
  tris <- matrix(c(n + 1L, n + 2L, n + 3L), 1)
  side <- turn(tris)

  for (i in seq_len(n)) {
    # Is the new point inside the circle of a triangle? The sign of a 3 x 3
    # determinant tells (it does not need the centre of the circle, so thin
    # triangles stay accurate); points on the circle do not count as inside.
    ax <- pts[tris[, 1], 1] - pts[i, 1]
    ay <- pts[tris[, 1], 2] - pts[i, 2]
    bx <- pts[tris[, 2], 1] - pts[i, 1]
    b_y <- pts[tris[, 2], 2] - pts[i, 2]
    cx <- pts[tris[, 3], 1] - pts[i, 1]
    cy <- pts[tris[, 3], 2] - pts[i, 2]
    a2 <- ax^2 + ay^2
    b2 <- bx^2 + b_y^2
    c2 <- cx^2 + cy^2
    inside_sign <- a2 * (bx * cy - cx * b_y) + b2 * (cx * ay - ax * cy) +
      c2 * (ax * b_y - bx * ay)
    size <- a2 * (abs(bx * cy) + abs(cx * b_y)) +
      b2 * (abs(cx * ay) + abs(ax * cy)) +
      c2 * (abs(ax * b_y) + abs(bx * ay))
    broken <- side * inside_sign > 1e-13 * size
    if (!any(broken)) next

    # The edge of the hole is every edge that belongs to exactly one of the
    # triangles that are cut away.
    gone  <- tris[broken, , drop = FALSE]
    edges <- rbind(gone[, c(1, 2)], gone[, c(2, 3)], gone[, c(1, 3)])
    edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
    key   <- edges[, 1] * (n + 4) + edges[, 2]
    rim   <- !(duplicated(key) | duplicated(key, fromLast = TRUE))
    added <- cbind(edges[rim, , drop = FALSE], i)
    tris  <- rbind(tris[!broken, , drop = FALSE], added)
    side  <- c(side[!broken], turn(added))
  }

  tris[rowSums(tris > n) == 0, , drop = FALSE]
}

# For every point, the numbers of the points it shares a triangle edge with
# (a list with one sorted vector per point).
.topo_vertex_neighbours <- function(mesh, n_points) {
  pairs <- rbind(mesh[, c(1, 2)], mesh[, c(2, 1)], mesh[, c(2, 3)],
                 mesh[, c(3, 2)], mesh[, c(1, 3)], mesh[, c(3, 1)])
  lapply(split(pairs[, 2], factor(pairs[, 1], levels = seq_len(n_points))),
         function(v) sort(unique(v)))
}

# For every triangle, the triangle across the edge opposite each corner (a
# matrix with one row per triangle and one column per corner; NA where the edge
# is on the outer boundary of the mesh).
.topo_triangle_neighbours <- function(mesh) {
  n_tri <- nrow(mesh)
  a   <- c(mesh[, 2], mesh[, 1], mesh[, 1])   # the edge opposite corner 1, 2, 3
  b   <- c(mesh[, 3], mesh[, 3], mesh[, 2])
  key <- pmin(a, b) * (max(mesh) + 1) + pmax(a, b)
  triangle <- rep(seq_len(n_tri), 3)
  corner   <- rep(1:3, each = n_tri)

  later    <- which(duplicated(key))          # the second time an edge shows up
  earlier  <- match(key[later], key)          # the first time
  neighbour <- matrix(NA_integer_, n_tri, 3)
  neighbour[cbind(triangle[later], corner[later])]     <- triangle[earlier]
  neighbour[cbind(triangle[earlier], corner[earlier])] <- triangle[later]
  neighbour
}

# For every triangle and every corner: the numbers that the Clough-Tocher
# scheme needs to choose how the surface bends across the edge opposite that
# corner (`num` and `den`; the direction is num / den). They are worked out
# from where the middle of the triangle across the edge lies, measured from the
# corners of this triangle. NA where there is no triangle across the edge.
#
# `den` is never zero: it is 3 * c - 1, where c is the weight of the corner
# opposite the edge for a point on the far side of that edge, so c is zero or
# negative and `den` is -1 or less.
.topo_edge_terms <- function(points, mesh, neighbour) {
  x1 <- points[mesh[, 1], 1]; y1 <- points[mesh[, 1], 2]
  x2 <- points[mesh[, 2], 1]; y2 <- points[mesh[, 2], 2]
  x3 <- points[mesh[, 3], 1]; y3 <- points[mesh[, 3], 2]
  area_2 <- (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)
  mid_x  <- (x1 + x2 + x3) / 3
  mid_y  <- (y1 + y2 + y3) / 3

  num <- den <- matrix(NA_real_, nrow(mesh), 3)
  for (k in 1:3) {
    has <- !is.na(neighbour[, k])
    if (!any(has)) next
    other <- neighbour[has, k]
    px <- mid_x[other]
    py <- mid_y[other]
    l1 <- ((y2[has] - y3[has]) * (px - x3[has]) +
             (x3[has] - x2[has]) * (py - y3[has])) / area_2[has]
    l2 <- ((y3[has] - y1[has]) * (px - x3[has]) +
             (x1[has] - x3[has]) * (py - y3[has])) / area_2[has]
    l3 <- 1 - l1 - l2
    if (k == 1) {
      num[has, k] <- 2 * l3 + l2 - 1
      den[has, k] <- 2 - 3 * l3 - 3 * l2
    } else if (k == 2) {
      num[has, k] <- 2 * l1 + l3 - 1
      den[has, k] <- 2 - 3 * l1 - 3 * l3
    } else {
      num[has, k] <- 2 * l2 + l1 - 1
      den[has, k] <- 2 - 3 * l2 - 3 * l1
    }
  }
  list(num = num, den = den)
}


# ============================================================================
#            THE CUBIC SURFACE (CLOUGH-TOCHER, PORTED FROM SCIPY)
# ============================================================================

# The slope (the change of the value per unit of x and of y) of the surface at
# every point. The slopes are chosen so that the surface bends as little as
# possible: each point's slope is set from the values and slopes around it,
# and the sweep over all points is repeated until nothing changes any more.
#
# Returns a matrix with one row per point and two columns (x slope, y slope).
.topo_slopes <- function(points, joined_to, f, max_sweeps = 400,
                         tolerance = 1e-6) {
  n <- nrow(points)
  slope <- matrix(0, n, 2)

  dx <- dy <- reach <- vector("list", n)
  for (i in seq_len(n)) {
    near     <- joined_to[[i]]
    dx[[i]]  <- points[near, 1] - points[i, 1]
    dy[[i]]  <- points[near, 2] - points[i, 2]
    reach[[i]] <- (dx[[i]]^2 + dy[[i]]^2)^1.5
  }

  for (sweep_no in seq_len(max_sweeps)) {
    biggest_change <- 0
    for (i in seq_len(n)) {
      near <- joined_to[[i]]
      ex <- dx[[i]]
      ey <- dy[[i]]
      l3 <- reach[[i]]

      along <- -ex * slope[near, 1] - ey * slope[near, 2]
      pull  <- (6 * (f[i] - f[near]) - 2 * along) / l3
      q0 <- sum(4 * ex * ex / l3)
      q1 <- sum(4 * ex * ey / l3)
      q3 <- sum(4 * ey * ey / l3)
      s0 <- sum(pull * ex)
      s1 <- sum(pull * ey)
      denominator <- q0 * q3 - q1 * q1
      if (denominator == 0) next
      r0 <- (q3 * s0 - q1 * s1) / denominator
      r1 <- (-q1 * s0 + q0 * s1) / denominator

      change <- max(abs(slope[i, 1] + r0), abs(slope[i, 2] + r1))
      slope[i, 1] <- -r0
      slope[i, 2] <- -r1
      biggest_change <- max(biggest_change, change / max(1, abs(r0), abs(r1)))
    }
    if (biggest_change < tolerance) return(slope)
  }

  warning("plot_topography(): the fit of the surface did not settle, so the ",
          "map may be slightly inaccurate.", call. = FALSE)
  slope
}

# The cubic surface over every triangle. The triangle is cut into three small
# triangles at its middle, and each piece is a cubic written in Bezier form,
# which is a list of 19 numbers per triangle (c3000, c2100, ...). They follow
# from the values and slopes at the corners and from the direction in which the
# surface bends across each edge.
#
# Returns a list of 19 vectors, each with one number per triangle.
.topo_pieces <- function(points, mesh, neighbour, f, slope) {
  v1 <- mesh[, 1]
  v2 <- mesh[, 2]
  v3 <- mesh[, 3]
  x1 <- points[v1, 1]; y1 <- points[v1, 2]
  x2 <- points[v2, 1]; y2 <- points[v2, 2]
  x3 <- points[v3, 1]; y3 <- points[v3, 2]
  e12x <- x2 - x1; e12y <- y2 - y1
  e23x <- x3 - x2; e23y <- y3 - y2
  e31x <- x1 - x3; e31y <- y1 - y3
  f1 <- f[v1]; f2 <- f[v2]; f3 <- f[v3]

  df12 <-   slope[v1, 1] * e12x + slope[v1, 2] * e12y
  df21 <- -(slope[v2, 1] * e12x + slope[v2, 2] * e12y)
  df23 <-   slope[v2, 1] * e23x + slope[v2, 2] * e23y
  df32 <- -(slope[v3, 1] * e23x + slope[v3, 2] * e23y)
  df31 <-   slope[v3, 1] * e31x + slope[v3, 2] * e31y
  df13 <- -(slope[v1, 1] * e31x + slope[v1, 2] * e31y)

  c3000 <- f1
  c2100 <- (df12 + 3 * c3000) / 3
  c2010 <- (df13 + 3 * c3000) / 3
  c0300 <- f2
  c1200 <- (df21 + 3 * c0300) / 3
  c0210 <- (df23 + 3 * c0300) / 3
  c0030 <- f3
  c1020 <- (df31 + 3 * c0030) / 3
  c0120 <- (df32 + 3 * c0030) / 3
  c2001 <- (c2100 + c2010 + c3000) / 3
  c0201 <- (c1200 + c0300 + c0210) / 3
  c0021 <- (c1020 + c0120 + c0030) / 3

  # The direction in which the surface bends across the edge opposite each
  # corner; -0.5 where there is no triangle across the edge.
  terms <- .topo_edge_terms(points, mesh, neighbour)
  g <- terms$num / terms$den
  g[is.na(g)] <- -0.5

  c0111 <- (g[, 1] * (-c0300 + 3 * c0210 - 3 * c0120 + c0030) +
              (-c0300 + 2 * c0210 - c0120 + c0021 + c0201)) / 2
  c1011 <- (g[, 2] * (-c0030 + 3 * c1020 - 3 * c2010 + c3000) +
              (-c0030 + 2 * c1020 - c2010 + c2001 + c0021)) / 2
  c1101 <- (g[, 3] * (-c3000 + 3 * c2100 - 3 * c1200 + c0300) +
              (-c3000 + 2 * c2100 - c1200 + c2001 + c0201)) / 2
  c1002 <- (c1101 + c1011 + c2001) / 3
  c0102 <- (c1101 + c0111 + c0201) / 3
  c0012 <- (c1011 + c0111 + c0021) / 3
  c0003 <- (c1002 + c0102 + c0012) / 3

  list(c3000 = c3000, c2100 = c2100, c2010 = c2010, c2001 = c2001,
       c1200 = c1200, c1101 = c1101, c1020 = c1020, c1011 = c1011,
       c1002 = c1002, c0300 = c0300, c0210 = c0210, c0201 = c0201,
       c0120 = c0120, c0111 = c0111, c0102 = c0102, c0030 = c0030,
       c0021 = c0021, c0012 = c0012, c0003 = c0003)
}

# Finds the triangle each query point lies in, and the point's position inside
# that triangle (its three barycentric weights). A point exactly on an edge
# belongs to either triangle; a point outside every triangle gets NA.
.topo_locate <- function(points, mesh, qx, qy) {
  tight <- 100 * .Machine$double.eps
  loose <- sqrt(.Machine$double.eps)
  pad   <- 1e-6 * max(diff(range(points[, 1])), diff(range(points[, 2])))

  triangle <- rep(NA_integer_, length(qx))
  weights  <- matrix(NA_real_, length(qx), 3)

  x1 <- points[mesh[, 1], 1]; y1 <- points[mesh[, 1], 2]
  x2 <- points[mesh[, 2], 1]; y2 <- points[mesh[, 2], 2]
  x3 <- points[mesh[, 3], 1]; y3 <- points[mesh[, 3], 2]
  area_2 <- (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)

  # First with a very strict test for "inside", then a looser one for the few
  # points that sit on an edge and were missed by rounding.
  for (slack in c(tight, loose)) {
    todo <- which(is.na(triangle))
    for (t in seq_len(nrow(mesh))) {
      if (length(todo) == 0) break
      near <- todo[qx[todo] >= min(x1[t], x2[t], x3[t]) - pad &
                     qx[todo] <= max(x1[t], x2[t], x3[t]) + pad &
                     qy[todo] >= min(y1[t], y2[t], y3[t]) - pad &
                     qy[todo] <= max(y1[t], y2[t], y3[t]) + pad]
      if (length(near) == 0) next
      l1 <- ((y2[t] - y3[t]) * (qx[near] - x3[t]) +
               (x3[t] - x2[t]) * (qy[near] - y3[t])) / area_2[t]
      l2 <- ((y3[t] - y1[t]) * (qx[near] - x3[t]) +
               (x1[t] - x3[t]) * (qy[near] - y3[t])) / area_2[t]
      l3 <- 1 - l1 - l2
      hit <- pmin(l1, l2, l3) >= -slack
      if (any(hit)) {
        found <- near[hit]
        triangle[found] <- t
        weights[found, ] <- cbind(l1[hit], l2[hit], l3[hit])
        todo <- setdiff(todo, found)
      }
    }
  }
  list(triangle = triangle, weights = weights)
}

# Reads the cubic pieces at the located points (NA for points outside the mesh).
.topo_evaluate <- function(pieces, located) {
  out  <- rep(NA_real_, length(located$triangle))
  has  <- !is.na(located$triangle)
  t    <- located$triangle[has]
  w    <- located$weights[has, , drop = FALSE]
  low  <- pmin(w[, 1], w[, 2], w[, 3])
  b1 <- w[, 1] - low
  b2 <- w[, 2] - low
  b3 <- w[, 3] - low
  b4 <- 3 * low
  cc <- function(name) pieces[[name]][t]

  out[has] <-
    b1^3 * cc("c3000") + 3 * b1^2 * b2 * cc("c2100") +
    3 * b1^2 * b3 * cc("c2010") + 3 * b1^2 * b4 * cc("c2001") +
    3 * b1 * b2^2 * cc("c1200") + 6 * b1 * b2 * b4 * cc("c1101") +
    3 * b1 * b3^2 * cc("c1020") + 6 * b1 * b3 * b4 * cc("c1011") +
    3 * b1 * b4^2 * cc("c1002") + b2^3 * cc("c0300") +
    3 * b2^2 * b3 * cc("c0210") + 3 * b2^2 * b4 * cc("c0201") +
    3 * b2 * b3^2 * cc("c0120") + 6 * b2 * b3 * b4 * cc("c0111") +
    3 * b2 * b4^2 * cc("c0102") + b3^3 * cc("c0030") +
    3 * b3^2 * b4 * cc("c0021") + 3 * b3 * b4^2 * cc("c0012") +
    b4^3 * cc("c0003")
  out
}
