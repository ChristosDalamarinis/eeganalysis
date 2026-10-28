#' ============================================================================
#'                     EEG Montage: Channel Scalp Positions
#' ============================================================================
#'
#' A montage stores the spatial layout of a set of channels on the scalp
#' (Cartesian x/y/z and spherical theta/phi/radius coordinates), built from
#' the built-in BioSemi 64-channel electrode database (see
#' get_electrode_database(), R/channel_info2.R) and attached to an eeg object
#' via set_montage() - mirroring MNE-Python's raw.set_montage() workflow.
#'
#' Author: Christos Dalamarinis
#' Date: July 2026
#' ============================================================================

#' Create a New Montage Object (low-level constructor)
#'
#' Low-level constructor for the 'montage' S3 class. Most users should call
#' \code{\link{create_montage}} instead, which builds a montage from the
#' package's built-in BioSemi electrode template. Use \code{new_montage()}
#' directly only when supplying custom channel positions of your own.
#'
#' @param channels Character vector of channel names.
#' @param positions Data frame with one row per channel (in the same order as
#'   \code{channels}) and columns \code{channel, x, y, z, theta, phi, radius}.
#' @param coord_frame Character string identifying the coordinate template or
#'   source this montage was built from (e.g. \code{"biosemi64"}).
#'
#' @return An object of class \code{'montage'}: a list with \code{channels},
#'   \code{positions}, and \code{coord_frame}.
#'
#' @export
new_montage <- function(channels, positions, coord_frame = "biosemi64") {
  structure(
    list(
      channels = as.character(channels),
      positions = positions,
      coord_frame = as.character(coord_frame)
    ),
    class = "montage"
  )
}

#' Create a Montage From a Built-in BioSemi Cap Layout
#'
#' Builds a \code{montage} object giving scalp positions (Cartesian
#' x/y/z and spherical theta/phi/radius) for a set of channel names, taken
#' from one of the package's built-in BioSemi cap layouts. Pick the layout
#' that matches the cap the data was recorded with.
#'
#' @param channels Character vector of channel names to include. If
#'   \code{NULL} (default), every electrode of the layout is used.
#' @param template Character string naming the cap layout: one of
#'   \code{"biosemi16"}, \code{"biosemi32"}, \code{"biosemi64"} (the
#'   default), \code{"biosemi128"}, \code{"biosemi160"} or
#'   \code{"biosemi256"}, the BioSemi caps with that many electrodes.
#'
#' @return An object of class \code{'montage'}. See \code{new_montage}.
#'
#' @details
#' Channel names are matched case-insensitively. The 16-, 32- and
#' 64-channel layouts use the standard 10-20/10-10 names (e.g.
#' \code{"Cz"}) and also accept BioSemi's socket names, which number the
#' electrodes in the order BioSemi lists them (\code{"A1"} to \code{"A16"} on
#' the 16-channel cap, \code{"A1"} to \code{"A32"} on the 32-channel cap, and
#' \code{"A1"} to \code{"B32"} on the 64-channel cap, where \code{"B16"} is
#' Cz). The 128-, 160- and 256-channel layouts have no 10-20 names, only
#' BioSemi's own: \code{"A1"} to \code{"D32"}, \code{"A1"} to \code{"E32"} and
#' \code{"A1"} to \code{"H32"} (32 sockets per letter).
#'
#' Positions are the electrode coordinates BioSemi publishes for each cap
#' (the workbook "Cap_coords_all.xls" on biosemi.com), with the head treated
#' as a sphere of radius 87.54 mm in every layout.
#'
#' Requested channels that have no scalp position in the layout - either
#' because the layout has no electrode of that name, or because they are
#' external channels (EOG/ECG/EMG/GSR/etc., which have no coordinates) - are
#' dropped with a warning naming exactly which channels were excluded. Use
#' \code{detect_electrode_naming_system} first if you are unsure which
#' naming convention a dataset uses.
#'
#' @examples
#' \dontrun{
#'   # All 64 standard electrodes
#'   m <- create_montage()
#'   print(m)
#'
#'   # Just a subset, mixing standard and BioSemi names
#'   m <- create_montage(c("Cz", "Fz", "Pz", "B1"))
#'
#'   # Every electrode of a 128-channel cap (named A1 to D32 in the data)
#'   m <- create_montage(template = "biosemi128")
#'
#'   # A 32-channel recording whose channels are named by BioSemi socket
#'   m <- create_montage(c("A1", "A2", "A3"), template = "biosemi32")
#' }
#'
#' @seealso \code{\link{set_montage}}, \code{\link{get_electrode_database}}
#'
#' @export
create_montage <- function(channels = NULL, template = "biosemi64") {

  available <- .montage_templates()
  if (!is.character(template) || length(template) != 1 ||
      !(template %in% available)) {
    stop("ERROR: Unknown montage template '", paste(template, collapse = ", "),
         "'. Available templates: ", paste(available, collapse = ", "), ".",
         call. = FALSE)
  }

  layout <- .montage_template_table(template)

  if (is.null(channels)) {
    channels <- layout$name
  }

  channels <- as.character(channels)
  lookup_keys <- tolower(channels)

  # A channel is found by its name or, failing that, by its alias (the BioSemi
  # socket name on the 16, 32 and 64 channel caps). A missing name never matches.
  row_of <- match(lookup_keys, tolower(layout$name))
  by_alias <- is.na(row_of) & !is.na(lookup_keys)
  row_of[by_alias] <- match(lookup_keys[by_alias], tolower(layout$alias))

  ok <- !is.na(row_of)
  if (any(!ok)) {
    warning("create_montage(): ", sum(!ok),
            " channel(s) have no scalp position in the '", template,
            "' template and will be dropped: ",
            paste(channels[!ok], collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }
  if (!any(ok)) {
    stop("ERROR: None of the requested channels have a scalp position in ",
         "the '", template, "' template.", call. = FALSE)
  }

  channels_ok <- channels[ok]
  entries <- layout[row_of[ok], , drop = FALSE]

  positions <- data.frame(
    channel = channels_ok,
    x      = entries$x,
    y      = entries$y,
    z      = entries$z,
    theta  = entries$theta,
    phi    = entries$phi,
    radius = entries$radius,
    stringsAsFactors = FALSE
  )

  # The lower-case lookup keys are the row names (when they are all different),
  # as they always were: matrices built from this table inherit them as dimnames.
  if (!anyDuplicated(lookup_keys[ok])) {
    rownames(positions) <- lookup_keys[ok]
  }

  new_montage(channels = channels_ok, positions = positions, coord_frame = template)
}

#' Attach a Montage to an EEG Object
#'
#' Attaches channel scalp positions to an \code{eeg} object, mirroring
#' MNE-Python's \code{raw.set_montage()}. Only channels classified as
#' \code{"eeg"} in \code{eeg_obj$channel_types} (see the internal
#' \code{classify_channels()}, R/eeg_class.R) receive positions - external
#' (EOG/ECG/EMG/GSR/etc.) and status channels are never spatial and are left
#' out of the attached montage.
#'
#' @param eeg_obj An object of class \code{'eeg'}.
#' @param montage An object of class \code{'montage'}, typically built with
#'   \code{\link{create_montage}}.
#'
#' @return The input \code{eeg} object with \code{$montage} set to a montage
#'   restricted to the channels actually present in \code{eeg_obj}, and an
#'   entry appended to \code{$preprocessing_history}.
#'
#' @details
#' Channels of type \code{"eeg"} in \code{eeg_obj} that have no matching
#' entry in \code{montage}, and channels in \code{montage} that are not
#' present in \code{eeg_obj}, both produce a warning (not an error) - such
#' channels simply have no scalp position to plot.
#'
#' @examples
#' \dontrun{
#'   eeg <- read_bdf_native("recording.bdf")
#'   eeg <- set_montage(eeg, create_montage())
#' }
#'
#' @seealso \code{\link{create_montage}}, \code{\link{plot_topography}}
#'
#' @export
set_montage <- function(eeg_obj, montage) {

  if (!inherits(eeg_obj, "eeg")) {
    stop("ERROR: 'eeg_obj' must be an object of class 'eeg'.", call. = FALSE)
  }
  if (!inherits(montage, "montage")) {
    stop("ERROR: 'montage' must be an object of class 'montage' ",
         "(see create_montage()).", call. = FALSE)
  }

  eeg_channels <- eeg_obj$channels[eeg_obj$channel_types == "eeg"]

  missing_in_montage <- setdiff(eeg_channels, montage$channels)
  if (length(missing_in_montage) > 0) {
    warning("set_montage(): ", length(missing_in_montage),
            " EEG channel(s) in the data have no position in this montage: ",
            paste(missing_in_montage, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  missing_in_data <- setdiff(montage$channels, eeg_obj$channels)
  if (length(missing_in_data) > 0) {
    warning("set_montage(): ", length(missing_in_data),
            " montage channel(s) are not present in the data and will be ",
            "ignored: ", paste(missing_in_data, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }

  kept <- intersect(montage$channels, eeg_channels)
  if (length(kept) == 0) {
    stop("ERROR: No overlap between montage channels and EEG channels in ",
         "this object.", call. = FALSE)
  }

  eeg_obj$montage <- new_montage(
    channels = kept,
    positions = montage$positions[montage$positions$channel %in% kept, , drop = FALSE],
    coord_frame = montage$coord_frame
  )

  eeg_obj$preprocessing_history <- c(
    eeg_obj$preprocessing_history,
    list(paste0("Montage attached: ", montage$coord_frame, " (",
                length(kept), " channel(s))"))
  )

  eeg_obj
}

#' Print Method for Montage Objects
#'
#' Displays a formatted summary of a montage: coordinate template, channel
#' count, and channel names. Mirrors the style of \code{\link{print.eeg}}.
#'
#' @param x An object of class \code{'montage'}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.montage <- function(x, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Montage Object Summary\n")
  cat(strrep("=", 70), "\n\n")

  cat("  Coordinate template: ", x$coord_frame, "\n")
  cat("  Channels:            ", length(x$channels), "\n")

  channel_display <- paste(x$channels[seq_len(min(5, length(x$channels)))],
                            collapse = ", ")
  if (length(x$channels) > 5) {
    channel_display <- paste0(channel_display, " ... (+",
                               length(x$channels) - 5, " more)")
  }
  cat("    List:              ", channel_display, "\n")

  cat("\n", strrep("=", 70), "\n\n", sep = "")

  invisible(x)
}











#' Read a Montage From a Positions File
#'
#' Reads electrode positions from a file and returns a \code{montage} object,
#' the same kind \code{\link{create_montage}} returns, so it can be attached
#' with \code{\link{set_montage}} and used by every function that needs
#' electrode positions. Use it for a cap that has no built-in layout, or for
#' positions measured on a participant with a 3D digitizer.
#'
#' @param file Path to the positions file. The type of file is taken from its
#'   extension: \code{.elc}, \code{.sfp}, \code{.csv}, \code{.tsv} or
#'   \code{.txt} (see Details).
#' @param unit Unit of the positions in the file: \code{"mm"}, \code{"cm"} or
#'   \code{"m"}. If \code{NULL} (default), the unit the file states is used
#'   (an \code{.elc} file states it). When the file does not state one, the
#'   unit is worked out from the size of the head and a message says which
#'   one was taken; if the size fits no unit, the function stops and asks for
#'   \code{unit}.
#' @param recenter Logical. If \code{TRUE} (default), the positions are
#'   shifted so that the centre of the sphere that best fits the electrodes
#'   sits at the origin. The functions that project or interpolate positions
#'   (for example \code{\link{plot_topography}} and
#'   \code{\link{interpolate_bads}}) assume the head is centred on the origin,
#'   and the origin of a digitized head is often somewhere else, for example
#'   between the ears. If \code{FALSE}, the positions are kept exactly as the
#'   file gives them, apart from the unit.
#'
#' @return An object of class \code{'montage'}, built like the one from
#'   \code{\link{create_montage}}: positions in millimetres, with the columns
#'   \code{channel, x, y, z, theta, phi, radius}. The coordinate frame is
#'   \code{"file: "} followed by the name of the file.
#'
#' @details
#' Supported files (one electrode per entry; lines that are blank or start with
#' \code{#} are ignored):
#' \itemize{
#'   \item \code{.elc}: an ASA electrode file. A header states the unit
#'     (\code{UnitPosition mm}), a \code{Positions} block holds one
#'     \code{x y z} line per electrode (or \code{E01 : x y z}), and a
#'     \code{Labels} block holds the names in the same order.
#'   \item \code{.sfp}: one electrode per line, written as name, x, y, z,
#'     separated by spaces or tabs.
#'   \item \code{.csv}, \code{.tsv}, \code{.txt}: one electrode per row, written
#'     as name, x, y, z, separated by commas (\code{.csv}), tabs
#'     (\code{.tsv}) or spaces (\code{.txt}). A first row with column names is
#'     allowed: columns called name, label, channel, electrode or site, and
#'     x, y and z, are found by name, otherwise the first four columns are used
#'     in that order. Quotes around names are removed, so names cannot contain
#'     the separator.
#' }
#'
#' Rows for the head landmarks (Nz, Nasion, LPA, RPA, FidNz, FidT9, FidT10) are
#' not scalp electrodes and are left out. T9 and T10 are kept, because they
#' are electrodes in the 10-10 and 10-05 systems.
#'
#' Positions are stored in millimetres, like the built-in layouts. The
#' \code{theta}, \code{phi} and \code{radius} columns are worked out from
#' \code{x}, \code{y} and \code{z} the way the built-in layouts give them:
#' \code{theta} is the angle from the top of the head (0 at the vertex, negative
#' on the left side, positive on the right), \code{phi} is the angle around the
#' ear line, and \code{radius} is the distance from the origin.
#'
#' \code{\link{set_montage}} matches names exactly, upper and lower case
#' included, so the names in the file must be spelled like the channels of the
#' recording (for example \code{"Fp1"}, not \code{"FP1"}).
#'
#' Recentring needs at least 5 electrodes spread over a head-like shape. When
#' that is not the case, the positions are kept as the file gives them and a
#' warning says so.
#'
#' @examples
#' \dontrun{
#'   # An .elc file from a digitizer or a cap maker
#'   m <- read_custom_montage("my_cap.elc")
#'   eeg <- set_montage(eeg, m)
#'
#'   # A table of name, x, y, z in centimetres, kept exactly as written
#'   m <- read_custom_montage("positions.csv", unit = "cm", recenter = FALSE)
#' }
#'
#' @seealso \code{\link{create_montage}}, \code{\link{set_montage}},
#'   \code{\link{new_montage}}
#'
#' @export
read_custom_montage <- function(file, unit = NULL, recenter = TRUE) {

  if (!is.character(file) || length(file) != 1 || is.na(file) ||
      !file.exists(file)) {
    stop("ERROR: 'file' must be the path of an existing file.", call. = FALSE)
  }
  if (!is.null(unit) && (!is.character(unit) || length(unit) != 1 ||
                         !(unit %in% c("mm", "cm", "m")))) {
    stop("ERROR: 'unit' must be \"mm\", \"cm\" or \"m\" (or NULL to take it ",
         "from the file).", call. = FALSE)
  }
  if (!is.logical(recenter) || length(recenter) != 1 || is.na(recenter)) {
    stop("ERROR: 'recenter' must be TRUE or FALSE.", call. = FALSE)
  }

  file_name <- basename(file)
  extension <- if (grepl(".", file_name, fixed = TRUE)) {
    tolower(sub("^.*\\.", "", file_name))
  } else {
    ""
  }
  supported <- c("elc", "sfp", "csv", "tsv", "txt")
  if (!(extension %in% supported)) {
    stop("ERROR: Unsupported file type '.", extension, "'. Supported types: ",
         paste0(".", supported, collapse = ", "), ".", call. = FALSE)
  }

  lines <- readLines(file, warn = FALSE)
  # Some programs put a byte-order mark at the very start of a text file
  lines <- sub("^\xef\xbb\xbf", "", lines, useBytes = TRUE)

  parsed <- switch(extension,
    elc = .montage_read_elc(lines),
    csv = .montage_read_table(lines, ","),
    tsv = .montage_read_table(lines, "\t"),
    .montage_read_table(lines, "[[:space:]]+")
  )
  rows <- parsed$table
  rows$name <- trimws(rows$name)

  if (any(!nzchar(rows$name))) {
    stop("ERROR: Some electrodes in the file have no name.", call. = FALSE)
  }
  not_numbers <- !(is.finite(rows$x) & is.finite(rows$y) & is.finite(rows$z))
  if (any(not_numbers)) {
    first_bad <- rows$name[not_numbers]
    stop("ERROR: The position of ",
         paste(first_bad[seq_len(min(5, length(first_bad)))], collapse = ", "),
         " is not made of three numbers.", call. = FALSE)
  }

  landmarks <- c("nz", "nasion", "fidnz", "fidt9", "fidt10", "lpa", "rpa")
  electrodes <- rows[!(tolower(rows$name) %in% landmarks), , drop = FALSE]
  if (nrow(electrodes) == 0) {
    stop("ERROR: The file holds no electrode positions (only head ",
         "landmarks).", call. = FALSE)
  }
  repeated <- unique(electrodes$name[duplicated(electrodes$name)])
  if (length(repeated) > 0) {
    stop("ERROR: These electrode names appear more than once in the file: ",
         paste(repeated, collapse = ", "), ".", call. = FALSE)
  }

  fit <- .montage_sphere_fit(electrodes$x, electrodes$y, electrodes$z)

  unit_used <- if (!is.null(unit)) unit else parsed$unit
  if (is.null(unit_used)) {
    size <- if (!is.null(fit)) {
      fit$radius
    } else {
      from_origin <- sort(sqrt(electrodes$x^2 + electrodes$y^2 + electrodes$z^2))
      from_origin[ceiling(length(from_origin) / 2)]
    }
    unit_used <- if (size >= 40 && size <= 250) {
      "mm"
    } else if (size >= 4 && size <= 25) {
      "cm"
    } else if (size >= 0.04 && size <= 0.25) {
      "m"
    } else {
      stop("ERROR: Cannot tell the unit of the positions (the head comes out ",
           "about ", signif(size, 3), " in the file's own unit). Give 'unit' ",
           "as \"mm\", \"cm\" or \"m\".", call. = FALSE)
    }
    message("read_custom_montage(): the file does not state its unit; the ",
            "head radius of about ", signif(size, 3), " looks like ",
            unit_used, ", so the positions are taken as ", unit_used, ".")
  }
  to_mm <- c(mm = 1, cm = 10, m = 1000)[[unit_used]]

  x <- electrodes$x * to_mm
  y <- electrodes$y * to_mm
  z <- electrodes$z * to_mm

  if (isTRUE(recenter)) {
    if (is.null(fit)) {
      warning("read_custom_montage(): the centre of the head could not be ",
              "found from these positions (at least 5 electrodes spread over ",
              "a head-like shape are needed), so they are kept as the file ",
              "gives them.", call. = FALSE, immediate. = TRUE)
    } else {
      shift <- fit$centre * to_mm
      x <- x - shift[1]
      y <- y - shift[2]
      z <- z - shift[3]
    }
  }

  angles <- .montage_angles(x, y, z)

  positions <- data.frame(
    channel = electrodes$name,
    x      = x,
    y      = y,
    z      = z,
    theta  = angles$theta,
    phi    = angles$phi,
    radius = angles$radius,
    stringsAsFactors = FALSE
  )

  # Same row names as the table from create_montage(): the lower-case names,
  # when they are all different.
  keys <- tolower(electrodes$name)
  if (!anyDuplicated(keys)) {
    rownames(positions) <- keys
  }

  new_montage(channels = electrodes$name, positions = positions,
              coord_frame = paste0("file: ", file_name))
}

# Reads the lines of an .elc file: a header that may state the unit, a
# "Positions" block (one "x y z" line per electrode, or "E01 : x y z") and a
# "Labels" block with the names in the same order.
.montage_read_elc <- function(lines) {
  lines <- trimws(lines)

  unit <- NULL
  unit_row <- grep("^UnitPosition([[:space:]]|$)", lines)
  if (length(unit_row) > 0) {
    unit <- tolower(strsplit(lines[unit_row[1]], "[[:space:]]+")[[1]][2])
    if (is.na(unit) || !(unit %in% c("mm", "cm", "m"))) {
      stop("ERROR: The .elc file gives an unknown unit for the positions ('",
           unit, "').", call. = FALSE)
    }
  }

  positions_row <- which(lines == "Positions")
  labels_row <- which(lines == "Labels")
  if (length(positions_row) != 1 || length(labels_row) != 1 ||
      labels_row < positions_row) {
    stop("ERROR: Not a valid .elc file: it needs one 'Positions' line ",
         "followed by one 'Labels' line.", call. = FALSE)
  }

  position_lines <- lines[seq_len(labels_row - positions_row - 1) + positions_row]
  position_lines <- position_lines[nzchar(position_lines)]
  if (length(position_lines) == 0) {
    stop("ERROR: The .elc file has no positions.", call. = FALSE)
  }
  new_style <- any(grepl(":", position_lines, fixed = TRUE))
  if (new_style) {
    position_lines <- sub("^[^:]*:", "", position_lines)
  }
  coords <- lapply(strsplit(trimws(position_lines), "[[:space:]]+"),
                   function(f) suppressWarnings(as.numeric(f)))
  if (any(lengths(coords) != 3)) {
    stop("ERROR: Every line of the 'Positions' block of the .elc file must ",
         "hold three numbers.", call. = FALSE)
  }
  coords <- do.call(rbind, coords)

  label_lines <- if (labels_row < length(lines)) {
    lines[(labels_row + 1):length(lines)]
  } else {
    character(0)
  }
  blank <- which(!nzchar(label_lines))
  if (length(blank) > 0) {
    label_lines <- label_lines[seq_len(blank[1] - 1)]
  }
  labels <- if (new_style) {
    unlist(strsplit(label_lines, "[[:space:]]+"))
  } else {
    label_lines
  }
  if (length(labels) != nrow(coords)) {
    stop("ERROR: The .elc file has ", nrow(coords), " positions but ",
         length(labels), " labels.", call. = FALSE)
  }

  list(table = data.frame(name = labels, x = coords[, 1], y = coords[, 2],
                          z = coords[, 3], stringsAsFactors = FALSE),
       unit = unit)
}

# Reads the lines of a table with one electrode per row (name, x, y, z),
# split at `sep`. A first row of column names is allowed.
.montage_read_table <- function(lines, sep) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) {
    stop("ERROR: The file holds no electrode positions.", call. = FALSE)
  }

  fields <- lapply(strsplit(lines, sep), function(f) {
    gsub("^[\"']|[\"']$", "", trimws(f))
  })
  too_short <- which(lengths(fields) < 4)
  if (length(too_short) > 0) {
    stop("ERROR: Line ", too_short[1], " of the file (not counting blank and ",
         "# lines) has fewer than four values (name, x, y, z).", call. = FALSE)
  }

  as_number <- function(v) suppressWarnings(as.numeric(v))
  has_header <- anyNA(as_number(fields[[1]][2:4]))

  cols <- 1:4
  rows <- fields
  if (has_header) {
    header <- tolower(fields[[1]])
    named <- c(match(TRUE, header %in% c("name", "label", "channel",
                                         "electrode", "site")),
               match("x", header), match("y", header), match("z", header))
    if (!anyNA(named)) {
      cols <- named
    }
    rows <- fields[-1]
    if (length(rows) == 0) {
      stop("ERROR: The file holds no electrode positions (only a header row).",
           call. = FALSE)
    }
  }

  cells <- do.call(rbind, lapply(rows, function(f) f[cols]))
  list(table = data.frame(name = cells[, 1], x = as_number(cells[, 2]),
                          y = as_number(cells[, 3]), z = as_number(cells[, 4]),
                          stringsAsFactors = FALSE),
       unit = NULL)
}

# Centre and radius of the sphere that best fits the electrodes (least
# squares). NULL when there are fewer than 5 electrodes, when they do not
# determine a sphere, when they are not close to one, or when they are close
# to flat (a sphere much larger than the electrodes themselves).
.montage_sphere_fit <- function(x, y, z) {
  if (length(x) < 5) {
    return(NULL)
  }
  design <- cbind(2 * x, 2 * y, 2 * z, 1)
  target <- x^2 + y^2 + z^2
  solution <- tryCatch(solve(crossprod(design), crossprod(design, target)),
                       error = function(e) NULL)
  if (is.null(solution)) {
    return(NULL)
  }
  centre <- as.numeric(solution[1:3])
  radius_squared <- solution[4] + sum(centre^2)
  if (!is.finite(radius_squared) || radius_squared <= 0) {
    return(NULL)
  }
  radius <- sqrt(radius_squared)
  distance <- sqrt((x - centre[1])^2 + (y - centre[2])^2 + (z - centre[3])^2)
  if (max(abs(distance - radius)) > 0.5 * radius) {
    return(NULL)
  }
  span <- max(diff(range(x)), diff(range(y)), diff(range(z)))
  if (radius > 2 * span) {
    return(NULL)
  }
  list(centre = centre, radius = radius)
}

# theta, phi and radius of each position, the way the built-in layouts give
# them: theta is the angle from the vertex (negative on the left side,
# positive on the right), phi is the angle around the ear line (-90 to 90),
# so that x = r sin(theta) cos(phi), y = r sin(theta) sin(phi) and
# z = r cos(theta). An electrode on the midline counts as being on the right
# side, even when rounding has put it a hair to the left of it. Straight above
# (or below) the origin the azimuth is undefined and phi is 0, as in the
# built-in layouts.
.montage_angles <- function(x, y, z) {
  radius <- sqrt(x^2 + y^2 + z^2)
  from_vertex <- ifelse(radius > 0,
                        acos(pmin(1, pmax(-1, z / radius))) * 180 / pi, 0)
  tiny <- 1e-6 * max(radius)
  right <- x > -tiny
  phi <- ifelse(right, atan2(y, x), atan2(-y, -x)) * 180 / pi
  phi[sqrt(x^2 + y^2) < tiny] <- 0
  list(theta  = ifelse(right, from_vertex, -from_vertex),
       phi    = phi,
       radius = radius)
}
