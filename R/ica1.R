#' ============================================================================
#'                  Independent Component Analysis (ICA) Functions
#' ============================================================================
#'
#' This module implements ICA-based artifact removal for EEG data. ICA
#' decomposes the multi-channel EEG signal into statistically independent
#' components, allowing identification and rejection of components that
#' capture non-neural artifacts (eye blinks, muscle activity, cardiac
#' signals, electrode noise).
#'
#' The decomposition follows the FastICA/Infomax approach commonly used in
#' EEG pipelines, including:
#'   - Whitening via PCA prior to ICA rotation
#'   - Component labelling helpers (correlation with EOG/ECG channels)
#'   - Reconstruction of clean data after component exclusion
#'
#' Design mirrors mne.preprocessing.ICA: the object returned by new_ica()
#' stores the requested *parameters* immediately, while data-derived
#' *fitted attributes* (named with a trailing underscore, e.g.
#' n_components_, mixing_matrix_) stay NULL until fit_ica() runs. The
#' current_fit field ("unfitted" | "raw" | "epochs") is the canonical flag
#' for whether the object has been fit.
#'
#' Depends on: eeg_class.R (for eeg object structure)
#'
#' Author: Christos Dalamarinis
#' Date: July 2026
#' Status: pre-whitening fixed to pool std by channel type (matches MNE's
#'         _compute_pre_whitener()). PCA whitening still to be ported
#'         natively (plain SVD, no external dependency needed). The FastICA
#'         rotation step itself remains the one open question - port
#'         sklearn's _fastica.py line-by-line, or call it via reticulate.
#' Tested: NaN
#' ============================================================================

#' Create a New ICA Object
#'
#' Creates an unfitted \code{eeg_ica} object that stores the configuration
#' for an Independent Component Analysis decomposition, following the design
#' of \code{mne.preprocessing.ICA}. Parameters are available immediately;
#' data-derived attributes (named with a trailing underscore, e.g.
#' \code{n_components_}, \code{mixing_matrix_}) stay \code{NULL} until the
#' object is fit to data (via a future \code{fit_ica()}).
#'
#' @param n_components \code{NULL}, an integer, or a float in (0, 1).
#'   Number of principal components (from the pre-whitening PCA step) passed
#'   to the ICA algorithm during fitting.
#'   \describe{
#'     \item{integer}{Must be greater than 1 and less than or equal to the
#'       number of channels.}
#'     \item{float in (0, 1)}{Selects the smallest number of components whose
#'       cumulative explained variance exceeds this threshold.}
#'     \item{NULL}{Default. \code{0.999999} will be used at fit time, to
#'       avoid numerical stability problems with rank-deficient data.}
#'   }
#'   The actual number resolved at fit time is stored in \code{n_components_}.
#' @param noise_cov \code{NULL} or a numeric channels x channels covariance
#'   matrix. Used for pre-whitening. If \code{NULL} (default), channels are
#'   scaled to unit variance ("z-standardized") prior to PCA whitening. This
#'   package does not yet implement a dedicated covariance-estimation
#'   workflow, so only \code{NULL} is functionally supported by the fitting
#'   step for now; a matrix may be supplied and is stored as-is for forward
#'   compatibility.
#' @param random_state \code{NULL} or a single integer. Seed for the random
#'   number generator used by the ICA algorithm. If \code{NULL} (default),
#'   results will generally differ between runs.
#' @param method Character string. The ICA method to use. Currently only
#'   \code{"fastica"} is implemented, backed by this package's own native
#'   port of scikit-learn's FastICA rotation (\code{.ica_par()} /
#'   \code{.ica_def()}) - not the \code{fastICA} package.
#' @param fit_params \code{NULL} or a named list of additional arguments
#'   controlling the FastICA rotation run by \code{fit_ica()}. Entries
#'   supplied here override these defaults:
#'   \describe{
#'     \item{algorithm}{\code{"parallel"} (default) or \code{"deflation"} -
#'       selects \code{.ica_par()} or \code{.ica_def()} at fit time.}
#'     \item{fun}{\code{"logcosh"} (default), \code{"exp"}, or \code{"cube"} -
#'       the neg-entropy G-function (\code{.g_logcosh()}, \code{.g_exp()},
#'       \code{.g_cube()}).}
#'     \item{fun_args}{Named list forwarded to \code{fun}, default
#'       \code{list(alpha = 1)}. Only \code{alpha} is read, and only by
#'       \code{"logcosh"}; must lie in \code{[1, 2]}.}
#'     \item{tol}{Numeric convergence tolerance, default \code{1e-04}.}
#'   }
#' @param max_iter Integer, or \code{"auto"} (default). If \code{"auto"}, it
#'   resolves immediately to \code{1000}. The actual number of iterations
#'   taken to fit will be stored in \code{n_iter_}.
#' @param verbose \code{NULL}, logical, or character. Reserved for future use
#'   controlling the verbosity of fitting/logging output.
#'
#' @return An object of class \code{eeg_ica}, a list containing:
#'  \describe{
#'    \item{n_components, noise_cov, random_state, method, fit_params,
#'      max_iter, verbose}{The (validated/resolved) parameters as supplied.}
#'    \item{current_fit}{\code{"unfitted"} until \code{fit_ica()} is run.}
#'    \item{ch_names}{\code{NULL} until fit; channel names used for fitting.}
#'    \item{n_components_, pre_whitener_, pca_components_, pca_mean_,
#'      pca_explained_variance_, mixing_matrix_, unmixing_matrix_,
#'      n_samples_, n_iter_}{\code{NULL} until fit; data-derived attributes.}
#'    \item{exclude}{Integer vector (initially empty) of component indices
#'      to exclude when reconstructing data. Present from creation and
#'      user/auto-detection editable — not a fitted attribute.}
#'    \item{labels_}{List (initially empty) of independent component
#'      indices grouped by type (e.g. \code{"eog"}, \code{"ecg"}). Populated
#'      by future automatic artifact-detection helpers.}
#'  }
#'
#' @details
#' Following scikit-learn/MNE convention, attribute names ending in a
#' trailing underscore (\code{n_components_}, \code{mixing_matrix_}, ...)
#' signify that the attribute is derived from data and only exists once the
#' object has been fit. Before fitting, these fields are present but
#' \code{NULL}. \code{current_fit} is the authoritative flag for whether an
#' \code{eeg_ica} object has been fit — check
#' \code{ica$current_fit != "unfitted"} rather than relying on individual
#' fields being non-\code{NULL}.
#'
#' @examples
#' \dontrun{
#'   # Default configuration (n_components resolved at fit time)
#'   ica <- new_ica()
#'
#'   # Request 15 components, custom max_iter
#'   ica <- new_ica(n_components = 15, max_iter = 500)
#'
#'   # Request components explaining 95% of variance
#'   ica <- new_ica(n_components = 0.95)
#' }
#'
#' @export
new_ica <- function(n_components = NULL,
                     noise_cov = NULL,
                     random_state = NULL,
                     method = "fastica",
                     fit_params = NULL,
                     max_iter = "auto",
                     verbose = NULL) {

  # ========== VALIDATE method ==========

  if (!is.character(method) || length(method) != 1) {
    stop("ERROR: 'method' must be a single character string.")
  }
  if (method != "fastica") {
    stop("ERROR: method = '", method, "' is not yet implemented. ",
         "Only 'fastica' is currently supported ('infomax' and 'picard' ",
         "are planned for a future release).")
  }

  # ========== VALIDATE n_components ==========

  if (!is.null(n_components)) {
    if (!is.numeric(n_components) || length(n_components) != 1) {
      stop("ERROR: 'n_components' must be a single number, or NULL.")
    }
    if (n_components != as.integer(n_components)) {
      # float case: must be in (0, 1) exclusive
      if (!(n_components > 0 && n_components < 1)) {
        stop("ERROR: Selecting ICA components by explained variance needs ",
             "'n_components' between 0.0 and 1.0 (exclusive), got ",
             n_components, ".")
      }
    } else if (n_components <= 1) {
      stop("ERROR: 'n_components' must be greater than 1 when selecting an ",
           "integer number of components, got ", n_components, ".")
    }
  }

  # ========== RESOLVE max_iter ==========

  if (identical(max_iter, "auto")) {
    max_iter <- 1000L
  } else {
    if (!is.numeric(max_iter) || length(max_iter) != 1 ||
        max_iter != as.integer(max_iter) || max_iter <= 0) {
      stop("ERROR: 'max_iter' must be a positive integer, or 'auto'.")
    }
    max_iter <- as.integer(max_iter)
  }

  # ========== RESOLVE fit_params ==========

  if (is.null(fit_params)) {
    fit_params <- list()
  } else if (!is.list(fit_params)) {
    stop("ERROR: 'fit_params' must be a named list, or NULL.")
  }

  default_fit_params <- list(
    algorithm = "parallel",
    fun       = "logcosh",
    fun_args  = list(alpha = 1),
    tol       = 1e-04
  )
  fit_params <- modifyList(default_fit_params, fit_params)

  if (!is.character(fit_params$algorithm) || length(fit_params$algorithm) != 1 ||
      !(fit_params$algorithm %in% c("parallel", "deflation"))) {
    stop("ERROR: fit_params$algorithm must be 'parallel' or 'deflation', got '",
         fit_params$algorithm, "'.")
  }
  if (!is.character(fit_params$fun) || length(fit_params$fun) != 1 ||
      !(fit_params$fun %in% c("logcosh", "exp", "cube"))) {
    stop("ERROR: fit_params$fun must be 'logcosh', 'exp', or 'cube', got '",
         fit_params$fun, "'.")
  }
  if (!is.list(fit_params$fun_args)) {
    stop("ERROR: fit_params$fun_args must be a named list.")
  }
  if (!is.numeric(fit_params$tol) || length(fit_params$tol) != 1 || fit_params$tol < 0) {
    stop("ERROR: fit_params$tol must be a single non-negative number.")
  }

  # ========== CONSTRUCT ICA OBJECT ==========

  ica_object <- structure(
    list(
      # ---- parameters ----
      n_components = n_components,
      noise_cov = noise_cov,
      random_state = random_state,
      method = method,
      fit_params = fit_params,
      max_iter = max_iter,
      verbose = verbose,

      # ---- fit state ----
      current_fit = "unfitted",
      ch_names = NULL,

      # ---- fitted attributes (trailing underscore, NULL until fit) ----
      n_components_ = NULL,
      pre_whitener_ = NULL,
      pca_components_ = NULL,
      pca_mean_ = NULL,
      pca_explained_variance_ = NULL,
      mixing_matrix_ = NULL,
      unmixing_matrix_ = NULL,
      n_samples_ = NULL,
      n_iter_ = NULL,

      # ---- present from creation, not fitted attributes ----
      exclude = integer(0),
      labels_ = list()
    ),
    class = "eeg_ica"
  )

  return(ica_object)
}

#' Print Method for ICA Objects
#'
#' Custom print method that displays a formatted summary of an \code{eeg_ica}
#' object: the configured parameters and, if the object has been fit, the
#' resolved fitted attributes.
#'
#' @param x An object of class \code{eeg_ica}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x} (standard R print method convention).
#'
#' @examples
#' \dontrun{
#'   ica <- new_ica(n_components = 15)
#'   print(ica)  # Calls this method automatically
#' }
#'
#' @export
print.eeg_ica <- function(x, ...) {

  # ========== HEADER ==========
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("ICA Object Summary\n")
  cat(strrep("=", 70), "\n\n")

  # ========== PARAMETERS ==========
  cat("PARAMETERS:\n")
  cat("  n_components:    ",
      if (is.null(x$n_components)) "None (resolved at fit time)" else x$n_components,
      "\n")
  cat("  method:          ", x$method, "\n")
  cat("  max_iter:        ", x$max_iter, "\n")
  cat("  noise_cov:       ",
      if (is.null(x$noise_cov)) "None (z-score standardization)" else "supplied",
      "\n")
  cat("  random_state:    ", if (is.null(x$random_state)) "None" else x$random_state, "\n")

  # ========== FIT STATUS ==========
  cat("\nFIT STATUS:\n")
  if (identical(x$current_fit, "unfitted")) {
    cat("  Not fitted. Call fit_ica() to run decomposition.\n")
  } else {
    cat("  Fitted on:       ", x$current_fit, "data\n")
    cat("  n_components_:   ", x$n_components_, "\n")
    cat("  n_samples_:      ", x$n_samples_, "\n")
    cat("  n_iter_:         ", x$n_iter_, "\n")
    cat("  Excluded ICs:    ", length(x$exclude), "\n")
    if (length(x$labels_) > 0) {
      cat("  Labels:          ",
          paste(names(x$labels_), "=", lengths(x$labels_), collapse = ", "), "\n")
    }
  }

  # ========== FOOTER ==========
  cat("\n", strrep("=", 70), "\n\n", sep = "")

  # Return invisibly (standard R convention)
  invisible(x)
}

#
# ============================================================================
#                    PRE-WHITENING HELPERS (private)
# ============================================================================
#
# The step run before PCA/ICA: every channel of a given *type* is scaled by
# one shared standard deviation (pooled across all channels of that type),
# so that channel types with naturally larger amplitude (e.g. EOG vs EEG)
# don't dominate the decomposition purely because of scale. Mirrors
# ICA._compute_pre_whitener() / ICA._pre_whiten() in python/ica/ica.py
# (python/ica/ica.py:841-857), but only the default (noise_cov = NULL)
# branch, which is self-contained in ica.py. The noise_cov branch there
# calls compute_whitener() from cov.py (a full noise-covariance
# eigendecomposition) - that path is not implemented here yet.
#
# ----------------------------------------------------------------------------
# .population_sd() - population standard deviation (numpy's ddof = 0)
# ----------------------------------------------------------------------------
#' Population standard deviation, matching numpy's default ddof
#'
#' R's \code{sd()} divides by \code{n - 1} (Bessel's correction). numpy's
#' \code{np.std()} - what \code{ICA._compute_pre_whitener()} actually calls
#' - divides by \code{n} (\code{ddof = 0}). Use this instead of \code{sd()}
#' anywhere exact numerical parity with the MNE reference matters.
#'
#' @param x Numeric vector or matrix (flattened before computing, matching
#'   numpy's default \code{axis = None} behavior).
#' @return A single numeric value.
#' @keywords internal
.population_sd <- function(x) {
  x <- as.double(x)
  m <- mean(x)
  sqrt(mean((x - m)^2))
}

# ----------------------------------------------------------------------------
# .compute_pre_whitener() - per-channel-*type* standard deviation
# ----------------------------------------------------------------------------
#' Compute the pre-whitening vector for ICA fitting
#'
#' Computes the pre-whitening scale used to z-standardize data before PCA
#' whitening. This is the default (\code{noise_cov = NULL}) pre-whitening
#' path from \code{mne.preprocessing.ICA}: MNE pools a single standard
#' deviation across *all* channels of a given type - e.g. every EEG channel
#' shares one scale factor - rather than scaling each channel independently
#' (see \code{ICA._compute_pre_whitener()}, \code{python/ica/ica.py:841-857}).
#'
#' @param data Numeric matrix, channels x time points (same layout as
#'   \code{eeg$data}). Must already be restricted to the channels intended
#'   for ICA fitting - e.g. the BioSemi status channel excluded - since that
#'   selection is a \code{fit_ica()}-level concern upstream of this helper,
#'   mirroring how MNE applies \code{picks} before \code{_fit()} ever sees
#'   the data.
#' @param channels Character vector of channel names, \code{length(channels)
#'   == nrow(data)}. Used only to name the returned vector.
#' @param channel_types Character vector, same length and order as
#'   \code{channels}, with values \code{"eeg"} / \code{"external"} /
#'   \code{"status"} - i.e. \code{eeg_obj$channel_types} as computed once by
#'   \code{classify_channels()} in \code{new_eeg()} (see R/eeg_class.R).
#'   Channels typed \code{"eeg"} form one pooled-std group and channels typed
#'   \code{"external"} form another; \code{"status"} entries are not
#'   assigned to either group and must already be excluded from \code{data}.
#' @param noise_cov \code{NULL} (default). A non-NULL value is not yet
#'   supported and raises an error.
#' @return Named numeric vector of length \code{nrow(data)}: one population
#'   standard deviation per channel, pooled within and broadcast across
#'   each type group, named by \code{channels}.
#' @keywords internal
.compute_pre_whitener <- function(data, channels, channel_types, noise_cov = NULL) {
  if (!is.null(noise_cov)) {
    stop("ERROR: noise_cov-based pre-whitening is not yet implemented. ",
         "Only the default z-score path (noise_cov = NULL) is supported.")
  }
  if (length(channels) != nrow(data)) {
    stop("ERROR: length(channels) (", length(channels),
         ") must match nrow(data) (", nrow(data), ").")
  }
  if (length(channel_types) != nrow(data)) {
    stop("ERROR: length(channel_types) (", length(channel_types),
         ") must match nrow(data) (", nrow(data), ").")
  }

  eeg_idx <- which(channel_types == "eeg")
  exg_idx <- which(channel_types == "external")

  pre_whitener <- numeric(length(channels))
  if (length(eeg_idx) > 0) {
    pre_whitener[eeg_idx] <- .population_sd(data[eeg_idx, , drop = FALSE])
  }
  if (length(exg_idx) > 0) {
    pre_whitener[exg_idx] <- .population_sd(data[exg_idx, , drop = FALSE])
  }

  names(pre_whitener) <- channels
  pre_whitener
}

# ----------------------------------------------------------------------------
# .pre_whiten() - apply the pre-whitening vector to data
# ----------------------------------------------------------------------------
#' Apply pre-whitening to data
#'
#' Divides each channel (row) of \code{data} by its corresponding entry in
#' \code{pre_whitener}. Used both when fitting ICA and later whenever new
#' data is projected through an already-fitted \code{eeg_ica} object (e.g.
#' \code{get_sources()}, \code{apply_ica()}), which is why it is kept as a
#' standalone helper rather than being inlined into the fitting step.
#'
#' @param data Numeric matrix, channels x time points.
#' @param pre_whitener Numeric vector of length \code{nrow(data)}, as
#'   returned by \code{.compute_pre_whitener()}.
#' @return Numeric matrix of the same shape as \code{data}, pre-whitened.
#' @keywords internal
.pre_whiten <- function(data, pre_whitener) {
  sweep(data, MARGIN = 1, STATS = pre_whitener, FUN = "/")
}

#
# ============================================================================
#                    PCA WHITENING HELPER (private)
# ============================================================================
#
# The full-rank PCA step run on the pre-whitened data before the ICA
# rotation itself. Mirrors the
#   pca = _PCA(n_components=self._max_pca_components, whiten=True)
#   data = pca.fit_transform(data.T)
# lines in ICA._fit() (python/ica/ica.py:900-901) - specifically sklearn's
# "covariance_eigh" solver path (eigendecompose the channel x channel
# covariance matrix), which is what PCA(svd_solver = "auto") resolves to
# whenever n_channels <= 1000 and n_samples >= 10 * n_channels - true for
# essentially any real EEG recording (tens of thousands of samples against
# a few dozen to a few hundred channels). The "full"/"randomized"/"arpack"
# solvers sklearn also offers, for small-sample or huge-channel-count edge
# cases, are not implemented here.
#
# ----------------------------------------------------------------------------
# .pca_whiten() - mean-center, eigendecompose, whiten to unit variance
# ----------------------------------------------------------------------------
#' PCA-whiten pre-whitened data before ICA rotation
#'
#' Mean-centers \code{data}, eigendecomposes its channel x channel covariance
#' matrix, and rescales each principal-component score to unit variance -
#' i.e. \code{whiten = TRUE} in \code{sklearn.decomposition.PCA}, as used by
#' \code{ICA._fit()} at \code{python/ica/ica.py:900}. Always returns the
#' *full* set of components; selecting how many of them (\code{n_components_})
#' go on to the ICA rotation step is a separate, later concern, mirroring how
#' MNE fits PCA at \code{max_pca_components} and only slices
#' \code{data[:, sel]} afterward. Because of that, trailing components whose
#' eigenvalue is (numerically) zero - e.g. the one rank-deficient dimension
#' introduced by average referencing - will legitimately contain \code{Inf}/
#' \code{NaN} in \code{data} after the unit-variance division; this matches
#' what MNE's own full PCA output looks like before its component-selection
#' step trims those dimensions away, and is why \code{new_ica()} defaults to
#' resolving \code{n_components = 0.999999} rather than requesting the full
#' rank.
#'
#' @param data Numeric matrix, channels x time points, already pre-whitened
#'   (see \code{.compute_pre_whitener()} / \code{.pre_whiten()}).
#' @return A list with:
#'  \describe{
#'    \item{data}{Numeric matrix, components x time points: the
#'      unit-variance-whitened principal component scores.}
#'    \item{pca_mean_}{Numeric vector, length \code{nrow(data)}: the
#'      per-channel mean removed before decomposition.}
#'    \item{pca_components_}{Numeric matrix, components x channels: the
#'      principal axes (covariance-matrix eigenvectors), sign-fixed so the
#'      largest-magnitude entry of each component is positive (matches
#'      sklearn's \code{svd_flip(..., u_based_decision = FALSE)}, making the
#'      result deterministic regardless of the underlying LAPACK
#'      eigensolver's arbitrary sign convention).}
#'    \item{pca_explained_variance_}{Numeric vector, one eigenvalue per
#'      component, decreasing order; tiny negative values (numerical noise)
#'      are clipped to 0.}
#'  }
#' @keywords internal
.pca_whiten <- function(data) {
  n_samples <- ncol(data)

  pca_mean_ <- rowMeans(data)
  data_centered <- sweep(data, MARGIN = 1, STATS = pca_mean_, FUN = "-")

  # channel x channel covariance, ddof = 1 (matches sklearn's
  # _cov(X, ddof=1, ...) in the covariance_eigh solver path)
  cov_matrix <- tcrossprod(data_centered) / (n_samples - 1)

  eig <- eigen(cov_matrix, symmetric = TRUE)
  eigenvalues  <- eig$values
  eigenvectors <- eig$vectors

  # Numerical noise can leave the smallest eigenvalue(s) slightly negative;
  # sklearn clips these the same way before storing explained_variance_.
  eigenvalues[eigenvalues < 0] <- 0

  # Deterministic sign convention (mirrors sklearn's svd_flip): the
  # largest-magnitude entry of each component becomes positive, so results
  # don't depend on the arbitrary sign LAPACK's eigensolver happens to
  # return.
  for (i in seq_len(ncol(eigenvectors))) {
    component <- eigenvectors[, i]
    peak <- which.max(abs(component))
    if (component[peak] < 0) {
      eigenvectors[, i] <- -component
    }
  }

  pca_components_ <- t(eigenvectors)

  whitened <- sweep(
    pca_components_ %*% data_centered,
    MARGIN = 1,
    STATS = sqrt(eigenvalues),
    FUN = "/"
  )

  list(
    data = whitened,
    pca_mean_ = pca_mean_,
    pca_components_ = pca_components_,
    pca_explained_variance_ = eigenvalues
  )
}


















































################################################################################
################################################################################
###                                                                          ###
###   >>>>>>>>>>>>  NEW CODE STARTS HERE — FASTICA ROTATION STEP  <<<<<<<<<< ###
###        (ported line-by-line from scikit-learn's _fastica.py)            ###
###                                                                          ###
################################################################################
################################################################################
#
# Everything below is a native R port of the FastICA rotation itself -
# the piece that turns the PCA-whitened data from .pca_whiten() into
# independent components. Ported from sklearn's python/sklearn/decomposition/
# _fastica.py (BSD-3-Clause), specifically lines 36-168: _gs_decorrelation(),
# _sym_decorrelation(), _ica_def(), _ica_par(), _logcosh(), _exp(), _cube().
# The rest of that file (the FastICA estimator class, its own whitening,
# parameter validation) is not ported - this package's whitening is already
# handled upstream by .compute_pre_whitener() / .pre_whiten() / .pca_whiten().
#
# All functions here operate on data laid out components x time points,
# matching the shape .pca_whiten() already returns.
#
# ============================================================================
#                    NEG-ENTROPY G-FUNCTIONS (private)
# ============================================================================
#
# Each g-function approximates the derivative of a neg-entropy contrast used
# to drive FastICA's fixed-point update. All three mirror sklearn exactly,
# including a numpy-broadcasting quirk that matters for parity: _logcosh()
# behaves differently depending on whether its input is a 1D vector (used by
# the deflation algorithm, one component's projection at a time) or a 2D
# matrix (used by the parallel algorithm, all components at once) - in the
# vector case it returns the elementwise derivative un-reduced, deferring the
# mean to the caller; in the matrix case it reduces to one value per row
# itself. _exp() and _cube() do not have that branch - they always reduce
# along the last axis regardless of input shape. `.axis_mean()` below
# reproduces that always-reduce behavior for both.
#
# ----------------------------------------------------------------------------
# .axis_mean() - numpy's `.mean(axis=-1)`, for either a vector or a matrix
# ----------------------------------------------------------------------------
#' Mean along the last axis, matching numpy's \code{.mean(axis=-1)}
#'
#' For a plain numeric vector, equivalent to \code{mean(x)} (numpy would
#' collapse a 1D array to a scalar). For a components x samples matrix,
#' equivalent to \code{rowMeans(x)} (numpy's last axis is the samples/columns
#' axis, so \code{axis=-1} reduces each row to one value). Used by
#' \code{.g_exp()} and \code{.g_cube()}, which - unlike \code{.g_logcosh()} -
#' always reduce this way regardless of input shape
#' (\code{_fastica.py:160-168}).
#'
#' @param x Numeric vector or matrix.
#' @return A single number (vector input) or a numeric vector of length
#'   \code{nrow(x)} (matrix input).
#' @keywords internal
.axis_mean <- function(x) {
  if (is.null(dim(x))) mean(x) else rowMeans(x)
}

# ----------------------------------------------------------------------------
# .g_logcosh() - logcosh neg-entropy approximation (sklearn's default `fun`)
# ----------------------------------------------------------------------------
#' logcosh G-function and its derivative
#'
#' Mirrors sklearn's \code{_logcosh()} (\code{_fastica.py:143-157}), the
#' default \code{fun} for both \code{.ica_par()} and \code{.ica_def()}. When
#' \code{x} is a plain vector (the deflation algorithm's one-component-at-a-
#' time projection), the derivative \code{g_x} is returned elementwise,
#' un-reduced - \code{.ica_def()} takes \code{mean(g_x)} itself afterward.
#' When \code{x} is a components x samples matrix (the parallel algorithm's
#' all-components-at-once projection), \code{g_x} is already reduced to one
#' value per row here, matching sklearn's per-row loop in the 2D branch.
#'
#' @param x Numeric vector (length n_samples) or matrix (components x
#'   n_samples): the current unmixing projection \code{w^T X} / \code{WX}.
#' @param fun_args Named list; only \code{alpha} is read (default \code{1.0}),
#'   matching sklearn's \code{fun_args.get("alpha", 1.0)}. Must lie in
#'   \code{[1, 2]}, per \code{FastICA}'s own parameter validation - not
#'   re-validated here.
#' @return A list with \code{gx} (same shape as \code{x}: elementwise
#'   \code{tanh(alpha * x)}) and \code{g_x} (elementwise derivative if
#'   \code{x} is a vector; one value per row if \code{x} is a matrix).
#' @keywords internal
.g_logcosh <- function(x, fun_args = list()) {
  alpha <- if (is.null(fun_args$alpha)) 1.0 else fun_args$alpha

  gx <- tanh(alpha * x)

  if (is.null(dim(gx))) {
    g_x <- alpha * (1 - gx^2)
  } else {
    g_x <- alpha * (1 - rowMeans(gx^2))
  }

  list(gx = gx, g_x = g_x)
}

# ----------------------------------------------------------------------------
# .g_exp() - exponential neg-entropy approximation
# ----------------------------------------------------------------------------
#' exp G-function and its derivative
#'
#' Mirrors sklearn's \code{_exp()} (\code{_fastica.py:160-164}). Unlike
#' \code{.g_logcosh()}, \code{g_x} is always reduced via \code{.axis_mean()}
#' regardless of whether \code{x} is a vector or a matrix - there is no
#' un-reduced branch, so \code{.ica_def()}'s own \code{mean(g_x)} is a no-op
#' when this function is used (matches numpy: \code{.mean()} of an
#' already-0-d/1-per-row result returns it unchanged).
#'
#' @inheritParams .g_logcosh
#' @return A list with \code{gx} (same shape as \code{x}) and \code{g_x}
#'   (already reduced: a single number for vector input, one value per row
#'   for matrix input).
#' @keywords internal
.g_exp <- function(x, fun_args = list()) {
  exp_term <- exp(-(x^2) / 2)
  gx <- x * exp_term
  g_x <- .axis_mean((1 - x^2) * exp_term)

  list(gx = gx, g_x = g_x)
}

# ----------------------------------------------------------------------------
# .g_cube() - cubic neg-entropy approximation
# ----------------------------------------------------------------------------
#' cube G-function and its derivative
#'
#' Mirrors sklearn's \code{_cube()} (\code{_fastica.py:167-168}). Same
#' always-reduced behavior as \code{.g_exp()} - see its docs.
#'
#' @inheritParams .g_logcosh
#' @return A list with \code{gx} (same shape as \code{x}) and \code{g_x}
#'   (already reduced: a single number for vector input, one value per row
#'   for matrix input).
#' @keywords internal
.g_cube <- function(x, fun_args = list()) {
  gx <- x^3
  g_x <- .axis_mean(3 * x^2)

  list(gx = gx, g_x = g_x)
}

#
# ============================================================================
#                    DECORRELATION HELPERS (private)
# ============================================================================
#
# ----------------------------------------------------------------------------
# .sym_decorrelation() - symmetric orthogonalization of the unmixing matrix
# ----------------------------------------------------------------------------
#' Symmetric decorrelation: \code{W <- (W W')^{-1/2} W}
#'
#' Mirrors sklearn's \code{_sym_decorrelation()} (\code{_fastica.py:61-73}).
#' Used both to initialize \code{.ica_par()}'s starting \code{W} and inside
#' its main loop after every update, so that all rows of \code{W} (all
#' estimated components) stay mutually orthonormal throughout - this is what
#' makes the "parallel" algorithm parallel, as opposed to \code{.ica_def()}'s
#' one-component-at-a-time Gram-Schmidt approach
#' (\code{.gs_decorrelation()}).
#'
#' \code{eigen()} returns eigenvalues in decreasing order versus numpy's
#' \code{eigh()} (increasing order), but that difference doesn't matter here:
#' the reconstruction \code{u \%*\% diag(f(s)) \%*\% t(u)} is invariant to the
#' order eigenvector/eigenvalue pairs are summed in, as long as each
#' eigenvector stays paired with its own eigenvalue - which \code{eigen()}
#' guarantees.
#'
#' @param W Numeric square matrix (n_components x n_components) to
#'   orthogonalize.
#' @return The orthogonalized matrix, same shape as \code{W}.
#' @keywords internal
.sym_decorrelation <- function(W) {
  eig <- eigen(W %*% t(W), symmetric = TRUE)
  s <- eig$values
  u <- eig$vectors

  # Avoid sqrt of a negative/zero value from rounding error - same clip
  # sklearn applies (np.finfo(W.dtype).tiny), which also prevents a
  # divide-by-zero next.
  s <- pmax(s, .Machine$double.xmin)

  u %*% diag(1 / sqrt(s), nrow = length(s)) %*% t(u) %*% W
}

# ----------------------------------------------------------------------------
# .gs_decorrelation() - Gram-Schmidt orthogonalization against prior rows
# ----------------------------------------------------------------------------
#' Orthonormalize a candidate row against already-estimated components
#'
#' Mirrors sklearn's \code{_gs_decorrelation()} (\code{_fastica.py:36-58}).
#' Used only by \code{.ica_def()}: after estimating a new candidate direction
#' \code{w} for component \code{j}, this projects out any overlap with the
#' \code{j - 1} components already stored in the earlier rows of \code{W},
#' so deflation extracts each component only once.
#'
#' @param w Numeric vector, length n_components: the candidate row to
#'   orthogonalize.
#' @param W Numeric matrix, n_components x n_components: rows
#'   \code{1:n_prior} hold the already-estimated (and already mutually
#'   orthonormal) components; remaining rows are not read.
#' @param n_prior Integer, number of already-estimated rows of \code{W} to
#'   orthogonalize \code{w} against - i.e. sklearn's \code{j} (their loop is
#'   0-based, so this is the same value; R callers pass \code{j - 1} from
#'   their own 1-based loop index). \code{0} is a no-op (first component).
#' @return \code{w}, orthogonalized against \code{W[1:n_prior, ]}. Not yet
#'   re-normalized to unit length - that is done by the caller afterward,
#'   matching sklearn.
#' @keywords internal
.gs_decorrelation <- function(w, W, n_prior) {
  if (n_prior > 0) {
    W_prior <- W[seq_len(n_prior), , drop = FALSE]
    w <- w - as.vector((w %*% t(W_prior)) %*% W_prior)
  }
  w
}

#
# ============================================================================
#                    FASTICA ROTATION ALGORITHMS (private)
# ============================================================================
#
# ----------------------------------------------------------------------------
# .ica_par() - parallel (symmetric) FastICA
# ----------------------------------------------------------------------------
#' Parallel FastICA rotation
#'
#' Mirrors sklearn's \code{_ica_par()} (\code{_fastica.py:111-140}) - the
#' main loop behind \code{algorithm = "parallel"}, sklearn's default. All
#' components are updated simultaneously each iteration, then re-orthogonalized
#' as a block via \code{.sym_decorrelation()}, rather than one at a time.
#'
#' @param X Numeric matrix, components x samples: PCA-whitened data, i.e.
#'   the \code{data} element returned by \code{.pca_whiten()} (sliced down to
#'   the \code{n_components_} columns actually going into the rotation).
#' @param w_init Numeric square matrix, n_components x n_components: starting
#'   unmixing matrix (sklearn draws this from a standard normal when not
#'   user-supplied; generating it is left to the caller, not this function).
#' @param g A G-function - one of \code{.g_logcosh}, \code{.g_exp},
#'   \code{.g_cube} - called as \code{g(x, fun_args)} and returning
#'   \code{list(gx, g_x)}.
#' @param fun_args Named list forwarded to \code{g} (e.g. \code{list(alpha =
#'   1)} for \code{.g_logcosh}).
#' @param max_iter Integer, maximum number of iterations.
#' @param tol Numeric, convergence tolerance: the loop stops once the largest
#'   per-row change in \code{W} (see \code{lim} below) drops under this.
#' @return A list with \code{W} (the estimated unmixing matrix, same shape as
#'   \code{w_init}) and \code{n_iter} (number of iterations actually run). A
#'   warning is issued - matching sklearn's \code{ConvergenceWarning} - if
#'   \code{max_iter} is reached without \code{lim < tol}.
#' @keywords internal
.ica_par <- function(X, w_init, g, fun_args, max_iter, tol) {
  W <- .sym_decorrelation(w_init)
  p_ <- ncol(X)
  converged <- FALSE

  for (ii in seq_len(max_iter)) {
    g_out <- g(W %*% X, fun_args)
    gwtx  <- g_out$gx
    g_wtx <- g_out$g_x

    W1 <- .sym_decorrelation((gwtx %*% t(X)) / p_ - g_wtx * W)

    # Row-wise dot product of W1 and W, i.e. diag(W1 %*% t(W)) without
    # materializing the full product - matches sklearn's
    # np.einsum("ij,ij->i", W1, W).
    lim <- max(abs(abs(rowSums(W1 * W)) - 1))
    W <- W1

    if (lim < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning("FastICA did not converge. Consider increasing tolerance or ",
            "the maximum number of iterations.")
  }

  list(W = W, n_iter = ii)
}

# ----------------------------------------------------------------------------
# .ica_def() - deflationary FastICA
# ----------------------------------------------------------------------------
#' Deflationary FastICA rotation
#'
#' Mirrors sklearn's \code{_ica_def()} (\code{_fastica.py:76-108}) - the main
#' loop behind \code{algorithm = "deflation"}. Components are estimated one
#' at a time: each new candidate direction is fit to convergence, orthogonalized
#' against every component already extracted via \code{.gs_decorrelation()},
#' and only then is the next component started.
#'
#' @inheritParams .ica_par
#' @return A list with \code{W} (the estimated unmixing matrix, same shape as
#'   \code{w_init}) and \code{n_iter} (the maximum number of iterations used
#'   by any single component, matching sklearn's \code{max(n_iter)} - unlike
#'   \code{.ica_par()}, no convergence warning is issued here, matching
#'   sklearn: \code{_ica_def()} silently accepts whatever \code{w} it has
#'   after \code{max_iter} iterations for a given component).
#' @keywords internal
.ica_def <- function(X, w_init, g, fun_args, max_iter, tol) {
  n_components <- nrow(w_init)
  W <- matrix(0, n_components, n_components)
  n_iter <- integer(n_components)

  for (j in seq_len(n_components)) {
    w <- w_init[j, ]
    w <- w / sqrt(sum(w^2))

    for (i in seq_len(max_iter)) {
      wtx <- as.vector(w %*% X)
      g_out <- g(wtx, fun_args)
      gwtx  <- g_out$gx
      g_wtx <- g_out$g_x

      w1 <- rowMeans(sweep(X, MARGIN = 2, STATS = gwtx, FUN = "*")) -
        mean(g_wtx) * w
      w1 <- .gs_decorrelation(w1, W, j - 1L)
      w1 <- w1 / sqrt(sum(w1^2))

      lim <- abs(abs(sum(w1 * w)) - 1)
      w <- w1

      if (lim < tol) break
    }

    n_iter[j] <- i
    W[j, ] <- w
  }

  list(W = W, n_iter = max(n_iter))
}

#
# ============================================================================
#                    N_COMPONENTS RESOLUTION HELPER (private)
# ============================================================================
#
# ----------------------------------------------------------------------------
# .resolve_n_components() - turn new_ica()'s n_components into a concrete int
# ----------------------------------------------------------------------------
#' Resolve \code{n_components} to a concrete integer at fit time
#'
#' \code{new_ica()} validates the *shape* of \code{n_components} (NULL / an
#' integer != 1 / a float in (0, 1)) but can't resolve what it actually means
#' until fit time, since that depends on data it doesn't have yet - either
#' the total channel count (integer case) or the PCA explained-variance
#' profile (NULL / float case, both selecting components by cumulative
#' variance). This is that resolution step, called from \code{fit_ica()}
#' right after \code{.pca_whiten()}.
#'
#' @param n_components The (already shape-validated) \code{new_ica()}
#'   parameter: \code{NULL}, an integer, or a float in (0, 1).
#' @param pca_explained_variance_ Numeric vector, decreasing order: the
#'   \code{pca_explained_variance_} element returned by \code{.pca_whiten()}
#'   for the *full* (unsliced) PCA - one eigenvalue per available channel.
#' @return A single integer, \code{>= 1} and \code{<= length(pca_explained_variance_)}.
#' @keywords internal
.resolve_n_components <- function(n_components, pca_explained_variance_) {
  n_available <- length(pca_explained_variance_)

  if (is.null(n_components) || (n_components > 0 && n_components < 1)) {
    threshold <- if (is.null(n_components)) 0.999999 else n_components
    cum_ratio <- cumsum(pca_explained_variance_) / sum(pca_explained_variance_)
    n_resolved <- which(cum_ratio >= threshold)[1]
    if (is.na(n_resolved)) n_resolved <- n_available
    return(as.integer(n_resolved))
  }

  n_components <- as.integer(n_components)
  if (n_components <= 1 || n_components > n_available) {
    stop("ERROR: 'n_components' (", n_components, ") must be greater than 1 ",
         "and less than or equal to the number of channels used for fitting ",
         "(", n_available, ").")
  }
  n_components
}

#
# ============================================================================
#                    FIT_ICA() - PUBLIC ORCHESTRATION
# ============================================================================
#
# This is the function that actually connects everything above to real data:
# it takes the *settings* stored on an unfitted eeg_ica object (new_ica())
# and an eeg object's data, runs them through pre-whitening -> PCA whitening
# -> n_components resolution -> the FastICA rotation (.ica_par()/.ica_def()),
# and writes the *results* back onto the eeg_ica object's fitted (trailing
# underscore) fields. Mirrors mne.preprocessing.ICA.fit(raw, picks=...), with
# one deliberate difference: R has no in-place mutation, so - unlike MNE,
# which mutates self and returns it - this returns a new, now-fitted
# eeg_ica object; the caller must reassign it (ica <- fit_ica(ica, eeg)).
#
# ----------------------------------------------------------------------------
# fit_ica() - fit an eeg_ica object to eeg data
# ----------------------------------------------------------------------------
#' Fit an ICA Decomposition to EEG Data
#'
#' Runs the full ICA fitting pipeline - pre-whitening, PCA whitening,
#' \code{n_components} resolution, and the FastICA rotation - on the data of
#' an \code{eeg} object, using the configuration stored on an unfitted
#' \code{eeg_ica} object (see \code{new_ica()}). The fitted attributes
#' (\code{mixing_matrix_}, \code{unmixing_matrix_}, \code{n_components_},
#' etc.) are written back onto that object, and \code{current_fit} changes
#' from \code{"unfitted"} to \code{"raw"}.
#'
#' @param ica An unfitted \code{eeg_ica} object, as returned by
#'   \code{new_ica()}.
#' @param eeg An \code{eeg} object (see \code{new_eeg()}) supplying the data
#'   to fit on.
#' @param picks \code{NULL} (default), a character vector of channel names,
#'   or a numeric vector of channel indices - which channels of \code{eeg}
#'   to fit ICA on. If \code{NULL}, every channel with
#'   \code{channel_types != "status"} is used (i.e. \code{"eeg"} and
#'   \code{"external"} channels; the BioSemi status channel is always
#'   excluded, matching \code{.compute_pre_whitener()}'s expectations).
#'
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with its
#'   fitted (trailing-underscore) fields populated and \code{current_fit} set
#'   to \code{"raw"}. Since R does not mutate arguments in place, the caller
#'   must reassign the result (\code{ica <- fit_ica(ica, eeg)}) - the input
#'   \code{ica} object itself is left untouched.
#'
#' @examples
#' \dontrun{
#'   eeg <- new_eeg(data, channels, sampling_rate = 512)
#'   ica <- new_ica(n_components = 0.95, random_state = 42)
#'   ica <- fit_ica(ica, eeg)
#'   print(ica)
#' }
#'
#' @export
fit_ica <- function(ica, eeg, picks = NULL) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).")
  }

  # ========== SELECT CHANNELS ==========

  if (is.null(picks)) {
    idx <- which(eeg$channel_types != "status")
  } else if (is.character(picks)) {
    idx <- match(picks, eeg$channels)
    if (anyNA(idx)) {
      stop("ERROR: 'picks' contains channel name(s) not found in eeg$channels: ",
           paste(picks[is.na(idx)], collapse = ", "))
    }
  } else {
    idx <- as.integer(picks)
  }

  if (length(idx) < 2) {
    stop("ERROR: at least 2 channels are required to fit ICA, got ",
         length(idx), ".")
  }

  data <- eeg$data[idx, , drop = FALSE]
  channels <- eeg$channels[idx]
  channel_types <- eeg$channel_types[idx]

  # ========== PRE-WHITEN ==========

  pre_whitener_ <- .compute_pre_whitener(data, channels, channel_types,
                                          noise_cov = ica$noise_cov)
  data_pw <- .pre_whiten(data, pre_whitener_)

  # ========== PCA WHITEN ==========

  pca <- .pca_whiten(data_pw)

  # ========== RESOLVE n_components_ ==========

  n_components_ <- .resolve_n_components(ica$n_components,
                                          pca$pca_explained_variance_)
  X <- pca$data[seq_len(n_components_), , drop = FALSE]

  # ========== GENERATE w_init ==========
  #
  # Scoped to this call: if random_state is supplied, the global RNG stream
  # is seeded, used, then restored to what it was beforehand - so a
  # reproducible fit doesn't silently reseed the user's session-wide random
  # stream (mirrors sklearn's check_random_state(), which never touches
  # numpy's global generator either).

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  if (!is.null(ica$random_state)) set.seed(ica$random_state)
  w_init <- matrix(rnorm(n_components_^2), n_components_, n_components_)
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  # ========== FASTICA ROTATION ==========

  g <- switch(ica$fit_params$fun,
              logcosh = .g_logcosh,
              exp     = .g_exp,
              cube    = .g_cube)

  rotation <- if (ica$fit_params$algorithm == "parallel") {
    .ica_par(X, w_init, g, ica$fit_params$fun_args, ica$max_iter,
             ica$fit_params$tol)
  } else {
    .ica_def(X, w_init, g, ica$fit_params$fun_args, ica$max_iter,
             ica$fit_params$tol)
  }

  # ========== STORE FITTED ATTRIBUTES ==========

  ica$ch_names <- channels
  ica$n_components_ <- n_components_
  ica$pre_whitener_ <- pre_whitener_
  ica$pca_components_ <- pca$pca_components_[seq_len(n_components_), , drop = FALSE]
  ica$pca_mean_ <- pca$pca_mean_
  ica$pca_explained_variance_ <- pca$pca_explained_variance_
  ica$unmixing_matrix_ <- rotation$W
  ica$mixing_matrix_ <- solve(rotation$W)
  ica$n_samples_ <- ncol(data)
  ica$n_iter_ <- rotation$n_iter
  ica$current_fit <- "raw"

  ica
}

#
# ============================================================================
#                    GET_SOURCES() - PUBLIC: COMPONENT TIME-COURSES
# ============================================================================
#
# Applies an already-fitted eeg_ica object's transform (pre-whiten -> PCA
# project -> unmix) to eeg data, producing the actual independent component
# time-courses. Mirrors mne.preprocessing.ICA.get_sources(inst): reuses the
# *fitted* pre_whitener_ / pca_mean_ / pca_components_ /
# pca_explained_variance_ / unmixing_matrix_ rather than recomputing anything
# from scratch, so this is always a transform of a fixed decomposition -
# never a new fit. Because of that, eeg does not have to be the same object
# fit_ica() was called on, as long as it has channels matching
# ica$ch_names - mirrors MNE, which lets get_sources() run on a different
# Raw than the one ICA was fit on.
#
# ----------------------------------------------------------------------------
# get_sources() - compute independent component time-courses
# ----------------------------------------------------------------------------
#' Compute Independent Component Time-Courses
#'
#' Applies a fitted \code{eeg_ica} object's transform to \code{eeg} data,
#' returning the actual independent component waveforms - the numbers behind
#' every other inspection tool (\code{plot_ica_sources()}, the component
#' summary table). Reuses the fitted pre-whitening scale, PCA projection, and
#' unmixing matrix stored on \code{ica} (from \code{fit_ica()}); does not
#' refit anything.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param eeg An \code{eeg} object supplying the data to transform. Does not
#'   have to be the same object \code{ica} was fit on, but must contain every
#'   channel named in \code{ica$ch_names}.
#'
#' @return A numeric matrix, \code{n_components_ x n_samples}, row-named
#'   \code{"IC1"}, \code{"IC2"}, ... in the same order as
#'   \code{ica$unmixing_matrix_}'s rows.
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   sources <- get_sources(ica, eeg)
#'   dim(sources)  # n_components_ x n_samples
#' }
#'
#' @export
get_sources <- function(ica, eeg) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.")
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).")
  }

  idx <- match(ica$ch_names, eeg$channels)
  if (anyNA(idx)) {
    stop("ERROR: 'eeg' is missing channel(s) that 'ica' was fit on: ",
         paste(ica$ch_names[is.na(idx)], collapse = ", "))
  }

  # ========== REPLAY THE FITTED TRANSFORM ==========

  data <- eeg$data[idx, , drop = FALSE]
  data_pw <- .pre_whiten(data, ica$pre_whitener_)

  centered <- sweep(data_pw, MARGIN = 1, STATS = ica$pca_mean_, FUN = "-")
  whitened <- sweep(
    ica$pca_components_ %*% centered,
    MARGIN = 1,
    STATS = sqrt(ica$pca_explained_variance_[seq_len(ica$n_components_)]),
    FUN = "/"
  )

  sources <- ica$unmixing_matrix_ %*% whitened
  rownames(sources) <- paste0("IC", seq_len(ica$n_components_))

  sources
}

#
# ============================================================================
#                    PLOT_ICA_SOURCES() - PUBLIC: STACKED COMPONENT VIEW
# ============================================================================
#
# Visual inspection tool: plots several (or all) independent components
# stacked on top of each other over time, like a mini multi-channel EEG
# viewer, one labeled trace per component. Built on get_sources() for the
# data and follows the same plotting conventions as plot_eeg_signal()
# (R/plot_signal.R) - tmin/tmax in ms, demean, event markers - so it looks
# and behaves consistently with the rest of the package, just for components
# instead of electrodes.
#
# ----------------------------------------------------------------------------
# plot_ica_sources() - stacked time-course view of independent components
# ----------------------------------------------------------------------------
#' Plot Independent Component Time-Courses (Stacked View)
#'
#' Plots several (or all) independent components of a fitted \code{eeg_ica}
#' object stacked on top of each other over a chosen time window - a
#' multi-channel-style view for visually scanning components before deciding
#' which to exclude (e.g. a component with sharp, occasional spikes is often
#' an eye blink; fast, jittery activity is often muscle). Each trace is
#' scaled to unit standard deviation purely for display, so all components
#' stack with comparable visual amplitude regardless of their true
#' (dimensionless) scale - this does not affect the underlying data, only how
#' it is drawn.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param eeg An \code{eeg} object, passed through to \code{get_sources()}
#'   (does not have to be the object \code{ica} was fit on - see its docs).
#' @param components \code{NULL} (default: every component,
#'   \code{1:ica$n_components_}) or an integer vector of component indices to
#'   plot, e.g. \code{c(1, 3, 5)}.
#' @param tmin,tmax Numeric, start/end of the time window to display, in
#'   milliseconds. \code{NULL} (default) plots from the start/to the end of
#'   the recording.
#' @param demean Logical, default \code{TRUE}. Subtracts each component's own
#'   mean over the plotted window before drawing (components are already
#'   near-zero-mean from the ICA fit; this only affects the plotted window).
#' @param show_events Logical, default \code{TRUE}. Draws a vertical dashed
#'   line (and label) at each \code{eeg$events} onset within the window, same
#'   as \code{plot_eeg_signal()}.
#' @param event_color Character string, color for event marker lines/labels.
#'   Default \code{"red"}.
#' @param line_color Character string, color of every component trace.
#'   Default \code{"black"}.
#' @param line_width Numeric, line width of the component traces. Default
#'   \code{1}.
#' @param title Character string, custom plot title. If \code{NULL} (default),
#'   a title is generated automatically.
#'
#' @return Invisibly returns a list with:
#'  \describe{
#'    \item{components}{Integer vector of component indices actually plotted.}
#'    \item{tmin_ms, tmax_ms}{Actual start/end of the plotted window, in ms.}
#'    \item{n_samples}{Number of samples in the plotted window.}
#'    \item{n_events}{Number of event markers drawn.}
#'  }
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   plot_ica_sources(ica, eeg)                       # all components, full recording
#'   plot_ica_sources(ica, eeg, components = c(1, 2))  # just IC1 and IC2
#'   plot_ica_sources(ica, eeg, tmin = 1000, tmax = 5000)
#' }
#'
#' @seealso \code{\link{get_sources}}, \code{\link{plot_eeg_signal}},
#'   \code{\link{fit_ica}}
#'
#' @importFrom graphics plot lines abline mtext axis par
#' @importFrom grDevices adjustcolor
#' @importFrom stats sd
#' @export
plot_ica_sources <- function(ica,
                              eeg,
                              components  = NULL,
                              tmin        = NULL,
                              tmax        = NULL,
                              demean      = TRUE,
                              show_events = TRUE,
                              event_color = "red",
                              line_color  = "black",
                              line_width  = 1,
                              title       = NULL) {

  # ========== INPUT VALIDATION ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).",
         call. = FALSE)
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.",
         call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  # ---- Resolve components ----

  if (is.null(components)) {
    components <- seq_len(ica$n_components_)
  } else {
    if (!is.numeric(components) || length(components) == 0) {
      stop("ERROR: 'components' must be a non-empty integer vector, or NULL.",
           call. = FALSE)
    }
    components <- as.integer(components)
    if (any(components < 1) || any(components > ica$n_components_)) {
      stop("ERROR: 'components' must be integer indices in [1, ",
           ica$n_components_, "], got: ",
           paste(components, collapse = ", "), ".", call. = FALSE)
    }
  }
  n_plot <- length(components)

  if (!is.logical(demean) || length(demean) != 1) {
    stop("ERROR: 'demean' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_events) || length(show_events) != 1) {
    stop("ERROR: 'show_events' must be TRUE or FALSE.", call. = FALSE)
  }

  # ========== COMPUTE SOURCES ==========

  sources      <- get_sources(ica, eeg)[components, , drop = FALSE]
  n_timepoints <- ncol(sources)
  times_ms     <- eeg$times * 1000

  # ---- Resolve time window (ms -> sample indices), matching plot_eeg_signal() ----

  rec_start_ms <- times_ms[1]
  rec_end_ms   <- times_ms[n_timepoints]

  if (is.null(tmin)) tmin <- rec_start_ms
  if (is.null(tmax)) tmax <- rec_end_ms

  if (!is.numeric(tmin) || length(tmin) != 1) {
    stop("ERROR: 'tmin' must be a single numeric value in milliseconds.",
         call. = FALSE)
  }
  if (!is.numeric(tmax) || length(tmax) != 1) {
    stop("ERROR: 'tmax' must be a single numeric value in milliseconds.",
         call. = FALSE)
  }
  if (tmin >= tmax) {
    stop("ERROR: 'tmin' (", tmin, " ms) must be less than 'tmax' (", tmax,
         " ms).", call. = FALSE)
  }
  if (tmin < rec_start_ms || tmax > rec_end_ms) {
    warning("Requested time window [", tmin, ", ", tmax, "] ms extends beyond ",
            "the recording [", round(rec_start_ms, 1), ", ",
            round(rec_end_ms, 1), "] ms. Clipping to recording bounds.",
            call. = FALSE, immediate. = TRUE)
    tmin <- max(tmin, rec_start_ms)
    tmax <- min(tmax, rec_end_ms)
  }

  samp_start <- which.min(abs(times_ms - tmin))
  samp_end   <- which.min(abs(times_ms - tmax))
  if (samp_start >= samp_end) {
    stop("ERROR: Time window too narrow - fewer than 2 samples selected. ",
         "Increase the range between 'tmin' and 'tmax'.", call. = FALSE)
  }

  # ========== EXTRACT, DEMEAN, NORMALIZE, STACK ==========

  sources_slice <- sources[, samp_start:samp_end, drop = FALSE]
  time_slice    <- times_ms[samp_start:samp_end]
  n_samples     <- length(time_slice)

  if (demean) {
    sources_slice <- sweep(sources_slice, MARGIN = 1,
                            STATS = rowMeans(sources_slice), FUN = "-")
  }

  # Unit-SD normalization for display only (does not affect returned data,
  # only vertical spacing on the plot).
  trace_sd <- apply(sources_slice, 1, sd)
  trace_sd[trace_sd == 0] <- 1
  sources_norm <- sources_slice / trace_sd

  # Spacing sized to the actual (post-normalization) peak amplitude present,
  # rather than a fixed constant - ICA components are often strongly
  # non-Gaussian (e.g. blink spikes far beyond +-1 SD), so a fixed spacing
  # tuned for Gaussian-ish traces would let spiky components collide with
  # their neighbors.
  max_abs <- max(abs(sources_norm), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  spacing <- max_abs * 2.2

  offsets <- (n_plot - seq_len(n_plot)) * spacing   # first component plotted at top
  stacked <- sweep(sources_norm, MARGIN = 1, STATS = offsets, FUN = "+")

  # ========== IDENTIFY EVENTS IN WINDOW ==========

  events_in_window <- data.frame(onset_ms = numeric(0), type = character(0),
                                  stringsAsFactors = FALSE)
  if (show_events && !is.null(eeg$events) && nrow(eeg$events) > 0) {
    ev        <- eeg$events
    ev_ms     <- ev$onset_time * 1000
    in_window <- ev_ms >= tmin & ev_ms <= tmax
    events_in_window <- data.frame(onset_ms = ev_ms[in_window],
                                    type = ev$type[in_window],
                                    stringsAsFactors = FALSE)
  }
  n_events <- nrow(events_in_window)

  # ========== BUILD TITLE ==========

  if (is.null(title)) {
    title <- paste0("ICA Sources  |  ", n_plot, " component(s)  |  ",
                     round(tmin), " - ", round(tmax), " ms")
  }

  # ========== DRAW PLOT ==========

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if (show_events && n_events > 0) {
    par(mar = c(4.5, 5, 4, 1))
  } else {
    par(mar = c(4.5, 5, 3, 1))
  }

  y_lim <- c(min(offsets) - spacing / 2, max(offsets) + spacing / 2)

  plot(
    x    = range(time_slice),
    y    = y_lim,
    type = "n",
    xlab = "Time (ms)",
    ylab = "",
    main = title,
    yaxt = "n",
    las  = 1
  )
  axis(2, at = offsets, labels = paste0("IC", components), las = 1, cex.axis = 0.8)
  abline(h = offsets, col = "gray85", lty = 3, lwd = 0.6)

  if (show_events && n_events > 0) {
    for (i in seq_len(n_events)) {
      abline(v   = events_in_window$onset_ms[i],
             col = adjustcolor(event_color, alpha.f = 0.7),
             lty = 2,
             lwd = 1.2)
    }
    mtext(
      text = events_in_window$type,
      side = 3,
      at   = events_in_window$onset_ms,
      col  = event_color,
      cex  = 0.65,
      line = 0.3
    )
  }

  for (i in seq_len(n_plot)) {
    lines(x = time_slice, y = stacked[i, ], col = line_color, lwd = line_width)
  }

  # ========== INVISIBLE RETURN ==========

  invisible(list(
    components = components,
    tmin_ms    = round(time_slice[1], 2),
    tmax_ms    = round(time_slice[n_samples], 2),
    n_samples  = n_samples,
    n_events   = n_events
  ))
}

#
# ============================================================================
#          GET_COMPONENT_TOPOGRAPHY() / PLOT_ICA_TOPOGRAPHY() - PUBLIC
# ============================================================================
#
# Spatial inspection tool: for one component, computes "how much does each
# electrode contribute to this component" - one number per channel - so it
# can be handed straight to the existing plot_topography() (R/topography.R)
# for a scalp heatmap, exactly like MNE's ica.plot_components().
#
# Derivation: fit_ica()/get_sources() computes, forward,
#   sources = unmixing_matrix_ %*% D %*% pca_components_ %*% (data / pre_whitener_ - pca_mean_)
# where D = diag(1 / sqrt(pca_explained_variance_[1:n_components_])). Dropping
# the constant offset (pca_mean_ - topographies describe the *linear*
# channel/component relationship, not a DC baseline), the topography is the
# pseudo-inverse of that combined linear operator, taken column-by-column
# (one component at a time) and inverted step by step in reverse order:
#   - pinv(unmixing_matrix_) = mixing_matrix_ (already stored: it's square
#     and invertible, so its ordinary inverse *is* its pseudo-inverse)
#   - pinv(D) = diag(sqrt(ev)) (D is diagonal and invertible)
#   - pinv(pca_components_) = t(pca_components_) (its rows are orthonormal -
#     pca_components_ %*% t(pca_components_) == I - which is exactly the
#     condition under which a matrix's transpose *is* its Moore-Penrose
#     pseudo-inverse)
#   - finally divide back through by pre_whitener_ to undo pre-whitening and
#     land back in original channel (e.g. microvolt) units, mirroring
#     sklearn's own mixing_ = pinv(components_) (components_ = W @ K in
#     _fastica.py) - same idea, just derived here step by step instead of
#     via a single generic pinv() call, since every step is already
#     invertible in closed form.
#
# ----------------------------------------------------------------------------
# get_component_topography() - one component's weight per channel
# ----------------------------------------------------------------------------
#' Compute an Independent Component's Scalp Topography
#'
#' For one fitted independent component, computes how strongly each channel
#' contributes to it - one number per channel, in the same units as the
#' original \code{eeg} data - by inverting the pre-whitening, PCA, and ICA
#' rotation steps \code{fit_ica()} applied going forward. The result is a
#' named numeric vector ready to pass directly as \code{values} to
#' \code{\link{plot_topography}} (see \code{\link{plot_ica_topography}} for a
#' one-call wrapper that does exactly that).
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param component Single integer, the component index (\code{1} to
#'   \code{ica$n_components_}) to compute the topography for.
#'
#' @return A named numeric vector, length \code{length(ica$ch_names)}, named
#'   by \code{ica$ch_names}. Channels with no scalp position (external/EOG/
#'   ECG/status channels) are included here - \code{plot_topography()} drops
#'   them itself, since it already only plots channels present in the
#'   attached montage.
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   topo <- get_component_topography(ica, component = 1)
#'   eeg  <- set_montage(eeg, create_montage())
#'   plot_topography(eeg, topo, main = "IC1 Topography")
#' }
#'
#' @seealso \code{\link{plot_ica_topography}}, \code{\link{plot_topography}},
#'   \code{\link{fit_ica}}
#'
#' @export
get_component_topography <- function(ica, component) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.")
  }
  if (!is.numeric(component) || length(component) != 1) {
    stop("ERROR: 'component' must be a single integer.")
  }
  component <- as.integer(component)
  if (component < 1 || component > ica$n_components_) {
    stop("ERROR: 'component' must be an integer in [1, ", ica$n_components_,
         "], got ", component, ".")
  }

  # ========== BACK-PROJECT ONTO CHANNELS (see derivation above) ==========

  ev_sel  <- ica$pca_explained_variance_[seq_len(ica$n_components_)]
  loading <- ica$mixing_matrix_[, component]

  topo_pw <- as.vector(t(ica$pca_components_) %*% (sqrt(ev_sel) * loading))
  topo    <- topo_pw * ica$pre_whitener_

  names(topo) <- ica$ch_names
  topo
}

# ----------------------------------------------------------------------------
# plot_ica_topography() - one-call scalp map for a single component
# ----------------------------------------------------------------------------
#' Plot an Independent Component's Scalp Topography
#'
#' One-call convenience wrapper: computes a component's topography with
#' \code{\link{get_component_topography}} and draws it with
#' \code{\link{plot_topography}}. Requires a montage attached to \code{eeg}
#' (or passed explicitly) - see \code{\link{set_montage}},
#' \code{\link{create_montage}}.
#'
#' @param ica A fitted \code{eeg_ica} object (see \code{fit_ica()}).
#' @param eeg An \code{eeg} object. Used to read \code{eeg$montage} when
#'   \code{montage} is not supplied directly, and passed through to
#'   \code{plot_topography()}.
#' @param component Single integer, the component index to plot.
#' @param montage An object of class \code{'montage'}. If \code{NULL}
#'   (default), \code{eeg$montage} is used.
#' @param ... Additional arguments forwarded to \code{\link{plot_topography}}
#'   (e.g. \code{interpolate_res}, \code{contour}, \code{col_palette}).
#'
#' @return Invisibly returns whatever \code{\link{plot_topography}} returns.
#'
#' @examples
#' \dontrun{
#'   eeg <- set_montage(eeg, create_montage())
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   plot_ica_topography(ica, eeg, component = 1)
#' }
#'
#' @seealso \code{\link{get_component_topography}}, \code{\link{plot_topography}}
#'
#' @export
plot_ica_topography <- function(ica, eeg, component, montage = NULL, ...) {
  topo <- get_component_topography(ica, component)
  plot_topography(eeg, values = topo, montage = montage,
                   main = paste0("IC", component, " Topography"), ...)
}

#
# ============================================================================
#                    ICA_COMPONENT_SUMMARY() - PUBLIC: NUMERIC SHORTCUTS
# ============================================================================
#
# Numeric inspection tool, complementing plot_ica_sources()/
# plot_ica_topography(): one row per component with two kinds of quick
# evidence for "is this an artifact?" - no eyeballing required.
#
#   - variance_explained: how much of the total *channel-space* variance
#     (across ica$ch_names, in eeg) this one component's contribution
#     accounts for on its own, using its topography and its own time-course
#     variance: sum(topo_k^2) * var(source_k) / total_data_variance. This is
#     not the same as ranking by whitened-space variance (which is close to
#     uniform across components by construction, since the ICA rotation is
#     orthonormal on unit-variance-whitened data) - it is the same idea MNE
#     uses to rank components by real-world importance. Because ICA sources
#     are (at least pairwise) decorrelated by construction, these ratios sum
#     to exactly 1 across all components whenever n_components_ == the
#     number of channels used to fit (verified empirically for both
#     "parallel" and "deflation"). When n_components_ is smaller (the usual
#     case), they sum to less than 1 - specifically, to the same fraction of
#     total variance the PCA step captured with only that many components
#     (sum(pca_explained_variance_[1:n_components_]) /
#     sum(pca_explained_variance_)) - since whatever variance lives outside
#     that PCA subspace was never available to any component to explain.
#   - corr_<channel>: for every "external" channel present in ica$ch_names
#     (EOG/ECG/EMG/GSR/etc. - see classify_channels(), R/eeg_class.R),
#     that channel's raw correlation with every component's time-course. A
#     component correlating strongly with an EOG channel is very likely an
#     eye-movement/blink artifact, etc. If there are no external channels in
#     ica$ch_names, no corr_* columns are added.
#
# ----------------------------------------------------------------------------
# ica_component_summary() - one row per component: variance + EOG/ECG corr
# ----------------------------------------------------------------------------
#' Summarize Independent Components for Exclusion Decisions
#'
#' Builds a table with one row per fitted independent component, giving two
#' numeric signals useful for deciding what to exclude before
#' \code{apply_ica()}: how much of the total channel-space variance each
#' component accounts for, and (when available) how strongly each component
#' correlates with any EOG/ECG/other "external" channel in the data.
#' Complements the visual tools \code{\link{plot_ica_sources}} and
#' \code{\link{plot_ica_topography}} with numbers that don't require
#' eyeballing every component.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param eeg An \code{eeg} object, passed through to \code{\link{get_sources}}
#'   (does not have to be the object \code{ica} was fit on - see its docs).
#'
#' @return A data frame, one row per component (in \code{ica}'s own component
#'   order - not resorted, so row \code{k} always corresponds to component
#'   \code{k} everywhere else in the package), with columns:
#'  \describe{
#'    \item{component}{Integer, \code{1:ica$n_components_}.}
#'    \item{variance_explained}{Numeric ratio of total channel-space variance
#'      this component alone accounts for. Sums to exactly 1 across all
#'      components when \code{n_components_} equals the number of channels
#'      fit on; otherwise sums to less than 1, by the same fraction of
#'      variance the PCA step captured with that many components (nothing
#'      outside that subspace was ever available for any component to
#'      explain).}
#'    \item{corr_<channel>}{One column per external (EOG/ECG/etc.) channel in
#'      \code{ica$ch_names}, if any - that channel's correlation with each
#'      component. Omitted entirely if there are no external channels.}
#'  }
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   summary_tbl <- ica_component_summary(ica, eeg)
#'   summary_tbl[order(-summary_tbl$variance_explained), ]
#' }
#'
#' @seealso \code{\link{get_sources}}, \code{\link{plot_ica_sources}},
#'   \code{\link{plot_ica_topography}}
#'
#' @importFrom stats cor var
#' @export
ica_component_summary <- function(ica, eeg) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.")
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).")
  }

  # ========== SOURCES + MATCHED CHANNEL DATA ==========

  sources <- get_sources(ica, eeg)
  idx <- match(ica$ch_names, eeg$channels)
  data <- eeg$data[idx, , drop = FALSE]
  channel_types <- eeg$channel_types[idx]

  # ========== VARIANCE EXPLAINED ==========

  total_var <- sum(apply(data, 1, var))
  variance_explained <- vapply(seq_len(ica$n_components_), function(k) {
    topo_k <- get_component_topography(ica, k)
    (sum(topo_k^2) * var(sources[k, ])) / total_var
  }, numeric(1))

  summary_df <- data.frame(
    component = seq_len(ica$n_components_),
    variance_explained = variance_explained
  )

  # ========== EXTERNAL CHANNEL CORRELATIONS ==========

  ext_idx <- which(channel_types == "external")
  if (length(ext_idx) > 0) {
    ext_names <- ica$ch_names[ext_idx]
    for (j in seq_along(ext_idx)) {
      ext_signal <- data[ext_idx[j], ]
      col_name <- paste0("corr_", ext_names[j])
      summary_df[[col_name]] <- apply(sources, 1, function(r) cor(r, ext_signal))
    }
  }

  summary_df
}

#
# ============================================================================
#                    SET_EXCLUDE() - PUBLIC: RECORD THE USER'S CHOICE
# ============================================================================
#
# The last step of "inspect, then decide": writes which component(s) the
# user picked (after plot_ica_sources() / plot_ica_topography() /
# ica_component_summary()) into ica$exclude, ready for apply_ica() to
# actually remove them. Replaces the exclude list rather than appending to
# it - matching the "set" in the name - so calling this again with an
# updated selection is always safe and predictable.
#
# ----------------------------------------------------------------------------
# set_exclude() - record which components to exclude
# ----------------------------------------------------------------------------
#' Mark Independent Components for Exclusion
#'
#' Records which component(s) should be dropped when reconstructing clean
#' data (\code{apply_ica()}), after inspecting them with
#' \code{\link{plot_ica_sources}}, \code{\link{plot_ica_topography}}, and/or
#' \code{\link{ica_component_summary}}. Replaces \code{ica$exclude}
#' entirely (it is not additive) - call it again with the full, updated set
#' any time the choice changes.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param components An integer vector of component indices to exclude
#'   (\code{1} to \code{ica$n_components_}), or \code{NULL} to clear the
#'   exclude list back to empty. Duplicates are silently removed and the
#'   result is sorted.
#'
#' @return The same \code{eeg_ica} object passed in as \code{ica}, with
#'   \code{exclude} updated. Since R does not mutate arguments in place, the
#'   caller must reassign the result (\code{ica <- set_exclude(ica, c(1, 3))}).
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- set_exclude(ica, c(1, 4))   # mark IC1 and IC4 as artifacts
#'   ica <- set_exclude(ica, NULL)      # change your mind, clear it
#' }
#'
#' @seealso \code{\link{ica_component_summary}}, \code{\link{plot_ica_sources}},
#'   \code{\link{plot_ica_topography}}
#'
#' @export
set_exclude <- function(ica, components) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.")
  }

  # ========== NULL CLEARS THE EXCLUDE LIST ==========

  if (is.null(components)) {
    ica$exclude <- integer(0)
    return(ica)
  }

  if (!is.numeric(components) || length(components) == 0) {
    stop("ERROR: 'components' must be a non-empty integer vector, or NULL ",
         "to clear the exclude list.")
  }
  components <- as.integer(components)
  if (any(components < 1) || any(components > ica$n_components_)) {
    stop("ERROR: 'components' must be integer indices in [1, ",
         ica$n_components_, "], got: ",
         paste(components, collapse = ", "), ".")
  }

  # ========== STORE (SORTED, DEDUPED) ==========

  ica$exclude <- sort(unique(components))
  ica
}

#
# ============================================================================
#                    APPLY_ICA() - PUBLIC: RECONSTRUCT CLEAN DATA
# ============================================================================
#
# The reconstruction step: zeroes out ica$exclude in the component
# time-courses, then re-mixes back through the fitted transform in reverse
# (undo unmixing -> undo PCA whitening/rotation -> undo centering -> undo
# pre-whitening) to land back in original channel units. This is exactly
# the reverse of get_sources()'s forward chain, and reuses the same
# pinv-by-construction logic verified for get_component_topography()'s
# round-trip test (reconstructing with nothing excluded reproduces the
# original data to floating-point precision when n_components_ equals the
# number of channels fit on). Mirrors mne.preprocessing.ICA.apply(raw).
#
# ----------------------------------------------------------------------------
# apply_ica() - reconstruct eeg data with ica$exclude removed
# ----------------------------------------------------------------------------
#' Reconstruct Clean EEG Data After Excluding ICA Components
#'
#' Removes the component(s) marked via \code{\link{set_exclude}} from
#' \code{eeg}'s data and reconstructs the cleaned signal: excluded
#' components are zeroed out in component space, then everything is mixed
#' back through the fitted transform in reverse, undoing the ICA rotation,
#' PCA whitening, and pre-whitening steps \code{fit_ica()} applied going
#' forward. Only the channels \code{ica} was fit on (\code{ica$ch_names})
#' are touched - any other channel in \code{eeg} (e.g. a status channel, or
#' one excluded via \code{fit_ica()}'s \code{picks}) passes through
#' unchanged.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}).
#' @param eeg An \code{eeg} object supplying the data to clean. Does not have
#'   to be the object \code{ica} was fit on, but must contain every channel
#'   named in \code{ica$ch_names} (see \code{\link{get_sources}}, used
#'   internally).
#'
#' @return A new \code{eeg} object (see \code{\link{new_eeg}}) with the
#'   \code{ica$ch_names} rows of \code{data} replaced by the reconstructed,
#'   cleaned signal, and a note appended to \code{preprocessing_history}.
#'   Since R does not mutate arguments in place, the caller must reassign the
#'   result (\code{eeg <- apply_ica(ica, eeg)}) - the input \code{eeg} object
#'   itself is left untouched. If \code{ica$exclude} is empty, a warning is
#'   issued and \code{eeg}'s data is returned effectively unchanged (a full
#'   forward-then-reverse round trip through the fitted transform, removing
#'   nothing).
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- set_exclude(ica, c(1, 4))
#'   eeg_clean <- apply_ica(ica, eeg)
#' }
#'
#' @seealso \code{\link{set_exclude}}, \code{\link{get_sources}},
#'   \code{\link{fit_ica}}
#'
#' @export
apply_ica <- function(ica, eeg) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).")
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.")
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).")
  }

  if (length(ica$exclude) == 0) {
    warning("ica$exclude is empty - no components marked for removal (see ",
            "set_exclude()). Returning a forward-then-reverse round trip ",
            "through the fitted transform, which removes nothing.")
  }

  # ========== ZERO OUT EXCLUDED COMPONENTS ==========

  sources <- get_sources(ica, eeg)
  sources[ica$exclude, ] <- 0

  # ========== REMIX BACK THROUGH THE FITTED TRANSFORM, IN REVERSE ==========

  ev_sel <- ica$pca_explained_variance_[seq_len(ica$n_components_)]

  whitened_clean <- ica$mixing_matrix_ %*% sources
  centered_clean <- t(ica$pca_components_) %*% (sqrt(ev_sel) * whitened_clean)
  data_pw_clean  <- sweep(centered_clean, MARGIN = 1,
                           STATS = ica$pca_mean_, FUN = "+")
  data_clean     <- sweep(data_pw_clean, MARGIN = 1,
                           STATS = ica$pre_whitener_, FUN = "*")

  # ========== WRITE BACK ONLY THE FITTED CHANNELS ==========

  idx <- match(ica$ch_names, eeg$channels)
  eeg_clean <- eeg
  eeg_clean$data[idx, ] <- data_clean

  eeg_clean$preprocessing_history <- c(
    eeg_clean$preprocessing_history,
    list(paste0("ICA applied: ", length(ica$exclude), " component(s) removed",
                if (length(ica$exclude) > 0) {
                  paste0(" (", paste(ica$exclude, collapse = ", "), ")")
                } else {
                  ""
                }))
  )

  eeg_clean
}

#
# ============================================================================
#                    PLOT_ICA_OVERLAY() - PUBLIC: BEFORE/AFTER PREVIEW
# ============================================================================
#
# Lets the user preview what apply_ica() would do to a channel *before*
# actually calling it (matching open_dvm's / MNE's ica.plot_overlay()
# workflow: pick components, see the effect, only then commit). Internally
# just calls apply_ica() and overlays the original and cleaned traces for
# one channel - no reconstruction logic is duplicated here.
#
# ----------------------------------------------------------------------------
# plot_ica_overlay() - preview the before/after effect of ica$exclude
# ----------------------------------------------------------------------------
#' Preview the Effect of Excluding ICA Components (Before/After Overlay)
#'
#' Overlays a channel's original signal against what it would look like after
#' \code{\link{apply_ica}} removes the components currently marked in
#' \code{ica$exclude} - a preview to check the artifact is actually gone (and
#' nothing else changed much) before committing to \code{apply_ica()} for
#' real. Follows the same time-window/demean/event-marker conventions as
#' \code{\link{plot_eeg_signal}} and \code{\link{plot_ica_sources}}.
#'
#' @param ica A fitted \code{eeg_ica} object (\code{current_fit != "unfitted"}
#'   - see \code{fit_ica()}). A warning is issued (propagated from
#'   \code{apply_ica()}) if \code{ica$exclude} is empty, since the two traces
#'   will then be effectively identical.
#' @param eeg An \code{eeg} object supplying the data.
#' @param channel Character string or integer: the channel to plot. Must be
#'   one of the channels \code{ica} was fit on (\code{ica$ch_names}), since
#'   only those have a cleaned version to compare against. \code{NULL}
#'   (default) uses the first fitted channel, \code{ica$ch_names[1]}.
#' @param tmin,tmax Numeric, start/end of the time window to display, in
#'   milliseconds. \code{NULL} (default) plots from the start/to the end of
#'   the recording.
#' @param demean Logical, default \code{TRUE}. Each trace is centered on its
#'   own mean over the plotted window (matching \code{plot_eeg_signal()}),
#'   so the comparison is about waveform shape, not absolute offset.
#' @param show_events Logical, default \code{TRUE}. Draws event markers from
#'   \code{eeg$events}, same as \code{plot_eeg_signal()}.
#' @param event_color Character string, color for event marker lines/labels.
#'   Default \code{"orange"}.
#' @param before_color Character string, color of the original ("before")
#'   trace. Default \code{"red"}.
#' @param after_color Character string, color of the cleaned ("after") trace.
#'   Default \code{"black"}.
#' @param line_width Numeric, line width for both traces. Default \code{1}.
#' @param title Character string, custom plot title. If \code{NULL} (default),
#'   a title is generated automatically.
#'
#' @return Invisibly returns a list with:
#'  \describe{
#'    \item{channel}{Channel name plotted.}
#'    \item{excluded}{Integer vector, \code{ica$exclude} at the time of the
#'      call (the components this overlay reflects).}
#'    \item{tmin_ms, tmax_ms}{Actual start/end of the plotted window, in ms.}
#'    \item{n_samples}{Number of samples in the plotted window.}
#'    \item{n_events}{Number of event markers drawn.}
#'  }
#'
#' @examples
#' \dontrun{
#'   ica <- fit_ica(new_ica(n_components = 0.95), eeg)
#'   ica <- set_exclude(ica, 1)
#'   plot_ica_overlay(ica, eeg, channel = "Fp1")   # preview
#'   eeg <- apply_ica(ica, eeg)                    # commit, once satisfied
#' }
#'
#' @seealso \code{\link{apply_ica}}, \code{\link{set_exclude}},
#'   \code{\link{plot_eeg_signal}}
#'
#' @importFrom graphics plot lines abline mtext legend par
#' @importFrom grDevices adjustcolor
#' @export
plot_ica_overlay <- function(ica,
                              eeg,
                              channel      = NULL,
                              tmin         = NULL,
                              tmax         = NULL,
                              demean       = TRUE,
                              show_events  = TRUE,
                              event_color  = "orange",
                              before_color = "red",
                              after_color  = "black",
                              line_width   = 1,
                              title        = NULL) {

  # ========== VALIDATE inputs ==========

  if (!inherits(ica, "eeg_ica")) {
    stop("ERROR: 'ica' must be an object of class 'eeg_ica' (see new_ica()).",
         call. = FALSE)
  }
  if (identical(ica$current_fit, "unfitted")) {
    stop("ERROR: 'ica' has not been fit yet. Call fit_ica() first.",
         call. = FALSE)
  }
  if (!inherits(eeg, "eeg")) {
    stop("ERROR: 'eeg' must be an object of class 'eeg' (see new_eeg()).",
         call. = FALSE)
  }

  # ---- Resolve channel ----

  if (is.null(channel)) channel <- ica$ch_names[1]

  if (is.character(channel)) {
    ch_idx <- match(channel, eeg$channels)
    if (is.na(ch_idx)) {
      stop("ERROR: Channel '", channel, "' not found in 'eeg'.", call. = FALSE)
    }
    ch_name <- channel
  } else if (is.numeric(channel)) {
    ch_idx <- as.integer(channel)
    if (ch_idx < 1 || ch_idx > length(eeg$channels)) {
      stop("ERROR: Channel index ", ch_idx, " is out of range [1, ",
           length(eeg$channels), "].", call. = FALSE)
    }
    ch_name <- eeg$channels[ch_idx]
  } else {
    stop("ERROR: 'channel' must be a channel name (character) or index ",
         "(integer), or NULL.", call. = FALSE)
  }

  if (!(ch_name %in% ica$ch_names)) {
    stop("ERROR: Channel '", ch_name, "' was not one of the channels 'ica' ",
         "was fit on (ica$ch_names) - no cleaned version exists to overlay.",
         call. = FALSE)
  }

  if (!is.logical(demean) || length(demean) != 1) {
    stop("ERROR: 'demean' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_events) || length(show_events) != 1) {
    stop("ERROR: 'show_events' must be TRUE or FALSE.", call. = FALSE)
  }

  # ========== COMPUTE CLEANED VERSION (reuses apply_ica()) ==========

  eeg_clean <- apply_ica(ica, eeg)

  # ========== RESOLVE TIME WINDOW (ms -> sample indices) ==========

  n_timepoints <- ncol(eeg$data)
  times_ms     <- eeg$times * 1000
  rec_start_ms <- times_ms[1]
  rec_end_ms   <- times_ms[n_timepoints]

  if (is.null(tmin)) tmin <- rec_start_ms
  if (is.null(tmax)) tmax <- rec_end_ms

  if (!is.numeric(tmin) || length(tmin) != 1) {
    stop("ERROR: 'tmin' must be a single numeric value in milliseconds.",
         call. = FALSE)
  }
  if (!is.numeric(tmax) || length(tmax) != 1) {
    stop("ERROR: 'tmax' must be a single numeric value in milliseconds.",
         call. = FALSE)
  }
  if (tmin >= tmax) {
    stop("ERROR: 'tmin' (", tmin, " ms) must be less than 'tmax' (", tmax,
         " ms).", call. = FALSE)
  }
  if (tmin < rec_start_ms || tmax > rec_end_ms) {
    warning("Requested time window [", tmin, ", ", tmax, "] ms extends beyond ",
            "the recording [", round(rec_start_ms, 1), ", ",
            round(rec_end_ms, 1), "] ms. Clipping to recording bounds.",
            call. = FALSE, immediate. = TRUE)
    tmin <- max(tmin, rec_start_ms)
    tmax <- min(tmax, rec_end_ms)
  }

  samp_start <- which.min(abs(times_ms - tmin))
  samp_end   <- which.min(abs(times_ms - tmax))
  if (samp_start >= samp_end) {
    stop("ERROR: Time window too narrow - fewer than 2 samples selected. ",
         "Increase the range between 'tmin' and 'tmax'.", call. = FALSE)
  }

  # ========== EXTRACT + DEMEAN SLICES ==========

  before_slice <- eeg$data[ch_idx, samp_start:samp_end]
  after_slice  <- eeg_clean$data[ch_idx, samp_start:samp_end]
  time_slice   <- times_ms[samp_start:samp_end]
  n_samples    <- length(time_slice)

  if (demean) {
    before_slice <- before_slice - mean(before_slice, na.rm = TRUE)
    after_slice  <- after_slice - mean(after_slice, na.rm = TRUE)
  }

  # ========== IDENTIFY EVENTS IN WINDOW ==========

  events_in_window <- data.frame(onset_ms = numeric(0), type = character(0),
                                  stringsAsFactors = FALSE)
  if (show_events && !is.null(eeg$events) && nrow(eeg$events) > 0) {
    ev        <- eeg$events
    ev_ms     <- ev$onset_time * 1000
    in_window <- ev_ms >= tmin & ev_ms <= tmax
    events_in_window <- data.frame(onset_ms = ev_ms[in_window],
                                    type = ev$type[in_window],
                                    stringsAsFactors = FALSE)
  }
  n_events <- nrow(events_in_window)

  # ========== BUILD TITLE ==========

  if (is.null(title)) {
    title <- paste0(
      "ICA Overlay  |  Channel: ", ch_name,
      "  |  ", round(tmin), " - ", round(tmax), " ms"
    )
  }

  # ========== DRAW PLOT ==========

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if (show_events && n_events > 0) {
    par(mar = c(4.5, 4.5, 4, 1))
  } else {
    par(mar = c(4.5, 4.5, 3, 1))
  }

  data_range <- range(c(before_slice, after_slice), na.rm = TRUE)
  padding    <- (data_range[2] - data_range[1]) * 0.10
  y_lim      <- c(data_range[1] - padding, data_range[2] + padding)

  plot(
    x    = range(time_slice),
    y    = y_lim,
    type = "n",
    xlab = "Time (ms)",
    ylab = "Amplitude (microV)",
    main = title,
    las  = 1
  )

  abline(h = 0, col = "gray70", lty = 2, lwd = 0.8)

  if (show_events && n_events > 0) {
    for (i in seq_len(n_events)) {
      abline(v   = events_in_window$onset_ms[i],
             col = adjustcolor(event_color, alpha.f = 0.7),
             lty = 2,
             lwd = 1.2)
    }
    mtext(
      text = events_in_window$type,
      side = 3,
      at   = events_in_window$onset_ms,
      col  = event_color,
      cex  = 0.65,
      line = 0.3
    )
  }

  lines(x = time_slice, y = before_slice, col = before_color, lwd = line_width)
  lines(x = time_slice, y = after_slice, col = after_color, lwd = line_width)

  legend(
    "topright",
    legend = c("Before (original)", "After (ICA cleaned)"),
    col    = c(before_color, after_color),
    lty    = 1,
    lwd    = line_width,
    bty    = "n",
    cex    = 0.8
  )

  # ========== INVISIBLE RETURN ==========

  invisible(list(
    channel   = ch_name,
    excluded  = ica$exclude,
    tmin_ms   = round(time_slice[1], 2),
    tmax_ms   = round(time_slice[n_samples], 2),
    n_samples = n_samples,
    n_events  = n_events
  ))
}
