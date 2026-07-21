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
#'   \code{"fastica"} is implemented (backed by the \code{fastICA} package).
#' @param fit_params \code{NULL} or a named list of additional arguments
#'   passed to the underlying ICA estimator (\code{fastICA::fastICA} when
#'   \code{method = "fastica"}). Entries supplied here override the package
#'   defaults (\code{alg.typ = "parallel"}, \code{fun = "logcosh"}).
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
    } else if (n_components == 1) {
      stop("ERROR: Selecting one component with n_components = 1 is not ",
           "supported.")
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

  default_fit_params <- list(alg.typ = "parallel", fun = "logcosh")
  fit_params <- modifyList(default_fit_params, fit_params)
  fit_params$maxit <- max_iter

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
# NOTE on duplication: the EEG-vs-external channel split below is a third
# copy of the same two-pass classification already duplicated between
# detect_external_channels() (R/setexchannels.R) and print.eeg()
# (R/eeg_class.R:203-214) - see the CLAUDE.md note on that drift risk. If
# the classification logic changes, all three call sites need to move
# together (or get factored into one shared helper).
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
#'   == nrow(data)}. Used to split channels into type groups via
#'   \code{detect_external_channels()} plus the renamed-channel regex
#'   fallback (the same EEG-vs-external classification \code{print.eeg()}
#'   uses): channels flagged as external (EOG/EMG/ECG/GSR/etc.) form one
#'   group, and every remaining channel forms the \code{"eeg"} group.
#' @param noise_cov \code{NULL} (default). A non-NULL value is not yet
#'   supported and raises an error.
#' @return Named numeric vector of length \code{nrow(data)}: one population
#'   standard deviation per channel, pooled within and broadcast across
#'   each type group, named by \code{channels}.
#' @keywords internal
.compute_pre_whitener <- function(data, channels, noise_cov = NULL) {
  if (!is.null(noise_cov)) {
    stop("ERROR: noise_cov-based pre-whitening is not yet implemented. ",
         "Only the default z-score path (noise_cov = NULL) is supported.")
  }
  if (length(channels) != nrow(data)) {
    stop("ERROR: length(channels) (", length(channels),
         ") must match nrow(data) (", nrow(data), ").")
  }

  # Same two-pass EEG/external split as print.eeg() (R/eeg_class.R:203-214):
  # pass 1 catches original BioSemi EXG-style names, pass 2 catches renamed
  # channels that kept the original name in parentheses (e.g. "MASTOID LEFT
  # (EXG5)").
  exg_pattern <- paste0(
    "\\((", "EXG[1-8]|GSR[12]|Plet|Temp|Resp|Erg[12]", ")\\)"
  )
  exg_pass1 <- detect_external_channels(channels)
  exg_pass2 <- channels[grepl(exg_pattern, channels, ignore.case = TRUE)]
  exg_idx   <- which(channels %in% unique(c(exg_pass1, exg_pass2)))
  eeg_idx   <- setdiff(seq_along(channels), exg_idx)

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
