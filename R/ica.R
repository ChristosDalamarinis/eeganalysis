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
#' Depends on: eeg_class.R (for eeg object structure)
#'
#' Author: Christos Dalamarinis
#' Date: July 2026
#' Status: In development
#' ============================================================================
