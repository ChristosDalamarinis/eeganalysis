# ============================================================================
#                       Test File for ica_detect.R
# ============================================================================
#
# Tests Phase 2 of the ICA pipeline (see notes/mne-ica-pipeline-reference.md
# §4 and R/ica_detect.R's own file header): the automatic detectors that
# score components against a reference and flag outliers, but never touch
# ica$exclude themselves ("detect, don't decide").
#
# Numerical correctness for .find_outliers() is established with hand-
# derivable deterministic vectors (documented inline). For the three public
# detectors, fixtures are built so the *intended* signal (a known blink/
# heartbeat/muscle-like source) is recovered with near-perfect correlation
# (checked directly via the internal .score_sources()) - which is what makes
# the subsequent outlier-flagging assertions meaningful rather than
# coincidental. One deliberately-sized fixture (16 components) exists
# because the iterative z-score rule has a hard mathematical ceiling: with
# n scores (n-1 near-zero + 1 dominant outlier), the maximum achievable
# |z| is sqrt(n-1) - so a small n_components (e.g. 4) can never cross the
# default threshold = 3.0 no matter how strong the true signal is. This
# isn't a bug; it's the same ceiling MNE's own algorithm has, and it's why
# real usage (10+ components) works fine while a too-small test fixture
# would not.
#
# Author: Christos Dalamarinis
# Date: Sep - 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                              SHARED HELPERS
# ============================================================================

# 16-component fixture: a blink-like source (wide Gaussian bumps, survives a
# 1-10Hz bandpass) and a heartbeat-like source (narrow periodic bumps,
# concentrated in 8-16Hz), each also present (independently, with its own
# noise) on a dedicated "EOG_L (EXG1)"/"ECG (EXG2)" channel exactly the way
# apply_external_labels() would name them - plus 14 brain-like sinusoids so
# n_components_ = 16 gives the z-score rule enough room (sqrt(15) ~= 3.87)
# to actually flag the true outlier under the default threshold = 3.0.
make_detect_fixture <- function(n_samples = 4000, sampling_rate = 256, seed = 42) {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate

  s_blink <- rep(0, n_samples)
  for (bt in seq(200, n_samples - 200, by = 400)) {
    win <- max(1, bt - 60):min(n_samples, bt + 60)
    s_blink[win] <- s_blink[win] + exp(-((win - bt)^2) / (2 * 30^2))
  }
  s_heart <- rep(0, n_samples)
  for (bt in seq(100, n_samples - 100, by = 150)) {
    win <- max(1, bt - 5):min(n_samples, bt + 5)
    s_heart[win] <- s_heart[win] + exp(-((win - bt)^2) / (2 * 2^2))
  }

  n_brain <- 14
  S_brain <- sapply(seq_len(n_brain), function(k) sin(2 * pi * (3 + k) * t + k))
  S <- t(cbind(s_blink, s_heart, S_brain))

  n_ch <- 20
  A <- matrix(stats::runif(n_ch * nrow(S), 0.3, 1), nrow = n_ch)
  eeg_data <- A %*% S + matrix(rnorm(n_ch * n_samples, sd = 0.01), nrow = n_ch)

  eog_ch <- s_blink + rnorm(n_samples, sd = 0.01)
  ecg_ch <- s_heart + rnorm(n_samples, sd = 0.01)

  data <- rbind(eeg_data, eog_ch, ecg_ch)
  channels <- c(paste0("Ch", seq_len(n_ch)), "EOG_L (EXG1)", "ECG (EXG2)")

  eeg <- new_eeg(data = data, channels = channels, sampling_rate = sampling_rate)
  list(eeg = eeg, n_components = nrow(S))
}

# 8-channel fixture with a real (biosemi64-template) montage: a broadband
# "muscle-like" source loading heavily on the two peripheral temporal
# channels (T7/T8), and a smooth sinusoidal "brain-like" source loading
# evenly across every channel - so one recovered component should be
# flagged as peripheral/focal (muscle) and the other should not.
make_muscle_fixture <- function(n_samples = 4000, sampling_rate = 256, seed = 7) {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate

  chans <- c("Fp1", "Fp2", "F7", "F8", "T7", "T8", "Cz", "Pz")
  mont  <- create_montage(chans)

  muscle_src <- rnorm(n_samples)
  brain_src  <- sin(2 * pi * 10 * t)

  A <- matrix(0.1, nrow = length(chans), ncol = 2)
  A[c(5, 6), 1] <- 3.0  # T7, T8 (peripheral) load heavily on the muscle-like source
  A[, 2] <- 1.0          # brain-like source loads broadly/evenly

  data <- A %*% rbind(muscle_src, brain_src) +
    matrix(rnorm(length(chans) * n_samples, sd = 0.05), nrow = length(chans))

  eeg <- new_eeg(data = data, channels = chans, sampling_rate = sampling_rate)
  eeg <- set_montage(eeg, mont)
  list(eeg = eeg)
}

# ============================================================================
#                    TEST SUITE 1: .find_outliers() - private helper
# ============================================================================

test_that(".find_outliers flags a lone outlier and survives the sd=0 edge case", {
  # After masking the outlier, the 20 remaining zeros have sd = 0 - z becomes
  # NaN, and this must resolve like numpy's `nan > threshold` (False), not
  # R's native NA (which would otherwise crash the `if (!any(local_bad))`
  # check on the second iteration).
  expect_equal(.find_outliers(c(rep(0, 20), 10), threshold = 3), 21L)
})

test_that(".find_outliers flags every point that stands out from a common baseline", {
  expect_equal(.find_outliers(c(rep(1, 20), -1, -1, 20), threshold = 3), c(21L, 22L, 23L))
})

test_that(".find_outliers respects the 'tail' parameter", {
  x <- c(rep(0, 20), 10, -10)
  expect_equal(.find_outliers(x, threshold = 3, tail = 0),  c(21L, 22L))
  expect_equal(.find_outliers(x, threshold = 3, tail = 1),  21L)
  expect_equal(.find_outliers(x, threshold = 3, tail = -1), 22L)
})

test_that(".find_outliers's max_iter controls how many masking passes run", {
  # A single pass: z(9) = 3.36 (>3, flagged) but z(8) = 2.95 (<3, not yet) -
  # masking the 9 and recomputing is what lets the second pass catch the 8.
  x <- c(rep(0, 20), 8, 9)
  expect_equal(.find_outliers(x, threshold = 3, max_iter = 1), 22L)
  expect_equal(.find_outliers(x, threshold = 3, max_iter = 2), c(21L, 22L))
})

test_that(".find_outliers returns nothing when there are no real outliers", {
  expect_equal(.find_outliers(c(-1, -0.5, 0, 0.5, 1), threshold = 3), integer(0))
})

# ============================================================================
#          TEST SUITE 2: .resolve_reference_channels() - private helper
# ============================================================================

test_that(".resolve_reference_channels auto-detects by name pattern among external channels", {
  fx <- make_detect_fixture()
  expect_equal(.resolve_reference_channels(fx$eeg, NULL, "EOG"), "EOG_L (EXG1)")
  expect_equal(.resolve_reference_channels(fx$eeg, NULL, "ECG"), "ECG (EXG2)")
})

test_that(".resolve_reference_channels errors clearly when nothing matches and no ch_name given", {
  fx <- make_detect_fixture()
  expect_error(.resolve_reference_channels(fx$eeg, NULL, "GSR"), "no channel name matching")
})

test_that(".resolve_reference_channels accepts and validates an explicit ch_name", {
  fx <- make_detect_fixture()
  expect_equal(.resolve_reference_channels(fx$eeg, "Ch1", "EOG"), "Ch1")
  expect_error(.resolve_reference_channels(fx$eeg, "Bogus", "EOG"), "not found in eeg\\$channels")
})

# ============================================================================
#                    TEST SUITE 3: find_bads_eog()
# ============================================================================

test_that("find_bads_eog recovers the true blink component and leaves ica$exclude untouched", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  eog_data <- fx$eeg$data[match("EOG_L (EXG1)", fx$eeg$channels), ]
  scores <- .score_sources(ica, fx$eeg, eog_data, l_freq = 1, h_freq = 10)
  expect_true(max(abs(scores)) > 0.99)  # the fixture's blink IC should be near-perfectly recovered

  ica2 <- find_bads_eog(ica, fx$eeg)
  expect_equal(ica2$labels_$eog, unname(which.max(abs(scores))))
  expect_length(ica2$exclude, 0)  # detect, don't decide
})

test_that("find_bads_eog's measure = 'correlation' mode flags the same component directly", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  ica_z <- find_bads_eog(ica, fx$eeg)
  ica_c <- find_bads_eog(ica, fx$eeg, threshold = 0.5, measure = "correlation")
  expect_equal(ica_c$labels_$eog, ica_z$labels_$eog)
})

test_that("find_bads_eog merges results across multiple EOG-like channels, sorted by strength", {
  fx <- make_detect_fixture()
  eog_l  <- fx$eeg$data[match("EOG_L (EXG1)", fx$eeg$channels), ]
  eog_r  <- eog_l * 0.7 + rnorm(length(eog_l), sd = 0.05)

  eeg2 <- new_eeg(data = rbind(fx$eeg$data, eog_r),
                   channels = c(fx$eeg$channels, "EOG_R (EXG3)"),
                   sampling_rate = fx$eeg$sampling_rate)

  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), eeg2))
  ica2 <- find_bads_eog(ica, eeg2)

  expect_true(all(c("eog/1/EOG_L (EXG1)", "eog/2/EOG_R (EXG3)", "eog") %in% names(ica2$labels_)))
  expect_true(length(ica2$labels_$eog) >= 1)
})

# ============================================================================
#                    TEST SUITE 4: find_bads_ecg()
# ============================================================================

test_that("find_bads_ecg recovers the true heartbeat component and leaves ica$exclude untouched", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  ecg_data <- fx$eeg$data[match("ECG (EXG2)", fx$eeg$channels), ]
  scores <- .score_sources(ica, fx$eeg, ecg_data, l_freq = 8, h_freq = 16)
  expect_true(max(abs(scores)) > 0.99)

  ica2 <- find_bads_ecg(ica, fx$eeg)
  expect_equal(ica2$labels_$ecg, unname(which.max(abs(scores))))
  expect_length(ica2$exclude, 0)
})

test_that("find_bads_ecg rejects an unknown method", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))
  expect_error(find_bads_ecg(ica, fx$eeg, method = "bogus"), "'method' must be")
})

test_that("find_bads_ecg validates 'measure' even when method = 'ctps' (matches MNE, unused but still checked)", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))
  expect_error(find_bads_ecg(ica, fx$eeg, method = "ctps", measure = "bogus"), "'measure' must be")
})

# ============================================================================
#         TEST SUITE 4b: find_bads_ecg(method = "ctps") - private helpers
# ============================================================================

test_that(".qrs_detector finds evenly-spaced synthetic heartbeats to within a couple of samples", {
  set.seed(1)
  sfreq <- 256
  n_samples <- 20 * sfreq
  true_beats <- seq(200, n_samples - 200, by = round(0.6 * sfreq))
  ecg <- rep(0, n_samples)
  for (bt in true_beats) {
    win <- max(1, bt - 5):min(n_samples, bt + 5)
    ecg[win] <- ecg[win] + exp(-((win - bt)^2) / (2 * 2^2)) * 5
  }
  ecg <- ecg + rnorm(n_samples, sd = 0.05)

  detected <- .qrs_detector(ecg, sfreq)
  expect_equal(length(detected), length(true_beats))
  errs <- vapply(true_beats, function(tb) min(abs(detected - tb)), numeric(1))
  expect_true(max(errs) <= 2)
})

test_that(".qrs_detector errors on less than 3 seconds of data", {
  expect_error(.qrs_detector(rnorm(256), sfreq = 256), "at least 3 seconds")
})

test_that(".hilbert_phase returns values in [0, 1)", {
  t  <- seq_len(512) / 256
  ph <- .hilbert_phase(sin(2 * pi * 5 * t))
  expect_true(all(ph >= 0 & ph <= 1))
})

test_that(".kuiper_test + .prob_kuiper score locked phase near 1 and random phase near 0", {
  n_trials <- 50
  locked <- matrix((rep(0.3, n_trials * 10) + rnorm(n_trials * 10, sd = 0.01)) %% 1,
                    nrow = n_trials, ncol = 10)
  pk_locked <- .prob_kuiper(.kuiper_test(locked), n_trials)
  expect_true(all(pk_locked > 0.8))

  set.seed(2)
  random <- matrix(runif(n_trials * 10), nrow = n_trials, ncol = 10)
  pk_random <- .prob_kuiper(.kuiper_test(random), n_trials)
  expect_true(all(pk_random < 0.1))
})

test_that(".get_ctps_threshold returns a plausible, sfreq-dependent cutoff in (0, 1)", {
  thr_128  <- .get_ctps_threshold(128)
  thr_256  <- .get_ctps_threshold(256)
  thr_1000 <- .get_ctps_threshold(1000)

  expect_true(all(c(thr_128, thr_256, thr_1000) > 0 & c(thr_128, thr_256, thr_1000) < 1))
  expect_true(thr_128 > thr_256)
  expect_true(thr_256 > thr_1000)
})

test_that("find_bads_ecg(method = 'ctps') recovers the true heartbeat component and agrees with correlation", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  ica_corr <- find_bads_ecg(ica, fx$eeg, method = "correlation")
  ica_ctps <- find_bads_ecg(ica, fx$eeg, method = "ctps")

  expect_equal(ica_ctps$labels_$ecg, ica_corr$labels_$ecg)
  expect_true("ecg/ECG (EXG2)" %in% names(ica_ctps$labels_))
  expect_length(ica_ctps$exclude, 0)
})

test_that("find_bads_ecg(method = 'ctps', threshold = 'auto') resolves via .get_ctps_threshold", {
  fx  <- make_detect_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  ica_auto     <- find_bads_ecg(ica, fx$eeg, method = "ctps")
  ica_explicit <- find_bads_ecg(ica, fx$eeg, method = "ctps",
                                 threshold = .get_ctps_threshold(fx$eeg$sampling_rate))

  expect_equal(ica_auto$labels_$ecg, ica_explicit$labels_$ecg)
})

test_that("find_bads_ecg(method = 'ctps') warns and uses only the first channel when several ECG-like channels exist", {
  fx <- make_detect_fixture()
  ecg_l <- fx$eeg$data[match("ECG (EXG2)", fx$eeg$channels), ]
  eeg2 <- new_eeg(data = rbind(fx$eeg$data, ecg_l),
                   channels = c(fx$eeg$channels, "ECG2 (EXG3)"),
                   sampling_rate = fx$eeg$sampling_rate)
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), eeg2))

  expect_warning(find_bads_ecg(ica, eeg2, method = "ctps"), "More than one ECG-like channel")
})

# ============================================================================
#                    TEST SUITE 5: find_bads_muscle()
# ============================================================================

test_that("find_bads_muscle flags the component loading on peripheral channels, not the central one", {
  fx  <- make_muscle_fixture()
  ica <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), fx$eeg))
  ica2 <- find_bads_muscle(ica, fx$eeg)

  topo1 <- abs(get_component_topography(ica, 1))
  topo2 <- abs(get_component_topography(ica, 2))
  peripheral_ic <- if (mean(topo1[c("T7", "T8")]) > mean(topo1[c("Fp1", "Fp2", "Cz")])) 1 else 2

  expect_true(peripheral_ic %in% ica2$labels_$muscle)
  expect_length(ica2$exclude, 0)
})

test_that("find_bads_muscle falls back to slope-only (with a warning) when there's no usable montage", {
  fx  <- make_detect_fixture()  # no set_montage() ever called on this fixture
  ica <- suppressWarnings(fit_ica(new_ica(n_components = fx$n_components, random_state = 1), fx$eeg))

  expect_warning(
    ica2 <- find_bads_muscle(ica, fx$eeg),
    "too few.*positions"
  )
  expect_true(is.list(ica2$labels_))
  expect_length(ica2$exclude, 0)
})

# ============================================================================
#                    TEST SUITE 6: corrmap()
# ============================================================================

# 6-channel, 3-source fixture with spatially DISTINCT topographies per
# source (frontal / central-parietal / occipital) - unlike make_detect_fixture()
# (whose two "brain" distractors share the same, spatially uniform loading,
# fine for score-against-a-channel detectors but not for a topography-
# matching detector like corrmap(), which would then have no reliable way
# to tell the distractors apart from the template purely by shape). One
# recognizable "blink-like" component (heavy Fp1/Fp2 loading) is shared in
# spirit across every subject, standing in for a manually-confirmed template.
# has_blink = FALSE fits only the two non-blink sources, for the "a subject
# genuinely has no match" case.
make_corrmap_subject <- function(seed, has_blink = TRUE, n_samples = 3000, sampling_rate = 256) {
  set.seed(seed)
  t <- seq_len(n_samples) / sampling_rate
  chans <- c("Fp1", "Fp2", "Cz", "Pz", "O1", "O2")

  s_brain1 <- sin(2 * pi * (5 + seed) * t)
  s_brain2 <- sign(sin(2 * pi * (9 + seed) * t))

  if (has_blink) {
    s_blink <- rep(0, n_samples)
    for (bt in seq(150, n_samples - 150, by = 300)) {
      win <- max(1, bt - 20):min(n_samples, bt + 20)
      s_blink[win] <- s_blink[win] + exp(-((win - bt)^2) / (2 * 8^2))
    }
    S <- rbind(s_blink, s_brain1, s_brain2)
    A <- matrix(0.05, nrow = length(chans), ncol = 3)
    A[c(1, 2), 1] <- 3.0  # blink: frontal
    A[c(3, 4), 2] <- 2.0  # brain1: central/parietal
    A[c(5, 6), 3] <- 2.0  # brain2: occipital
    n_comp <- 3
  } else {
    S <- rbind(s_brain1, s_brain2)
    A <- matrix(0.05, nrow = length(chans), ncol = 2)
    A[c(3, 4), 1] <- 2.0
    A[c(5, 6), 2] <- 2.0
    n_comp <- 2
  }

  data <- A %*% S + matrix(rnorm(length(chans) * n_samples, sd = 0.02), nrow = length(chans))
  eeg  <- new_eeg(data = data, channels = chans, sampling_rate = sampling_rate)
  suppressWarnings(fit_ica(new_ica(n_components = n_comp, random_state = 1), eeg))
}

# Programmatically finds whichever component loads most heavily on Fp1/Fp2
# relative to the rest - used to get the "blink" template without hardcoding
# a component index that FastICA's convergence order could shuffle.
find_frontal_component <- function(ica) {
  scores <- vapply(seq_len(ica$n_components_), function(k) {
    topo <- abs(get_component_topography(ica, k))
    mean(topo[c("Fp1", "Fp2")]) - mean(topo[c("Cz", "Pz", "O1", "O2")])
  }, numeric(1))
  which.max(scores)
}

test_that(".find_max_corrs matches sign-flipped topographies and RMS-normalizes the average", {
  target <- c(1, 1, -1, -1)
  m1 <- rbind(c(0.1, -0.2, 0.3, 0.1), c(0.9, 0.95, -0.9, -0.85), c(-0.1, 0.2, 0.1, -0.3))
  m2 <- rbind(c(-0.9, -0.95, 0.9, 0.85), c(0.2, -0.1, 0.3, 0.2), c(0.1, 0.1, -0.2, 0.1))

  res <- .find_max_corrs(list(m1, m2), target, threshold = 0.8)
  expect_equal(res$subj_idx, list(2L, 1L))
  # m2[1, ] == -m1[2, ] exactly, so after sign-correction both normalized
  # topographies are identical and their average equals either one alone.
  expected <- m1[2, ] / sqrt(sum(m1[2, ]^2))
  expect_equal(res$newtarget, expected, tolerance = 1e-8)
})

test_that(".find_max_corrs returns a NULL newtarget when nothing clears the threshold", {
  target <- c(1, 1, -1, -1)
  m1 <- rbind(c(0.1, -0.2, 0.3, 0.1), c(0.9, 0.95, -0.9, -0.85))
  res <- .find_max_corrs(list(m1), target, threshold = 0.9999)
  expect_null(res$newtarget)
  expect_equal(res$median_corr, 0)
})

test_that("corrmap recovers the shared component across every subject", {
  icas <- list(make_corrmap_subject(1), make_corrmap_subject(2), make_corrmap_subject(3))
  template <- get_component_topography(icas[[1]], find_frontal_component(icas[[1]]))

  icas2 <- corrmap(icas, template, label = "blink")

  for (i in seq_along(icas2)) {
    expect_equal(icas2[[i]]$labels_$blink, find_frontal_component(icas[[i]]))
    expect_length(icas2[[i]]$exclude, 0)  # detect, don't decide
  }
})

test_that("corrmap gives a subject with no matching component an explicit empty entry", {
  icas <- list(make_corrmap_subject(1, TRUE), make_corrmap_subject(2, TRUE),
               make_corrmap_subject(4, FALSE))
  template <- get_component_topography(icas[[1]], find_frontal_component(icas[[1]]))

  icas2 <- corrmap(icas, template, label = "blink")

  expect_length(icas2[[3]]$labels_$blink, 0)
  expect_false(is.null(icas2[[3]]$labels_$blink))  # present but empty, not absent
})

test_that("corrmap with label = NULL is a dry run - labels_ stays untouched", {
  icas <- list(make_corrmap_subject(1), make_corrmap_subject(2))
  template <- get_component_topography(icas[[1]], find_frontal_component(icas[[1]]))

  icas2 <- corrmap(icas, template)
  expect_true(all(vapply(icas2, function(ic) length(ic$labels_) == 0, logical(1))))
})

test_that("corrmap appends to (not overwrites) an existing label on a second call", {
  icas <- list(make_corrmap_subject(1), make_corrmap_subject(2))
  template <- get_component_topography(icas[[1]], find_frontal_component(icas[[1]]))

  icas2 <- corrmap(icas, template, label = "blink")
  icas3 <- corrmap(icas2, template, label = "blink")

  expect_equal(icas3[[1]]$labels_$blink, icas2[[1]]$labels_$blink)
})

test_that("corrmap errors when the icas don't share the same ch_names", {
  icas <- list(make_corrmap_subject(1))
  eeg_other <- new_eeg(data = matrix(rnorm(3 * 200), nrow = 3),
                        channels = c("A", "B", "C"), sampling_rate = 256)
  ica_other <- suppressWarnings(fit_ica(new_ica(n_components = 2, random_state = 1), eeg_other))

  template <- get_component_topography(icas[[1]], find_frontal_component(icas[[1]]))
  expect_error(corrmap(list(icas[[1]], ica_other), template), "same ch_names")
})

test_that("corrmap validates 'template' length against the fitted channel count", {
  icas <- list(make_corrmap_subject(1))
  expect_error(corrmap(icas, template = c(1, 2, 3)), "topography vector of length")
})
