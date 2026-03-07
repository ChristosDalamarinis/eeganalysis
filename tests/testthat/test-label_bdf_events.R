# ============================================================================
#               Test File for label_bdf_events.R
# ============================================================================
#
# Functions tested:
#   1. label_bdf_events()       - Interactive trigger labeling (non-interactive
#                                 guard tested; label application mocked)
#   2. apply_trigger_labels()   - Programmatic label application
#
# Testing strategy for label_bdf_events():
#   Because the function calls readline() it cannot be driven end-to-end in
#   an automated suite. The tests therefore cover:
#     (a) all input-validation and guard branches (no readline() reached)
#     (b) the non-interactive guard itself
#     (c) apply_trigger_labels(), which executes the full post-prompt logic
#         and is fully testable without user interaction.
#
# Author: Christos Dalamarinis
# Date: March 2026
# ============================================================================

library(testthat)
library(eeganalysis)

# ============================================================================
#                      SHARED HELPER
# ============================================================================

#' Build a minimal mock events data frame
make_mock_events <- function(types      = c("1", "2", "1", "2", "3"),
                             onset      = c(100L, 200L, 300L, 400L, 500L),
                             onset_time = c(0.5, 1.0, 1.5, 2.0, 2.5)) {
  df <- data.frame(
    onset       = onset,
    onset_time  = onset_time,
    type        = types,
    description = paste0("Trigger: ", types),
    stringsAsFactors = FALSE
  )
  attr(df, "sample_rate") <- 200L
  attr(df, "bdf_file")    <- "mock.bdf"
  df
}


# ============================================================================
#   TEST SUITE 1: label_bdf_events() — input validation & guard branches
# ============================================================================

# ----------------------------------------------------------------------------
# Test 1.1: Rejects non-data-frame input
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The very first guard in label_bdf_events() that calls
# stop() when 'events' is not a data frame.
test_that("label_bdf_events() rejects non-data-frame input", {
  expect_error(
    label_bdf_events(list(type = "1", onset_time = 0.5)),
    "must be a data frame"
  )
  expect_error(
    label_bdf_events("not_a_df"),
    "must be a data frame"
  )
  expect_error(
    label_bdf_events(NULL),
    "must be a data frame"
  )
})

# ----------------------------------------------------------------------------
# Test 1.2: Rejects events data frame missing required columns
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The column-check guard that produces a targeted error
# message listing whichever required columns are absent.
test_that("label_bdf_events() rejects missing required columns", {
  # Missing 'onset_time'
  bad1 <- data.frame(type = c("1", "2"), stringsAsFactors = FALSE)
  expect_error(label_bdf_events(bad1), "onset_time")
  
  # Missing 'type'
  bad2 <- data.frame(onset_time = c(0.5, 1.0), stringsAsFactors = FALSE)
  expect_error(label_bdf_events(bad2), "type")
  
  # Missing both
  bad3 <- data.frame(onset = c(100L, 200L), stringsAsFactors = FALSE)
  expect_error(label_bdf_events(bad3), "type")
})

# ----------------------------------------------------------------------------
# Test 1.3: Returns early with warning for zero-row events
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The early-return branch triggered by nrow(events) == 0,
# which should warn rather than error and return the frame with a 'label' col.
test_that("label_bdf_events() warns and returns early on empty events", {
  empty_df <- data.frame(
    onset      = integer(0),
    onset_time = numeric(0),
    type       = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- label_bdf_events(empty_df),
    "zero rows"
  )
  expect_true("label" %in% names(result))
  expect_equal(nrow(result), 0L)
})

# ----------------------------------------------------------------------------
# Test 1.4: Non-interactive guard raises an informative error
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The interactive() check that prevents readline() being
# called during automated runs (R CMD check, CI, testthat).
# NOTE: When running tests from RStudio, interactive() returns TRUE, so we
# must mock it explicitly with local_mocked_bindings() rather than relying
# on the environment being non-interactive.
test_that("label_bdf_events() errors in non-interactive context", {
  # This guard is only reachable in non-interactive sessions (R CMD check, CI).
  # In RStudio interactive() == TRUE so we skip — the guard is tested naturally
  # during automated runs.
  skip_if(interactive(), "Skipping in interactive session — guard runs in CI/R CMD check")
  events <- make_mock_events()
  expect_error(
    label_bdf_events(events),
    "interactive R session"
  )
})


# ============================================================================
#   TEST SUITE 2: apply_trigger_labels() — full programmatic path
# ============================================================================

# ----------------------------------------------------------------------------
# Test 2.1: Adds 'label' column matching the supplied mapping
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: Core mapping logic — every row in events receives the
# correct label according to its 'type' value.
test_that("apply_trigger_labels() adds correct label column", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard", "2" = "target", "3" = "response")
  
  result <- apply_trigger_labels(events, mapping, verbose = FALSE)
  
  expect_true("label" %in% names(result))
  expect_equal(result$label[result$type == "1"], rep("standard", 2))
  expect_equal(result$label[result$type == "2"], rep("target",   2))
  expect_equal(result$label[result$type == "3"], rep("response", 1))
})

# ----------------------------------------------------------------------------
# Test 2.2: Attaches trigger_labels attribute identical to mapping
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: That attr(result, "trigger_labels") is set and matches the
# input mapping (when all codes are present in the mapping).
test_that("apply_trigger_labels() attaches trigger_labels attribute", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard", "2" = "target", "3" = "response")
  
  result  <- apply_trigger_labels(events, mapping, verbose = FALSE)
  
  attached <- attr(result, "trigger_labels")
  expect_false(is.null(attached))
  expect_equal(attached[["1"]], "standard")
  expect_equal(attached[["2"]], "target")
  expect_equal(attached[["3"]], "response")
})

# ----------------------------------------------------------------------------
# Test 2.3: Unmapped codes receive the unmapped_label default ("unknown")
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The setdiff() branch that detects codes absent from the
# mapping and assigns the fallback label.
test_that("apply_trigger_labels() assigns 'unknown' to unmapped codes", {
  events  <- make_mock_events()  # types: "1", "2", "3"
  mapping <- c("1" = "standard")  # "2" and "3" are deliberately absent
  
  result <- apply_trigger_labels(events, mapping,
                                 unmapped_label = "unknown",
                                 verbose = FALSE)
  
  expect_equal(unique(result$label[result$type == "2"]), "unknown")
  expect_equal(unique(result$label[result$type == "3"]), "unknown")
})

# ----------------------------------------------------------------------------
# Test 2.4: Custom unmapped_label is honoured
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: That the unmapped_label parameter correctly overrides the
# default "unknown" string.
test_that("apply_trigger_labels() respects custom unmapped_label", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard")
  
  result <- apply_trigger_labels(events, mapping,
                                 unmapped_label = "NEEDS_LABEL",
                                 verbose = FALSE)
  
  expect_true(all(result$label[result$type %in% c("2", "3")] == "NEEDS_LABEL"))
})

# ----------------------------------------------------------------------------
# Test 2.5: Verbose = TRUE reports unmapped codes to the console
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: That the cat() block inside the unmapped-codes branch
# actually produces output when verbose is TRUE.
test_that("apply_trigger_labels() prints warning for unmapped codes when verbose", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard")
  
  output <- capture.output(
    apply_trigger_labels(events, mapping, verbose = TRUE)
  )
  
  expect_true(any(grepl("WARNING", output, fixed = TRUE)))
})

# ----------------------------------------------------------------------------
# Test 2.6: No output when all codes are mapped and verbose = TRUE
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: That no warning message is emitted when every code in the
# events frame has a corresponding entry in the mapping.
test_that("apply_trigger_labels() is silent when all codes mapped", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard", "2" = "target", "3" = "response")
  
  output <- capture.output(
    apply_trigger_labels(events, mapping, verbose = TRUE)
  )
  
  expect_equal(length(output), 0L)
})

# ----------------------------------------------------------------------------
# Test 2.7: Rejects non-data-frame input
# ----------------------------------------------------------------------------
test_that("apply_trigger_labels() rejects non-data-frame events", {
  expect_error(
    apply_trigger_labels(list(), c("1" = "x")),
    "must be a data frame"
  )
})

# ----------------------------------------------------------------------------
# Test 2.8: Rejects events missing the 'type' column
# ----------------------------------------------------------------------------
test_that("apply_trigger_labels() rejects events without 'type' column", {
  bad <- data.frame(onset_time = c(0.5, 1.0), stringsAsFactors = FALSE)
  expect_error(
    apply_trigger_labels(bad, c("1" = "standard")),
    "'type' column"
  )
})

# ----------------------------------------------------------------------------
# Test 2.9: Rejects unnamed mapping vectors
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The guard that checks is.character(mapping) && !is.null(names(mapping)).
test_that("apply_trigger_labels() rejects unnamed mapping", {
  events <- make_mock_events()
  expect_error(
    apply_trigger_labels(events, c("standard", "target")),
    "named character vector"
  )
  expect_error(
    apply_trigger_labels(events, 1:3),
    "named character vector"
  )
})

# ----------------------------------------------------------------------------
# Test 2.10: Returns early with warning for zero-row events
# ----------------------------------------------------------------------------
test_that("apply_trigger_labels() warns and returns early on empty events", {
  empty_df <- data.frame(
    onset      = integer(0),
    onset_time = numeric(0),
    type       = character(0),
    stringsAsFactors = FALSE
  )
  mapping <- c("1" = "standard")
  
  expect_warning(
    result <- apply_trigger_labels(empty_df, mapping, verbose = FALSE),
    "zero rows"
  )
  expect_true("label" %in% names(result))
  expect_equal(nrow(result), 0L)
})

# ----------------------------------------------------------------------------
# Test 2.11: Original columns are preserved unchanged
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: That apply_trigger_labels() is non-destructive — existing
# columns like 'onset', 'onset_time', 'type', 'description' are untouched.
test_that("apply_trigger_labels() does not modify existing columns", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard", "2" = "target", "3" = "response")
  
  result  <- apply_trigger_labels(events, mapping, verbose = FALSE)
  
  expect_equal(result$onset,      events$onset)
  expect_equal(result$onset_time, events$onset_time)
  expect_equal(result$type,       events$type)
  expect_equal(result$description, events$description)
})

# ----------------------------------------------------------------------------
# Test 2.12: Returns invisibly
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The invisible() wrapper on the return value, consistent
# with the conventions used across extract_bdf_events.R.
test_that("apply_trigger_labels() returns invisibly", {
  events  <- make_mock_events()
  mapping <- c("1" = "standard", "2" = "target", "3" = "response")
  
  # withVisible() detects whether the value is returned visibly or not
  vis <- withVisible(apply_trigger_labels(events, mapping, verbose = FALSE))
  expect_false(vis$visible)
})

# ----------------------------------------------------------------------------
# Test 2.13: Integration — mapping round-trip from label_bdf_events attribute
# ----------------------------------------------------------------------------
# WHAT THIS TESTS: The intended multi-participant workflow: extract the
# trigger_labels attribute written by label_bdf_events (simulated here by
# setting it manually) and pass it to apply_trigger_labels().
test_that("apply_trigger_labels() round-trips a trigger_labels attribute", {
  # Simulate events_p1 that was already labeled interactively
  events_p1 <- make_mock_events()
  attr(events_p1, "trigger_labels") <- c(
    "1" = "standard",
    "2" = "target",
    "3" = "response"
  )
  
  # Apply to a second participant with the same trigger codes
  events_p2  <- make_mock_events()
  mapping    <- attr(events_p1, "trigger_labels")
  result_p2  <- apply_trigger_labels(events_p2, mapping, verbose = FALSE)
  
  # apply_trigger_labels() strips names when assigning into a data frame column,
  # but mapping[events_p2$type] returns a named vector — unname() both sides.
  expect_equal(result_p2$label, unname(mapping[events_p2$type]))
  expect_equal(attr(result_p2, "trigger_labels"), mapping)
})

# ============================================================================
# SUMMARY
# ============================================================================
#
# label_bdf_events() — 4 tests:
#   1.1  Non-data-frame input rejected
#   1.2  Missing required columns rejected
#   1.3  Zero-row events: warning + early return
#   1.4  Non-interactive guard raised
#
# apply_trigger_labels() — 13 tests:
#   2.1  Correct label column values
#   2.2  trigger_labels attribute set correctly
#   2.3  Unmapped codes → "unknown"
#   2.4  Custom unmapped_label respected
#   2.5  Verbose prints warning for unmapped codes
#   2.6  Silent when all codes mapped
#   2.7  Non-data-frame events rejected
#   2.8  Missing 'type' column rejected
#   2.9  Unnamed mapping rejected
#   2.10 Zero-row events: warning + early return
#   2.11 Existing columns preserved
#   2.12 Returns invisibly
#   2.13 Round-trip integration with trigger_labels attribute
# ============================================================================