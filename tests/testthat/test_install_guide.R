test_that("guide detects Phase 1 for a fresh wf_state", {
  wf <- new_wf_state(mode = "learn", stage = "explore")
  out <- capture.output(guide(wf))
  expect_true(any(grepl("Phase 1", out)))
  expect_true(any(grepl("phase1_exploration", out)))
})

test_that("guide detects diagnostic block when diagnostics failed and unacknowledged", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$fit_timestamp            <- Sys.time()
  wf$diagnostics$passed       <- FALSE
  wf$diagnostics$acknowledged <- FALSE
  out <- capture.output(guide(wf))
  expect_true(any(grepl("diagnose", out)))
})

test_that("guide detects exit stage", {
  wf <- new_wf_state(mode = "practice", stage = "exit")
  out <- capture.output(guide(wf))
  expect_true(any(grepl("exit", tolower(out))))
})
