# tests/testthat/test_display_contract.R
# UT-2: print.wf_state() display contract tests.

test_that("print.wf_state withholds coefficients on failed unacknowledged diagnostics", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$diagnostics$passed       <- FALSE
  wf$diagnostics$acknowledged <- FALSE
  wf$diagnostics$failed_criteria <- c("rhat_max")
  out <- capture.output(print(wf))
  expect_false(any(grepl("Estimate", out)))
  expect_true(any(grepl("diagnose", out)))
})

test_that("print.wf_state shows coefficient placeholder after acknowledgment in practice mode", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$diagnostics$passed       <- FALSE
  wf$diagnostics$acknowledged <- TRUE
  wf$diagnostics$failed_criteria <- c("rhat_max")
  out <- capture.output(print(wf))
  expect_true(any(grepl("Coefficient table", out)))
  expect_true(any(grepl("WARNING", out)))
})
