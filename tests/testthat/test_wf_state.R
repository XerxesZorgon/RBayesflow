# tests/testthat/test_wf_state.R
# UT-1: new_wf_state() constructor validation tests.

source(file.path(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), "R", "wf_state.R"))

test_that("new_wf_state validates mode and stage", {
  expect_error(new_wf_state(mode = "invalid", stage = "explore"))
  expect_error(new_wf_state(mode = "learn",   stage = "invalid"))
  wf <- new_wf_state(mode = "learn", stage = "explore")
  expect_s3_class(wf, "wf_state")
  expect_identical(wf$mode, "learn")
  expect_false(wf$diagnostics$acknowledged)
})
