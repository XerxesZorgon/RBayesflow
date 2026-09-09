# tests/testthat/test_context.R
# UT-5: export_context() JSON schema tests.

root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
source(file.path(root, "R", "context.R"))

test_that("export_context produces valid JSON with required fields", {
  wf <- new_wf_state(mode = "learn", stage = "explore")
  withr::with_tempdir({
    path <- export_context(wf, path = "wf_context.json")
    json <- jsonlite::fromJSON(path)
    expect_equal(json$mode, "learn")
    expect_equal(json$stage, "explore")
    expect_false(is.null(json$rbayesflow_version))
  })
})
