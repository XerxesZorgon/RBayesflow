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

test_that("export_context writes to data/<name>/ when called from an analysis subfolder [v0.2.0]", {
  withr::with_tempdir({
    # Simulate an RBayesflow project root with a DESCRIPTION file
    file.create("DESCRIPTION")
    dir.create(file.path("data", "my_analysis"), recursive = TRUE)
    file.create(file.path("data", "my_analysis", ".Rprofile"))
    withr::with_dir(file.path("data", "my_analysis"), {
      wf <- new_wf_state(mode = "learn", stage = "explore")
      path <- export_context(wf)   # no explicit path — use new default
      expect_true(file.exists("wf_context.json"))
      expect_equal(normalizePath(path), normalizePath("wf_context.json"))
    })
  })
})
