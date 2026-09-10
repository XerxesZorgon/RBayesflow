test_that("rbf_new creates correct subfolder structure", {
  withr::with_tempdir({
    dir.create("data")
    file.create("DESCRIPTION")
    rbf_new("my_analysis")
    expect_true(dir.exists(file.path("data", "my_analysis")))
    expect_true(file.exists(file.path("data", "my_analysis", ".Rprofile")))
    expect_true(file.exists(file.path("data", "my_analysis", "wf_context.json")))
    expect_true(file.exists(file.path("data", "my_analysis", "README.md")))
  })
})

test_that("rbf_new errors if analysis folder already exists", {
  withr::with_tempdir({
    file.create("DESCRIPTION")
    dir.create(file.path("data", "existing"), recursive = TRUE)
    expect_error(rbf_new("existing"), regexp = "already exists")
  })
})
