test_that("rbf_install returns a named logical vector with 7 entries", {
  with_mocked_bindings(
    check_cmdstan_toolchain = function(...) invisible(NULL),
    cmdstan_version         = function(...) "2.35.0",
    .package = "cmdstanr",
    {
      result <- rbf_install(dry_run = TRUE)
      expect_type(result, "logical")
      expect_length(result, 7)
      expect_named(result)
    }
  )
})
