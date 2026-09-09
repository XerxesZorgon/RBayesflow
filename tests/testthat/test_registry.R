# tests/testthat/test_registry.R
# UT-4: DIAGNOSTIC_REGISTRY dispatch and family_key() tests.

root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
source(file.path(root, "R", "context.R"))
source(file.path(root, "R", "diagnostic_registry.R"))

test_that("DIAGNOSTIC_REGISTRY bernoulli entry is a function", {
  expect_true(is.function(DIAGNOSTIC_REGISTRY[["bernoulli"]]))
})

test_that("family_key emits warning and returns 'unknown' for unrecognised family", {
  fit_custom <- structure(
    list(
      data    = data.frame(y = 1:5),
      formula = structure(list(), class = "formula")
    ),
    class = "brmsfit"
  )
  expect_warning(
    key <- family_key(fit_custom),
    regexp = "No family-specific"
  )
  expect_identical(key, "unknown")
})
