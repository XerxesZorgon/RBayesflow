# tests/testthat/test_offramps.R
# UT-6: assess_offramps() output contract tests.

test_that("assess_offramps returns at least 2 alternatives for binary rare event", {
  dat <- data.frame(
    y = c(rep(0L, 145L), rep(1L, 5L)),
    x = rnorm(150)
  )
  result <- assess_offramps(
    data         = dat,
    outcome_var  = "y",
    outcome_type = "binary",
    goal         = "coefficient estimation",
    n            = 150
  )
  expect_gte(length(result$alternatives), 2)
  methods <- sapply(result$alternatives, function(a) a$method)
  expect_true("full_stan" %in% methods)
})
