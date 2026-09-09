# tests/testthat/helpers.R
# Shared test helpers loaded automatically by testthat.

root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
source(file.path(root, "R", "wf_state.R"))
source(file.path(root, "R", "offramps.R"))

# mock_brmsfit(): minimal brmsfit stub with a fixef() method.
mock_brmsfit <- function() {
  m <- list(
    data    = data.frame(y = rnorm(10), x = rnorm(10)),
    formula = y ~ x
  )
  class(m) <- "brmsfit"
  m
}

# fixef.brmsfit stub so cat_coefficient_table() can call it
fixef.brmsfit <- function(object, ...) {
  mat <- matrix(c(0.5, 0.1, 5.0, 0.0003), nrow = 1,
                dimnames = list("Intercept",
                                c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")))
  mat
}

# mock_fit(): minimal stub for family dispatch tests (used in Task 022)
mock_fit <- function(family = "gaussian") {
  m <- list(family = list(family = family))
  class(m) <- "brmsfit"
  m
}
