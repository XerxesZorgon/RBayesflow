# tests/testthat/test_exit_workflow.R
# UT-3: exit_workflow() YAML schema compliance test.

test_that("exit_workflow produces a valid YAML log with all required fields", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$declared_goal  <- "coefficient estimation"
  wf$n_observations <- 150L
  wf$event_rate     <- 0.03

  withr::with_tempdir({
    local_mocked_bindings(readline = function(...) "yes", .package = "base")
    log_path <- exit_workflow(
      wf           = wf,
      method       = "firth_logistic",
      alternatives = c("logistic_bootstrap", "full_stan"),
      path         = "exit.yaml"
    )
    log <- yaml::read_yaml("exit.yaml")
    expect_true(log$user_confirmed)
    expect_equal(log$method, "firth_logistic")
    expect_false(is.null(log$justification))
    expect_false(is.null(log$evidence_inspected$sample_size))
    expect_true("full_stan" %in% log$alternatives_considered)
  })
})
