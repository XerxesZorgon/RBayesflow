# R/exit_workflow.R
# exit_workflow(): writes machine-generated YAML exit log and closes stage.
# DESIGN.md §7; SC-3; ADR-008 (readline for confirmation).
# Justification drawn from wf_state fields, not user self-report.

exit_workflow <- function(wf, method, alternatives,
                          path = "exit.yaml") {
  stopifnot(inherits(wf, "wf_state"))
  stopifnot(is.character(method),      length(method)       == 1)
  stopifnot(is.character(alternatives), length(alternatives) >= 1)

  # Rendering guard — consistent with ADR-008
  if (isTRUE(getOption("knitr.in.progress"))) {
    stop(
      "exit_workflow() must be called interactively, not inside a rendered document."
    )
  }

  # Prompt for user confirmation (mocked in tests via withr)
  ack <- readline(sprintf(
    "Confirm exit to '%s' as your analysis method? (yes/no): ", method
  ))
  if (tolower(trimws(ack)) != "yes") {
    message("Exit not confirmed. Workflow remains open.")
    return(invisible(wf))
  }

  # Build justification from wf_state fields — not user-typed
  justification <- paste0(
    "Method selected based on: ",
    if (!is.null(wf$n_observations))
      paste0("n = ", wf$n_observations, "; ")
    else "n = unspecified; ",
    if (!is.null(wf$event_rate) && !is.na(wf$event_rate))
      paste0("event rate = ", sprintf("%.1f%%", wf$event_rate * 100), "; ")
    else "",
    if (!is.null(wf$declared_goal))
      paste0("declared goal = '", wf$declared_goal, "'.")
    else "declared goal = unspecified."
  )

  # evidence_inspected: drawn from wf_state, not user input
  evidence <- list(
    sample_size    = wf$n_observations,
    event_rate     = wf$event_rate,
    declared_goal  = wf$declared_goal,
    diagnostics_run = !is.na(wf$diagnostics$passed)
  )

  # Build exit log
  exit_log <- list(
    method                = method,
    justification         = justification,
    evidence_inspected    = evidence,
    alternatives_considered = as.list(alternatives),
    user_confirmed        = TRUE,
    timestamp             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    rbayesflow_version    = wf$rbayesflow_version
  )

  # Write YAML
  yaml::write_yaml(exit_log, path)

  # Update wf_state
  wf$stage    <- "exit"
  wf$exit_log <- exit_log
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase     = 0,
    action    = "exit_workflow",
    timestamp = Sys.time(),
    method    = method,
    path      = path
  )))

  export_context(wf)

  message("Exit logged to: ", path)
  invisible(wf)
}
