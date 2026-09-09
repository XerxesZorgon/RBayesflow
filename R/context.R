# R/context.R
# Exports wf_state to wf_context.json for Posit Assistant integration.
# ADR-006: file-based context; Posit Assistant reads this file from
# the project root at the start of each session.

export_context <- function(wf, path = "wf_context.json") {
  stopifnot(inherits(wf, "wf_state"))

  # Build a serialisable summary — strip non-serialisable fields
  ctx <- list(
    rbayesflow_version    = wf$rbayesflow_version,
    mode                  = wf$mode,
    stage                 = wf$stage,
    fit_timestamp         = if (!is.null(wf$fit_timestamp))
                              format(wf$fit_timestamp, "%Y-%m-%dT%H:%M:%S")
                            else NULL,
    parameterization      = wf$parameterization,
    diagnostics_passed    = wf$diagnostics$passed,
    diagnostics_acknowledged = wf$diagnostics$acknowledged,
    rhat_max              = wf$diagnostics$rhat_max,
    bulk_ess_min          = wf$diagnostics$bulk_ess_min,
    tail_ess_min          = wf$diagnostics$tail_ess_min,
    n_divergences         = wf$diagnostics$n_divergences,
    failed_criteria       = wf$diagnostics$failed_criteria,
    ppc_complete          = wf$ppc_complete,
    loo_complete          = wf$loo_complete,
    declared_goal         = wf$declared_goal,
    n_observations        = wf$n_observations,
    event_rate            = wf$event_rate,
    audit_trail_length    = length(wf$audit_trail),
    stage_is_exit         = identical(wf$stage, "exit")
  )

  jsonlite::write_json(ctx, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}
