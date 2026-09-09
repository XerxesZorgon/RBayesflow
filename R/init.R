# R/init.R
# init_workflow(): entry point for every RBayesflow session.
# Creates a fresh wf_state, exports context, prints startup message.
# ADR-011: mode = "expert" accepted and fully supported.

init_workflow <- function(mode = "learn", stage = "explore") {
  # new_wf_state() validates mode and stage — will stop() on invalid values
  wf <- new_wf_state(mode = mode, stage = stage)

  # Write initial context file for Posit Assistant
  ctx_path <- export_context(wf, path = "wf_context.json")

  # Startup message
  cat("=== RBayesflow", wf$rbayesflow_version, "===\n")
  cat("Mode :", wf$mode,  "\n")
  cat("Stage:", wf$stage, "\n")
  cat("Posit Assistant context written to:", ctx_path, "\n")
  cat("Run source(\"R/source_all.R\") if not already loaded.\n")

  invisible(wf)
}
