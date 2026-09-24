# R/source_all.R
# Master source script — loads all RBayesflow workflow components
# in dependency order. Source this at the top of every phase template.
#
# Usage (from project root):
#   source("R/source_all.R")
#
# Order matters:
#   1. wf_state.R     — S3 object, display contract stubs, core methods
#   2. context.R      — export_context() depends on wf_state
#   3. init.R         — init_workflow() depends on wf_state + context
#   4. diagnostic_registry.R — DIAGNOSTIC_REGISTRY global, family_key()
#   5. diagnostics.R  — run_diagnostics(), detect_parameterization(),
#                        refit_noncentered() — depends on registry
#   6. offramps.R     — assess_offramps() — no upstream dependency
#   7. exit_workflow.R— exit_workflow() depends on wf_state + context
#   8. phase1_exploration.R — depends on offramps
#   9. phase2_priors.R      — depends on wf_state + context
#  10. phase3_fit.R         — depends on wf_state + diagnostics
#  11. display.R     — re-defines display stubs from wf_state.R with
#                       fuller implementations (sourced last to win)

root <- tryCatch(
  rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
  error = function(e) {
    # Strategy 2: walk up from the script's own location
    here <- normalizePath(
      Sys.getenv("RENV_PROJECT",
        unset = file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..")),
      mustWork = FALSE
    )
    if (file.exists(file.path(here, "DESCRIPTION"))) return(here)
    # Strategy 3: walk up from getwd()
    d <- getwd()
    for (i in 1:4) {
      if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
      d <- dirname(d)
    }
    stop("Cannot locate RBayesflow project root (no DESCRIPTION found).")
  }
)

source(file.path(root, "R", "wf_state.R"))
source(file.path(root, "R", "context.R"))
source(file.path(root, "R", "init.R"))
source(file.path(root, "R", "diagnostic_registry.R"))
source(file.path(root, "R", "diagnostics.R"))
source(file.path(root, "R", "diagnose_cmdstan.R"))
source(file.path(root, "R", "offramps.R"))
source(file.path(root, "R", "exit_workflow.R"))
source(file.path(root, "R", "phase1_exploration.R"))
source(file.path(root, "R", "phase2_priors.R"))
source(file.path(root, "R", "phase3_fit.R"))
source(file.path(root, "R", "display.R"))

cat("RBayesflow workflow loaded.\n")
cat("Run init_workflow(mode = 'learn') to begin.\n")
source(file.path(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), "R", "install.R"))

