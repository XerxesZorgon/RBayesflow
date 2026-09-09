# R/wf_state.R
# Core wf_state S3 object: constructor, display contract, methods.
# Schema defined in DESIGN.md §2. ADR-004 (S3), ADR-009 (separate object).

RBAYESFLOW_VERSION <- "0.1.0"

new_wf_state <- function(mode = "learn", stage = "explore") {
  valid_modes  <- c("learn", "practice", "expert")
  valid_stages <- c("explore", "confirm", "exit")

  if (!mode  %in% valid_modes)  stop("mode must be one of: ", paste(valid_modes,  collapse = ", "))
  if (!stage %in% valid_stages) stop("stage must be one of: ", paste(valid_stages, collapse = ", "))

  structure(
    list(
      # Identity
      rbayesflow_version = RBAYESFLOW_VERSION,
      mode               = mode,
      stage              = stage,

      # Model specification
      formula      = NULL,
      family       = NULL,
      data_hash    = NULL,

      # Phase 2: Priors
      priors_text      = NULL,
      priors_objects   = NULL,
      prior_pred_draws = NULL,

      # Phase 3: Fit
      fit_hash         = NULL,
      fit_timestamp    = NULL,
      stan_backend     = NULL,
      parameterization = NULL,

      # Phase 4: Diagnostics
      diagnostics = list(
        passed            = NA,
        acknowledged      = FALSE,
        rhat_max          = NA,
        bulk_ess_min      = NA,
        tail_ess_min      = NA,
        n_divergences     = NA,
        bfmi              = NA,
        max_treedepth_hit = NA,
        family_checks     = list(),
        failed_criteria   = character()
      ),

      # Phase 5: PPC
      ppc_complete = FALSE,
      ppc_summary  = NULL,

      # Phase 6: LOO
      loo_complete = FALSE,
      loo_table    = NULL,

      # Workflow metadata
      declared_goal  = NULL,
      n_observations = NULL,
      event_rate     = NULL,

      # Audit trail and exit log
      audit_trail = list(),
      exit_log    = NULL
    ),
    class = c("wf_state", "list")
  )
}

# Stubs — fully implemented in later tasks
print.wf_state   <- function(wf, ...) cat("[wf_state stub] mode:", wf$mode, "| stage:", wf$stage, "\n")
diagnose.wf_state <- function(wf, ...) invisible(wf)
