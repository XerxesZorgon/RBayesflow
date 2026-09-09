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

# --- Display helper stubs (fully implemented in R/display.R, Task 035) ---

cat_wf_header <- function(wf) {
  cat("=== RBayesflow workflow state ===\n")
  cat("Mode:", wf$mode, "| Stage:", wf$stage, "\n")
  if (!is.null(wf$formula))
    cat("Formula:", deparse(wf$formula), "\n")
  cat("Diagnostics: not yet run\n")
}

cat_health_summary_one_line <- function(wf) {
  cat("Diagnostics: PASSED (Rhat_max =", wf$diagnostics$rhat_max,
      "| ESS_bulk_min =", wf$diagnostics$bulk_ess_min,
      "| divergences =", wf$diagnostics$n_divergences, ")\n")
}

cat_coefficient_table <- function(wf) {
  cat("[Coefficient table: source fit object with fixef() or brms::fixef()]\n")
}

cat_diagnostic_failure_message <- function(wf) {
  cat("!!! DIAGNOSTIC FAILURE !!!\n")
  cat("Failed criteria:", paste(wf$diagnostics$failed_criteria, collapse = ", "), "\n")
}

plot_prior_posterior_overlay <- function(wf) {
  cat("[Prior-vs-posterior overlay: implemented in R/display.R]\n")
  invisible(NULL)
}

# --- Full display contract (DESIGN.md §3) ---

print.wf_state <- function(wf, ...) {
  stopifnot(inherits(wf, "wf_state"))

  if (is.na(wf$diagnostics$passed)) {
    # Fit not yet run
    cat_wf_header(wf)

  } else if (isTRUE(wf$diagnostics$passed)) {
    # Clean diagnostics
    if (wf$mode == "learn") {
      plot_prior_posterior_overlay(wf)
      cat_health_summary_one_line(wf)
    } else {
      cat_coefficient_table(wf)
      cat_health_summary_one_line(wf)
    }

  } else {
    # Failed diagnostics
    if (!isTRUE(wf$diagnostics$acknowledged)) {
      cat_diagnostic_failure_message(wf)
      cat("-> Call wf$diagnose() to review the failing checks.\n")
    } else {
      if (wf$mode == "learn") {
        cat_diagnostic_failure_message(wf)
        cat("(Diagnostics reviewed and acknowledged.)\n")
      } else {
        cat_coefficient_table(wf)
        cat("WARNING: Diagnostics failed (acknowledged). Coefficients shown with caveat.\n")
      }
    }
  }

  invisible(wf)
}

# diagnose stub — fully implemented in Task 026
diagnose.wf_state <- function(wf, ...) invisible(wf)
