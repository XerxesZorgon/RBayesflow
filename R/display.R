# R/display.R
# Mode-aware display helpers for wf_state print contract.
# Replaces stubs in R/wf_state.R (Task 010). DESIGN.md §3.

# --- cat_wf_header() ---
cat_wf_header <- function(wf) {
  cat("=== RBayesflow workflow state ===\n")
  cat("Mode :", wf$mode, "| Stage:", wf$stage,
      "| Version:", wf$rbayesflow_version, "\n")
  if (!is.null(wf$formula))
    cat("Formula:", deparse(wf$formula), "\n")
  ts <- if (!is.null(wf$fit_timestamp))
    format(wf$fit_timestamp) else "not yet fitted"
  cat("Fit timestamp:", ts, "\n")
  cat("Diagnostics: not yet run\n")
}

# --- cat_health_summary_one_line() ---
cat_health_summary_one_line <- function(wf) {
  cat(sprintf(
    "Diagnostics: PASSED (Rhat_max = %.4f | ESS_bulk_min = %d | divergences = %d)\n",
    wf$diagnostics$rhat_max,
    as.integer(wf$diagnostics$bulk_ess_min),
    as.integer(wf$diagnostics$n_divergences)
  ))
}

# --- cat_coefficient_table() ---
# Shows fixef() summary from brms (or placeholder if fit not attached).
cat_coefficient_table <- function(wf) {
  # wf itself doesn't hold the fit object; coefficient display is
  # context-dependent. In interactive use the fit is in the user's
  # environment. Display a reminder unless wf carries a fit reference.
  cat("[Coefficient table: call fixef(fit) or brms::fixef(fit) directly.]\n")
  cat("[Tip: summary(fit) shows the full posterior summary.]\n")
}

# --- cat_diagnostic_failure_message() ---
cat_diagnostic_failure_message <- function(wf) {
  cat("!!! DIAGNOSTIC FAILURE !!!\n")
  if (length(wf$diagnostics$failed_criteria) > 0) {
    cat("Failed criteria:\n")
    for (criterion in wf$diagnostics$failed_criteria) {
      val <- switch(criterion,
        rhat_max      = sprintf("  Rhat_max = %.4f (threshold < 1.01)",
                                wf$diagnostics$rhat_max),
        bulk_ess_min  = sprintf("  bulk_ESS_min = %d (threshold > 400)",
                                as.integer(wf$diagnostics$bulk_ess_min)),
        tail_ess_min  = sprintf("  tail_ESS_min = %d (threshold > 400)",
                                as.integer(wf$diagnostics$tail_ess_min)),
        n_divergences = sprintf("  Divergent transitions = %d (threshold = 0)",
                                as.integer(wf$diagnostics$n_divergences)),
        bfmi          = sprintf("  BFMI = %s (threshold > 0.3)",
                                paste(round(wf$diagnostics$bfmi, 3),
                                      collapse = ", ")),
        max_treedepth = "  Max treedepth hit (increase max_treedepth)",
        criterion     # fallback: print the name
      )
      cat(val, "\n")
    }
  }
  # Family-specific warnings
  fw <- wf$diagnostics$family_checks$warnings
  if (!is.null(fw) && length(fw) > 0) {
    cat("Family-specific warnings:\n")
    for (w in fw) cat(" *", w, "\n")
  }
}

# --- plot_prior_posterior_overlay() ---
# learn mode: overlay prior predictive draws with posterior draws.
plot_prior_posterior_overlay <- function(wf) {
  if (is.null(wf$prior_pred_draws)) {
    cat("[Prior-vs-posterior overlay: prior predictive draws not yet computed.",
        "Run Phase 2 first.]\n")
    return(invisible(NULL))
  }
  cat("[Prior-vs-posterior overlay plot]\n")
  # In a real session the fit object must be in scope; this helper
  # prints a placeholder when called without a fit argument.
  # Full overlay: bayesplot::mcmc_areas() + geom_density(prior_pred_draws)
  # See templates/phase4_diagnostics.qmd for the interactive version.
  invisible(NULL)
}

# --- show_diagnostic_plots() ---
# Called from diagnose.wf_state(); displays bayesplot panels.
show_diagnostic_plots <- function(wf) {
  if (length(wf$diagnostics$failed_criteria) == 0 &&
      !is.na(wf$diagnostics$passed)) {
    cat("All diagnostics passed. No panels to display.\n")
    return(invisible(NULL))
  }
  cat("--- Diagnostic Plot Panels ---\n")
  cat("[mcmc_trace():   run bayesplot::mcmc_trace(fit) to see chain mixing]\n")
  cat("[mcmc_rhat():    run bayesplot::mcmc_rhat(brms::rhat(fit)) for Rhat distribution]\n")
  cat("[mcmc_neff():    run bayesplot::mcmc_neff(brms::neff_ratio(fit)) for ESS ratios]\n")
  if (!is.null(wf$diagnostics$n_divergences) &&
      !is.na(wf$diagnostics$n_divergences) &&
      wf$diagnostics$n_divergences > 0) {
    cat("[mcmc_pairs():   run bayesplot::mcmc_pairs(fit, ...) to see divergence geometry]\n")
  }
  invisible(NULL)
}
