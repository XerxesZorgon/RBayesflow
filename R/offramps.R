# R/offramps.R
# assess_offramps(): non-Bayesian off-ramp assessment (DESIGN.md §4).
# ADR-010: event rate computed from raw data, not prior predictive.
# All alternatives presented with equal visual weight; none marked recommended.

assess_offramps <- function(data, outcome_var, outcome_type, goal,
                            n = nrow(data)) {
  stopifnot(is.data.frame(data))
  stopifnot(outcome_type %in% c("binary", "count", "continuous"))
  stopifnot(is.character(goal), length(goal) == 1)

  # ADR-010: compute event rate from raw data before any fitting
  event_rate <- NA_real_
  if (outcome_type == "binary") {
    y_raw <- data[[outcome_var]]
    # Handle factor, character, or 0/1 coding
    if (is.factor(y_raw) || is.character(y_raw)) {
      y_int <- as.integer(factor(y_raw)) - 1L
    } else {
      y_int <- as.integer(y_raw)
    }
    if (!all(y_int %in% c(0L, 1L))) {
      warning("outcome_var '", outcome_var,
              "' does not appear to be binary (0/1). Event rate may be unreliable.")
    }
    event_rate <- mean(y_int, na.rm = TRUE)
  }

  # Detect hierarchical structure from formula if present
  is_hierarchical <- FALSE  # Phase 1 only — formula not yet specified

  # Build alternatives via decision matrix (DESIGN.md §4.1)
  alternatives <- list()

  if (outcome_type == "binary" && n < 200 && !is_hierarchical) {
    rationale_base <- sprintf(
      "Binary outcome, n = %d (< 200)", n)
    er_note <- if (!is.na(event_rate) && event_rate < 0.05)
      sprintf(", rare event rate = %.1f%%", event_rate * 100) else ""

    alternatives <- list(
      list(method    = "logistic_bootstrap",
           label     = "Logistic regression + bootstrap CI",
           rationale = paste0(rationale_base, er_note,
                              ". Fast, interpretable, no MCMC required.")),
      list(method    = "firth_logistic",
           label     = "Firth penalized logistic regression",
           rationale = paste0(rationale_base, er_note,
                              ". Recommended when event rate < 5% to reduce separation bias.")),
      list(method    = "full_stan",
           label     = "Full Bayesian (Stan via brms)",
           rationale = paste0(rationale_base,
                              ". Propagates uncertainty fully; appropriate when priors are meaningful."))
    )

  } else if (outcome_type == "binary" && n >= 200 && !is_hierarchical) {
    alternatives <- list(
      list(method    = "logistic_bootstrap",
           label     = "Logistic regression + bootstrap CI",
           rationale = sprintf(
             "Binary outcome, n = %d. Adequate power for standard MLE; bootstrap CI is reliable.", n)),
      list(method    = "full_stan",
           label     = "Full Bayesian (Stan via brms)",
           rationale = sprintf(
             "Binary outcome, n = %d. Bayesian approach justified when priors carry substantive information.", n))
    )

  } else if (outcome_type == "count") {
    alternatives <- list(
      list(method    = "poisson_quasi",
           label     = "Poisson GLM + quasi-likelihood",
           rationale = "Count outcome. Quasi-likelihood corrects for overdispersion without MCMC."),
      list(method    = "negbinomial_glm",
           label     = "Negative binomial GLM",
           rationale = "Count outcome. Directly models overdispersion via NegBin family."),
      list(method    = "full_stan",
           label     = "Full Bayesian (Stan via brms)",
           rationale = "Count outcome. Bayesian approach allows full prior specification and uncertainty propagation.")
    )

  } else if (outcome_type == "continuous") {
    alternatives <- list(
      list(method    = "ttest_cohend",
           label     = "t-test + Cohen's d",
           rationale = "Continuous outcome, simple comparison goal. Efficient; well understood."),
      list(method    = "lm_bootstrap",
           label     = "Linear model + bootstrap CI",
           rationale = "Continuous outcome. Parametric flexibility with non-parametric CI."),
      list(method    = "full_stan",
           label     = "Full Bayesian (Stan via brms)",
           rationale = "Continuous outcome. Bayesian approach for full posterior inference.")
    )

  } else {
    # Hierarchical or unmatched — offer GLMM + Full Stan
    alternatives <- list(
      list(method    = "glmm_reml",
           label     = "GLMM with REML (lme4 / glmmTMB)",
           rationale = sprintf(
             "Hierarchical structure detected, n = %d. REML is efficient for variance component estimation.", n)),
      list(method    = "full_stan",
           label     = "Full Bayesian (Stan via brms)",
           rationale = "Hierarchical structure. Bayesian approach propagates group-level uncertainty fully.")
    )
    if (n < 50) {
      warning(sprintf(
        "Small n (%d) with hierarchical structure detected. Power for group-level variance estimation may be low.",
        n))
    }
  }

  list(
    n            = n,
    outcome_type = outcome_type,
    event_rate   = event_rate,
    goal         = goal,
    alternatives = alternatives
  )
}

# Pretty-print assess_offramps() result with equal visual weight
print_offramps <- function(offramp_result) {
  cat("=== Off-Ramp Assessment ===\n")
  cat("Outcome type:", offramp_result$outcome_type, "\n")
  cat("n =", offramp_result$n, "\n")
  if (!is.na(offramp_result$event_rate))
    cat("Event rate:", sprintf("%.1f%%", offramp_result$event_rate * 100), "\n")
  cat("Goal:", offramp_result$goal, "\n\n")
  cat("Available analysis paths (no preference implied):\n\n")
  for (i in seq_along(offramp_result$alternatives)) {
    alt <- offramp_result$alternatives[[i]]
    cat(sprintf("  [%d] %s\n", i, alt$label))
    cat(sprintf("      %s\n\n", alt$rationale))
  }
  invisible(offramp_result)
}
