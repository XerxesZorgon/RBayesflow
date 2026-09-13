# R/phase2_priors.R
# Phase 2: prior specification and prior predictive simulation.
# Called interactively from templates/phase2_priors.qmd.
# DESIGN.md §§2, 5.

run_phase2 <- function(wf, formula, family, priors, data, seed = 42) {
  stopifnot(inherits(wf, "wf_state"))
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(data))
  stopifnot(is.list(priors) || inherits(priors, "brmsprior"))

  cat("=== Phase 2: Prior Specification ===\n")

  # --- Store priors on wf ---
  wf$priors_objects <- if (inherits(priors, "brmsprior")) list(priors) else priors
  wf$priors_text    <- paste(deparse(priors), collapse = "\n")
  wf$formula        <- formula
  wf$family         <- family

  cat("Formula:", deparse(formula), "\n")
  cat("Family: ", family$family, "(", family$link, ")\n")
  cat("Priors:\n")
  prior_df <- if (inherits(priors, "brmsprior")) priors else do.call(c, priors)
  for (i in seq_len(nrow(prior_df))) {
    row <- prior_df[i, ]
    label <- if (nchar(trimws(row$coef)) > 0)
      paste0(row$class, "[", row$coef, "]")
    else
      row$class
    cat(" ", label, "~", row$prior, "\n")
  }
  cat("\n")

  # --- Prior predictive simulation ---
  cat("Running prior predictive simulation...\n")
  fit_prior <- brms::brm(
    formula       = formula,
    data          = data,
    family        = family,
    prior         = if (inherits(priors, "brmsprior")) priors
                    else do.call(c, priors),
    sample_prior  = "only",
    backend       = "cmdstanr",
    seed          = seed,
    chains        = 2,
    iter          = 1000,
    warmup        = 500,
    refresh       = 0
  )

  # Store prior predictive draws
  wf$prior_pred_draws <- brms::posterior_predict(fit_prior, ndraws = 200)

  cat("Prior predictive draws stored (", nrow(wf$prior_pred_draws),
      " draws x", ncol(wf$prior_pred_draws), " observations).\n\n")

  # --- Mode-specific plot ---
  if (wf$mode == "learn") {
    y_vals <- as.numeric(data[[as.character(formula[[2]])]])
    y_lo   <- round(min(y_vals) - sd(y_vals, na.rm = TRUE), 1)
    y_hi   <- round(max(y_vals) + sd(y_vals, na.rm = TRUE), 1)
    cat("[Learn mode — Prior Predictive Check]\n")
    cat("The plot below shows 50 simulated datasets drawn from your priors,\n")
    cat("overlaid on the observed data distribution (dark line).\n")
    cat("Ask yourself: do the prior draws cover the plausible range of outcomes?\n")
    cat(sprintf("  Observed outcome range: %.1f to %.1f\n", min(y_vals), max(y_vals)))
    cat(sprintf("  Prior draws should stay roughly within: %.1f to %.1f\n", y_lo, y_hi))
    cat("If prior draws are wildly outside this range, tighten your priors before fitting.\n\n")
    outcome_name <- as.character(formula[[2]])
    p <- bayesplot::ppc_dens_overlay(
      y    = as.numeric(data[[outcome_name]]),
      yrep = wf$prior_pred_draws[seq_len(min(50, nrow(wf$prior_pred_draws))), ]
    ) +
      ggplot2::labs(x = outcome_name) +
      ggplot2::coord_cartesian(xlim = c(y_lo - 1, y_hi + 1))
    print(p)
  }

  # --- Audit trail ---
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase     = 2,
    action    = "priors_set",
    timestamp = Sys.time(),
    n_priors  = length(wf$priors_objects),
    seed      = seed
  )))

  export_context(wf)
  invisible(wf)
}
