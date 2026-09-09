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
  cat("Family: ", deparse(family), "\n")
  cat("Priors:\n")
  cat(wf$priors_text, "\n\n")

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
    cat("[Phase 2 learn mode: prior predictive distribution plot]\n")
    p <- bayesplot::ppc_dens_overlay(
      y   = as.numeric(data[[as.character(formula[[2]])]]),
      yrep = wf$prior_pred_draws[seq_len(min(50, nrow(wf$prior_pred_draws))), ]
    )
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
