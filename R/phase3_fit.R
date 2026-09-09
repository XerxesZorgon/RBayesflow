# R/phase3_fit.R
# Phase 3: model fitting via brms/cmdstanr.
# Called interactively from templates/phase3_fit.qmd.
# DESIGN.md §§2, 6; ADR-002 (cmdstanr backend); ADR-009 (hash linkage).

run_phase3 <- function(wf, formula, data, family, priors,
                       seed      = 42,
                       chains    = 4,
                       iter      = 2000,
                       warmup    = 1000,
                       refresh   = 100) {
  stopifnot(inherits(wf, "wf_state"))
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(data))

  cat("=== Phase 3: Model Fitting ===\n")
  cat("Formula:", deparse(formula), "\n")
  cat("Family: ", deparse(family), "\n")
  cat("Backend: cmdstanr | seed:", seed, "| chains:", chains,
      "| iter:", iter, "| warmup:", warmup, "\n\n")
  cat("Fitting model — this may take several minutes...\n")

  fit <- brms::brm(
    formula = formula,
    data    = data,
    family  = family,
    prior   = if (inherits(priors, "brmsprior")) priors
               else do.call(c, priors),
    backend = "cmdstanr",
    seed    = seed,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    refresh = refresh
  )

  cat("\nFit complete.\n")

  # --- Hash linkage (ADR-009) ---
  wf <- record_fit(wf, fit)

  # --- Parameterization detection (DESIGN.md §6) ---
  wf$parameterization <- detect_parameterization(fit)
  param_note <- if (is.na(wf$parameterization)) {
    "Non-hierarchical model; parameterization not applicable."
  } else {
    paste0("Parameterization: ", wf$parameterization)
  }
  cat(param_note, "\n")

  # --- Warn if centered + hierarchical ---
  re_terms <- reformulas::findbars(formula)
  if (identical(wf$parameterization, "centered") &&
      !is.null(re_terms) && length(re_terms) > 0) {
    cat("Note: Centered parameterization detected for hierarchical model.\n")
    cat("If divergences occur, call refit_noncentered(wf, fit).\n")
  }

  # --- Audit trail ---
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase            = 3,
    action           = "parameterization_logged",
    timestamp        = Sys.time(),
    parameterization = wf$parameterization,
    seed             = seed,
    chains           = chains,
    iter             = iter,
    warmup           = warmup
  )))

  export_context(wf)

  cat("\nSave your fit with: saveRDS(fit, 'fit.rds')\n")
  list(fit = fit, wf = wf)
}
