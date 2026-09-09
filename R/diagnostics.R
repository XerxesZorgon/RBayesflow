# R/diagnostics.R
# run_diagnostics(): generic MCMC health checks + family registry dispatch.
# detect_parameterization(): centered vs non-centered detection.
# refit_noncentered(): one-call reparameterisation switch.
# DESIGN.md §6, §8; ADR-007.

run_diagnostics <- function(fit, wf) {
  stopifnot(inherits(fit,  "brmsfit"))
  stopifnot(inherits(wf,   "wf_state"))

  check_fit_hash(wf, fit)

  # --- Generic MCMC checks via posterior ---
  draws_summary <- posterior::summarise_draws(
    fit,
    "mean", "sd", "rhat",
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )

  rhat_max     <- max(draws_summary$rhat,     na.rm = TRUE)
  bulk_ess_min <- min(draws_summary$ess_bulk, na.rm = TRUE)
  tail_ess_min <- min(draws_summary$ess_tail, na.rm = TRUE)

  # --- Divergences, BFMI, treedepth from cmdstanr diagnostics ---
  diag_info <- tryCatch(fit$diagnostic_summary(quiet = TRUE),
                        error = function(e) NULL)

  n_divergences     <- if (!is.null(diag_info)) sum(diag_info$num_divergent)   else NA_integer_
  bfmi              <- if (!is.null(diag_info)) diag_info$ebfmi                else NA_real_
  max_treedepth_hit <- if (!is.null(diag_info)) any(diag_info$num_max_treedepth > 0) else NA

  # --- Evaluate each criterion ---
  failed <- character()
  if (!is.na(rhat_max)     && rhat_max     >= 1.01)   failed <- c(failed, "rhat_max")
  if (!is.na(bulk_ess_min) && bulk_ess_min <= 400)     failed <- c(failed, "bulk_ess_min")
  if (!is.na(tail_ess_min) && tail_ess_min <= 400)     failed <- c(failed, "tail_ess_min")
  if (!is.na(n_divergences) && n_divergences > 0)      failed <- c(failed, "n_divergences")
  if (!all(is.na(bfmi))    && any(bfmi < 0.3, na.rm = TRUE)) failed <- c(failed, "bfmi")
  if (isTRUE(max_treedepth_hit))                       failed <- c(failed, "max_treedepth")

  passed <- length(failed) == 0

  # --- Store generic results in wf ---
  wf$diagnostics$rhat_max          <- rhat_max
  wf$diagnostics$bulk_ess_min      <- bulk_ess_min
  wf$diagnostics$tail_ess_min      <- tail_ess_min
  wf$diagnostics$n_divergences     <- n_divergences
  wf$diagnostics$bfmi              <- bfmi
  wf$diagnostics$max_treedepth_hit <- max_treedepth_hit
  wf$diagnostics$failed_criteria   <- failed
  wf$diagnostics$passed            <- passed

  # --- Family-specific registry dispatch ---
  fkey         <- family_key(fit)
  registry_fn  <- DIAGNOSTIC_REGISTRY[[fkey]]
  family_result <- tryCatch(
    registry_fn(fit, wf),
    error = function(e) {
      warning("Family-specific diagnostic failed for '", fkey, "': ", e$message)
      list(checks = list(), warnings = character(), plots = list())
    }
  )
  wf$diagnostics$family_checks <- family_result

  # --- Append divergence + parameterization note to warnings if needed ---
  if (!is.na(n_divergences) && n_divergences > 0 &&
      !is.null(wf$parameterization) &&
      identical(wf$parameterization, "centered")) {
    p_warn <- "Centered parameterization with divergences detected. Consider calling refit_noncentered(wf)."
    wf$diagnostics$family_checks$warnings <-
      c(wf$diagnostics$family_checks$warnings, p_warn)
  }

  # --- Audit trail ---
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase           = 4,
    action          = "diagnostics_run",
    timestamp       = Sys.time(),
    passed          = passed,
    failed_criteria = failed,
    family_key      = fkey
  )))

  export_context(wf)
  invisible(wf)
}

# --- detect_parameterization() — DESIGN.md §6 ---

detect_parameterization <- function(fit) {
  stopifnot(inherits(fit, "brmsfit"))

  # Only relevant for hierarchical models
  re_terms <- reformulas::findbars(formula(fit))
  if (is.null(re_terms) || length(re_terms) == 0) {
    return(NA_character_)
  }

  # Inspect generated Stan code for parameterization indicators
  stan_code <- tryCatch(
    brms::stancode(fit),
    error = function(e) ""
  )

  # Non-centered: z_ prefix variables present (e.g., z_1, z_2)
  # Centered: r_ prefix variables present (e.g., r_Subject)
  has_z_prefix <- grepl("\\bz_[0-9]", stan_code)
  has_r_prefix <- grepl("\\br_[A-Za-z]", stan_code)

  param <- if (has_z_prefix && !has_r_prefix) {
    "non-centered"
  } else if (has_r_prefix && !has_z_prefix) {
    "centered"
  } else if (has_z_prefix && has_r_prefix) {
    # Mixed — report non-centered (dominant pattern)
    "non-centered"
  } else {
    # Cannot determine — fall back gracefully
    NA_character_
  }

  param
}

# --- refit_noncentered() — DESIGN.md §6 ---
# Switches a centered-parameterization hierarchical model to non-centered
# by rebuilding the brms formula with (0 + Intercept | group) syntax.
# v1.0: handles simple (1 | group) terms only.
# Complex random-effects structures require manual formula revision.

refit_noncentered <- function(wf, fit) {
  stopifnot(inherits(wf,  "wf_state"))
  stopifnot(inherits(fit, "brmsfit"))

  if (!identical(wf$parameterization, "centered")) {
    stop(
      "refit_noncentered() requires wf$parameterization == 'centered'.\n",
      "Current parameterization: ", wf$parameterization
    )
  }

  # Rebuild formula: replace (1 | group) with (0 + Intercept | group)
  original_formula <- deparse(formula(fit))
  nc_formula_str   <- gsub(
    pattern     = "\\(1\\s*\\|\\s*([^)]+)\\)",
    replacement = "(0 + Intercept | \\1)",
    x           = original_formula
  )

  if (identical(nc_formula_str, original_formula)) {
    warning(
      "No simple (1 | group) terms found in formula. ",
      "The formula was not modified. ",
      "For complex random-effects structures, revise the formula manually."
    )
  }

  nc_formula <- as.formula(nc_formula_str)

  # Refit with non-centered parameterization
  message("Refitting with non-centered parameterization: ", nc_formula_str)
  fit_nc <- brms::brm(
    formula  = nc_formula,
    data     = fit$data,
    family   = brms::family(fit),
    prior    = brms::prior_summary(fit),
    backend  = "cmdstanr",
    seed     = 42,
    refresh  = 100
  )

  # Update wf with new fit
  wf_nc                  <- wf
  wf_nc$parameterization <- "non-centered"
  wf_nc                  <- record_fit(wf_nc, fit_nc)
  wf_nc$audit_trail      <- append(wf_nc$audit_trail, list(list(
    phase     = 4,
    action    = "refit_noncentered",
    timestamp = Sys.time(),
    formula   = nc_formula_str
  )))

  list(fit_new = fit_nc, wf_new = wf_nc)
}


