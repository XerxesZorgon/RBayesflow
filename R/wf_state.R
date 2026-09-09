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

# --- summary and format methods ---

summary.wf_state <- function(object, ...) {
  wf <- object
  cat("=== RBayesflow Workflow Summary ===\n")
  cat("Mode: ", wf$mode, "\n", sep = "")
  cat("Stage:", wf$stage, "\n")
  cat("Version:", wf$rbayesflow_version, "\n")
  cat("Fit timestamp:",
      if (!is.null(wf$fit_timestamp)) format(wf$fit_timestamp) else "not yet fitted",
      "\n")
  diag_status <- if (is.na(wf$diagnostics$passed)) {
    "not yet run"
  } else if (isTRUE(wf$diagnostics$passed)) {
    "PASSED"
  } else {
    paste0("FAILED (", paste(wf$diagnostics$failed_criteria, collapse = ", "), ")")
  }
  cat("Diagnostics:", diag_status, "\n")
  cat("Acknowledged:", wf$diagnostics$acknowledged, "\n")
  cat("Audit trail entries:", length(wf$audit_trail), "\n")
  cat("PPC complete:", wf$ppc_complete, "\n")
  cat("LOO complete:", wf$loo_complete, "\n")
  invisible(wf)
}

format.wf_state <- function(x, ...) {
  wf <- x
  diag_status <- if (is.na(wf$diagnostics$passed)) {
    "diagnostics_not_run"
  } else if (isTRUE(wf$diagnostics$passed)) {
    "diagnostics_passed"
  } else {
    "diagnostics_failed"
  }
  paste0(
    "RBayesflow v", wf$rbayesflow_version,
    " | mode=", wf$mode,
    " | stage=", wf$stage,
    " | ", diag_status,
    " | audit_entries=", length(wf$audit_trail)
  )
}

# --- record_fit() and check_fit_hash() — ADR-009 hash linkage ---

record_fit <- function(wf, fit) {
  stopifnot(inherits(wf, "wf_state"))
  stopifnot(inherits(fit, "brmsfit"))

  ts <- Sys.time()

  # Compute SHA-256 hash linking this wf to this specific fit
  wf$fit_hash <- digest::digest(
    list(
      formula   = as.character(formula(fit)),
      data_hash = digest::digest(fit$data, algo = "sha256"),
      timestamp = ts
    ),
    algo = "sha256"
  )
  wf$fit_timestamp <- ts
  wf$stan_backend  <- "cmdstanr"

  # Append to audit trail
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase     = 3,
    action    = "fit_recorded",
    timestamp = wf$fit_timestamp,
    fit_hash  = wf$fit_hash
  )))

  export_context(wf)
  invisible(wf)
}

check_fit_hash <- function(wf, fit) {
  stopifnot(inherits(wf, "wf_state"))
  stopifnot(inherits(fit, "brmsfit"))

  if (is.null(wf$fit_hash)) {
    stop("wf has no recorded fit. Call record_fit(wf, fit) first.")
  }

  current_hash <- digest::digest(
    list(
      formula   = as.character(formula(fit)),
      data_hash = digest::digest(fit$data, algo = "sha256"),
      timestamp = wf$fit_timestamp
    ),
    algo = "sha256"
  )

  if (!identical(current_hash, wf$fit_hash)) {
    stop(
      "Fit hash mismatch: the supplied fit object does not match the ",
      "model recorded in this wf_state.\n",
      "Did you pass the wrong fit object, or re-run brm() without ",
      "calling record_fit() again?"
    )
  }

  invisible(TRUE)
}

