# R/diagnose_cmdstan.R
#
# diagnose_cmdstan(): populate wf$diagnostics from a CmdStanFit object.
#
# Used when a case study bypasses brms and fits a hand-written Stan program
# directly with cmdstanr (e.g. Ch 19 disease-prevalence, Ch 20 mixture
# classification). The brms family registry is not dispatched; the display
# contract and diagnostic gate in wf_state operate normally.
#
# Promoted from inline copies in:
#   data/covid_ch19/covid_ch19_analysis.R
#   data/classification_ch20/classification_ch20_analysis.R
# after both chapters used the identical function. A third cmdstanr-native
# case study triggered promotion per the ADR note in those scripts.
#
# Arguments:
#   fit    -- a CmdStanFit object (from cmdstanr::cmdstan_model()$sample())
#   wf     -- a wf_state object (from new_wf_state() / init_workflow())
#   params -- character vector of sampled parameter names to compute
#             Rhat/ESS on. Should exclude generated quantities variables
#             (e.g. "p[n,k]", "log_lik") because those have no convergence
#             meaning. Defaults to NULL, which uses all variables present
#             in the draws; pass an explicit vector for large models.
#
# Returns: updated wf_state with diagnostics fields populated and
#   wf$diagnostics$passed set to TRUE or FALSE. Does NOT call export_context()
#   or trigger the interactive diagnose.wf_state() gate; callers do that.
#
# Thresholds used (matching run_diagnostics() in R/diagnostics.R):
#   Rhat_max       > 1.01  => fail
#   bulk_ESS_min   < 400   => fail
#   tail_ESS_min   < 400   => fail
#   n_divergences  > 0     => fail
#   any BFMI       < 0.3   => fail
#   max treedepth hit      => fail

diagnose_cmdstan <- function(fit, wf, params = NULL) {
  stopifnot(inherits(fit, "CmdStanFit"))
  stopifnot(inherits(wf,  "wf_state"))

  draws <- fit$draws(variables = params)
  smry  <- posterior::summarise_draws(
    draws,
    posterior::default_convergence_measures()
  )
  diag  <- fit$diagnostic_summary(quiet = TRUE)

  wf$fit_timestamp    <- Sys.time()
  wf$stan_backend     <- "cmdstanr"
  wf$parameterization <- NA_character_   # no brms random-effects term to detect

  wf$diagnostics$rhat_max          <- max(smry$rhat,     na.rm = TRUE)
  wf$diagnostics$bulk_ess_min      <- as.integer(min(smry$ess_bulk, na.rm = TRUE))
  tail_ess_raw <- min(smry$ess_tail, na.rm = TRUE)
  wf$diagnostics$tail_ess_min <- if (is.finite(tail_ess_raw)) as.integer(tail_ess_raw) else NA_integer_
  wf$diagnostics$n_divergences     <- as.integer(sum(diag$num_divergent))
  wf$diagnostics$bfmi              <- diag$ebfmi
  wf$diagnostics$max_treedepth_hit <- any(diag$num_max_treedepth > 0)
  wf$diagnostics$family_checks     <- list(
    note = paste0(
      "Hand-written Stan program via cmdstanr; ",
      "brms family registry bypassed. ",
      "Stan file: ", tryCatch(fit$metadata()$stan_file, error = function(e) "unknown")
    )
  )

  failed <- character()
  rhat_max     <- wf$diagnostics$rhat_max
  bulk_ess_min <- wf$diagnostics$bulk_ess_min
  tail_ess_min <- wf$diagnostics$tail_ess_min
  bfmi_vals    <- wf$diagnostics$bfmi

  if (!is.na(rhat_max)     && is.finite(rhat_max)     && rhat_max     > 1.01) failed <- c(failed, "rhat")
  if (!is.na(bulk_ess_min) && is.finite(bulk_ess_min) && bulk_ess_min < 400)  failed <- c(failed, "bulk_ess")
  if (!is.na(tail_ess_min) && is.finite(tail_ess_min) && tail_ess_min < 400)  failed <- c(failed, "tail_ess")
  if (is.na(rhat_max) || is.na(bulk_ess_min) || is.na(tail_ess_min))          failed <- c(failed, "ess_na")
  if (wf$diagnostics$n_divergences > 0)                                        failed <- c(failed, "divergences")
  if (length(bfmi_vals) > 0 && any(!is.na(bfmi_vals) & bfmi_vals < 0.3))     failed <- c(failed, "bfmi")
  if (isTRUE(wf$diagnostics$max_treedepth_hit))                                failed <- c(failed, "treedepth")

  wf$diagnostics$failed_criteria <- failed
  wf$diagnostics$passed          <- (length(failed) == 0L)

  wf$fit_hash <- digest::digest(
    list(
      stan_file    = tryCatch(fit$metadata()$stan_file, error = function(e) NA_character_),
      diag_summary = diag,
      timestamp    = wf$fit_timestamp
    ),
    algo = "sha256"
  )

  wf
}

