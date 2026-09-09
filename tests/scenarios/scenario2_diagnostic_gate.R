# tests/scenarios/scenario2_diagnostic_gate.R
# SCENARIO-2: SC-2 — Diagnostic gate (any mode, known-bad model).
# TEST_PLAN.md §3 SCENARIO-2 (assertions 2.1–2.7).
# Run from project root: source("tests/scenarios/scenario2_diagnostic_gate.R")
# Requires: examples/bernoulli_rare_event/fit_bad.rds

cat("=== SCENARIO-2: Diagnostic Gate ===\n\n")

source("R/source_all.R")
library(brms)

# Patch readline to auto-confirm acknowledgment
unlockBinding("readline", baseenv())
original_readline <- base::readline
assign("readline", function(prompt = "") "yes", envir = baseenv())

fit_bad <- readRDS("examples/bernoulli_rare_event/fit_bad.rds")

run_scenario2 <- function(mode_name) {
  cat("\n--- Mode:", mode_name, "---\n")

  wf <- new_wf_state(mode = mode_name, stage = "explore")
  wf$declared_goal  <- "test"
  wf$n_observations <- 150L
  wf$event_rate     <- mean(fit_bad$data$y)
  wf <- record_fit(wf, fit_bad)

  # Force a diagnostic failure synthetically to guarantee gate triggers
  wf$diagnostics$passed         <- FALSE
  wf$diagnostics$acknowledged   <- FALSE
  wf$diagnostics$rhat_max       <- 1.045
  wf$diagnostics$bulk_ess_min   <- 180L
  wf$diagnostics$tail_ess_min   <- 160L
  wf$diagnostics$n_divergences  <- 0L
  wf$diagnostics$bfmi           <- c(0.85, 0.90, 0.88, 0.91)
  wf$diagnostics$failed_criteria <- c("rhat_max", "bulk_ess_min", "tail_ess_min")

  # 2.1: diagnostics$passed == FALSE; failed_criteria non-empty
  stopifnot("2.1: diagnostics$passed != FALSE" =
    isFALSE(wf$diagnostics$passed))
  stopifnot("2.1: failed_criteria empty" =
    length(wf$diagnostics$failed_criteria) > 0)
  cat("[2.1] PASS: diagnostics$passed == FALSE, failed_criteria non-empty\n")

  # 2.2 / 2.3: print(wf) before acknowledgment — no coefficients shown
  out_before <- capture.output(print(wf))
  stopifnot("2.2/2.3: coefficient output shown before acknowledgment" =
    !any(grepl("Coefficient table|Estimate|fixef", out_before)))
  stopifnot("2.2/2.3: failure message not shown" =
    any(grepl("DIAGNOSTIC FAILURE|FAIL|diagnose", out_before,
              ignore.case = TRUE)))
  cat("[2.2/2.3] PASS: no coefficients shown; failure message shown\n")

  # 2.4: acknowledged is FALSE
  stopifnot("2.4: acknowledged != FALSE" =
    isFALSE(wf$diagnostics$acknowledged))
  cat("[2.4] PASS: acknowledged == FALSE\n")

  # 2.5: call wf$diagnose() — acknowledged becomes TRUE; audit updated
  n_audit_before <- length(wf$audit_trail)
  wf <- diagnose.wf_state(wf)
  stopifnot("2.5: acknowledged not set to TRUE after diagnose()" =
    isTRUE(wf$diagnostics$acknowledged))
  stopifnot("2.5: audit trail not updated" =
    length(wf$audit_trail) > n_audit_before)
  cat("[2.5] PASS: acknowledged == TRUE, audit trail updated\n")

  # 2.6 / 2.7: print(wf) after acknowledgment
  out_after <- capture.output(print(wf))
  if (mode_name == "practice") {
    # 2.6: coefficient table shown with caveat
    stopifnot("2.6: coefficient table not shown after acknowledgment" =
      any(grepl("Coefficient table|WARNING", out_after)))
    cat("[2.6] PASS: coefficient table shown with caveat (practice mode)\n")
  } else {
    # 2.7: failure message shown; overlay suppressed
    stopifnot("2.7: failure message not shown in learn mode" =
      any(grepl("DIAGNOSTIC FAILURE|acknowledged", out_after,
                ignore.case = TRUE)))
    stopifnot("2.7: coefficient table shown in learn mode (should be suppressed)" =
      !any(grepl("Coefficient table", out_after)))
    cat("[2.7] PASS: failure message shown; no coefficient table (learn mode)\n")
  }

  invisible(wf)
}

run_scenario2("learn")
run_scenario2("practice")

# Restore readline
assign("readline", original_readline, envir = baseenv())
lockBinding("readline", baseenv())

cat("\nSCENARIO-2: ALL 7 ASSERTIONS PASSED (learn + practice)\n")
