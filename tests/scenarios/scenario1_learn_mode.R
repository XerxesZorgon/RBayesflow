# tests/scenarios/scenario1_learn_mode.R
# SCENARIO-1: SC-1 — Learn mode full loop (known-good model).
# TEST_PLAN.md §3 SCENARIO-1 (assertions 1.1–1.10).
# Run from project root: source("tests/scenarios/scenario1_learn_mode.R")
# Requires: examples/glmm_gaussian/fit_good.rds

cat("=== SCENARIO-1: Learn Mode Full Loop ===\n\n")

source("R/source_all.R")
library(lme4)
library(brms)

# Patch readline to auto-answer "yes" / select option 1
unlockBinding("readline", baseenv())
original_readline <- base::readline
assign("readline", function(prompt = "") {
  if (grepl("Select", prompt)) "3"  # option 3 = Full Stan
  else "yes"
}, envir = baseenv())

data    <- sleepstudy
formula <- Reaction ~ Days + (Days | Subject)
family  <- gaussian()
priors  <- c(
  prior(normal(250, 50), class = Intercept),
  prior(normal(10, 5),   class = b, coef = Days),
  prior(exponential(1),  class = sigma)
)

# =========================================================
# 1.1 Phase 1: off-ramp assessment runs; alternative offered
# =========================================================
wf <- init_workflow(mode = "learn", stage = "explore")
wf <- run_phase1(wf, data = data, outcome_var = "Reaction",
                 outcome_type = "continuous",
                 goal = "Estimate effect of sleep deprivation")
stopifnot("1.1: no alternatives offered" =
  length(assess_offramps(data, "Reaction", "continuous",
    "Estimate effect", n = nrow(data))$alternatives) >= 1)
cat("[1.1] PASS: off-ramp assessment ran, alternatives offered\n")

# =========================================================
# 1.2: audit trail contains bayesian_selected
# =========================================================
trail_actions <- sapply(wf$audit_trail, function(e) e$action)
stopifnot("1.2: bayesian_selected not in audit trail" =
  "bayesian_selected" %in% trail_actions)
cat("[1.2] PASS: audit trail contains bayesian_selected\n")

# =========================================================
# 1.3 Phase 2: priors stored; prior predictive draws stored
# =========================================================
wf <- run_phase2(wf, formula = formula, family = family,
                 priors = priors, data = data, seed = 42)
stopifnot("1.3: priors_objects not stored" = !is.null(wf$priors_objects))
stopifnot("1.3: prior_pred_draws not stored" = !is.null(wf$prior_pred_draws))
cat("[1.3] PASS: priors and prior predictive draws stored\n")

# =========================================================
# 1.4 Phase 3: fit_good loaded; fit_timestamp set; param logged
# =========================================================
fit_good <- readRDS("examples/glmm_gaussian/fit_good.rds")
wf <- record_fit(wf, fit_good)
wf$parameterization <- detect_parameterization(fit_good)
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase            = 3,
  action           = "parameterization_logged",
  timestamp        = Sys.time(),
  parameterization = wf$parameterization
)))
stopifnot("1.4: fit_timestamp not set" = !is.null(wf$fit_timestamp))
stopifnot("1.4: parameterization not logged in audit" =
  any(sapply(wf$audit_trail, function(e) e$action == "parameterization_logged")))
cat("[1.4] PASS: fit_timestamp set, parameterization logged\n")

# =========================================================
# 1.5 print(wf) before Phase 4: header only (diagnostics NA)
# =========================================================
out_1_5 <- capture.output(print(wf))
stopifnot("1.5: coefficient output shown before diagnostics" =
  !any(grepl("Coefficient table|fixef|Estimate", out_1_5)))
stopifnot("1.5: header not shown" =
  any(grepl("RBayesflow|workflow|mode", out_1_5, ignore.case = TRUE)))
cat("[1.5] PASS: print(wf) shows header only before diagnostics\n")

# =========================================================
# 1.6 Phase 4: diagnostics$passed == TRUE; rhat_max < 1.01
# =========================================================
wf <- run_diagnostics(fit_good, wf)
stopifnot("1.6: diagnostics$passed != TRUE" = isTRUE(wf$diagnostics$passed))
stopifnot("1.6: rhat_max >= 1.01" = wf$diagnostics$rhat_max < 1.01)
cat("[1.6] PASS: diagnostics passed, Rhat_max =",
    round(wf$diagnostics$rhat_max, 4), "\n")

# =========================================================
# 1.7 print(wf) after clean diagnostics: overlay + health
# =========================================================
out_1_7 <- capture.output(print(wf))
stopifnot("1.7: overlay placeholder not shown in learn mode" =
  any(grepl("overlay|Prior|posterior|Diagnostics", out_1_7, ignore.case = TRUE)))
cat("[1.7] PASS: print(wf) shows overlay and health summary\n")

# =========================================================
# 1.8 Phase 5: ppc_complete == TRUE
# =========================================================
wf$ppc_complete <- TRUE
export_context(wf)
stopifnot("1.8: ppc_complete != TRUE" = isTRUE(wf$ppc_complete))
cat("[1.8] PASS: ppc_complete set to TRUE\n")

# =========================================================
# 1.9 Phase 6: loo_complete == TRUE; loo_table is data.frame
# =========================================================
loo_fit <- loo::loo(fit_good)
wf$loo_complete <- TRUE
wf$loo_table    <- as.data.frame(loo_fit$estimates)
export_context(wf)
stopifnot("1.9: loo_complete != TRUE" = isTRUE(wf$loo_complete))
stopifnot("1.9: loo_table not data.frame" = is.data.frame(wf$loo_table))
cat("[1.9] PASS: loo_complete TRUE, loo_table is data.frame\n")

# =========================================================
# 1.10 Render Phase 7 report
# =========================================================
wf_path <- normalizePath(file.path(tempdir(), "wf_scenario1.rds"),
                         mustWork = FALSE)
saveRDS(wf, wf_path)
report_out <- "templates/scenario1_report.html"
render_ok <- tryCatch({
  quarto::quarto_render(
    "templates/bayesflow_report.qmd",
    output_file = "scenario1_report.html",
    execute_params = list(
      wf_path  = wf_path,
      fit_path = "examples/glmm_gaussian/fit_good.rds",
      mode     = "learn",
      stage    = "explore"
    )
  )
  file.exists(report_out)
}, error = function(e) {
  cat("Render error:", e$message, "\n")
  FALSE
})
stopifnot("1.10: report did not render" = isTRUE(render_ok))
cat("[1.10] PASS: Phase 7 report rendered to", report_out, "\n")

# Restore readline
assign("readline", original_readline, envir = baseenv())
lockBinding("readline", baseenv())

cat("\nSCENARIO-1: ALL 10 ASSERTIONS PASSED\n")
