# tests/scenarios/check_templates.R
# M-4 gate: validate all six phase templates with quarto_inspect.
# Run from project root: source("tests/scenarios/check_templates.R")
# Requires: QUARTO_PATH set in .Renviron (set during Task 039).

cat("=== M-4 Template Inspection Check ===\n\n")

templates <- c(
  "templates/phase1_exploration.qmd",
  "templates/phase2_priors.qmd",
  "templates/phase3_fit.qmd",
  "templates/phase4_diagnostics.qmd",
  "templates/phase5_ppc.qmd",
  "templates/phase6_loo.qmd"
)

results <- logical(length(templates))
for (i in seq_along(templates)) {
  tpl <- templates[i]
  ok <- tryCatch({
    quarto::quarto_inspect(tpl)
    TRUE
  }, error = function(e) {
    cat("FAIL:", tpl, "\n     Error:", e$message, "\n")
    FALSE
  })
  results[i] <- ok
  if (ok) cat("PASS:", tpl, "\n")
}

cat("\n--- Summary ---\n")
cat("Passed:", sum(results), "/", length(results), "\n")
cat("Failed:", sum(!results), "/", length(results), "\n")

# Report template (bayesflow_report.qmd) not yet written — skip with note
cat("\nnote: bayesflow_report.qmd deferred to M-5 (Task 049).\n")

stopifnot("One or more phase templates failed quarto_inspect" = all(results))
cat("\nM-4 TEMPLATE CHECK: ALL 6 PHASE TEMPLATES PASSED\n")
