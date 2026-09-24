---
{
  "id": "file_f6o12ir3",
  "filetype": "document",
  "filename": "04-workflow-phases",
  "created_at": "2026-09-24T19:08:19.938Z",
  "updated_at": "2026-09-24T19:08:19.938Z",
  "meta": {
    "location": "/",
    "tags": [],
    "categories": [],
    "description": "",
    "source": "markdown"
  }
}
---
# RBayesflow — Workflow Phases

Work through these phases in order. Each phase has a corresponding Quarto template in `templates/`. Open the template for that phase, run the cells, and move on when the acceptance note says you are ready.

---

## Phase 1: Goal Declaration and Data Inspection

**Purpose:** Understand your data and decide whether a full Bayesian fit is warranted before committing to one.

**Run:** Open `templates/phase1_exploration.qmd`. Call `run_phase1(wf, data, outcome_var, outcome_type, goal)`.

**What to expect:** An exploratory plot of your outcome variable, followed by a numbered list of equally-weighted alternatives to a full Bayesian fit (the off-ramp assessment). Choose one; your choice is logged to `wf$audit_trail`.

**Move on when:** You have chosen a path — either a non-Bayesian alternative (which triggers `exit_workflow()`) or full Stan.

---

## Phase 2: Prior Specification and Prior Predictive Simulation

**Purpose:** Set priors and verify they produce plausible data before seeing the real data.

**Run:** Open `templates/phase2_priors.qmd`. Call `run_phase2(wf, formula, family, priors, data)`.

**What to expect:** A prior predictive density overlay plot showing draws from your priors mapped through the likelihood. Prior draws that produce impossible values (e.g. negative counts, probabilities outside 0–1) signal priors that need revision.

**Move on when:** Prior predictive draws are plausible given your domain knowledge.

---

## Phase 3: Model Fitting

**Purpose:** Fit the model and record the fit for reproducibility.

**Run:** Open `templates/phase3_fit.qmd`. Call `run_phase3(wf, formula, data, family, priors)`.

**What to expect:** CmdStan sampling progress, then a brief fit summary. The parameterization (centered or non-centered for hierarchical models) is logged to `wf$parameterization`.

**Move on when:** The fit completes without error and `wf$fit_timestamp` is set.

---

## Phase 4: MCMC Diagnostics

**Purpose:** Verify the sampler explored the posterior reliably before reading any results.

**Run:** Open `templates/phase4_diagnostics.qmd`. Call `run_diagnostics(fit, wf)`.

**What to expect:** Rhat values (target < 1.01), bulk and tail ESS (target > 400), divergence count (target 0), BFMI, and treedepth summary. Coefficient output is gated until this phase passes.

> **Important:** If diagnostics fail, you must call `wf <- wf$diagnose()` **at the R console**, not from inside the Quarto template. This is by design — the acknowledgment prompt uses `readline()`, which cannot run in a rendered document. The template contains a clearly marked cell reminding you of this.

**Move on when:** `wf$diagnostics$passed == TRUE`, or diagnostics have failed, you have called `wf$diagnose()`, and you have decided to proceed with caveats.

---

## Phase 5: Posterior Predictive Checks

**Purpose:** Verify the fitted model reproduces the structure of your observed data.

**Run:** Open `templates/phase5_ppc.qmd`. Call `brms::pp_check(fit)`.

**What to expect:** PPC density overlay and test statistic plots. In `mode = "learn"`, a Socratic prompt asks what discrepancies you see between the model and data. In `mode = "practice"`, standard PPC output is shown.

**Move on when:** `wf$ppc_complete == TRUE`.

---

## Phase 6: Model Comparison

**Purpose:** Compare candidate models using leave-one-out cross-validation (LOO-CV).

**Run:** Open `templates/phase6_loo.qmd`. Call `loo(fit)` and, if comparing models, `loo_compare(loo1, loo2)`.

**What to expect:** A LOO comparison table showing expected log-predictive density (ELPD) differences and standard errors between models.

**Move on when:** `wf$loo_complete == TRUE`.

---

## Phase 7: Reporting

**Purpose:** Render a single reproducible document in which every claim traces to a saved model object or `wf_state` field.

**Run:** Open `templates/bayesflow_report.qmd`. Render it with `quarto::quarto_render()` or the Render button in your IDE, supplying `wf_path` and `fit_path` as parameters.

**What to expect:** A complete analysis report. No models are re-fitted at render time — all summaries, plots, and tables come from the saved `wf_state` and `brmsfit` objects.

**Move on when:** The report renders without error and all claims are traceable to stored objects.
