# RBayesflow — Project Plan (PLAN.md)

**Status:** Draft  
**Last updated:** 2026-09-08  
**Owner:** John Peach  
**Spec Kit role:** The **WHEN** document — phased implementation strategy, milestones, and sequencing constraints. Feeds `tasks.md` (atomic work items). Each phase produces a testable deliverable.

---

## Dependency Order

```
SDD → DESIGN → PLAN → tasks.md → implementation
ADR-001 through ADR-007 (all accepted before implementation begins)
```

ODQ-1 and ODQ-2 (Open Design Questions in DESIGN.md) must be resolved as ADRs before the phases that depend on them begin. ODQ-1 blocks Phase 4 template design. ODQ-2 blocks Phase 3 script design.

---

## Phase 0: Environment Setup and Repository

**ID:** PHASE-0 | **Duration:** ~1 week | **Dependencies:** None

### Description
Establish the project folder structure, pin all dependencies, and verify that Stan is operational. This phase produces a runnable R session that can install all required libraries. No workflow logic is implemented.

### Milestones
- [ ] `RBayesflow/` folder initialized with the directory structure defined in DESIGN.md §1
- [ ] `DESCRIPTION` file listing all required packages (brms, cmdstanr, bayesplot, posterior, tidybayes, loo, ggplot2, esquisse, jsonlite, yaml, testthat, withr)
- [ ] `renv.lock` generated with pinned versions of all dependencies
- [ ] `README.md` contains a working "Getting Started" section
- [ ] Stan compiles a simple test model (`stan_demo.R`) without error
- [ ] All seven ADRs finalized and accepted

### Open Decisions Required Before This Phase
- ODQ-4 (ADR candidate): Does `mode = "expert"` appear as a stub in v1.0 or not?

---

## Phase 1: Core Infrastructure — `wf_state` and Display Contract

**ID:** PHASE-1 | **Duration:** ~2 weeks | **Dependencies:** PHASE-0 complete; ODQ-2 resolved

### Description
Implement the foundational S3 object, its constructor, and the display contract. This is the most architecturally critical phase — every subsequent phase builds on the `wf_state` schema and `print.wf_state()` correctness.

### Milestones
- [ ] `R/wf_state.R`: `new_wf_state()`, `print.wf_state()`, `summary.wf_state()`, `format.wf_state()`
- [ ] `R/context.R`: `export_context()` producing a valid `wf_context.json`
- [ ] `R/init.R`: `init_workflow(mode, stage)` — constructs initial `wf_state`, writes context
- [ ] Unit tests UT-1, UT-5 pass (no Stan required)
- [ ] Display contract unit tests UT-2 pass for all 4 branches (using mock fit objects)

### Acceptance Criterion
`print(wf)` on a `wf_state` object with `diagnostics$passed = FALSE` and `acknowledged = FALSE` shows the failure message and `wf$diagnose()` prompt, and **no coefficient output**. Verified by UT-2.

---

## Phase 2: Diagnostic Infrastructure

**ID:** PHASE-2 | **Duration:** ~2 weeks | **Dependencies:** PHASE-1 complete; ODQ-1 resolved

### Description
Implement the diagnostic registry, the generic MCMC diagnostic runner, the `diagnose.wf_state()` method, and the four built-in family-specific registry entries.

### Milestones
- [ ] `R/diagnostic_registry.R`: `DIAGNOSTIC_REGISTRY` with entries for bernoulli, poisson, negbinomial, gaussian-hierarchical, time-series, unknown
- [ ] `R/diagnostics.R`: `run_diagnostics(fit, wf_state)` — generic Rhat/ESS/divergences/BFMI/treedepth checks + registry dispatch
- [ ] `diagnose.wf_state()` method: displays bayesplot panels, prompts acknowledgment, updates `wf_state$diagnostics$acknowledged`, appends to audit trail
- [ ] `detect_parameterization(fit)` helper: returns "centered" or "non-centered"
- [ ] `refit_noncentered(wf)` one-call switch (stub for Phase 4 use)
- [ ] Unit tests UT-3 (smoke), UT-4 pass

### Acceptance Criterion
`run_diagnostics(fit_good, wf)` returns `wf_state$diagnostics$passed == TRUE` for the known-good model. `run_diagnostics(fit_bad, wf)` returns `passed == FALSE` with `failed_criteria` non-empty and the Bernoulli rare-event warning present.

---

## Phase 3: Phase Scripts — Fitting and Off-Ramps

**ID:** PHASE-3 | **Duration:** ~2 weeks | **Dependencies:** PHASE-2 complete; ODQ-2 resolved; ODQ-3 clarified

### Description
Implement the R scripts for Phases 1–3 of the Bayesian workflow (goal declaration, prior specification, fitting). Implement the off-ramp assessment and exit workflow.

### Milestones
- [ ] `R/offramps.R`: `assess_offramps(outcome_type, n, event_rate, goal)` — decision matrix from DESIGN.md §4
- [ ] `R/exit_workflow.R`: `exit_workflow(wf, method, alternatives)` — writes YAML exit log
- [ ] Phase 1 script: data inspection, `assess_offramps()` call, choice logging
- [ ] Phase 2 script: `set_prior()` / `prior()` interface; prior predictive simulation; stores in `wf_state`
- [ ] Phase 3 script: `brm()` call with cmdstanr backend; `detect_parameterization()`; stores results in `wf_state`; calls `export_context()`
- [ ] Unit tests UT-6 (off-ramp contract), UT-3 (exit log schema) pass

### Acceptance Criterion
Running Phase 1 on `dat_bad` (binary, n=150, event rate 3%) produces an off-ramp list with ≥ 2 alternatives including "Firth penalized logistic" and "full_stan". Running `exit_workflow()` produces a valid YAML file. SCENARIO-3 assertions 3.1–3.6 pass.

---

## Phase 4: Quarto Templates — Phases 1–6

**ID:** PHASE-4 | **Duration:** ~3 weeks | **Dependencies:** PHASE-3 complete; ODQ-1 resolved

### Description
Implement all six phase Quarto templates (`.qmd` files). Each template sources the appropriate Phase script, uses the Posit Assistant context block, and renders correctly in both `mode = "learn"` and `mode = "practice"`.

### Milestones
- [ ] `templates/phase1_exploration.qmd` — data inspection, off-ramp table, choice cell
- [ ] `templates/phase2_priors.qmd` — prior specification UI, prior predictive plot cell
- [ ] `templates/phase3_fit.qmd` — `brm()` call cell, parameterization note
- [ ] `templates/phase4_diagnostics.qmd` — `run_diagnostics()` call; mode-aware display (ODQ-1 resolution applied)
- [ ] `templates/phase5_ppc.qmd` — bayesplot PPC panels; mode-aware plot selection
- [ ] `templates/phase6_loo.qmd` — `loo()` + `loo_compare()` calls; LOO table
- [ ] Each template renders without error in both `learn` and `practice` modes

### Acceptance Criterion
Running Phase 4 template on `fit_good` with `mode = "learn"`: renders without error; shows prior-vs-posterior overlay in the output; no coefficient table visible. Running with `mode = "practice"` and clean diagnostics: shows coefficient table.

---

## Phase 5: Phase 7 Report Template and Integration

**ID:** PHASE-5 | **Duration:** ~2 weeks | **Dependencies:** PHASE-4 complete

### Description
Implement the final Quarto report template that integrates all phases into a single reproducible document. Every claim in the report must be traceable to a stored model object or `wf_state` field.

### Milestones
- [ ] `templates/bayesflow_report.qmd` — parameterized template (`params: [mode, stage, model_path]`) that `readRDS(params$model_path)` and renders all summaries from `wf_state`
- [ ] Report renders from a saved `brmsfit` + `wf_state` without re-fitting
- [ ] All diagnostic summaries, coefficient tables, PPC plots, and LOO tables in the report are sourced from `wf_state` fields, not re-computed
- [ ] SCENARIO-1 step 1.10 passes

### Acceptance Criterion
SCENARIO-1 passes end-to-end (all 10 assertions).

---

## Phase 6: Acceptance Testing and Documentation Polish

**ID:** PHASE-6 | **Duration:** ~1 week | **Dependencies:** PHASE-5 complete

### Description
Run all three success-criterion scenarios, close any failures, and polish the user-facing documentation.

### Milestones
- [ ] SCENARIO-1 (SC-1) passes in `mode = "learn"`, `stage = "explore"`
- [ ] SCENARIO-2 (SC-2) passes in both `mode = "learn"` and `mode = "practice"`
- [ ] SCENARIO-3 (SC-3) passes (Phase 1 exit and Phase 4 exit)
- [ ] SCENARIO-4 passes (off-ramp equal weighting)
- [ ] SCENARIO-5 passes (parameterization detection)
- [ ] All unit tests (UT-1 through UT-6) pass
- [ ] README "Getting Started" section verified by a test user (student persona)
- [ ] CHANGELOG.md updated with v0.1.0 entry
- [ ] All ODQs resolved or explicitly deferred with rationale

### Acceptance Criterion
All three success criteria (SC-1, SC-2, SC-3) from SDD §6 are verified. v0.1.0 tag is applied.

---

## Summary Timeline

| Phase | Duration | Cumulative | Key Deliverable |
|---|---|---|---|
| PHASE-0 | ~1 week | Week 1 | Repo + pinned dependencies + Stan smoke test |
| PHASE-1 | ~2 weeks | Week 3 | `wf_state` S3 object + display contract |
| PHASE-2 | ~2 weeks | Week 5 | Diagnostic registry + `run_diagnostics()` + `diagnose()` |
| PHASE-3 | ~2 weeks | Week 7 | Phase scripts 1–3 + off-ramps + exit workflow |
| PHASE-4 | ~3 weeks | Week 10 | All six Quarto phase templates |
| PHASE-5 | ~2 weeks | Week 12 | Report template + SCENARIO-1 end-to-end |
| PHASE-6 | ~1 week | Week 13 | All scenarios passing; v0.1.0 |

**Total estimated duration:** ~13 weeks (part-time, single developer).

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| ODQ-1 (diagnose acknowledgment mechanism) blocks Phase 4 | High | High | Resolve as ADR before PHASE-3 ends; `readline()` is acceptable for v1.0 as a pragmatic choice |
| Posit Assistant API changes before v1.0 | Medium | Medium | File-based `wf_context.json` (Option Y) is the primary path; not dependent on any private API |
| brms version bump changes `stancode()` output format | Low | Medium | Pin brms version in `renv.lock`; test `detect_parameterization()` against pinned version |
| Stan compilation fails on target platform | Low | High | Smoke-test in PHASE-0; document `cmdstanr::install_cmdstan()` in README |

---

*End of PLAN.md*
