# RBayesflow — Project Plan (PLAN.md)

**Status:** Active
**Last updated:** 2026-09-10
**Owner:** John Peach
**Spec Kit role:** The **WHEN** document — phased implementation strategy, milestones, and sequencing constraints. Feeds `tasks.md` (atomic work items). Each phase produces a testable deliverable.

---

## Dependency Order

```
SDD → DESIGN → PLAN → tasks.md → implementation
ADR-001 through ADR-011 (all accepted before v0.1.0 implementation)
ADR-012 accepted before PHASE-7 implementation [v0.2.0]
```

ODQ-1 and ODQ-2 (Open Design Questions in DESIGN.md) were resolved as ADRs before their dependent phases began. ODQ-1 blocked Phase 4 template design; ODQ-2 blocked Phase 3 script design. Both are settled in v0.1.0.

---

## Phase 0: Environment Setup and Repository

**ID:** PHASE-0 | **Duration:** ~1 week | **Dependencies:** None

### Description
Establish the project folder structure, pin all dependencies, and verify that Stan is operational. This phase produces a runnable R session that can install all required libraries. No workflow logic is implemented.

### Milestones
- [x] `RBayesflow/` folder initialized with the directory structure defined in DESIGN.md §1
- [x] `DESCRIPTION` file listing all required packages (brms, cmdstanr, bayesplot, posterior, tidybayes, loo, ggplot2, esquisse, jsonlite, yaml, testthat, withr)
- [x] `renv.lock` generated with pinned versions of all dependencies
- [x] `README.md` contains a working "Getting Started" section
- [x] Stan compiles a simple test model (`stan_demo.R`) without error
- [x] All seven ADRs finalized and accepted

### Open Decisions Required Before This Phase
- ODQ-4 (ADR candidate): Does `mode = "expert"` appear as a stub in v1.0 or not? — resolved.

---

## Phase 1: Core Infrastructure — `wf_state` and Display Contract

**ID:** PHASE-1 | **Duration:** ~2 weeks | **Dependencies:** PHASE-0 complete; ODQ-2 resolved

### Description
Implement the foundational S3 object, its constructor, and the display contract. This is the most architecturally critical phase — every subsequent phase builds on the `wf_state` schema and `print.wf_state()` correctness.

### Milestones
- [x] `R/wf_state.R`: `new_wf_state()`, `print.wf_state()`, `summary.wf_state()`, `format.wf_state()`
- [x] `R/context.R`: `export_context()` producing a valid `wf_context.json`
- [x] `R/init.R`: `init_workflow(mode, stage)` — constructs initial `wf_state`, writes context
- [x] Unit tests UT-1, UT-5 pass (no Stan required)
- [x] Display contract unit tests UT-2 pass for all 4 branches (using mock fit objects)

### Acceptance Criterion
`print(wf)` on a `wf_state` object with `diagnostics$passed = FALSE` and `acknowledged = FALSE` shows the failure message and `wf$diagnose()` prompt, and **no coefficient output**. Verified by UT-2.

---

## Phase 2: Diagnostic Infrastructure

**ID:** PHASE-2 | **Duration:** ~2 weeks | **Dependencies:** PHASE-1 complete; ODQ-1 resolved

### Description
Implement the diagnostic registry, the generic MCMC diagnostic runner, the `diagnose.wf_state()` method, and the four built-in family-specific registry entries.

### Milestones
- [x] `R/diagnostic_registry.R`: `DIAGNOSTIC_REGISTRY` with entries for bernoulli, poisson, negbinomial, gaussian-hierarchical, time-series, unknown
- [x] `R/diagnostics.R`: `run_diagnostics(fit, wf_state)` — generic Rhat/ESS/divergences/BFMI/treedepth checks + registry dispatch
- [x] `diagnose.wf_state()` method: displays bayesplot panels, prompts acknowledgment, updates `wf_state$diagnostics$acknowledged`, appends to audit trail
- [x] `detect_parameterization(fit)` helper: returns "centered" or "non-centered"
- [x] `refit_noncentered(wf)` one-call switch (stub for Phase 4 use)
- [x] Unit tests UT-3 (smoke), UT-4 pass

### Acceptance Criterion
`run_diagnostics(fit_good, wf)` returns `wf_state$diagnostics$passed == TRUE` for the known-good model. `run_diagnostics(fit_bad, wf)` returns `passed == FALSE` with `failed_criteria` non-empty and the Bernoulli rare-event warning present.

---

## Phase 3: Phase Scripts — Fitting and Off-Ramps

**ID:** PHASE-3 | **Duration:** ~2 weeks | **Dependencies:** PHASE-2 complete; ODQ-2 resolved; ODQ-3 clarified

### Description
Implement the R scripts for Phases 1–3 of the Bayesian workflow (goal declaration, prior specification, fitting). Implement the off-ramp assessment and exit workflow.

### Milestones
- [x] `R/offramps.R`: `assess_offramps(outcome_type, n, event_rate, goal)` — decision matrix from DESIGN.md §4
- [x] `R/exit_workflow.R`: `exit_workflow(wf, method, alternatives)` — writes YAML exit log
- [x] Phase 1 script: data inspection, `assess_offramps()` call, choice logging
- [x] Phase 2 script: `set_prior()` / `prior()` interface; prior predictive simulation; stores in `wf_state`
- [x] Phase 3 script: `brm()` call with cmdstanr backend; `detect_parameterization()`; stores results in `wf_state`; calls `export_context()`
- [x] Unit tests UT-6 (off-ramp contract), UT-3 (exit log schema) pass

### Acceptance Criterion
Running Phase 1 on `dat_bad` (binary, n=150, event rate 3%) produces an off-ramp list with ≥ 2 alternatives including "Firth penalized logistic" and "full_stan". Running `exit_workflow()` produces a valid YAML file. SCENARIO-3 assertions 3.1–3.6 pass.

---

## Phase 4: Quarto Templates — Phases 1–6

**ID:** PHASE-4 | **Duration:** ~3 weeks | **Dependencies:** PHASE-3 complete; ODQ-1 resolved

### Description
Implement all six phase Quarto templates (`.qmd` files). Each template sources the appropriate Phase script, uses the Posit Assistant context block, and renders correctly in both `mode = "learn"` and `mode = "practice"`.

### Milestones
- [x] `templates/phase1_exploration.qmd` — data inspection, off-ramp table, choice cell
- [x] `templates/phase2_priors.qmd` — prior specification UI, prior predictive plot cell
- [x] `templates/phase3_fit.qmd` — `brm()` call cell, parameterization note
- [x] `templates/phase4_diagnostics.qmd` — `run_diagnostics()` call; mode-aware display (ODQ-1 resolution applied)
- [x] `templates/phase5_ppc.qmd` — bayesplot PPC panels; mode-aware plot selection
- [x] `templates/phase6_loo.qmd` — `loo()` + `loo_compare()` calls; LOO table
- [x] Each template renders without error in both `learn` and `practice` modes

### Acceptance Criterion
Running Phase 4 template on `fit_good` with `mode = "learn"`: renders without error; shows prior-vs-posterior overlay in the output; no coefficient table visible. Running with `mode = "practice"` and clean diagnostics: shows coefficient table.

---

## Phase 5: Phase 7 Report Template and Integration

**ID:** PHASE-5 | **Duration:** ~2 weeks | **Dependencies:** PHASE-4 complete

### Description
Implement the final Quarto report template that integrates all phases into a single reproducible document. Every claim in the report must be traceable to a stored model object or `wf_state` field.

### Milestones
- [x] `templates/bayesflow_report.qmd` — parameterized template (`params: [mode, stage, model_path]`) that `readRDS(params$model_path)` and renders all summaries from `wf_state`
- [x] Report renders from a saved `brmsfit` + `wf_state` without re-fitting
- [x] All diagnostic summaries, coefficient tables, PPC plots, and LOO tables in the report are sourced from `wf_state` fields, not re-computed
- [x] SCENARIO-1 step 1.10 passes

### Acceptance Criterion
SCENARIO-1 passes end-to-end (all 10 assertions).

---

## Phase 6: Acceptance Testing and Documentation Polish

**ID:** PHASE-6 | **Duration:** ~1 week | **Dependencies:** PHASE-5 complete

### Description
Run all three success-criterion scenarios, close any failures, and polish the user-facing documentation.

### Milestones
- [x] SCENARIO-1 (SC-1) passes in `mode = "learn"`, `stage = "explore"`
- [x] SCENARIO-2 (SC-2) passes in both `mode = "learn"` and `mode = "practice"`
- [x] SCENARIO-3 (SC-3) passes (Phase 1 exit and Phase 4 exit)
- [x] SCENARIO-4 passes (off-ramp equal weighting)
- [x] SCENARIO-5 passes (parameterization detection)
- [x] All unit tests (UT-1 through UT-6) pass
- [x] README "Getting Started" section verified by a test user (student persona)
- [x] CHANGELOG.md updated with v0.1.0 entry
- [x] All ODQs resolved or explicitly deferred with rationale

### Acceptance Criterion
All three success criteria (SC-1, SC-2, SC-3) from SDD §6 are verified. **v0.1.0 tag has been applied.**

---

## PHASE-7: User Guidance Layer **[v0.2.0]**

**ID:** PHASE-7 | **Duration:** ~2 weeks | **Dependencies:** PHASE-6 complete (v0.1.0 green); ADR-012 accepted

### Description
Implement the three-tier user guidance layer. No changes to core workflow logic. All new code goes in `R/install.R` and `docs/user-guide/`. The one code change to existing files is the `export_context()` path default (DESIGN.md §10.5).

### Milestones
- [ ] ADR-012 accepted (analysis subfolder layout and `export_context()` path change)
- [ ] `R/install.R`: `rbf_install()`, `rbf_new()`, `guide()`, `rbf_analysis_path()` — all function contracts from DESIGN.md §10
- [ ] `export_context()` default path updated; UT-5 updated to test subfolder path
- [ ] `docs/user-guide/00-overview.md` written
- [ ] `docs/user-guide/01-installation.md` written — covers R, RTools, the one manual `cmdstanr` install step, then `rbf_install()`
- [ ] `docs/user-guide/02-posit-assistant-setup.md` written — covers RStudio and Positron as parallel paths, OpenRouter API key, model selection guidance
- [ ] `docs/user-guide/03-starting-an-analysis.md` written — covers `rbf_new()`, `init_workflow()`, modes, and `guide()`
- [ ] `docs/user-guide/04-workflow-phases.md` written — one section per phase; what to run, what to expect, when to move on
- [ ] `docs/user-guide/05-plotting-reference.md` written — see plotting reference spec in the v0.2.0 brief
- [ ] Unit tests UT-7 through UT-9 pass (see TEST_PLAN additions)
- [ ] SC-1 through SC-5 verified

### Acceptance Criterion
`rbf_install()` runs from a fresh R session on a machine with R ≥ 4.3 and RTools, prints ✓ for all 7 steps, and requires no manual package installation beyond the documented `cmdstanr` prerequisite. `rbf_new("test_analysis")` creates the correct subfolder structure. `guide(wf)` returns the correct next step for all 10 phase-detection branches. All five user-guide documents render as valid Markdown.

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
| PHASE-7 **[v0.2.0]** | ~2 weeks | Week 15 | Guidance layer + user guide |

**Total estimated duration through v0.2.0:** ~15 weeks (part-time, single developer).

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| ODQ-1 (diagnose acknowledgment mechanism) blocks Phase 4 | High | High | Resolved as ADR before PHASE-3 ended; `readline()` accepted for v1.0 as a pragmatic choice |
| Posit Assistant API changes before v1.0 | Medium | Medium | File-based `wf_context.json` (Option Y) is the primary path; not dependent on any private API |
| brms version bump changes `stancode()` output format | Low | Medium | Pin brms version in `renv.lock`; test `detect_parameterization()` against pinned version |
| Stan compilation fails on target platform | Low | High | Smoke-test in PHASE-0; document `cmdstanr::install_cmdstan()` in README |
| **[v0.2.0]** OpenRouter free-tier model retirement leaves the recommended model unavailable | Medium | Low | `02-posit-assistant-setup.md` names two free models and states the selection criterion (context window ≥ `wf_context.json` size), so the user can substitute another OpenRouter model without re-writing the doc |
| **[v0.2.0]** `.Rprofile` in the analysis subfolder conflicts with a user's global `.Rprofile` | Low | Low | The generated `.Rprofile` sources `source_all.R` only; documented in `03-starting-an-analysis.md` so a user with conflicting global settings can adjust |

---

*End of PLAN.md*
