---
{
  "id": "file_s4114syw",
  "filetype": "document",
  "filename": "CHANGELOG",
  "created_at": "2026-09-24T19:07:36.619Z",
  "updated_at": "2026-09-24T19:07:36.619Z",
  "meta": {
    "location": "/",
    "tags": [],
    "categories": [],
    "description": "",
    "source": "markdown"
  }
}
---
# Changelog

All notable changes to RBayesflow will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).  
Versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.2.0] — 2026-09-10

### Added
- `R/install.R` — four new functions: `rbf_install()`,
  `rbf_new()`, `guide()`, `rbf_analysis_path()`
- `docs/user-guide/` — five user-guide documents: 00-overview,
  01-installation, 02-posit-assistant-setup,
  03-starting-an-analysis, 04-workflow-phases,
  05-plotting-reference
- `docs/adr/ADR-012-analysis-subfolder-layout.md` — accepted

### Changed
- `R/context.R` — `export_context()` path default changed to
  `file.path(rbf_analysis_path(), "wf_context.json")` so multiple
  concurrent analyses do not overwrite each other's context file
- `R/source_all.R` — `R/install.R` added to load order
- `docs/SDD.md`, `docs/DESIGN.md`, `docs/PLAN.md`,
  `docs/TEST_PLAN.md` — all updated with `[v0.2.0]` sections

### Verified
- SC-4: all 6 user-guide files exist and render as valid Markdown
- SC-5: plotting reference contains esquisse, mcmc_trace,
  ppc_dens_overlay, add_epred_draws; 10 "What to look for" sections
- UT-7, UT-8, UT-9: all pass (10 tests, 36 expectations, 0 failures)

## [0.1.0] — 2026-09-08

### Added
- Seven-phase Bayesian workflow (goal declaration → prior specification
  → fitting → diagnostics → posterior predictive checks → model
  comparison → reporting) implemented as interactive Quarto templates
- Three user modes: `learn`, `practice`, `expert`
- `stage` dimension per fit: `explore`, `confirm`, `exit`
- Evidence gate: coefficient output withheld until diagnostic failures
  are acknowledged via `wf$diagnose()`
- Non-Bayesian off-ramps at Phase 1 with equal visual weight;
  machine-generated YAML exit log
- Diagnostic registry with family-specific checks: Bernoulli (rare
  event), Poisson/NegBinomial (overdispersion), Gaussian hierarchical
  (few groups / parameterization), time-series (ACF lag-1)
- Parameterization detection and `refit_noncentered()` helper
- Posit Assistant integration via `wf_context.json`
- `bayesflow_report.qmd` — the only rendered template; reads from
  saved `wf_state` and `brmsfit` objects
- Full audit trail logged to `wf_state` and exported to JSON
- 11 Architecture Decision Records (ADR-001 through ADR-011)

### Verified
- SC-1: Learn mode full loop — SCENARIO-1 all 10 assertions passed
- SC-2: Diagnostic gate — SCENARIO-2 all 7 assertions passed (learn +
  practice modes)
- SC-3: Exit log — SCENARIO-3 all 7 assertions passed
- SCENARIO-4: Off-ramp equal weighting — 4 assertions passed
- SCENARIO-5: Parameterization transparency — 3 assertions passed
- Unit tests: 8 tests, 22 expectations, 0 failures

### Dependencies
- R ≥ 4.3, brms ≥ 2.21, CmdStan ≥ 2.33 (2.39.0 tested)
- Note: CmdStan 2.39.0 requires `array[N] real` syntax in raw Stan
  code (old `real y[N]` syntax rejected)
- Note: `formula(brmsfit)` returns `brmsformula`, not base `formula`;
  wrap with `as.formula()` before using `[[2]]` or `findbars()`
- SC-4 verified: all 6 user-guide files exist (TRUE). SC-5 verified: esquisse, mcmc_trace, ppc_dens_overlay, add_epred_draws present; "What to look for" sections: 10.
