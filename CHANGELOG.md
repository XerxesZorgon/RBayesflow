# Changelog

All notable changes to RBayesflow will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).  
Versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- Project planning documents: SDD.md, DESIGN.md, PLAN.md, TEST_PLAN.md
- Architecture Decision Records: ADR-001 through ADR-007
- README.md with Getting Started section

### Open Design Questions (must be resolved before v0.1.0)
- ODQ-1: `diagnose()` acknowledgment mechanism — `readline()` vs. Quarto checkbox
- ODQ-2: `wf_state` storage — attribute of `brmsfit` vs. separate object
- ODQ-3: Minimum draws for `assess_offramps()` event rate estimation
- ODQ-4: `mode = "expert"` — stub in v1.0 or reserved only

---

## [0.1.0] — *Target: ~13 weeks from project start*

### Planned for this release
- `R/wf_state.R`: `new_wf_state()`, `print.wf_state()`, `summary.wf_state()`, `diagnose.wf_state()`
- `R/diagnostic_registry.R`: extensible family-specific diagnostic registry with four built-in entries (bernoulli, poisson/negbinomial, gaussian-hierarchical, time-series)
- `R/diagnostics.R`: `run_diagnostics()` — generic Rhat/ESS/divergences/BFMI/treedepth checks + registry dispatch
- `R/offramps.R`: `assess_offramps()` — equal-weight off-ramp assessment at Phase 1
- `R/exit_workflow.R`: `exit_workflow()` — YAML exit log generation
- `R/context.R`: `export_context()` — `wf_state` → `wf_context.json` for Posit Assistant
- `R/init.R`: `init_workflow(mode, stage)`
- All six phase Quarto templates (`phase1_exploration.qmd` through `phase6_loo.qmd`)
- Final report template (`bayesflow_report.qmd`)
- Examples: `glmm_gaussian/` (known-good model), `bernoulli_rare_event/` (known-bad model)
- `renv.lock` with pinned dependencies
- All three success criteria verified: SC-1, SC-2, SC-3

### Not in this release
- `mode = "expert"` (architecturally reserved; not implemented)
- Custom Stan programs
- CRAN submission
