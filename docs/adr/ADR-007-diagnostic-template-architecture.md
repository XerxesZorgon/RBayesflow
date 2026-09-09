# ADR-007: Family-Specific Diagnostic Template Architecture — Hardcoded vs. Extensible Registry

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §5 Phase 4 (MCMC Diagnostics), DESIGN.md §6 (Diagnostic Registry)

---

## Context

RBayesflow's Phase 4 (MCMC Diagnostics) and Phase 5 (Posterior Predictive Checks) require family-specific checks beyond generic MCMC diagnostics. The brief identifies four known failure modes that generic tools miss:

1. **Bernoulli with rare events** — check proportion of zeros vs. predicted
2. **Count models (Poisson / NegBinomial)** — check variance-to-mean ratio (overdispersion)
3. **Hierarchical with few groups** — warn about centered parameterization inefficiency; suggest non-centered
4. **Time-series outcomes** — check autocorrelated residuals

Two architectural approaches:

**Option A — Hardcoded family-specific branches**  
- A single large `run_diagnostics(fit, wf_state)` function with `if (family == "bernoulli") {...} else if (family == "poisson") {...}` branches
- Simple to implement; trivial to read
- Extending to a new family requires editing core code
- Coupling: adding a new check for an existing family also requires editing core code
- Suitable if the set of families is fixed and small

**Option B — Extensible registry (named list of diagnostic functions)**  
- A global `DIAGNOSTIC_REGISTRY` object: a named list where each key is a family name (or pattern) and each value is a function `function(fit, wf_state) -> list(checks = ..., warnings = ...)`
- `run_diagnostics(fit, wf_state)` looks up the family in the registry, runs the generic checks (Rhat, ESS, divergences, BFMI, treedepth), then calls the family-specific function if one is registered
- Extending to a new family = adding one entry to the registry; no core code changes
- Users (or advanced plugins) can register custom diagnostic functions without forking the codebase

---

## Decision

**Option B: extensible registry.**

`DIAGNOSTIC_REGISTRY` is a named list defined in `R/diagnostic_registry.R`. Generic diagnostics run unconditionally. Family-specific diagnostics are looked up by family name (e.g., `"bernoulli"`, `"poisson"`, `"negbinomial"`, `"gaussian"`) and called if found. An `"unknown"` fallback entry runs a generic PPC check with a warning that no family-specific template is registered.

---

## Rationale

1. **The set of relevant families is not fixed.** brms supports ≥ 40 families. The brief identifies four critical ones, but a working scientist fitting Beta-regression, hurdle models, or ordinal outcomes needs the same framework. A hardcoded switch can only grow by editing core functions.
2. **Low implementation cost.** The registry is a named list of functions — three to five lines per family entry. The dispatch loop is ten lines. The complexity cost of Option B over Option A is negligible.
3. **Separation of concerns.** Each family-specific function owns its own checks, warnings, and bayesplot calls. The generic runner owns iteration and result aggregation. These concerns should not be interleaved in a single large function.
4. **Extensibility enables future Wild Peaches articles.** New family templates can be added as the workflow is used in practice, documented, and published — without requiring a fork or a version bump of core code.

---

## Registry Entry Contract

Each registry entry is a function with the following signature and return contract:

```r
# Signature
function(fit, wf_state) {
  list(
    checks = list(
      # Named logical list of per-check pass/fail results
      # e.g., list(zero_inflation_ok = TRUE, overdispersion_ratio = 1.3)
    ),
    warnings = character(),   # Zero or more warning strings
    plots    = list()         # Named list of ggplot/bayesplot objects (may be empty)
  )
}
```

The runner merges `checks` into `wf_state$diagnostics` and appends `warnings` to the audit trail. Plots are returned for display according to the mode/stage display contract.

---

## Built-in Registry Entries (v1.0)

| Family | Key checks | Specific warning |
|---|---|---|
| `bernoulli` | Zero proportion vs. predicted; calibration | "Rare event: consider Firth logistic or penalized regression" if event rate < 5% |
| `poisson` | Variance-to-mean ratio (overdispersion) | "Overdispersion detected (ratio > 2): consider `negbinomial`" |
| `negbinomial` | Overdispersion shape parameter; zero inflation | "Zero inflation suspected: consider hurdle model" if excess zeros |
| `gaussian` (hierarchical) | Centered vs. non-centered flag; divergence check | "Hierarchical model with few groups (n < 5): centered parameterization may be inefficient" |
| `(time-series)` | ACF of residuals; Durbin-Watson proxy | "Autocorrelated residuals detected: consider AR(1) error structure" |
| `unknown` | Generic PPC only | "No family-specific diagnostic template registered for family `{family}`" |

---

## Consequences

- **Positive:** Adding a new family requires editing only `R/diagnostic_registry.R`, not any core runner function.
- **Positive:** Family-specific checks are independently testable as pure functions.
- **Positive:** The registry is inspectable at the R console: `names(DIAGNOSTIC_REGISTRY)` shows all registered families.
- **Negative:** If a user fits a model with `family = negbinomial()` (an object, not a string), the registry lookup requires normalizing family objects to string keys. **Mitigation:** a `family_key(fit)` helper function extracts the canonical family name string from a `brmsfit` object.
- **Negative:** Registry state is global. If multiple projects with different registry customizations are loaded in one R session, entries may conflict. **Mitigation:** document that the registry is project-scoped; `init_workflow()` resets it to defaults.
