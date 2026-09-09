# ADR-009: Workflow State Storage (Attribute on brmsfit vs. Separate Top-Level Object)

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** DESIGN.md §2 (wf_state schema), DESIGN.md §9 (ODQ-2), ADR-004

---

## Context

`wf_state` is the central workflow state object (an S3 list; see ADR-004). It
must be accessible alongside the fitted `brmsfit` object throughout the
analysis. Two storage strategies are available:

**Option A — Attribute on `brmsfit`**
- `attr(fit, "wf_state") <- wf`
- Both objects travel together: the user manages one object (`fit`) and accesses
  workflow state via `attr(fit, "wf_state")` or a helper `get_wf_state(fit)`.
- `brmsfit` objects are frequently large (100–500 MB for complex hierarchical
  models with many posterior draws). `saveRDS(fit)` serializes everything,
  including the attached `wf_state`, bloating saved files.
- If the user saves `fit` without `wf_state` (e.g., earlier code that
  predates RBayesflow), they lose the workflow record silently.

**Option B — Separate top-level object**
- `wf` is a standalone named S3 object in the user's environment.
- The user manages two objects (`fit` and `wf`) explicitly.
- `wf_state` references the fit object by a stored hash (SHA-256 of the data,
  formula string, and fit timestamp) rather than by direct pointer.
- `saveRDS(wf, "my_analysis_wf.rds")` is small (kilobytes) and independent.
- The audit trail, diagnostic results, and exit log remain accessible even if
  the large `fit` object is not loaded.

---

## Decision

**Option B: `wf_state` is a separate top-level object (`wf`) in the user's
R environment.**

The user works with two objects: `fit` (the `brmsfit` result) and `wf` (the
`wf_state` S3 object). They are linked by a stored hash, not by R reference.
The workflow initializes both:

```r
wf  <- init_workflow(mode = "learn", stage = "explore")
fit <- brm(formula, data, family, prior, backend = "cmdstanr", seed = 42)
wf  <- record_fit(wf, fit)   # stores hash, timestamp, parameterization
```

The `print.wf_state()` display contract operates on `wf` directly. The Posit
Assistant reads `wf` (via `wf_context.json`) independently of whether `fit` is
loaded in the session.

---

## Rationale

1. **`brmsfit` objects are large.** Complex hierarchical models routinely
   produce `brmsfit` objects exceeding 100 MB due to stored posterior draws.
   Attaching `wf_state` to every such object adds to serialization cost for no
   computational benefit — `wf_state` itself is kilobytes.
2. **Separation of concerns.** `brmsfit` is a computational artifact (Stan
   samples, model code, data). `wf_state` is an epistemic record (decisions
   made, diagnostics run, audit trail). Conflating them inside one object
   obscures this distinction.
3. **Posit Assistant independence.** `wf_context.json` (exported from `wf`) can
   be read by Posit Assistant without loading the large `fit` object. In
   `mode = "learn"`, a student may want to review their workflow log at the
   start of a new session before reloading the fit — this is natural with
   separate storage, awkward with attribute storage.
4. **Explicit beats implicit for the target audience.** Students learning
   Bayesian workflow benefit from seeing `wf` as a distinct object they can
   inspect (`wf$diagnostics`, `wf$audit_trail`) without needing to know about
   R's attribute system.
5. **Robustness to partial saves.** A user who saves only `fit` (the large
   object) to free memory retains their full workflow record in the small `wf`
   file. A user who saves only `wf` can re-run diagnostics on a reloaded `fit`
   and the hash check confirms it is the same model.

---

## Hash Linkage Schema

`wf_state` stores a `fit_hash` field computed at `record_fit()` time:

```r
fit_hash <- digest::digest(
  list(
    formula   = as.character(formula(fit)),
    data_hash = digest::digest(fit$data),
    timestamp = wf$fit_timestamp
  ),
  algo = "sha256"
)
```

When `run_diagnostics(fit, wf)` is called, it verifies `fit_hash` matches
before writing diagnostic results to `wf`. This prevents silently applying
diagnostics from one model to another model's workflow record.

---

## Consequences

- **Positive:** `saveRDS(wf, "wf.rds")` produces a small, portable file that
  is the complete epistemic record of the analysis.
- **Positive:** Posit Assistant can read `wf_context.json` without the large
  `fit` object being loaded.
- **Positive:** Students see `wf` as a first-class object they can inspect and
  reason about independently.
- **Negative:** The user manages two objects. Mitigated by: (a) clear naming
  convention (`fit` and `wf` throughout all templates), (b) `record_fit()`
  making the linkage explicit, and (c) `check_fit_hash()` catching mismatches.
- **Negative:** If a user accidentally passes the wrong `fit` to
  `run_diagnostics()`, the hash check raises an informative error rather than
  silently corrupting the workflow record.
- **Watch:** If users find two-object management confusing in practice,
  consider a thin wrapper `wf_fit <- list(fit = fit, wf = wf)` with class
  dispatch — but only if evidence from user testing warrants it.
