# ADR-004: Workflow State as Named List vs. R6 / S3 Object

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** DESIGN.md §2 (wf_state schema), ADR-005

---

## Context

The `wf_state` workflow state object is the central data structure of RBayesflow. It must:
- Cache mode, stage, priors, prior predictive samples, MCMC diagnostic summary, PPC results, LOO table, audit trail, and exit log
- Be attached to every fitted model object
- Be readable by Posit Assistant from the live R session
- Support a display contract where `print()` behavior depends on mode and diagnostic state
- Survive serialization (for logging, audit, and Posit Assistant context injection)

Three R object paradigms are candidates:

**Option A — Plain named list**  
- `wf_state <- list(mode = ..., stage = ..., diagnostics = list(...), ...)`
- No class system; no methods; any code can read or modify any field
- Simplest to understand and inspect at the R console
- No dispatch mechanism; display contract must be implemented separately

**Option B — S3 object (named list with a class attribute)**  
- `structure(list(...), class = "wf_state")`
- Enables `print.wf_state()`, `summary.wf_state()`, `format.wf_state()` dispatch
- No reference semantics (copy-on-modify); consistent with R's functional style
- Familiar to any R user who has encountered `lm()` output
- Lightweight: no external dependencies

**Option C — R6 class**  
- Reference semantics (mutate in place)
- Requires the `R6` package
- `$` method syntax differs from S3 dispatch
- Appropriate when many methods mutate state and reference semantics matter
- Overkill for a workflow state object that is primarily read, displayed, and serialized

---

## Decision

**Option B: S3 object (named list with `class = "wf_state"`).**

`wf_state` is constructed by a `new_wf_state()` constructor function that returns a named list with `class = c("wf_state", "list")`. S3 methods (`print`, `summary`, `format`, `to_json`) dispatch on this class.

---

## Rationale

1. **Display contract requires dispatch.** The brief's display contract (what `print(wf)` shows depends on mode and diagnostic state) is the most important architectural decision. S3 `print.wf_state()` is the idiomatic R way to implement this; a plain list has no dispatch mechanism.
2. **No reference semantics needed.** The workflow advances through stages by constructing and replacing `wf_state` objects, not by mutating a shared instance. Copy-on-modify is the right behavior.
3. **No external dependencies.** A plain named list + class attribute requires zero additional packages. R6 requires the `R6` package (widely available but still an added dependency; in a project folder without package infrastructure, dependencies must be managed manually).
4. **Inspectability.** An S3 object is still a list; users can type `wf_state$diagnostics` at the console and see the raw values. R6 objects are less transparent at the console.
5. **Familiar idiom.** The output of `lm()`, `glm()`, and `brm()` are all S3 objects. Users fitting Bayesian models are already comfortable with this pattern.

---

## Consequences

- **Positive:** `print(fit)` dispatches to `print.wf_state()` and implements the full display contract (see ADR-005) without any wrapper infrastructure.
- **Positive:** `to_json(wf_state)` can be a simple S3 method that serializes the named list to JSON for Posit Assistant context injection (see ADR-006).
- **Positive:** Users can inspect any field directly: `fit$wf_state$diagnostics$rhat_max`.
- **Negative:** No reference semantics; if multiple R objects share a `wf_state`, they do not automatically stay in sync. This is not a workflow pattern RBayesflow supports, so it is not a practical limitation.
- **Negative:** The constructor must validate all required fields explicitly (no R6 `initialize()` with active bindings). Mitigated by a thorough `new_wf_state()` function that validates inputs on construction.
- **Watch:** If the diagnostic and audit trail grows significantly more complex (e.g., interactive drill-down), revisit R6. That complexity is not in scope for v1.0.
