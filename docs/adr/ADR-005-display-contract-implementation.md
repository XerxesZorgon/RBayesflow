# ADR-005: Display Contract Implementation — S3 Print Method vs. Wrapper Function

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §4 (Modes and Stages), DESIGN.md §3 (Display Contract), ADR-004

---

## Context

The display contract is the most important architectural decision in RBayesflow:

> What `print(wf)` shows depends on mode and diagnostic state. In any failure case, the coefficient summary is withheld until the user calls `wf$diagnose()` and acknowledges the finding.

| State | `mode = "learn"` | `mode = "practice"` |
|---|---|---|
| Diagnostics clean | Prior-vs-posterior overlay + one-line health summary | Coefficients + credible intervals + one-line health summary |
| Diagnostics failed | Overlay suppressed; failure message + `wf$diagnose()` prompt | Coefficients suppressed; failure message + `wf$diagnose()` prompt |

Two implementation strategies:

**Option A — S3 `print` method on `wf_state` class**  
- Define `print.wf_state()` which reads `wf_state$mode`, `wf_state$diagnostics$passed`, and `wf_state$diagnostics$acknowledged`
- Dispatched automatically whenever user types a `wf_state` object at the console or calls `print(wf)`
- Standard R idiom; zero extra machinery
- Requires ADR-004 decision (S3 object) as prerequisite

**Option B — Wrapper function `show_wf(wf_state)`**  
- An explicit function users must call instead of relying on print dispatch
- Does not intercept accidental `print()` calls
- Breaks the natural R workflow of typing an object name at the console
- Could be used as a fallback if the `wf_state` is not a formal S3 object (plain list)

---

## Decision

**Option A: S3 `print.wf_state()` method, plus `wf$diagnose()` as a method that mutates `wf_state$diagnostics$acknowledged <- TRUE` and returns the updated object.**

The display contract is enforced by S3 dispatch. Users interact with the workflow by typing the object name at the console — exactly the same as they would with `lm()` or `brm()` output.

---

## Rationale

1. **The gate must be automatic.** If users can bypass the gate by calling any function other than the standard R idiom, the gate fails its purpose. S3 dispatch on `print()` is the only mechanism that fires automatically when a user types the object name.
2. **Consistency with the R ecosystem.** Every major R modeling library (lm, glm, brm, stanfit) uses S3 print dispatch. Deviating from this idiom creates a confusing double standard.
3. **`wf$diagnose()` as acknowledgment gate.** The method reads the current diagnostic state, displays the relevant plots (bayesplot panels), prompts the user to confirm they have reviewed the finding, sets `acknowledged = TRUE`, and returns the updated `wf_state`. The print method checks `acknowledged` before displaying coefficients.
4. **Immutable audit trail.** `wf$diagnose()` appends to `wf_state$audit_trail` when acknowledgment occurs, recording the timestamp and the failed criterion. This feeds the exit log schema.

---

## Consequences

- **Positive:** Typing a `wf_state` object at the console always invokes the display contract. No user can accidentally bypass the coefficient gate.
- **Positive:** `wf$diagnose()` provides a natural interaction pattern: inspect → acknowledge → proceed.
- **Positive:** The audit trail entry from `wf$diagnose()` is available for the exit log at any stage.
- **Negative:** `print.wf_state()` must be careful not to conflict with brm's own print method if a `wf_state` is attached to a `brmsfit` object. **Resolution:** `wf_state` is stored as a separate object (or as `attr(fit, "wf_state")`); print dispatch is on `wf_state` objects, not on `brmsfit` objects. Users type `wf` (the `wf_state`), not `fit` (the `brmsfit`), to see the workflow-governed output.
- **Negative:** Users who call `summary(fit)` directly on a `brmsfit` bypass the gate. **Resolution:** documented as a known bypass in the README. The gate is a workflow guardrail, not a security control. Advanced users who call `summary(fit)` directly have intentionally stepped outside the workflow.
- **Watch:** Evaluate whether attaching `wf_state` as an attribute of the `brmsfit` object (via `attr()`) and defining `print.brmsfit_wf()` on a subclass is cleaner. Deferred to implementation phase.
