# ADR-011: Expert Mode Inclusion and Scope in v1.0

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §3 (Users), SDD §4 (Modes), DESIGN.md §3 (Display Contract),
DESIGN.md §9 (ODQ-4)

---

## Context

RBayesflow defines three user modes: `"learn"`, `"practice"`, and `"expert"`.
The primary design target is `"learn"` and `"practice"` (students and working
scientists). The science council explicitly scoped out expert Bayesian
statisticians as the *primary* audience on the grounds that they already own
their workflow.

However, two arguments support including `"expert"` in v1.0:

1. **Low implementation cost.** Expert mode is the *simplest* display branch —
   it removes guardrails rather than adding them. The architecture already
   carries `mode` as a field on `wf_state`; adding `else if (mode == "expert")`
   branches to the display contract and LLM behavior logic is approximately
   10–15% additional implementation effort.
2. **Architectural coherence.** A `"not yet implemented"` error for a documented
   mode is confusing; its absence from the interface entirely is a silent
   constraint. Including it cleanly avoids both problems and leaves the door
   open for expert users who find the tool useful for reproducibility or
   teaching, even if they would not use it for their own analysis.

The question is whether to include it fully, stub it, or omit it.

---

## Decision

**Include `mode = "expert"` in v1.0 as a fully implemented, deliberately minimal
mode.**

Expert mode is not a stripped-down version of `"practice"` — it is a distinct
design with a specific purpose: providing a reproducible audit bundle with
maximum transparency and minimal interference. It is appropriate for:
- Experts who want `wf_state`'s audit trail and exit log for reproducibility
  without the guardrails designed for non-experts.
- Instructors running demonstrations where students observe expert workflow.
- Users who have graduated from `"practice"` and no longer need the
  acknowledgment gate.

---

## Expert Mode Specification

### Display contract

`print.wf_state(wf)` when `wf$mode == "expert"`:

- Always shows the full coefficient summary (`fixef(fit)` or equivalent) with
  credible intervals, regardless of diagnostic state.
- Shows a one-line diagnostic health summary beneath coefficients.
- If diagnostics failed: shows a `⚠` warning with the names of failed criteria
  but does **not** suppress coefficients and does **not** require
  `wf$diagnose()` acknowledgment.
- No prior-vs-posterior overlay (available on explicit request via
  `plot_prior_posterior(wf)` but not shown by default).

### `wf$diagnose()` behavior

- Displays the full `bayesplot` diagnostic panel suite.
- Does **not** prompt for acknowledgment — calling `wf$diagnose()` itself
  constitutes review for an expert user.
- Appends to `audit_trail` automatically on call (no readline prompt).

### LLM (Posit Assistant) behavior

- Silent by default. Posit Assistant does not generate unsolicited commentary
  on diagnostic results or next steps.
- Responds to explicit queries (`wf$ask("Why did this diverge?")`) with full
  technical detail, no Socratic scaffolding.
- Code proposals are still shown in non-executable preview by default
  (consistent with ADR-006; this is a Posit Assistant configuration setting,
  not mode-specific).

### Off-ramps

- Offered at Phase 1 (goal declaration) with equal visual weight, consistent
  with all modes.
- Not offered proactively during later phases (the expert is assumed capable of
  recognizing when to exit without prompting).

### Audit trail and exit log

- Identical to `"practice"` mode — full audit trail is always written.
- `exit_workflow()` produces the same machine-generated YAML exit log.
- The audit bundle (`wf_context.json`) is the primary value proposition for
  expert users: a complete, serialized record of every decision made during
  the analysis.

---

## Implementation Notes

The following functions require an `"expert"` branch:

| Function | Expert behavior |
|---|---|
| `print.wf_state()` | Show full coefficients regardless of diagnostic state |
| `diagnose.wf_state()` | No `readline()` prompt; auto-append to audit trail |
| `init_workflow()` | Accept `mode = "expert"` without error |
| `export_context()` | Unchanged — full `wf_state` JSON export |
| `cat_wf_header()` | Display `[expert mode]` tag in header |

No new R files are required. Expert branches are added inline in existing
functions.

---

## Rationale

1. **Low cost, meaningful value.** The display contract for expert mode is the
   simplest of the three (show everything, no gating). Implementing it costs
   less than the scaffolding for learn or practice mode.
2. **Avoids confusing stubs.** A `stop("not yet implemented")` error in a
   documented, named mode creates a worse user experience than either full
   implementation or complete omission. Since full implementation is cheap,
   there is no reason to stub it.
3. **Teaching use case.** An instructor demonstrating Bayesian workflow to
   students benefits from running in `mode = "expert"` (no guardrails
   interrupting the demonstration) while the students follow along in
   `mode = "learn"`. This use case is realistic and served at zero additional
   cost.
4. **Reproducibility value for experts.** Even a Bayesian expert benefits from
   `wf_state`'s audit trail, parameterization logging, and exit log — these
   are reproducibility tools, not scaffolding. Expert mode makes the tool
   useful to this audience without imposing non-expert guardrails.

---

## Consequences

- **Positive:** `mode = "expert"` is fully usable in v1.0 with no stubs or
  error messages.
- **Positive:** Instructors can run demonstrations in expert mode alongside
  student learn-mode sessions.
- **Positive:** The audit bundle remains available to all users regardless of
  mode — reproducibility is a universal feature, not a beginner feature.
- **Negative:** Expert mode requires approximately 10–15% additional
  implementation effort above the `"learn"` + `"practice"` baseline.
  Acceptable given the value.
- **Negative:** Expert users who do not want any scaffolding may still find
  `wf_state` management (two objects: `fit` and `wf`) unfamiliar. Mitigated
  by clear README documentation and the `record_fit()` helper that handles
  linkage in one call.
- **Watch:** If expert users request raw `cmdstanr` integration (bypassing
  brms entirely), revisit ADR-002 and consider adding a direct `cmdstanr` path
  gated to `mode = "expert"`. This is out of scope for v1.0.
