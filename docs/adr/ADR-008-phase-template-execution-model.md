# ADR-008: Phase Template Execution Model (Interactive Script vs. Rendered Document)

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** DESIGN.md §9 (ODQ-1), PLAN.md Phase 4, TEST_PLAN.md §1

---

## Context

RBayesflow delivers its workflow as seven Quarto (`.qmd`) phase templates plus a
final report template. A key architectural question is how users are expected to
execute these templates: as fully rendered documents (HTML/PDF output produced by
hitting "Render") or as interactive scripts (chunks run one at a time in RStudio
with Ctrl+Enter, producing inline output in the editor).

This decision directly controls how `wf$diagnose()` can request user
acknowledgment of failed diagnostics, because `readline()` — the natural R
mechanism for interactive prompts — is silently swallowed when Quarto renders a
document in batch mode (it returns `""` without blocking).

**Three execution models are candidates:**

**Model A — Interactive script (chunks run interactively; only Phase 7 rendered)**
- Phases 1–6 are `.qmd` files used as structured, annotated scripts. The user
  runs chunks interactively in RStudio; inline output appears in the editor.
- Only the Phase 7 report template is rendered to HTML/PDF.
- `readline()` works correctly throughout Phases 1–6 because they are never
  batch-rendered.
- Consistent with how most R users work with `.Rmd` and `.qmd` files day-to-day.

**Model B — Per-phase render**
- The user renders each phase `.qmd` at the end of that phase to produce a
  per-phase HTML record.
- `readline()` fails silently on render; requires either a rendering guard that
  raises an explicit error or a document-native acknowledgment mechanism.
- Adds complexity to the acknowledgment step; requires the user to understand
  why the console step must precede rendering.

**Model C — Single rendered document**
- All phases live in one `.qmd`; the user renders the entire workflow at once.
- `readline()` fails everywhere in the final render.
- Incompatible with the iterative nature of the Gelman workflow (each phase
  depends on decisions made in the previous one).

---

## Decision

**Model A: phase templates are interactive scripts; only Phase 7 is rendered.**

Phases 1–6 `.qmd` templates are interactive scripts. Users run chunks in RStudio
using Ctrl+Enter (or the Run button). The Phase 7 report template is the only
file intended to be rendered, and it reads exclusively from the completed
`wf_state` object — it contains no calls to `wf$diagnose()` or `readline()`.

A rendering guard is added to `wf$diagnose()` as a safety net:

```r
if (isTRUE(getOption("knitr.in.progress"))) {
  stop(
    "wf$diagnose() must be called interactively, not inside a rendered document.\n",
    "Run this chunk in the RStudio console, then render your report."
  )
}
```

This guard is defensive infrastructure, not a primary design path.

---

## Rationale

1. **Matches the Gelman workflow's iterative nature.** Bayesian workflow is
   deliberately non-linear: a failed PPC in Phase 5 may require returning to
   Phase 2 to revise priors. A rendered document implies a linear, completed
   analysis. Interactive chunks match the actual working pattern.
2. **`readline()` works correctly in interactive mode.** The acknowledgment
   mechanism in `wf$diagnose()` requires no additional machinery when chunks
   are run interactively. The prompt appears in the RStudio Console pane; the
   user types `yes` and proceeds.
3. **Eliminates Quarto-checkbox complexity.** A document-native acknowledgment
   widget (Shiny checkbox, Observable JS input) would require either a Shiny
   server dependency or custom JavaScript — both violate the "use existing
   libraries, minimize new code" constraint.
4. **Familiar to the target audience.** Students and working scientists who use
   RStudio typically run `.Rmd` / `.qmd` chunks interactively and render only
   for final output. This is the standard RStudio workflow.
5. **Phase 7 report is a clean read-only render.** Because all interactive work
   completes before rendering, the report template reads from `wf_state` without
   any side effects, making renders deterministic and repeatable.

---

## Analogy to Pluto (Julia)

Quarto interactive chunks in RStudio are analogous to Pluto notebooks in that
both show inline output after each cell. The critical difference: Pluto is
*reactive* (changing one cell reruns all downstream cells automatically), while
Quarto is *sequential and non-reactive* (cells are independent; re-running one
does not trigger others). RBayesflow's phase templates exploit Quarto's
sequential model intentionally — each chunk advances state explicitly, matching
the analyst's deliberate step-by-step decisions.

---

## Consequences

- **Positive:** `readline()` works correctly throughout Phases 1–6 with no
  additional infrastructure beyond the rendering guard.
- **Positive:** No Shiny or JavaScript dependency introduced.
- **Positive:** Phase 7 report renders deterministically from completed
  `wf_state`; no interactive elements to fail.
- **Positive:** The execution model matches how RStudio users already work.
- **Negative:** Users cannot produce a per-phase HTML record by rendering
  individual phase templates. Mitigated by: (a) inline chunk output in RStudio
  is visible during the session, and (b) the Phase 7 report captures the full
  audit trail from `wf_state`.
- **Negative:** Users unfamiliar with interactive Quarto may attempt to render
  Phase 1–6 templates and encounter the rendering guard error. Mitigated by
  clear template headers and README instructions.
- **Watch:** If demand for per-phase rendered records emerges (e.g., for
  teaching portfolios), consider a `render_phase_summary(wf, phase = 4)` helper
  that generates a read-only HTML snapshot from `wf_state` without re-running
  any interactive logic.
