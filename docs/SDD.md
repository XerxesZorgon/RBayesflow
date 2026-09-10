# RBayesflow — Software Description Document (SDD.md)

**Status:** Active
**Last updated:** 2026-09-10
**Owner:** John Peach (john.x.peach@gmail.com)
**Spec Kit role:** The **WHY** document — defines what RBayesflow is, who it is for, and the non-negotiable success criteria against which all design decisions are evaluated. Feeds `DESIGN.md` (HOW) and `PLAN.md` (WHEN).

---

## 1. Purpose

RBayesflow fills the gap between a Bayesian textbook and a working analysis. Students and working scientists using R and Stan have access to excellent individual libraries (brms, cmdstanr, bayesplot, posterior, tidybayes, loo) but no tool that sequences them into a reproducible, decision-aware pipeline while keeping every intermediate step visible and understandable. RBayesflow is that sequencer: a set of R scripts, Quarto templates, and a workflow-state object that guides users through the iterative Bayesian workflow described in Gelman et al. (2020), with Posit Assistant reading live session state to interpret diagnostics and suggest next steps.

No new statistical methods are implemented. No new R package is built or submitted to CRAN.

**Primary reference:** Gelman, A. et al. (2020). *Bayesian Workflow.* arXiv:2011.01808.
<https://users.aalto.fi/~ave/Bayesian-Workflow.pdf>

**[v0.2.0]** v0.2.0 adds a user guidance layer comprising automated environment verification (`rbf_install()`), analysis scaffolding (`rbf_new()`), in-session step guidance (`guide()`), and a complete written user guide in `docs/user-guide/`. No new statistical methods. No package submission. The guidance layer operates on top of the v0.1.0 core; core workflow logic is unchanged. One small behaviour change: `export_context()` writes `wf_context.json` into each analysis subfolder (`data/<name>/`) rather than the project root when called from inside such a folder, so multiple concurrent analyses do not overwrite each other's context — see ADR-012.

---

## 2. Principles

1. **No statistical black boxes.** Every computation delegates to an established library. RBayesflow is glue, sequencing, and display logic only.
2. **The workflow enforces epistemic discipline.** A user cannot read coefficient output from a failed fit. The workflow state object withholds the summary display until the diagnostic failure is acknowledged and logged.
3. **Mode-appropriate transparency.** What the workflow shows depends on the user's declared learning mode. A learner sees overlays and Socratic prompts; a practitioner sees standard output.
4. **Non-Bayesian paths are first-class.** At goal-declaration time, equally weighted alternatives to a full Bayesian fit are offered. Exiting to a simpler method is a documented, logged outcome, not a fallback.
5. **All code proposals require explicit user action to run.** LLM suggestions (from Posit Assistant) must appear in a non-executable preview by default.
6. **Reproducibility over convenience.** Every fitted object carries a `wf_state` list recording the priors, diagnostic results, and audit trail that produced it.

---

## 3. Users

### 3.1 Student
Knows R and frequentist statistics; has read some Bayesian material; has installed Stan but only run tutorial examples. Needs hand-holding at the conceptual level, not at the R syntax level. Will use `mode = "learn"`. **[v0.2.0]** The guidance layer (`rbf_install()`, `rbf_new()`, `guide()`, and the written user guide) is designed for this persona.

### 3.2 Working Scientist
Uses brms or rstanarm occasionally; competent in their domain; not a statistician; values reproducibility over deep understanding. Will use `mode = "practice"`. **[v0.2.0]** The guidance layer is designed for this persona.

### 3.3 Expert (out of scope for primary design; architecture must not preclude)
Already owns a Bayesian workflow. May use `mode = "expert"` for audit-bundle access and silent LLM. Unaffected by the v0.2.0 guidance layer.

---

## 4. Modes and Stages

### 4.1 User Modes (set once per project)

| Mode | Print output | LLM behavior | Coefficient gate |
|---|---|---|---|
| `"learn"` | Prior-vs-posterior overlay plots | Source pointer + Socratic question | Withheld until user articulates finding |
| `"practice"` | Smart print methods; standard coefficients | Interprets diagnostics on request | Withheld until `wf$diagnose()` acknowledged |
| `"expert"` | Standard output; full audit bundle | Silent unless queried | No gate |

### 4.2 Stage Dimension (per fit)

| Stage | Prior policy | Verification protocol | LLM role | Off-ramps |
|---|---|---|---|---|
| `"explore"` | Weakly informative | Basic numerical health checks | Accelerator | Offered proactively at goal declaration |
| `"confirm"` | Must justify priors | Full: prior sensitivity, PPC, artifact audit | Skeptical auditor | Offered, but not proactively |
| `"exit"` | N/A | Machine-generated YAML exit log | Silent | This IS the off-ramp |

---

## 5. Workflow Phases

Phases execute in order; each is a distinct R script and Quarto section:

1. **Goal declaration and data inspection** — exploratory graphics via ggplot2 / esquisse; non-Bayesian off-ramp assessment
2. **Prior specification and prior predictive simulation** — `brms::prior()` objects + predictive draws
3. **Model fitting** — via brms or cmdstanr formula interfaces
4. **MCMC diagnostics** — Rhat, bulk/tail ESS, divergences, BFMI, treedepth
5. **Posterior predictive checks** — bayesplot, tidybayes
6. **Model comparison** — loo (LOO-CV)
7. **Reporting** — Quarto template; all claims traced to saved model objects

---

## 6. Success Criteria

All success criteria must be satisfied before the corresponding release is declared.

**v0.1.0 (verified; complete):**

**SC-1 (Learn mode completeness):** A user in `mode = "learn"` can complete a full prior predictive → fit → diagnostics → posterior predictive check → model comparison → reporting loop on a standard GLMM using only the workflow scripts and existing libraries, without writing code outside the templates.

**SC-2 (Diagnostic gate):** A user in any mode cannot read coefficient output from a failed fit — the workflow state object withholds the summary display until the diagnostic failure is acknowledged and logged.

**SC-3 (Exit log):** A user can exit the Bayesian path at any stage via `exit_workflow()`, which writes a machine-generated YAML log (method chosen, justification drawn from session objects, evidence inspected, user confirmation) before closing the stage.

**[v0.2.0]:**

**SC-4 (Written user guide):** `docs/user-guide/` contains a complete written guide covering installation, Posit Assistant setup for both RStudio and Positron with OpenRouter, analysis initialization, all seven workflow phases, and a plotting reference with interpretation notes for each phase's plots.

**SC-5 (Plotting reference correctness):** The plotting reference correctly distinguishes the Phase 1 role of `esquisse` from the Phases 4–6 diagnostic role of `bayesplot` and `tidybayes`, with a one-paragraph interpretation guide for each distinct plot type.

Additional acceptance targets for v0.2.0 (verified via unit tests UT-7 through UT-9 and manual acceptance):

- A user with R ≥ 4.3 and RTools installed can run `rbf_install()` and see a ✓/✗ result for each prerequisite; no manual package installation steps are required afterward beyond the one documented `cmdstanr` prerequisite.
- `rbf_new("name")` creates `data/name/` with the required starter files, and `init_workflow()` called from that analysis writes `wf_context.json` to `data/name/`, not to the project root.
- `guide(wf)` reads the current `wf_state` and prints — in plain language — the current phase, what the user should expect to see, and the exact next step, with no output requiring interpretation.

---

## 7. Non-Goals (Out of Scope)

- Building a new R package or submitting to CRAN
- Writing custom Stan programs (workflow uses brms / cmdstanr formula interfaces only)
- Python, Julia, or any non-R interface
- Production deployment or model serving
- Building new graphics libraries (integrates bayesplot, tidybayes, ggplot2, esquisse)
- Expert Bayesian statisticians who already own their workflow

**[v0.2.0] additions to out-of-scope:**

- Any GUI or interactive widget. `rbf_install()`, `rbf_new()`, and `guide()` are console functions. No Shiny, no widgets, no menus.
- Changes to core workflow logic (`wf_state`, diagnostic gate, display contract, off-ramps).
- Programmatic Posit Assistant configuration. Setup is documented for the user to perform through the IDE settings UI.
- The Wild Peaches article on RBayesflow (separate future project).

---

## 8. Constraints

- **Technology:** R (≥ 4.3), Stan (via brms ≥ 2.21 or cmdstanr ≥ 0.7), Quarto (≥ 1.4), RStudio (≥ 2026.04 with Posit Assistant 0.7.7+) or Positron (current release)
- **No new statistical code** — every computation delegates to an existing library
- **No CRAN submission** — deliverable is a project folder (R scripts + Quarto templates + configuration)
- **Parameterization transparency** — when brms chooses a parameterization (e.g., centered vs. non-centered for hierarchical models), the choice must be logged to `wf_state` with an explicit note
- **Code proposals from Posit Assistant** must appear in a non-executable preview by default, requiring explicit user action to run

**[v0.2.0]** `rbf_install()` is a standalone R script with no dependency on `source_all.R`; it must be runnable from a fresh R session before any other RBayesflow code is loaded. The one exception: `cmdstanr` must already be installed from r-universe before `rbf_install()` can call it. This one manual step is documented in `docs/user-guide/01-installation.md`.

---

## 9. Dependencies

All existing libraries; no new statistical implementations required:

| Library | Role |
|---|---|
| brms (≥ 2.21) | Primary formula interface for model fitting |
| cmdstanr (≥ 0.7) | Alternative backend; required by brms or standalone |
| bayesplot | Diagnostic and PPC plots |
| posterior | Draws manipulation, Rhat/ESS computation |
| tidybayes | Tidy extraction of posterior draws |
| loo | LOO-CV model comparison |
| ggplot2 | Base graphics layer |
| esquisse | Data-to-viz decision aid (Phase 1) |
| Quarto (≥ 1.4) | Report templates |
| Posit Assistant | LLM integration (reads live `wf_state`) |
| **[v0.2.0]** rprojroot | Root detection used by `source_all.R` and `rbf_analysis_path()` — already in use, now explicit |

---

## 10. Open Questions at SDD Stage

- Exact serialization format for `wf_state` to Posit Assistant context (JSON vs. YAML vs. `dput()`) — resolved in ADR-006 (JSON)
- Whether `mode = "expert"` ships in v1.0 or is architecturally reserved — resolved in ADR-011 (fully implemented)
- Minimum R and Stan version requirements to pin — resolved at environment setup

**[v0.2.0]** All new design questions for the guidance layer are resolved in ADR-012 (analysis subfolder layout and `export_context()` path).

---

*End of SDD.md*
