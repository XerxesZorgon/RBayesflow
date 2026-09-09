# RBayesflow — Software Description Document (SDD.md)

**Status:** Draft  
**Last updated:** 2026-09-08  
**Owner:** John Peach (john.x.peach@gmail.com)  
**Spec Kit role:** The **WHY** document — defines what RBayesflow is, who it is for, and the non-negotiable success criteria against which all design decisions are evaluated. Feeds `DESIGN.md` (HOW) and `PLAN.md` (WHEN).

---

## 1. Purpose

RBayesflow fills the gap between a Bayesian textbook and a working analysis. Students and working scientists using R and Stan have access to excellent individual libraries (brms, cmdstanr, bayesplot, posterior, tidybayes, loo) but no tool that sequences them into a reproducible, decision-aware pipeline while keeping every intermediate step visible and understandable. RBayesflow is that sequencer: a set of R scripts, Quarto templates, and a workflow-state object that guides users through the iterative Bayesian workflow described in Gelman et al. (2020), with Posit Assistant reading live session state to interpret diagnostics and suggest next steps.

No new statistical methods are implemented. No new R package is built or submitted to CRAN.

**Primary reference:** Gelman, A. et al. (2020). *Bayesian Workflow.* arXiv:2011.01808.  
<https://users.aalto.fi/~ave/Bayesian-Workflow.pdf>

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
Knows R and frequentist statistics; has read some Bayesian material; has installed Stan but only run tutorial examples. Needs hand-holding at the conceptual level, not at the R syntax level. Will use `mode = "learn"`.

### 3.2 Working Scientist
Uses brms or rstanarm occasionally; competent in their domain; not a statistician; values reproducibility over deep understanding. Will use `mode = "practice"`.

### 3.3 Expert (out of scope for primary design; architecture must not preclude)
Already owns a Bayesian workflow. May use `mode = "expert"` for audit-bundle access and silent LLM.

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

All three must be satisfied before v1.0 is declared:

**SC-1 (Learn mode completeness):** A user in `mode = "learn"` can complete a full prior predictive → fit → diagnostics → posterior predictive check → model comparison → reporting loop on a standard GLMM using only the workflow scripts and existing libraries, without writing code outside the templates.

**SC-2 (Diagnostic gate):** A user in any mode cannot read coefficient output from a failed fit — the workflow state object withholds the summary display until the diagnostic failure is acknowledged and logged.

**SC-3 (Exit log):** A user can exit the Bayesian path at any stage via `exit_workflow()`, which writes a machine-generated YAML log (method chosen, justification drawn from session objects, evidence inspected, user confirmation) before closing the stage.

---

## 7. Non-Goals (Out of Scope)

- Building a new R package or submitting to CRAN
- Writing custom Stan programs (workflow uses brms / cmdstanr formula interfaces only)
- Python, Julia, or any non-R interface
- Production deployment or model serving
- Building new graphics libraries (integrates bayesplot, tidybayes, ggplot2, esquisse)
- Expert Bayesian statisticians who already own their workflow

---

## 8. Constraints

- **Technology:** R (≥ 4.3), Stan (via brms ≥ 2.21 or cmdstanr ≥ 0.7), Quarto (≥ 1.4), RStudio (≥ 2024.04) with Posit Assistant enabled
- **No new statistical code** — every computation delegates to an existing library
- **No CRAN submission** — deliverable is a project folder (R scripts + Quarto templates + configuration)
- **Parameterization transparency** — when brms chooses a parameterization (e.g., centered vs. non-centered for hierarchical models), the choice must be logged to `wf_state` with an explicit note
- **Code proposals from Posit Assistant** must appear in a non-executable preview by default, requiring explicit user action to run

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

---

## 10. Open Questions at SDD Stage

- Exact serialization format for `wf_state` to Posit Assistant context (JSON vs. YAML vs. `dput()`) — deferred to ADR-006
- Whether `mode = "expert"` ships in v1.0 or is architecturally reserved — deferred to PLAN.md
- Minimum R and Stan version requirements to pin — resolved at environment setup

---

*End of SDD.md*
