# RBayesflow

A structured RStudio workflow for the iterative Bayesian analysis process described in Gelman et al. (2020). RBayesflow sequences the existing R/Stan ecosystem (brms, bayesplot, posterior, tidybayes, loo) into a reproducible, mode-aware pipeline, with Posit Assistant reading live session state to interpret diagnostics and guide the next step.

**No new statistical methods are implemented. No new R package is built.** RBayesflow is glue, sequencing, and display logic.

---

## What It Does

- Guides you through the seven phases of the Bayesian workflow: goal declaration → prior specification → fitting → diagnostics → posterior predictive checks → model comparison → reporting
- Enforces an evidence gate: you cannot read coefficient output from a failed fit until you acknowledge the diagnostic finding
- Offers equally-weighted non-Bayesian alternatives at the start of each analysis, before any fitting occurs
- Logs every decision (method choice, diagnostic acknowledgment, off-ramp selection) to a machine-readable audit trail
- Generates a Posit Assistant context file so the LLM can read your current workflow state and provide stage-appropriate guidance

---

## Modes

Set once per project via `init_workflow(mode = ...)`:

- `mode = "learn"` — prior-vs-posterior overlays; Posit Assistant asks Socratic questions; coefficients withheld until you articulate what you see in the diagnostic plot
- `mode = "practice"` — standard output; diagnostics suppressed on failure until acknowledged; Posit Assistant interprets on request
- `mode = "expert"` — standard output; full audit bundle; Posit Assistant silent unless queried *(architecturally reserved; not fully implemented in v0.1.0)*

---

## Prerequisites

- R ≥ 4.3
- RStudio ≥ 2024.04 (with Quarto bundled) and Posit Assistant enabled
- CmdStan ≥ 2.33 (installed via `cmdstanr::install_cmdstan()`)
- brms ≥ 2.21
- Internet access for first-time Stan compilation

---

## Getting Started

### 1. Clone or download

```bash
git clone https://github.com/<your-username>/RBayesflow.git
```

### 2. Install dependencies

Open `RBayesflow.Rproj` in RStudio, then:

```r
# Restore pinned environment
renv::restore()

# Verify Stan compiles
source("R/stan_demo.R")
```

### 3. Initialize a workflow session

```r
source("R/init.R")
wf <- init_workflow(mode = "learn", stage = "explore")
```

### 4. Work through the phase templates

Open the templates in order from the `templates/` folder:

1. `phase1_exploration.qmd` — declare your goal; inspect your data; consider off-ramps
2. `phase2_priors.qmd` — specify and visualize priors
3. `phase3_fit.qmd` — fit the model
4. `phase4_diagnostics.qmd` — review MCMC diagnostics
5. `phase5_ppc.qmd` — posterior predictive checks
6. `phase6_loo.qmd` — model comparison
7. `bayesflow_report.qmd` — render the final reproducible report

### 5. Check Posit Assistant context

After each phase, call:

```r
export_context(wf)
```

This writes `wf_context.json` to the project root. Posit Assistant reads this file for context-aware guidance.

---

## Key Functions

| Function | Description |
|---|---|
| `init_workflow(mode, stage)` | Create a `wf_state` object and initialize the session |
| `print(wf)` | Display mode- and diagnostic-appropriate output |
| `wf$diagnose()` | Review failing diagnostic plots and acknowledge the finding |
| `assess_offramps(...)` | Get equally-weighted non-Bayesian alternatives |
| `exit_workflow(wf, method)` | Log the exit decision and write `exit.yaml` |
| `export_context(wf)` | Serialize `wf_state` to `wf_context.json` for Posit Assistant |
| `refit_noncentered(wf)` | Switch a hierarchical model to non-centered parameterization |

---

## The Evidence Gate

The most important behavioral guarantee: **a user in any mode cannot read coefficient output from a failed fit**.

```r
# After fitting, if diagnostics fail:
print(wf)
# ✗ MCMC diagnostics failed: Rhat > 1.05 (max: 1.23), 47 divergences
# → Call wf$diagnose() to review the failing checks.

wf <- wf$diagnose()   # Shows bayesplot panels; prompts acknowledgment
# Have you reviewed the diagnostic plots? (yes/no): yes
# Acknowledged. You may now access coefficient output via print(wf).

print(wf)   # Now shows coefficients (with caveat in learn/practice mode)
```

Calling `summary(fit)` directly on the underlying `brmsfit` object bypasses the gate. This is a documented limitation — the gate is a workflow guardrail, not a security control.

---

## Non-Bayesian Off-Ramps

At Phase 1, before any fitting, `assess_offramps()` inspects your data and offers alternatives with equal weight:

```
Your data: binary outcome, n = 150, observed event rate: 3%

Alternative methods for "coefficient estimation":
  [1] Logistic regression + bootstrap CI
      Stakes: Standard; valid for n > ~50, but unreliable at this event rate
  [2] Firth penalized logistic regression
      Stakes: Designed for rare events; recommended when event rate < 5%
  [3] Full Bayesian (Stan via brms)
      Stakes: Full uncertainty propagation; requires prior specification

Choose an alternative (1–3), or press Enter to continue with Bayesian:
```

Your choice is logged to the audit trail. If you choose an alternative, `exit_workflow()` writes a machine-generated YAML exit log.

---

## Reference

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B., Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020). *Bayesian Workflow.* arXiv:2011.01808.  
<https://users.aalto.fi/~ave/Bayesian-Workflow.pdf>

---

## Project Documentation

| Document | Contents |
|---|---|
| `docs/SDD.md` | What RBayesflow is and why it exists |
| `docs/DESIGN.md` | How it works: schemas, contracts, logic |
| `docs/PLAN.md` | Phased implementation plan |
| `docs/TEST_PLAN.md` | Acceptance tests and unit tests |
| `docs/adr/` | Architecture Decision Records (ADR-001 through ADR-007) |
| `CHANGELOG.md` | Version history |

---

## License

MIT. See `LICENSE`.
