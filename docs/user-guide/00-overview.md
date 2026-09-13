# RBayesflow — Overview

RBayesflow is a sequencer for Bayesian analysis in R. It is not an R package and implements no new statistical methods. It is a set of R scripts, Quarto templates, and a workflow-state object (`wf_state`) that guides you through the seven-phase iterative Bayesian workflow described in Gelman et al. (2020), using established libraries — brms, cmdstanr, bayesplot, posterior, tidybayes, loo, ggplot2, and esquisse — at every step.

---

## The Three-Tier Help System

RBayesflow provides three layers of help:

1. **Written documentation** (`docs/user-guide/`) — what to do at each step, when to do it, and why. Start here.
2. **Software guidance** (`guide(wf)`) — call this function at any point in an R session and it reads your current workflow state and tells you exactly what phase you are in and what to run next.
3. **AI assistance** (Posit Assistant) — once configured with an OpenRouter API key, Posit Assistant reads `wf_context.json` and can interpret diagnostics, explain results, and suggest next steps in plain language. It never runs code for you.

---

## User Modes

RBayesflow has three modes, set once at the start of each analysis by calling `init_workflow(mode = ...)`.

**learn** — designed for students and newcomers to Bayesian analysis. In this mode the workflow shows prior-vs-posterior overlay plots instead of coefficient tables, and Posit Assistant asks Socratic questions pointing to the relevant section of Gelman et al. (2020) after each diagnostic display. Coefficient output is withheld until you have articulated what the diagnostics show.

**practice** — designed for working scientists who use Bayesian methods occasionally and value reproducibility. Standard coefficient output is shown once diagnostics pass (or are explicitly acknowledged). Posit Assistant interprets diagnostics on request only.

**expert** — for users who already own a Bayesian workflow and want the audit-bundle access and silent LLM that RBayesflow provides. No guardrails; no diagnostic gate.

---

## The Seven Workflow Phases

1. **Goal declaration and data inspection** — explore your data with esquisse and decide whether a full Bayesian fit is warranted.
2. **Prior specification and prior predictive simulation** — set priors and verify they produce plausible data before fitting.
3. **Model fitting** — fit via brms or cmdstanr formula interfaces.
4. **MCMC diagnostics** — check Rhat, ESS, divergences, BFMI, and treedepth. Coefficient output is gated until this step passes.
5. **Posterior predictive checks** — verify the model reproduces the structure of your data.
6. **Model comparison** — compare candidate models with LOO-CV.
7. **Reporting** — render a Quarto report in which every claim traces to a saved model object or `wf_state` field.

---

## Quick Start

To set up the environment and start a new analysis, run these three commands from the RBayesflow project root:

```r
source("R/install.R")
rbf_install()
rbf_new("my_analysis")
```

Then open `data/my_analysis/` as your working directory in RStudio or Positron and call:

```r
wf <- init_workflow(mode = "learn")
guide(wf)
```

See `docs/user-guide/01-installation.md` for full installation instructions.

---

*Primary reference: Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808. <https://users.aalto.fi/~ave/Bayesian-Workflow.pdf>*
