# ADR-002: brms vs. cmdstanr as Primary Fitting Backend (or Both)

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §5 (Workflow Phases), SDD §9 (Dependencies)

---

## Context

RBayesflow must support model fitting via existing R-to-Stan interfaces. Two mature options exist:

**Option A — brms only**  
- Formula interface (`brm(y ~ x + (1 | group), data, family, prior)`) familiar to lme4 users
- Handles prior specification, data transformation, and Stan code generation automatically
- Wide family support: Gaussian, Bernoulli, Binomial, Poisson, NegBinomial, student-t, Beta, and many more
- Generates readable Stan code (inspectable via `stancode()`)
- Active development; well-documented
- Limitation: opinionated parameterization (centered hierarchical by default); limited control over sampler settings beyond brms-exposed arguments

**Option B — cmdstanr only**  
- Direct interface to CmdStan; maximum flexibility
- Requires user-written Stan programs → immediately out of scope per the brief
- Not usable for target audience (students and working scientists with no Stan experience)

**Option C — Both: brms as primary, cmdstanr as backend**  
- brms uses cmdstanr (or rstan) as its computational backend
- brms formula interface serves as the primary user-facing API
- cmdstanr is available as the execution engine (recommended over rstan for performance and maintenance)
- Users never write raw Stan; brms generates it
- Advanced users can inspect and extend the generated Stan code if needed

---

## Decision

**Option C: brms as primary formula interface; cmdstanr as the computational backend.**

brms is the only fitting interface users interact with. cmdstanr is required as brms's backend (preferred over rstan for CmdStan 2.x compatibility and future-proofing). No raw Stan code is written by users or by RBayesflow.

---

## Rationale

1. **Target audience constraint.** The brief explicitly states "workflow uses brms / cmdstanr formula interfaces only" and "no custom Stan programs." brms is the only library that provides a formula interface over Stan.
2. **Prior specification in R.** brms's `set_prior()` / `prior()` interface lets users specify and inspect priors as R objects, which can be stored in `wf_state` directly. This is essential for the audit trail.
3. **Parameterization transparency.** brms's centered-vs-non-centered choice is inspectable (`stancode()`), can be overridden (`control = list(...)`), and can be logged to `wf_state`.
4. **cmdstanr over rstan.** cmdstanr is the current recommendation from the Stan team: faster compilation, better CmdStan compatibility, and active maintenance. rstan 2.x has lagged on CmdStan updates.

---

## Consequences

- **Positive:** Users fit models with a lme4-style formula; no Stan experience required.
- **Positive:** Prior objects are R-native and serializable; `wf_state` can store them directly.
- **Positive:** `stancode()` makes the generated Stan program inspectable without requiring users to write it.
- **Negative:** brms's parameterization choices (centered hierarchical by default) are not always optimal. **Mitigation:** when divergences appear and a hierarchical model is in use, RBayesflow must recommend non-centered parameterization and make switching one function call (see DESIGN.md §6).
- **Negative:** brms occasionally lags behind the latest CmdStan features. Acceptable for this use case.
- **Watch:** If a future phase introduces custom Stan code (e.g., for user-defined likelihoods), revisit this ADR and add a `cmdstanr`-direct path gated behind `mode = "expert"`.
