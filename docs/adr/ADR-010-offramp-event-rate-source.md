# ADR-010: Off-Ramp Event Rate Source (Raw Data vs. Prior Predictive Simulation)

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** DESIGN.md §4 (Off-Ramp Logic), DESIGN.md §9 (ODQ-3), PLAN.md Phase 3

---

## Context

`assess_offramps()` is called at the end of Phase 1 (goal declaration) to
evaluate non-Bayesian alternatives before any model is fitted. For binary
outcomes, one of the key inputs is the observed event rate — the proportion of
outcome = 1 in the data. This is used to trigger warnings (e.g., rare event
rate < 5% → offer Firth penalized logistic) and to construct the one-sentence
stakes rationale shown to the user alongside each alternative.

The question is where this event rate comes from:

**Option A — Raw data (computed before any fitting)**
- `event_rate <- mean(data[[outcome_var]])`
- Available immediately at Phase 1, before any Stan compilation or sampling.
- Requires only that the user has loaded their data and declared their outcome
  variable.
- No additional computation; no sampling required.

**Option B — Prior predictive simulation**
- Event rate estimated from draws from the prior predictive distribution:
  `mean(posterior_predict(fit_prior_only))`.
- Requires a prior predictive fit (Phase 2) to have already completed.
- Provides an estimate of what the model *expects* rather than what the data
  show — useful for model criticism but not for initial off-ramp assessment,
  which is specifically a pre-fitting decision.

---

## Decision

**Option A: raw data is sufficient. `assess_offramps()` computes event rate
directly from the outcome variable before any fitting.**

```r
assess_offramps <- function(data, outcome_var, outcome_type, goal, n = nrow(data)) {
  event_rate <- if (outcome_type == "binary") mean(data[[outcome_var]]) else NA_real_
  # ... decision matrix logic ...
}
```

---

## Rationale

1. **Off-ramps are a pre-fitting decision.** The purpose of `assess_offramps()`
   is to present equally weighted alternatives *before* the user commits to a
   full Bayesian fit. Requiring a prior predictive fit first defeats this
   purpose — the user would have already begun the Bayesian path before being
   offered a way off it.
2. **Raw event rate is the correct input for the decision matrix.** The
   decision to offer Firth penalized logistic depends on whether the *data*
   have a rare event problem, not on whether the *prior* predicts rare events.
   A prior can be misspecified in ways that obscure a genuine rare-event
   problem; the raw data cannot.
3. **No computational cost.** `mean(data[[outcome_var]])` is instantaneous.
   Prior predictive simulation requires Stan compilation (30–90 seconds on
   first run) and sampling. Blocking Phase 1 on sampling is a friction the
   target audience will not tolerate for a step that is supposed to help them
   decide whether to fit at all.
4. **Edge case (completely unlabeled outcome) is handled.** If the outcome
   variable is not binary-coded (e.g., factors "yes"/"no" rather than 1/0),
   `assess_offramps()` converts via `as.integer(factor(data[[outcome_var]])) - 1L`
   with a warning. This covers the realistic range of user data without
   requiring prior simulation.

---

## Consequences

- **Positive:** Phase 1 runs instantly; no Stan compilation required before
  the user sees their off-ramp options.
- **Positive:** The event rate in the exit log is traceable to the raw data,
  making it auditable and reproducible without rerunning any sampling.
- **Positive:** No ordering dependency between Phase 1 and Phase 2; the user
  can assess off-ramps before specifying priors.
- **Negative:** The raw event rate does not account for imbalance introduced
  by stratified sampling or complex survey designs. For such cases,
  `assess_offramps()` emits a note: "Event rate computed from raw outcome
  variable; verify this reflects your sampling design." This is a warning,
  not a blocker.
- **Watch:** If a future phase introduces design-weighted analyses
  (e.g., complex survey data), revisit this ADR to add a `weights` argument
  to `assess_offramps()`.
