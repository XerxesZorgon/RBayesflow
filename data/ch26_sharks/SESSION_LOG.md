# Session Log
**Updated:** 2026-09-22
**Active skill:** software-project (case study series — Ch 26)
**Last confirmed state:** Ch 26 complete — article written, all figures saved, all three models converged — no pending tasks

## What happened this session

Wrote and ran the full Ch 26 (white sharks HMM) analysis from scratch: three Stan files, the R script, and all 12 figures. The first run hit Rhat 1.73 / ESS 6 on Models 1 and 2 (label-switching); fixed by switching from a pooled to state-specific prior on `mu` and tightening `sigma ~ exponential(5)`. Two hotfixes went through Antigravity to correct a `$pathfinder()` init-argument bug. After those fixes, all three models converged cleanly (Rhat ≤ 1.002, ESS ≥ 407). The MDX article was then written from the actual SVG figure data and saved.

## Decisions made (not yet in an ADR)

- **State-specific mu priors for HMMs:** For any HMM with `positive_ordered` means, a pooled prior on `mu` is insufficient to prevent label-switching across chains. Use state-specific Normal priors anchored to known biological or domain values. Add this as a standing note to the Stan HMM template.
- **Pathfinder init with simplex parameters:** `$pathfinder()` does not accept per-path init lists (`init = list_of_lists`). Pass a single zero-argument function. More importantly, simplex parameters declared as `array[K] simplex[K]` in Stan (e.g. `tpm_raw`) cause Pathfinder to fail silently when supplied as `list(c(...), c(...))` in the init. Either omit simplex parameters from the init or use `init = 0`. The `tryCatch` + NULL-guard fallback pattern is now the standard wrapper for all `$pathfinder()` calls.
- **Tighter sigma prior for HMMs:** `student_t(3, 0, 0.5)` on step-length SD allows heavy tails that widen the Gamma distributions until states overlap. Use `exponential(5)` (mean 0.2 km) for animal movement data at this spatial scale.

## Blocked on / open question

- LOO not computed for any of the three models. The per-observation log-likelihood from the forward algorithm is not the individual-observation contribution needed for PSIS-LOO. A correct LOO would require leave-one-observation-out re-marginalisation at O(T²) cost. Noted in the article; left unresolved by design.
- Model 3 ESS_bulk_min = 407 (single chain). Adequate for the ribbon visualisations but worth revisiting if the article is extended to include formal uncertainty intervals on covariate effects.
- The article has not yet gone through the /humanizer pass or the editor checklist.

## Next action

Run `/humanizer` on `ch26_sharks_article.mdx`, then run the editor checklist skill against it before publication. File is at:
`C:\Users\johnx\Documents\WildPeaches\Projects\RBayesflow\data\ch26_sharks\ch26_sharks_article.mdx`

## Lessons learned this session (Ch 26 — to add to the code rubric)

**Lesson 6 — State-specific priors for ordered HMM means:**
`positive_ordered` on `mu` is a necessary but not sufficient condition for avoiding label-switching. When the two Gamma distributions overlap substantially (which happens when `sigma` is large relative to the gap between means), chains can assign different parameter values to the two states even though `mu[1] < mu[2]` holds within each chain. The fix is state-specific Normal priors anchored well apart: `mu[1] ~ normal(0.08, 0.05)`, `mu[2] ~ normal(0.35, 0.12)`. Rhat went from 1.73 to 1.002.

**Lesson 7 — Pathfinder's `init` argument does not accept simplex parameters as lists of vectors:**
`array[K] simplex[K]` parameters (e.g. `tpm_raw`) passed as `list(c(0.95, 0.05), c(0.05, 0.95))` cause `$pathfinder()` to fail silently and return NULL. Standard pattern going forward: wrap every `$pathfinder()` call in `tryCatch(..., error = function(e) NULL)` and guard the downstream draw-extraction block with `if (!is.null(pf_result))`. Pass simplex parameters via a zero-arg function that either omits them or uses `init = 0`.

**Lesson 8 — Sigma prior width controls state separability:**
A heavy-tailed prior on step-length SD (e.g. `student_t(3, 0, 0.5)`) can produce posterior samples where both Gamma distributions are so wide they fully overlap, destroying the state contrast the ordered-means constraint was meant to enforce. Use `exponential(λ)` with mean `1/λ` chosen so the 95th percentile is well below the gap between the two state means.
