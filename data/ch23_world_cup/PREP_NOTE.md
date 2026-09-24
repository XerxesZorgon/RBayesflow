# Ch 23 Prep Note — "Debugging a model: World Cup football"

**Chapter:** 23
**Book title (exact):** Debugging a model: World Cup football
**Book URL:** https://avehtari.github.io/Bayesian-Workflow/world_cup/world_cup.html
**GitHub source:** https://github.com/avehtari/Bayesian-Workflow/tree/master/world_cup/
**Target script:** `data/ch23_world_cup/world_cup_analysis.R`

---

## What this chapter is about

Ch 23 is the debugging chapter. The central theme is what to do when a model
produces wrong or implausible results — not wrong in the "diagnostics failed"
sense, but wrong in the "the posterior doesn't match what we know about
football" sense. The chapter walks through a model of 2014 World Cup match
outcomes, shows it produces absurd team-strength estimates, and traces the
failure to a prior that is too vague. The fix is prior tightening, not
structural model change.

The key concepts are:
- **Prior sensitivity as a debugging tool**, not just a robustness check.
- **Domain knowledge as a prior constraint**: football teams can't be 50×
  stronger than average; a prior that allows this is wrong.
- **Posterior retrodiction** (checking whether the fitted model can reproduce
  the observed data, not just predict future data) as the diagnostic that
  reveals the bug.
- **The difference between "the sampler ran fine" and "the model is right"**:
  Rhat and ESS are clean throughout, so the usual diagnostic gate doesn't
  catch this failure.

---

## What to expect from the models

The chapter builds two or three versions of a Poisson model for goals scored
in each match, with team-strength parameters for attack and defence. Expect:

- Outcome variable: goals scored by each team in each match (count data,
  likely Poisson).
- Team-strength parameters: one attack and one defence parameter per team,
  or a single net-strength parameter. The exact parameterisation depends on
  the chapter; check the Stan file.
- The broken model uses a vague prior (e.g., `normal(0, 10)`) on team
  strengths; the fixed model uses a tighter prior informed by the realistic
  range of match scores.
- The comparison between models is via posterior retrodiction plots, not
  LOO-CV — this is a debugging chapter, not a model-selection chapter.

---

## Key differences from Ch 22 to anticipate

| Aspect | Ch 22 (cat adoptions) | Ch 23 (World Cup) |
|---|---|---|
| Stan model type | Geometric survival, imputation | Poisson goals model |
| brms usable? | No | Possibly yes, or Stan-native |
| LOO-CV | Not run | Probably not (debugging focus) |
| Primary diagnostic | K-M posterior overlay | Posterior retrodiction of match scores |
| "Bug" being demonstrated | Censoring ignored → biased p | Vague prior → absurd team strengths |
| Fix | Add censoring likelihood term | Tighten prior using domain knowledge |

---

## Before writing any code — download checklist

Per the code rubric §2, fetch the chapter page and download all assets before
writing `world_cup_analysis.R`.

- [ ] Fetch `https://avehtari.github.io/Bayesian-Workflow/world_cup/world_cup.html`
- [ ] Identify every `.stan` file shown via `print_stan_file()` calls
- [ ] Identify the data file (likely a CSV of 2014 World Cup results)
- [ ] Write each Stan file verbatim to `data/ch23_world_cup/`
- [ ] Cache the data file locally (download-on-first-run block in the script,
      or write directly if the file is small)
- [ ] List the folder and confirm every file referenced by the script is present

---

## Watch-outs from Ch 22

1. **Figure interpretation rule**: read the rendered figure before writing
   the caption. Ch 22 had a wrong caption on Fig 22.6 (posteriors shifted
   left, not sitting on the true-value lines). For this chapter, the key
   comparison figures are the two retrodiction plots — check that the broken
   model's predictions are genuinely wider/more extreme than the data before
   writing "the model predicts absurd scores."

2. **LOO absence**: if LOO is not run, log the reason explicitly with
   `log_result("LOO-CV: not performed. Reason: ...")` before setting
   `wf$loo_complete <- FALSE`.

3. **Prior shrinkage note pattern**: if any recovery check shows a posterior
   shifted from a true value, add a NOTE in `results.txt` explaining whether
   the shift is expected or a bug.

4. **Verify `adoptions_varying.stan`**: the Ch 22 M5 Stan code was
   reconstructed from the chapter page rather than copied from GitHub.
   Before starting Ch 23, open that file and confirm it looks right.

---

## Style note for the article

The article's opening hook should anchor on a human question, not a
statistical one. Something like: "How do you know when your model is
lying to you?" The chapter theme — discovering that a model with clean
diagnostics can still be completely wrong — is the strongest pedagogical
point in Part 4. The article should make that tension vivid before
explaining the prior-tightening fix.

The "How the Models Evolved" table for this chapter tracks prior evolution,
not model structure: the structural model is the same throughout; what
changes is the prior width on team strengths and the evidence that each
prior choice produces.
