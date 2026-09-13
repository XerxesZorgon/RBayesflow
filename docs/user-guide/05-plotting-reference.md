# RBayesflow — Plotting Reference

This reference covers every plot produced in the RBayesflow workflow: which phase it belongs to, the exact R call to produce it, what a good result looks like, what a bad result looks like, and what to do about it. Plots are listed in phase order.

---

## Phase 1 — Data Exploration

### esquisse::esquisser()

**Purpose:** Interactive drag-and-drop interface for building exploratory ggplot2 plots of your raw data before any modelling.

**Call:**
```r
esquisse::esquisser(data)
```

**What to look for:** Use this to understand your outcome variable's distribution (skew, outliers, zero-inflation), the relationship between predictors and outcome, and whether any data-quality issues are visible (implausible values, unexpected clusters). When you find a plot you want to keep, click **Export code** — the generated ggplot2 code can be pasted directly into the Phase 1 template cell.

**Learn-mode note:** No difference between modes — esquisse is always interactive and always produces the same interface.

---

## Phase 2 — Prior Predictive Simulation

### bayesplot::ppc_dens_overlay() — prior predictive

**Purpose:** Overlay density curves of prior predictive draws on top of the observed outcome distribution, verifying that priors produce plausible data before fitting.

**Call:**
```r
bayesplot::ppc_dens_overlay(y, prior_pred_draws[1:50, ])
```

**What to look for:** Good result — the prior predictive curves (light blue) roughly bracket the observed data density (dark line) without being wildly wider or narrower. Bad result — curves extend far outside the observed data range (e.g. predicting negative counts or probabilities above 1), or are so tight they exclude the observed data entirely. Either signals priors that need revision before fitting. Return to Phase 2 and adjust `brms::prior()` objects.

**Learn-mode note:** In `mode = "learn"` a callout prompt asks you to describe in one sentence why the prior draws are or are not plausible before moving on.

---

## Phase 4 — MCMC Diagnostics

### bayesplot::mcmc_trace()

**Purpose:** Display sampler chain trajectories over iterations, revealing whether all chains mixed well and explored the same region of parameter space.

**Call:**
```r
bayesplot::mcmc_trace(fit)
```

**What to look for:** Good result — all chains (different colours) overlap throughout the plot with no systematic separation, producing a "fuzzy caterpillar" appearance. Bad result — chains drift apart, one chain occupies a different region, or a chain gets stuck (flat horizontal segments). Any of these indicate poor mixing; the corresponding Rhat value will be above 1.01. Try `refit_noncentered(wf)` for hierarchical models, or revise priors.

**Learn-mode note:** In `mode = "learn"` a Socratic prompt asks what the trace plots reveal about chain mixing before the Rhat summary is shown.

### bayesplot::mcmc_rhat()

**Purpose:** Display the distribution of Rhat (potential scale reduction factor) values across all model parameters. Rhat measures how well the chains agree; values close to 1.00 indicate convergence.

**Call:**
```r
bayesplot::mcmc_rhat(brms::rhat(fit))
```

**What to look for:** Good result — all bars fall to the left of the 1.01 threshold line. Bad result — one or more bars exceed 1.01 (shown in a warning colour). Parameters with Rhat > 1.01 have not converged; do not interpret their posteriors. Increase iterations, reparameterize, or revise the model structure.

**Learn-mode note:** In `mode = "learn"` a Socratic prompt asks what an Rhat of 1.00 versus 1.05 implies about the reliability of posterior summaries.

### bayesplot::mcmc_neff()

**Purpose:** Display the distribution of effective sample size (ESS) ratios across all parameters. ESS measures how many independent draws the chains are equivalent to after accounting for autocorrelation.

**Call:**
```r
bayesplot::mcmc_neff(brms::neff_ratio(fit))
```

**What to look for:** Good result — all bars are above the 0.1 threshold (ESS ratio ≥ 10% of total draws). Bad result — bars below 0.1 indicate high autocorrelation; posterior summaries for those parameters are unreliable. Increase `iter_sampling`, use `adapt_delta`, or reparameterize.

**Learn-mode note:** In `mode = "learn"` a Socratic prompt asks how a low ESS ratio would affect the width of a credible interval.

### bayesplot::mcmc_pairs()

**Purpose:** Display pairwise scatter plots of parameter samples, with divergent transitions highlighted in red. Only relevant when `wf$diagnostics$n_divergences > 0`.

**Call:**
```r
bayesplot::mcmc_pairs(fit, np = nuts_params(fit))
```

**What to look for:** Good result (no divergences) — funnel shapes or banana curves are normal for correlated parameters; no red dots. Bad result — red dots (divergent transitions) cluster in a narrow region of one or more pairwise plots, indicating the sampler struggled with the posterior geometry in that region. For hierarchical models, try `refit_noncentered(wf)`. For other models, increase `adapt_delta` toward 0.99.

**Learn-mode note:** In `mode = "learn"` this plot is shown only when divergences are present, with a note explaining what red dots represent.

---

## Phase 5 — Posterior Predictive Checks

### bayesplot::ppc_dens_overlay() — posterior predictive

**Purpose:** Overlay density curves of posterior predictive draws on top of the observed outcome distribution, checking whether the fitted model reproduces the data-generating structure.

**Call:**
```r
bayesplot::ppc_dens_overlay(y, posterior_predict(fit)[1:50, ])
```

**What to look for:** Good result — the posterior predictive curves (light blue) closely track the observed data density (dark line). Bad result — systematic discrepancies: the model consistently misses the mode, fails to capture skew, or cannot reproduce zero-inflation. These indicate model misspecification; consider changing the likelihood family or adding a predictor.

**Learn-mode note:** In `mode = "learn"` a Socratic prompt asks you to describe any visible discrepancy between the curves before showing the test statistic results.

### bayesplot::ppc_stat()

**Purpose:** Compare a test statistic (e.g. mean, standard deviation, proportion of zeros) computed on observed data to the distribution of that statistic across posterior predictive draws.

**Call:**
```r
bayesplot::ppc_stat(y, posterior_predict(fit), stat = "mean")
```

**What to look for:** Good result — the observed statistic (vertical line) falls near the centre of the posterior predictive distribution (histogram). Bad result — the observed statistic falls in the tail (p-value < 0.05 or > 0.95), indicating the model systematically over- or under-predicts that aspect of the data. Try `stat = "sd"`, `stat = "max"`, or a custom function to check other aspects.

**Learn-mode note:** No difference between modes for this plot.

### tidybayes::add_epred_draws()

**Purpose:** Compute and plot the expected value of the posterior predictive distribution (the fitted mean) as a ribbon over your predictor range, giving a tidy alternative to `pp_check()` for continuous outcomes.

**Call:**
```r
library(tidybayes)
library(ggplot2)
data |>
  tidybayes::add_epred_draws(fit) |>
  ggplot(aes(x = predictor, y = outcome)) +
  tidybayes::stat_lineribbon(aes(y = .epred), alpha = 0.25) +
  geom_point(data = data)
```

**What to look for:** Good result — the posterior ribbon covers the observed data points without being implausibly wide or narrow. Use this plot instead of `ppc_dens_overlay()` when you want to visualise the fitted relationship across a continuous predictor rather than the marginal outcome distribution.

**Learn-mode note:** In `mode = "learn"` a callout note explains the difference between `.epred` (expected prediction, averaging over parameter uncertainty) and `.prediction` (individual observation prediction, including residual variance).

---

## Phase 6 — Model Comparison

### loo::loo_compare()

**Purpose:** Compare two or more fitted models on expected log-predictive density (ELPD), estimated by leave-one-out cross-validation (LOO-CV). Higher ELPD is better.

**Call:**
```r
loo1 <- loo(fit1)
loo2 <- loo(fit2)
loo::loo_compare(loo1, loo2)
```

**What to look for:** The output table ranks models from best to worst ELPD. The `elpd_diff` column shows the difference relative to the best model; `se_diff` is the standard error of that difference. A difference is meaningful when `|elpd_diff| > 2 * se_diff`. If the best model's advantage is less than twice its standard error, the models are not distinguishable on predictive accuracy and you should prefer the simpler one. Check `wf$diagnostics$passed == TRUE` for both models before comparing — LOO estimates from poorly-mixing chains are unreliable.

**Learn-mode note:** In `mode = "learn"` a Socratic prompt asks you to interpret the ELPD difference before the table is displayed.

---

## Phase 7 — Reporting

All plots in the Phase 7 report (`templates/bayesflow_report.qmd`) are reproduced from saved objects — no new plots are generated at render time. The report calls `readRDS(params$fit_path)` and reads all summaries and plot data from `wf_state` fields. If a plot looks different from what you saw during the analysis, check that `params$fit_path` and `params$wf_path` point to the correct saved files.
