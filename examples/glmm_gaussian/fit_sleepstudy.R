# Known-good model for RBayesflow integration testing.
# Fits a Gaussian GLMM on the sleepstudy dataset (lme4).
# Expected outcome: Rhat < 1.01, ESS > 400, 0 divergences.
# Run from the project root: source("examples/glmm_gaussian/fit_sleepstudy.R")

library(brms)
library(lme4)

fit_good <- brm(
  formula  = Reaction ~ Days + (Days | Subject),
  data     = lme4::sleepstudy,
  family   = gaussian(),
  prior    = c(
    prior(normal(250, 50), class = Intercept),
    prior(normal(10, 5),   class = b, coef = Days),
    prior(exponential(1),  class = sigma)
  ),
  backend  = "cmdstanr",
  seed     = 42,
  chains   = 4,
  iter     = 2000,
  warmup   = 1000,
  refresh  = 100
)

saveRDS(fit_good, "examples/glmm_gaussian/fit_good.rds")
cat("fit_good saved to examples/glmm_gaussian/fit_good.rds\n")
