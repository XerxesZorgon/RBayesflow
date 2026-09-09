# Known-bad model for RBayesflow integration testing.
# Simulates a rare-event binary outcome (3% base rate) with a
# deliberately diffuse prior to force diagnostic failures.
# Expected outcome: Bernoulli rare-event warning; possible Rhat issues.
# Run from the project root: source("examples/bernoulli_rare_event/fit_rare_event.R")

library(brms)

set.seed(99)
n <- 150
x <- rnorm(n)
p <- plogis(-4 + 0.5 * x)   # ~3% baseline event rate
y <- rbinom(n, 1, p)
dat_bad <- data.frame(y = y, x = x)

fit_bad <- brm(
  formula = y ~ x,
  data    = dat_bad,
  family  = bernoulli(),
  prior   = prior(normal(0, 10), class = Intercept),
  backend = "cmdstanr",
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  refresh = 100
)

saveRDS(fit_bad, "examples/bernoulli_rare_event/fit_bad.rds")
cat("fit_bad saved to examples/bernoulli_rare_event/fit_bad.rds\n")
