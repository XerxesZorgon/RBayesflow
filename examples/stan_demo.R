# Stan smoke test — verifies CmdStan compiles and samples correctly.
# Run with: source("examples/stan_demo.R")

library(cmdstanr)

# Inline Stan model: normal-normal conjugate (estimate mu, known sigma)
stan_code <- "
data {
  int<lower=0> N;
  array[N] real y;
  real<lower=0> sigma;
}
parameters {
  real mu;
}
model {
  mu ~ normal(0, 10);       // prior
  y ~ normal(mu, sigma);    // likelihood
}
"

# Write model to temp file and compile
model_file <- tempfile(fileext = ".stan")
writeLines(stan_code, model_file)
mod <- cmdstan_model(model_file)

# Sample
fit <- mod$sample(
  data    = list(N = 8, y = c(2.1, 1.9, 2.4, 1.8, 2.3, 2.0, 2.2, 1.7),
                 sigma = 0.5),
  chains        = 2,
  iter_warmup   = 200,
  iter_sampling = 200,
  refresh       = 0,
  seed          = 42
)

# Print summary
print(fit$summary())

# Verify Rhat for mu
rhat_mu <- fit$summary("mu")$rhat
cat("mu rhat:", rhat_mu, "\n")
stopifnot("Stan smoke test failed: rhat >= 1.05" = rhat_mu < 1.05)
cat("Stan smoke test PASSED\n")
