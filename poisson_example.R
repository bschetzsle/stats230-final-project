source("bsl.R")

# symmetric proposal
rand_proposal <- function(old_theta) {
  bandwidth <- 3
  old_lambda <- old_theta$lambda

  list(lambda = abs(runif(1, old_lambda - bandwidth, old_lambda + bandwidth)))
}

prior_log_density <- function(theta)
  dgamma(theta$lambda, shape = 0.001, rate = 0.001, log = T)



true_lambda <- 30

set.seed(100)
y <- rpois(100, true_lambda)

result <- mcmc_bsl(
    y                 = y,
    prior_log_density = prior_log_density,
    rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
    rand_proposal     = rand_proposal,
    statistic         = function(y) mean(y),
    initial_theta     = list(lambda = 6),
    iterations        = 1000
)

result <- mcmc_bsl(
  y                 = y,
  prior_log_density = prior_log_density,
  rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
  rand_proposal     = rand_proposal,
  statistic         = function(y) max(y),
  initial_theta     = list(lambda = 6),
  iterations        = 1000
)

result <- mcmc_bsl(
  y                 = y,
  prior_log_density = prior_log_density,
  rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
  rand_proposal     = rand_proposal,
  statistic         = function(y) lm(sort(y) ~ poly(seq_along(y), degree = 3))$coef,
  initial_theta     = list(lambda = 6),
  iterations        = 1000
)