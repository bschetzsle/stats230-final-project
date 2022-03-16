source("bsl.R")

statistic <- function(y) {
  mean(y)
  # max(y)
  # c(min(y), max(y))
}

# symmetric proposal
rand_proposal <- function(old_theta) {
  bandwidth <- 3
  old_lambda <- old_theta$lambda

  list(lambda = abs(runif(1, old_lambda - bandwidth, old_lambda + bandwidth)))
}

prior_log_density <- function(theta)
  dgamma(theta$lambda, shape = 15, rate = 0.5, log = T)



true_lambda <- 30

result <- mcmc_bsl(
    y                 = rpois(100, true_lambda),
    prior_log_density = prior_log_density,
    rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
    rand_proposal     = rand_proposal,
    statistic         = statistic,
    initial_theta     = list(lambda = 6.29)
)

plot_param(result, lambda, true_lambda) / plot_log_prob(result)
