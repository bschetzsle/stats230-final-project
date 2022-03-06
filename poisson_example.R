source("bsl.R")

statistic <- function(X) {
  return(mean(X))
  #return(max(X))
  #return(c(min(X), max(X)))
}

rand_proposal <- function(old_theta) {
  bandwidth <- 3

  old_lambda <- old_theta$lambda
  lambda <- abs(runif(1, old_lambda - bandwidth, old_lambda + bandwidth))

  list(lambda = lambda)
}

model <- list(lambda = 30)

result <- MCMC_BSL(
    y                = rpois(100, model$lambda),
    prior_density    = function(theta) dgamma(theta$lambda, shape = 15, rate = 0.5),
    proposal_density = function(old_theta, new_theta) 1, # it's symmetric
    rand_proposal    = rand_proposal,
    rand_model       = function(n, theta) rpois(n, theta$lambda),
    statistic        = statistic,
    initial_theta    = list(lambda = 6.29),
    iterations       = 1000
)

view_results(result, model)
