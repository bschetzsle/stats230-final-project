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
  dgamma(theta$lambda, shape = 0.001, rate = 0.001, log = T)



true_lambda <- 30

set.seed(100)
y = rpois(100, true_lambda)

result <- mcmc_bsl(
    y                 = y,
    prior_log_density = prior_log_density,
    rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
    rand_proposal     = rand_proposal,
    statistic         = function(y) mean(y),
    initial_theta     = list(lambda = 6),
    iterations        = 1000
)

p = plot_param(result, lambda, true_lambda) / plot_log_prob(result)
ggsave("./plots/Poisson_mean.png", plot = p, width = 30, height = 20, units = "cm")

#I need help with this plot
plot(density(result$lambda[100:1000]), main="Posterior Density", xlab="Lambda")
curve(dgamma(x,shape=0.001+sum(y), rate=100.001), from=25, to=35, col="red", add=T)
legend(31,0.8,legend=c("Estimated","True"),col=c("black","red"),lty=1)


temp = density(result$lambda)
p1 <- ggplot(cbind(temp$x, temp$y), aes(x, y)) +
  geom_line() +
  geom_hline(yintercept = true_value, color = "dodgerblue4")



result <- mcmc_bsl(
  y                 = y,
  prior_log_density = prior_log_density,
  rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
  rand_proposal     = rand_proposal,
  statistic         = function(y) max(y),
  initial_theta     = list(lambda = 6),
  iterations = 1000
)

p = plot_param(result, lambda, true_lambda) / plot_log_prob(result)
ggsave("./plots/Poisson_max.png", plot = p, width = 30, height = 20, units = "cm")

#why are the estimates of lambda off? it turns out y had an unusually large max
plot(density(replicate(1000, max(rpois(100,30)))))
abline(v=max(y))


result <- mcmc_bsl(
  y                 = y,
  prior_log_density = prior_log_density,
  rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
  rand_proposal     = rand_proposal,
  statistic         = function(y) c(min(y), max(y)),
  initial_theta     = list(lambda = 6),
  iterations = 1000
)

p = plot_param(result, lambda, true_lambda) / plot_log_prob(result)
ggsave("./plots/Poisson_min_max.png", plot = p, width = 30, height = 20, units = "cm")

plot(density(replicate(1000, min(rpois(100,30)))))
abline(v=min(y))
plot(density(replicate(1000, max(rpois(100,30)))))
abline(v=max(y))


result <- mcmc_bsl(
  y                 = y,
  prior_log_density = prior_log_density,
  rand_likelihood   = function(n, theta) rpois(n, theta$lambda),
  rand_proposal     = rand_proposal,
  statistic         = function(y) lm(sort(y) ~ poly(1:length(y), degree=3))$coef,
  initial_theta     = list(lambda = 6),
  iterations        = 1000
)

ggsave("./plots/Poisson_sorted_Y.png", plot(sort(y), ylab="Count", main="Sorted Observations"),
       width = 30, height = 20, units = "cm")

p = plot_param(result, lambda, true_lambda) / plot_log_prob(result)
ggsave("./plots/Poisson_poly_regress.png", plot = p, width = 30, height = 20, units = "cm")

plot(density(replicate(1000, min(rpois(100,30)))))

temp = replicate(1000, lm(sort(rpois(100,30)) ~ poly(1:length(y), degree=3))$coef)

plot(density(temp[,1]))
abline(v=lm(sort(y) ~ poly(1:length(y), degree=3))$coef[1])
plot(density(temp[,2]))
abline(v=lm(sort(y) ~ poly(1:length(y), degree=3))$coef[2])
plot(density(temp[,3]))
abline(v=lm(sort(y) ~ poly(1:length(y), degree=3))$coef[3])
plot(density(temp[,4]))
abline(v=lm(sort(y) ~ poly(1:length(y), degree=3))$coef[4])
