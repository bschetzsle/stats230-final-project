# Recreating https://www.nature.com/articles/nature09319#ref-CR2

source("bsl.R")

rand_ricker <- function(n, theta) {
    rand_N_t <- function(prev_N)
        rlnorm(1, log(prev_N) + theta$log_r - prev_N, theta$sigma)

    N <- Reduce(
        function(prev_N, ...) rand_N_t(prev_N),
        c(theta$N_0, seq_len(n)),
        accumulate = T
    )

    y <- rpois(length(N), theta$phi * N)
    names(y) <- paste0("t", 0:n)

    y
}

prior_log_density <- function(theta) {
    dunif(theta$log_r, -20, 20, log = T) +
    dunif(theta$phi,     0, 50, log = T) +
    dunif(theta$sigma,   0,  3, log = T)

}

# Symmetric
rand_proposal <- function(old_theta) {
    #noise <- function(bandwidth)
    #    runif(1, -bandwidth, bandwidth)

    #list(
    #    log_r = old_theta$log_r + noise(1),
    #    phi   = abs(old_theta$phi + noise(0.5)),
    #    sigma = abs(old_theta$sigma + noise(0.05)),
    #    N_0   = old_theta$N_0
    #)

    mu <- with(old_theta, c(log_r, phi, sigma))

    cov_mat <- matrix(
        c(2.25, -5.25, -1.35, -5.25, 25, 3, -1.35, 3, 2.25),
        nrow = 3
    ) * 10^(-2)

    new_theta <- rmvnorm(1, mu, cov_mat)[1, ]

    # They automatically use the same theta as before if sigma comes out negative... hm.
    if (new_theta[3] <= 0)
        return(old_theta)

    list(
        log_r = new_theta[1],
        phi   = new_theta[2],
        sigma = new_theta[3],
        N_0   = old_theta$N_0
    )
}

statistic <- function(y) {
    # To avoid the error "degree must be lower than number of unique values"
    # when fitting the autoregression, we add some noise to the covariates
    noise <- rnorm(length(y) - 1, 0, 0.001)

    c(
        mean(y),
        sum(y == 0),
        lm(diff(y) ~ poly(y[-1] + noise, degree = 3))$coef,
        lm(I(y[-1]^0.3) ~ I(y[-50]^0.3) + I(y[-50]^0.6))$coef,
        acf(y, lag.max = 5, plot = F, type = "covariance")$acf[, 1, 1]
    )
}





true_theta <- list(
    log_r = 3.8,
    sigma = 0.3,
    phi = 10,
    N_0 = 1
)

y <- rand_ricker(50, true_theta)

result <- mcmc_bsl(
    y                 = y,
    prior_log_density = prior_log_density,
    rand_likelihood   = rand_ricker,
    rand_proposal     = rand_proposal,
    statistic         = statistic,
    initial_theta     = true_theta, # this is what they do
    simulations       = 50,
    iterations        = 500000
)

saveRDS(
    list(y = y, result = result),
    "mcmc_bsl.rds"
)



plot_param(result, log_r, true_theta$log_r) /
plot_param(result, phi, true_theta$phi) /
plot_param(result, sigma, true_theta$sigma) /
plot_log_prob(result)










iterations = length(result$theta)
new_simulation = rand_ricker(n = 50, model=ricker_model(
    log_r = result$theta[[iterations]]$log_r,
    sigma = result$theta[[iterations]]$sigma,
    phi = result$theta[[iterations]]$phi,
    N_0 = result$theta[[iterations]]$N_0
))
plot(simulation$y, type="l")
lines(new_simulation$y, col="red")
