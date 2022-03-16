
rand_ricker <- function(n, theta) {
    e <- rnorm(n, 0, theta$sigma^2)
    names(e) <- 0:(n - 1) # Note N_t uses e_{t-1} and N_{t-1}

    N <- ricker_N(e, theta)
    names(N) <- 0:n

    y <- rpois(length(N), theta$phi * N)
    names(y) <- names(N)

    list(
        time = as.numeric(names(y)),
        e = e,
        y = y,
        log_lik = ricker_log_lik(model, y, e)
    )
}

ricker_log_lik <- function(model, y, e) {
    # Recreate what the N's would be, given these e's and parameters
    N <- ricker_N(e, model)

    # The likelihoods are so small, we can't log them ourselves...
    # must use hardcoded log-lik with `log = T`
    sum(dpois(y, model$phi * N, log = T)) +
        sum(dnorm(e, 0, model$sigma, log = T))
}

# Uses N_0 and e_0 to e_{n-1} to generate N_1 to N_n.
# Returns N_0, N_1, ..., N_n
ricker_N <- function(e, model) {
    Reduce(
        \(N_t, e_t) N_t * exp(model$log_r - N_t + e_t),
        c(model$N_0, e),
        accumulate = T
    )
}


model <- list(
    log_r = 3.8,
    sigma = 0.3,
    phi = 10,
    N_0 = 1
)

simulation <- rand_ricker(n = 50, model)

# This is our observed data
fig_1a <- ggplot(mapping = aes(simulation$time, simulation$y)) +
    geom_line() +
    geom_point()

# We only observe y, not e. e just is noise that adds variation to the y.
# To do inference on this model, we need a likelihood marginalized over e,
# ie. we'd have to integrate out e_1, ..., e_n. That's hard analytically, and apparently
# numerically too. If we wanted to marginalize out e_1 alone, take a look at the function
# we'd have to integrate over:

e1 <- seq(-2, 2, 0.005)
log_lik <- sapply(e1,
    \(e1) ricker_log_lik(
        model,
        simulation$y,
        c(e1, simulation$e[-1])
    )
)

fig_1b <- ggplot(mapping = aes(e1, log_lik)) +
    geom_line()

# Compare this to, say, integrating out sigma in a N(mu, sigma^2) model:
data <- rnorm(50, 0, 2)
sigma <- seq(1, 5, 0.005)
log_lik <- sapply(sigma,
    \(sigma) sum(dnorm(data, 0, sigma, log = T))
)

ggplot(mapping = aes(sigma, log_lik)) +
    geom_line()

# Smooth. Easy.
# noise bad


# Something something sample an r from this, bayes is hard too
log_r <- seq(3, 5, 0.005)
log_lik <- sapply(log_r,
    \(log_r) ricker_log_lik(
        modifyList(model, list(log_r = log_r)),
        simulation$y,
        simulation$e
    )
)

fig_1c <- ggplot(mapping = aes(log_r, log_lik)) +
    geom_line()