ricker_model <- function(log_r, phi, sigma, N_0)
    list(log_r = log_r, phi = phi, sigma = sigma, N_0 = N_0)

rand_ricker <- function(n, model) {
    e <- rnorm(n, 0, model$sigma^2)
    names(e) <- 0:(n - 1) # Note N_t uses e_{t-1} and N_{t-1}

    N <- ricker_N(e, model)
    names(N) <- 0:n

    y <- rpois(length(N), model$phi * N)
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