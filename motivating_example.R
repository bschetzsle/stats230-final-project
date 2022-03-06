# Recreating https://www.nature.com/articles/nature09319#ref-CR2

library(tidyverse)
library(patchwork)
source("ricker.R")

model <- ricker_model(
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
