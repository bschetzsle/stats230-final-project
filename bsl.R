library(mvtnorm)
library(patchwork)
library(tidyverse)

#this function will generate samples from the posterior P(theta|Y) without computing P(Y|theta)
MCMC_BSL <- function(y,
                     prior_density,
                     proposal_density,
                     rand_proposal,
                     rand_model,
                     statistic,
                     initial_theta,
                     iterations = 1000,
                     simulations = 100) {
    #simulations is n in the papers

    s_y <- statistic(y)

    #generate a structure to keep track of the samples
    result <- list(theta = list(), log_joint_prob = seq_len(iterations))

    #simulate x from starting theta0 and calculate statistic
    theta <- initial_theta
    s <- Reduce(rbind, replicate(simulations, statistic(rand_model(length(y), theta)), simplify = F))
    mu <- apply(s, 2, mean)
    Sigma <- t(s - mu) %*% (s - mu) / (simulations - 1)

    for (i in 1:iterations) {
        new_theta <- rand_proposal(theta)

        s <- Reduce(rbind, replicate(simulations, statistic(rand_model(length(y), new_theta)), simplify = F))
        new_mu <- apply(s, 2, mean)
        new_Sigma <- t(s - new_mu) %*% (s - new_mu) / (simulations - 1)

        log_MH_ratio <- (dmvnorm(s_y, new_mu, new_Sigma, log = T) +
                             prior_density(new_theta, log = T) +
                             log(proposal_density(new_theta, theta))) -
            (dmvnorm(s_y, mu, Sigma, log = T) +
                 prior_density(theta, log = T) +
                 log(proposal_density(theta, new_theta)))
        r <- min(1, exp(log_MH_ratio))


        if (is.na(r) || runif(1) < r) {
            theta <- result$theta[[i]] <- new_theta
            mu    <- new_mu
            Sigma <- new_Sigma
        }

        else {
            result$theta[[i]] <- theta
        }

        result$log_joint_prob[i] <- dmvnorm(s_y, mu, Sigma, log = T) + prior_density(theta, log = T)
    }

    result
}


view_results <- function(result, model) {
    plot_param <- function(param) {
        samples <- lapply(result$theta, \(theta) theta[[param]]) %>% unlist()

        p1 <- ggplot(mapping = aes(seq_along(samples), samples)) +
            geom_line() +
            geom_hline(yintercept = model[[param]], color = "dodgerblue4") +
            ylab(param)

        p2 <- ggplot(mapping = aes(samples)) +
            geom_histogram() +
            geom_vline(xintercept = model[[param]], color = "dodgerblue4") +
            xlab(param)

        p1 + p2
    }

    Map(plot_param, names(model)) %>%
    Reduce(`/`, .) / (
        (ggplot(mapping = aes(seq_along(result$log_joint_prob), result$log_joint_prob)) + geom_line()) + ggplot()
    )
}
