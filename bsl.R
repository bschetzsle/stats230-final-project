library(mvtnorm)
library(tidyverse)
library(patchwork)

# Implements MCMC BSL as written in Appendix B.
# Note that `simulations` here is n in the paper
mcmc_bsl <- function(y,
                     prior_log_density,
                     rand_likelihood,
                     rand_proposal,
                     statistic,
                     initial_theta,
                     iterations = 1000,
                     simulations = 100) {
    # Save timestamp for checkpoint filenames
    start_time <- as.numeric(Sys.time())
    # First, some helper functions

    # Monte Carlo approximation of the synthetic likelihood mu and sigma
    estimate_synth_lik_params <- function(theta) {
      s <- replicate(
        simulations,
        statistic(rand_likelihood(length(y), theta)),
        simplify = F
      )

      s <- Reduce(rbind, s)

      mu <- apply(s, 2, mean)


      mu_matrix <- matrix(mu, nrow = simulations, ncol = length(mu), byrow = T)
      sigma <- t(s - mu_matrix) %*% (s - mu_matrix) / (simulations - 1)

      list(mu = mu, sigma = sigma)
    }

    # Evaluate the log of the synthetic likelihood
    synthetic_log_lik <- function(s, params) {
      dmvnorm(s, params$mu, params$sigma, log = T)
    }

    # End helper functions

    s_y <- statistic(y)

    result <- list(
      theta          = list(),
      log_joint_prob = seq_len(iterations)
    )

    theta <- initial_theta
    lik_params <- estimate_synth_lik_params(theta)

    for (i in 1:iterations) {
        if (i %% 1000 == 0) {
          message("Saving checkpoint at iteration ", i)
          filename <- paste0("mcmc_bsl_", start_time, "_iteration_", i, ".rds")
          saveRDS(result, filename)
        }

        new_theta <- rand_proposal(theta)
        new_lik_params <- estimate_synth_lik_params(new_theta)

        log_mh_ratio <- synthetic_log_lik(s_y, new_lik_params) +
                        prior_log_density(new_theta) -
                        synthetic_log_lik(s_y, lik_params) -
                        prior_log_density(theta)

        r <- min(1, exp(log_mh_ratio))


        if (is.na(r) || runif(1) < r) {
            theta <- new_theta
            lik_params <- new_lik_params
        }

        result$theta[[i]] <- theta
        result$log_joint_prob[i] <- synthetic_log_lik(s_y, lik_params)
                                    prior_log_density(theta)
    }

    final_result <- as_tibble(result) %>%
      mutate(iteration = seq_along(theta)) %>%
      unnest_wider(theta)

    filename <- paste0("mcmc_bsl_", start_time, "_final", ".rds")
    saveRDS(final_result, filename)

    final_result
}


plot_param_traceplot <- function(result, param, true_value, ...) {
  ggplot(result, aes(iteration, {{param}})) +
    geom_line() +
    geom_hline(
      yintercept = true_value,
      color = "dodgerblue4"
    ) +
    ggtitle("Traceplot (samples)") +
    theme_gray(...)
}

plot_param_hist <- function(result, param, true_value, ...) {
  ggplot(result, aes({{param}})) +
    geom_histogram() +
    geom_vline(
      xintercept = true_value,
      color = "dodgerblue4",
      size = 2,
      alpha = 0.5
    ) +
    ggtitle("Posterior histogram") +
    theme_gray(...)
}

plot_log_prob_traceplot <- function(result, ...) {
  ggplot(result, aes(iteration, log_joint_prob)) +
    geom_line() +
    ggtitle("Traceplot (log joint density)") +
    ylab("log joint density") +
    theme_gray(...)
}
