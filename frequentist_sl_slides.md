# Motivating the synthetic likelihood

Synthetic likelihood originates in a *frequentist setting*, when the likelihood function is too irregular to easily maximize (analytically or numerically).
As long as we can generate samples, the synthetic likelihood facilitates parameter inference.

Motivating example: Ricker population model

# Motivating the synthetic likelihood

In the Ricker model, $N_t$ models the time course of a population, where

\begin{align*}
N_0 &= 1 \\
N_t &= r N_{t-1} \mathrm{e}^{- N_{t-1} + e_{t-1}} \\
e_t &\overset{\text{iid}}{\sim} \mathcal{N}(0, \sigma^2)
\end{align*}

In other words, $N_t | N_{t-1} \sim LogNormal(\log N_{t-1} + \log r - N_{t-1}, \sigma^2)$.

Additionally, suppose $N_t$ is not actually observed, but a sample from the population $Y_t | N_t \sim Poisson(\phi N_t)$ is observed.

Can we obtain the likelihood $f(\bm{Y}; r, \sigma, \phi)$ to do inference?

# Motivating the synthetic likelihood

Let $\bm{\theta} = (r, \sigma, \phi)$.
The joint likelihood is easier to work out:

\begin{align*}
f(\bm{Y}, \bm{N}; \bm{\theta})
&= f(\bm{N}; \bm{\theta}) f(\bm{Y} | \bm{N}; \bm{\theta}) \\
&= f(N_1; \bm{\theta}) \prod_{t = 2}^n f(N_t | N_{t-1}; \bm{\theta}) \prod_{t = 1}^n f(Y_t | N_t; \bm{\theta}) \\
&= LogNormal(N_1; ...) \prod_{t = 2}^n LogNormal(N_t; ...) \prod_{t = 1}^n Poisson(Y_t; ...)
\end{align*}

But, to do frequentist inference, we would need to marginalize over the $N_t$:

$$
f(\bm{Y}; r, \sigma, \phi) =
\int_{N_1} ... \int_{N_n} f(\bm{Y}, \bm{N}; \bm{\theta}) \, dN_1 ... dN_n
$$

This is expensive, in part because the number of integrals grows with the number of data points collected.