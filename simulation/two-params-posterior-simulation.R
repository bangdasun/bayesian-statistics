# Two parameters model simulation
# 
# (alpha, beta) are two parameters with prior
# prior could be informative: Normal (here I use this one)
# or noninformative: proper prior
#
# data is in raw data
# simulate the posterior based on Bayesian inference and draw random sample from posterior


## ------------------------------------------------------------------------
library("reshape2")
library("dplyr")
library("ggplot2")

## ------------------------------------------------------------------------
raw_data <- data.frame(
  x = c(-.86, -.3, -.05, .73),
  n = c(5, 5, 5, 5),
  y = c(0, 1, 3, 5)
)
S <- 1000

## ------------------------------------------------------------------------
log_reg <- glm(cbind(y, n - y) ~ x, data = raw_data, family = binomial(link = 'logit'))
coef(log_reg)

## ------------------------------------------------------------------------
Alpha <- seq(-4, 8, by = .05)
Beta  <- seq(-5, 45, by = .2)
grid_df <- expand.grid(alpha = Alpha, beta = Beta)

## ------------------------------------------------------------------------
bivariate_normal <- function(x, y, mean = c(0, 10), sd = c(2, 10), rho = 0.5) {
  const <- 1 / (2 * pi * sd[1] * sd[2] * sqrt(1 - rho^2))
  expo  <- exp(-((x - mean[1])^2/(sd[1]^2) + (y - mean[2])^2/(sd[2]^2) - 
      (2 * rho * (x - mean[1]) * (y - mean[2]))/(sd[1] * sd[2])) /
      (2 - 2 * rho^2))
  return(const * expo)
}
bivariate_normal(-5, -10)
grid_df$prior <- mapply(bivariate_normal, x = grid_df$alpha, y = grid_df$beta)

## ------------------------------------------------------------------------
logit_inv <- function(x) return(exp(x) / (1 + exp(x)))

lik <- function(Alpha, Beta, n, x, y) {
  return(logit_inv(Alpha + Beta * x)^y * (1 - logit_inv(Alpha + Beta * x))^(n - y))
}

prod_lik <- function(Alpha, Beta) {
  val <- rep(NA, length(Alpha))
  for (i in seq(nrow(raw_data))) {
    val[i] <- lik(Alpha, Beta, raw_data$n[i], raw_data$x[i], raw_data$y[i])
  }
  return(prod(val))
}
grid_df$lik <- mapply(prod_lik, Alpha = grid_df$alpha, Beta = grid_df$beta)
grid_df$posterior <- grid_df$prior * grid_df$lik

## ------------------------------------------------------------------------
# Comparison of prior and posterior
grid_df %>%
  melt(id.vars = c("alpha", "beta", "lik")) %>%
  ggplot() + 
  geom_contour(aes(x = alpha, y = beta, z = value, col = variable), size = .6)

## ------------------------------------------------------------------------
delta <- (Alpha[2] - Alpha[1]) / 2
epsilon <- (Beta[2] - Beta[1]) / 2

## ------------------------------------------------------------------------
normalizer <- 1 / (0.01 * sum(grid_df$posterior))
grid_df$posterior <- grid_df$posterior * normalizer
(prior_mode <- max(grid_df$prior))
(posterior_mode <- max(grid_df$posterior))

## ------------------------------------------------------------------------
Alpha_marginal_posterior <- grid_df %>%
  group_by(alpha) %>%
  summarise(marginal_posterior = sum(posterior))

U <- runif(S)
sample_Alpha <- rep(NA, length(U))
for (i in seq(length(U))) {
  s <- sample(length(Alpha), 1, prob = Alpha_marginal_posterior$marginal_posterior)
  sample_Alpha[i] <- Alpha_marginal_posterior$alpha[s]
}

## ------------------------------------------------------------------------
sample_Beta <- rep(NA, length(U))
for (i in seq(length(U))) {
  sample_range <- grid_df[grid_df$alpha == sample_Alpha[i], c("beta", "posterior")]
  s <- sample(length(Beta), 1, prob = sample_range$posterior)
  sample_Beta[i] <- sample_range$beta[s]
}

## ------------------------------------------------------------------------
sample_Alpha <- sample_Alpha + runif(length(sample_Alpha), -delta, delta)
sample_Beta <- sample_Beta + runif(length(sample_Beta), -epsilon, epsilon)
sample_Alpha_Beta <- data.frame(alpha = sample_Alpha, beta = sample_Beta)

# Random sample on posterior
ggplot() + 
  geom_contour(data = grid_df, aes(x = alpha, y = beta, z = posterior), col = "red") + 
  geom_point(data = sample_Alpha_Beta, aes(x = alpha, y = beta), alpha = I(1/5))
