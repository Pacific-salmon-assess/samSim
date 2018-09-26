rtnorm <- function(n, mean = 0, sd = 1, min = 0, max = 1) {
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  qnorm(u, mean, sd)
}

mu = 1.2; sig = 0.1
set.seed(123)
dd <- rtnorm(10000, mean = mu, sd = sig)
set.seed(123)
loc <- log(mu^2 / sqrt(sig^2 + mu^2))
shape <- sqrt(log(1 + (sig^2 / mu^2)))
dd2 <- rlnorm(10000, loc, shape)
exp(rnorm(5, mean = log(mu), sd = log(sig)))
