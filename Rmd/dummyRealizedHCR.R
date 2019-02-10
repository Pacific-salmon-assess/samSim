f <- function(mu) {
  location <- pmax(0.001,
                   mu^2 * (((1 - mu) / sigma^2) - (1 / mu)))
  shape <- pmax(0.001,
                location * (1 / mu - 1))
  h <- rbeta(n = length(mu), shape1 = location, shape2 = shape)
  # set.seed(999)
  return(h)
}

x1 <- runif(10, 0, 0.99)
x2 <- runif(10, 0.5, 0.99)

set.seed(123)
f(x1)
runif(1)

set.seed(123)
f(x2)
runif(1)
