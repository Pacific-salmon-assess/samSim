realizedHarvestRate <- function(targetHarvest, sigmaHarvest = 0.1,
                                location = NULL, shape = NULL, normal = FALSE) {
  if (is.null(location)) {
    location <- pmax(0.001,
                     targetHarvest^2 * (((1 - targetHarvest) / sigmaHarvest^2)
                                        - (1 / targetHarvest)))
  }
  if (is.null(shape)) {
    shape <- pmax(0.001,
                  location * (1 / targetHarvest - 1))
  }
  if (normal == TRUE) {
    h <- rnorm(n = length(targetHarvest), targetHarvest, sigmaHarvest)
  } else {
    h <- rbeta(n = length(targetHarvest), shape1 = location, shape2 = shape)
  }
  return(h)
}

x1 <- runif(10, 0, 0.99)
x2 <- runif(10, 0.5, 0.99)

# default behavior
set.seed(123)
f(x1)
runif(1)
set.seed(123)
f(x2)
runif(1)

# see it works with normal!
set.seed(123)
f(x1, normal = TRUE)
runif(1)
set.seed(123)
f(x2, normal = TRUE)
runif(1)


doMatch <- numeric(1000)
for(i in 1:1000){
  # default behavior
  set.seed(123)
  realizedHarvestRate(targetHarvest = x1[i], sigmaHarvest = 0.1)
  y1 <- sample.int(.Machine$integer.max, 1)
  set.seed(123)
  realizedHarvestRate(targetHarvest = x2[i], sigmaHarvest = 0.1)
  y2 <- sample.int(.Machine$integer.max, 1)

  doMatch[i] <- ((y1-y2) == 0)
}

plot(x1, x2, col=doMatch+1)
text(0.5, 0.5, "MATCH", col=2, font=2)
text(0.5, 0.99, "DON'T MATCH", col=1, font=2)
