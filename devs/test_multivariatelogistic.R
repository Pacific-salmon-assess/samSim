#=========================================
#test the multivariate logistic routines for age comps
#=====================================
#library(samSim)
estTaus<-rep(0,100)
for(u in 1:100){
  #simulate soem data
ppnMat<-matrix(0,nrow=40,ncol=4)

for(i in 1:40){
   ppnMat[i,]<-ppnAgeErr(c(0.1,0.3,0.4,0.2), 0.25,
        error = runif(4, 0.0001, 0.9999)) 
}

#estimate and recover tau


esttau<-getTau(ppnMat, plotTaus = TRUE)

estTaus[u]<-esttau$bestTau
  
}

boxplot(estTaus)
abline(h=0.25, lty=2, col="darkred", lwd=2)


#============================
#test code from Paul

A <- 4
T <- 40 
sigma <- 0.25
mu <- c(0.1,0.3,0.4,0.2)
mu <- mu/sum(mu)
epsilon <- matrix(rnorm(A*T), nrow = A, ncol = T)
xobs <- matrix(0, nrow = A, ncol = T)
pobs <- matrix(0, nrow = A, ncol = T)
pobs2 <- pobs
for( i in 1:T ) {
  xobs[,i] <- log(mu) + epsilon[,i] - 1/A*sum(log(mu) - epsilon[,i])
  pobs[,i] <- exp(xobs[,i])/sum(exp(xobs[,i]))
  pobs2[,i] <- exp(log(mu) + epsilon[,i])/sum(exp(log(mu) + epsilon[,i]))
}
pobs2-pobs




mu <- rnorm(A-1, 0, 3)
ptrue <- exp(c(mu, 0))/sum(exp(c(mu, 0)))
epsilon <- matrix(rnorm((A-1)*T, 0, sd = sigma), nrow = (A-1), ncol = T)
xobs <- matrix(0, nrow = A, ncol = T)
pobs <- matrix(0, nrow = A, ncol = T)
for( i in 1:T ) {
  xobs[,i] <- c(mu + epsilon[,i], 0)
  pobs[,i] <- exp(xobs[,i])/sum(exp(xobs[,i]))
}
 



dmvlogit <- function(x, mean, sd, log = FALSE){
  np <- length(x)
  ans <- sum( dnorm(log( x[-np]/x[np] ), mean, sd, log = TRUE) ) - sum(log(x))
  # ans <- -(np-1)/2*log(2*pi*sd^2) - 0.5*sum((log( x[-np]/x[np] ) - mean)^2)/sd^2 - sum(log(x)) 
  if(log) return(ans)
  else return(exp(ans))
}

negll <- function(pars, obs){
  np <- length(pars)
  mu <- pars[1:(np-1)]
  sigma <- exp(pars[np])

  nll <- 0
  for( i in 1:ncol(obs)){
    nll <- nll - dmvlogit(obs[,i], mu, sigma, log = TRUE)
  }
  return(nll)
}



trans <- function (alpha){
 expa <- exp(alpha)
 c(expa , 1+numeric(1))/(1+sum(expa))
}

negll_schnute <- function(pars, props){
  np <- length(pars)
  mu <- trans(pars[1:(np-1)])
  sigma <- exp(pars[np])
  
  nll <- 0
  for( i in 1:ncol(props)){
    y <- log(props[,i]) - mean(log(props[,i]))
    mean <- log(mu) - mean(log(mu))
    nll <- nll - sum(dnorm(y, mean, sigma, log = TRUE)) 
  }
  return(nll)
} 


#test with schunute simulation and Paul's estimation
est_tau_grid<-rep(0,100)
est_tau_nll<-rep(0,100)
est_tau_nll_paul<-rep(0,100)

est_pa_nll<-matrix(0,nrow=100,ncol=4)
est_pa_nll_paul<-matrix(0,nrow=100,ncol=4)

for(u in 1:100){

  for(i in 1:40){
    ppnMat[i,]<-ppnAgeErr(c(0.1,0.3,0.4,0.2), 0.25,
        error = runif(4, 0.0001, 0.9999)) 
  }

   #estimate and recover tau

   #grid
   estgrid<-getTau(ppnMat, plotTaus = TRUE)
   est_tau_grid[u]<-estgrid$bestTau

   #paul
   fit <- nlminb(numeric(A), negll, obs = t(ppnMat))
   phat <- c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
   est_tau_nll_paul[u]<-exp(fit$par[A])
   est_pa_nll_paul[u]<-c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
   
   #schnute
   fit_schnute <- nlminb(numeric(A), negll_schnute, props = t(ppnMat))
   est_tau_nll[u]<-exp(fit_schnute$par[A])
   est_pa_nll[u]<-trans(fit$par[1:(A-1)])
}


par(mfrow=c(1,3))
boxplot(estTaus)
abline(h=0.25, lty=2, col="darkred", lwd=2)






fit <- nlminb(numeric(A), negll, obs = pobs)
fit$par[1:(A-1)] - mu
phat <- c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
exp(fit$par[A])
sigma
 
library(RTMB)
 
obs <- pobs
pars <- list(mu = numeric(A-1), logsigma = 0)
negll <- function(pars){
  getAll(pars)
  sigma <- exp(logsigma)
 
  nll <- 0
  for( i in 1:ncol(obs)){
    nll <- nll - dmvlogit(obs[,i], mu, sigma, log = TRUE)
  }
  phat <- c(exp(mu), 1)
  phat <- phat/sum(phat)
  REPORT(phat)
  return(nll)
}
 
sigmai <- seq(0, 5, by = 0.1)
dll <- NULL
for( sigj in sigmai ) dll <- c(dll, dmvlogit(obs[,i], mu, sigj))
plot(sigmai, dll)
 
obj <- MakeADFun(negll, pars, silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit$par[1:(A-1)] - mu
phat <- c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
exp(fit$par[A])
sigma
 
sdrep <- sdreport(obj)
plot(obj$report()$phat, ptrue)
abline(a=0, b=1)