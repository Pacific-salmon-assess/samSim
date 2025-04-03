#=========================================
#test the multivariate logistic routines for age comps
#=====================================
#library(samSim)
library(ggplot2)

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

ppnMat<-matrix(0,nrow=40,ncol=4)


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
   est_pa_nll_paul[u,]<-c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
   
   #schnute
   fit_schnute <- nlminb(numeric(A), negll_schnute, props = t(ppnMat))
   est_tau_nll[u]<-exp(fit_schnute$par[A])
   est_pa_nll[u,]<-trans(alpha=fit$par[1:(A-1)])
}



taudf<-data.frame(estmod=rep(c("grid","paul","schnute"),100),
                  value=c(est_tau_grid,est_tau_nll_paul,est_tau_nll))


ggplot(taudf)+
geom_boxplot(aes(x=estmod,y=value))+
theme_bw(15)


padf<-data.frame(estmod=rep(c("paul","schnute"),100*A),

                  value=c(est_pa_nll_paul,est_pa_nll),
                  age=rep(rep(1:4, each=100),2))

ggplot()+
geom_boxplot(data=padf,aes(x=as.factor(age),y=value,fill=estmod))+
geom_line(aes(x=1:4,y=c(0.1,0.3,0.4,0.2)))+
geom_point(aes(x=1:4,y=c(0.1,0.3,0.4,0.2)))+
theme_bw(15)


#==============
#try with different numbers



#test with schunute simulation and Paul's estimation
est_tau_grid<-rep(0,100)
est_tau_nll<-rep(0,100)
est_tau_nll_paul<-rep(0,100)

est_pa_nll<-matrix(0,nrow=100,ncol=4)
est_pa_nll_paul<-matrix(0,nrow=100,ncol=4)


truepa<-c(0.07052628, 0.32951836, 0.48814346, 0.11181190)
truetau<-0.25


for(u in 1:100){

  for(i in 1:40){
    ppnMat[i,]<-ppnAgeErr(truepa, truetau,
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
   est_pa_nll_paul[u,]<-c(exp(fit$par[1:(A-1)]),1)/(1+sum(exp(fit$par[1:(A-1)])))
   
   #schnute
   ans<-mvLogisticLL(agePropData=ppnMat)
   est_tau_nll[u]<-ans$tau
   est_pa_nll[u,]<-ans$pa

}



taudf<-data.frame(estmod=rep(c("grid","paul","schnute"),100),
                  value=c(est_tau_grid,est_tau_nll_paul,est_tau_nll))


ggplot(taudf)+
geom_boxplot(aes(x=estmod,y=value))+
theme_bw(15)


padf<-data.frame(estmod=rep(c("paul","schnute"),100*A),

                  value=c(est_pa_nll_paul,est_pa_nll),
                  age=rep(rep(1:4, each=100),2))

ggplot()+
geom_boxplot(data=padf,aes(x=as.factor(age),y=value,fill=estmod))+
geom_line(aes(x=1:4,y=c(0.1,0.3,0.4,0.2)))+
geom_point(aes(x=1:4,y=c(0.1,0.3,0.4,0.2)))+
theme_bw(15)





#==============
#test and write example

data(widedatprop)

as.matrix(widedatprop[,-1])

mvLogisticLL(agePropData=as.matrix(widedatprop[,-1]))


#=====================================
#what to do when we have zeroes



 #eps<-min(agePropData[agePropData>0])/2

 #  agePropData<-agePropData+eps
 # agePropData<-agePropData/apply(agePropData,1,sum)

est_tau_nll<-rep(0,500)
est_tau_nll_zero<-rep(0,500)

est_pa_nll<-matrix(0,nrow=500,ncol=4)
est_pa_nll_zero<-matrix(0,nrow=500,ncol=4)

truepa<-c(0.07052628, 0.32951836, 0.48814346, 0.11181190)
truetau<-0.25

for(u in 1:500){

  ppnMat<-matrix(0,nrow=40,ncol=4)

  for(i in 1:40){
    ppnMat[i,]<-ppnAgeErr(truepa, truetau,
        error = runif(4, 0.0001, 0.9999)) 
  }
  ppnMat0<-ppnMat

  #add in zeros artificially
  ppnMat0[min(ppnMat[,1]),1]<-0
  ppnMat0[order(ppnMat[,4])<4,4]<-0
  ppnMat0<-ppnMat0/apply(ppnMat0,1,sum)
   #estimate and recover tau
  
   #no zero data
   ans<-mvLogisticLL(agePropData=ppnMat)
   est_tau_nll[u]<-ans$tau
   est_pa_nll[u,]<-ans$pa


  # zero data
   ans<-mvLogisticLL(agePropData=ppnMat0)
   est_tau_nll_zero[u]<-ans$tau
   est_pa_nll_zero[u,]<-ans$pa
}


taudf<-data.frame(estmod=rep(c("schnute_full","schnute_zero"),500),
                  value=c(est_tau_nll,est_tau_nll_zero))


ggplot()+
geom_boxplot(data=taudf,aes(x=estmod,y=value,fill=estmod))+
geom_hline(aes(yintercept=truetau))+
theme_bw(15)


padf<-data.frame(estmod=rep(c("schnute_full","schnute_zero"),500*A),
                  value=c(est_pa_nll,est_pa_nll_zero),
                  age=rep(rep(1:4, each=500),2))

ggplot()+
geom_boxplot(data=padf,aes(x=as.factor(age),y=value,fill=estmod))+
geom_line(aes(x=1:4,y=truepa))+
geom_point(aes(x=1:4,y=truepa))+
theme_bw(15)

