#=========================================
#test the multivariate logistic routines for age comps
#=====================================

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
