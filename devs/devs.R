




devtools::document()

devtools::load_all()




#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)


#debug step HCR

lowerObsBM=100
upperObsBM=300

foreRecRY=seq(10,1200,20)
trendCanER.iter<-rep(0,length(foreRecRY))
trendCanER=.5

amberStatusER<-.25
redStatusER<-.05


if(is.na(amberStatusER)){amberStatusER<-max(trendCanER*bmERAdj,0.05,na.rm=T)}
if(is.na(redStatusER)){redStatusER<-max(trendCanER*bmERAdj,0.05,na.rm=T)}

for(y in seq_along(foreRecRY)){
   
  if((foreRecRY[y]*(1-amberStatusER))<lowerObsBM){
    #red status
    trendCanER.iter[y] <- redStatusER
    
  }else if((foreRecRY[y]*(1-trendCanER)>=lowerObsBM&
           foreRecRY[y]*(1-trendCanER)<upperObsBM)|
           (foreRecRY[y]*(1-amberStatusER)>=lowerObsBM&
             foreRecRY[y]*(1-amberStatusER)<upperObsBM)){
    
    trendCanER.iter[y]<- amberStatusER
  
   
  }else{
    trendCanER.iter[y] <- max(trendCanER,0.05,na.rm=T)
  }
}



plot(foreRecRY*(1-trendCanER.iter),trendCanER.iter)
abline(v=lowerObsBM, col="darkred")
abline(v=upperObsBM, col="darkred")

cbind(foreRecRY*(1-trendCanER.iter),trendCanER.iter)

plot(foreRecRY*(1-trendCanER),trendCanER.iter)
abline(v=lowerObsBM, col="darkred")
abline(v=upperObsBM, col="darkred")



#second:

for(y in seq_along(foreRecRY)){
   
  if((foreRecRY[y]*(1-amberStatusER))<lowerObsBM){
    #red status
    trendCanER.iter[y] <- redStatusER
    
  }else if(foreRecRY[y]*(1-trendCanER)>lowerObsBM&
           foreRecRY[y]*(1-trendCanER)<upperObsBM|
           foreRecRY[y]*(1-amberStatusER)>lowerObsBM&
             foreRecRY[y]*(1-amberStatusER)<upperObsBM){
    
    trendCanER.iter[y]<- amberStatusER
  
   
  }else{
    trendCanER.iter[y] <- max(trendCanER,0.05,na.rm=T)
  }
}



plot(foreRecRY*(1-trendCanER.iter),trendCanER.iter)
abline(v=lowerObsBM, col="darkred")
abline(v=upperObsBM, col="darkred")