ertarget=.4
erlow=.1

usr=90
lrp=30


forecast<-seq(0,200,by=5)
er<-NULL



for(i in seq_along(forecast)){
    if(forecast[i]*(1-ertarget)>usr){
      er[i]<-ertarget
     
    }else if(forecast[i]*(1-erlow)>lrp){

      er[i]<- max(ertarget*forecast[i]*(1-ertarget)/(usr),erlow)
      
    }else{
        er[i]<-erlow
    }
}

catch<-forecast*er
spawners<-forecast*(1-er)

par(mfrow=c(2,2))
plot(spawners,er,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,er,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(spawners,catch,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,catch,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
###############################33
#three steps
ermed=.2
erstep<-NULL

for(i in seq_along(forecast)){
    if(forecast[i]*(1-ertarget)>usr){
      erstep[i]<-ertarget
     
    }else if(forecast[i]*(1-ertarget)<=usr&
             forecast[i]*(1-ermed)>lrp){
      
      erstep[i]<- max(ermed,(forecast[i]-usr)/forecast[i])
      
    }else{
        erstep[i]<-erlow
    }
}



catchstep<-forecast*erstep
spawnersstep<-forecast*(1-erstep)

par(mfrow=c(2,2))
plot(spawnersstep,erstep,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,erstep,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(spawnersstep,catchstep,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,catchstep,type="b")
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")


#combine


ermed=.2

ertarget=0.7215355
erlow=.1

usr=104268.1*.8
lrp=15734.37


forecast<-seq(0,144508.7*5,by=2000)
er<-NULL
erstep<-NULL
erretro<-NULL
scaler<-NULL
scaler2<-NULL
 retrospawner<-seq(0,144508.7*5,by=2000)
 adj<-NULL
forecastupperbm<-usr/(1-ertarget)
        forecastlowbm<-lrp/(1-erlow)
        
for(i in seq_along(forecast)){
    if(forecast[i]*(1-ertarget)>usr){
      er[i]<-ertarget
      erstep[i]<-ertarget
      scaler[i]<-1
      scaler2[i]<-1
     
    }else if(forecast[i]*(1-erlow)<=lrp){
      er[i]<-erlow
      erstep[i]<-erlow 
      scaler[i]<-0 
      scaler2[i]<-0 
    }else{
      if(forecast[i]*(1-ermed)<=lrp){
        erstep[i]<-max(erlow,(forecast[i]-lrp)/forecast[i])  
      }else{
        erstep[i]<-max(ermed,(forecast[i]-usr)/forecast[i])
      }
      scaler[i]<-max(forecast[i]*(1-ertarget)-lrp,0)/(usr-lrp)
      scaler2[i]<-forecast[i]*(1-ertarget)/(usr) 
        ertmp<-ertarget
      #if(forecast[i]*(1-ertarget)<=lrp){
      #  er[i]<-max(erlow,(forecast[i]-lrp)/lrp)  
      #}else{
        
        er[i]<-erlow+(ertarget-erlow)*(forecast[i]-forecastlowbm)/(forecastupperbm-forecastlowbm)

        #plot(forecast,erlow+(ertarget-erlow)*(forecast-forecastlowbm)/(forecastupperbm-forecastlowbm))   

        
        #er[i]<-erlow+(ertarget-erlow)*((forecast[i]*(1-ertarget)-lrp)/(usr-lrp))
        #for(o in 1:120){
        #  er[i]<-erlow+(ertarget-erlow)*((forecast[i]*(1-er[i])-lrp)/(usr-lrp))
        #}
        

        adj[i]<-(forecast[i]*(1-ertarget)-lrp)/(usr-lrp)

      #}

        
      #er[i]<-max(min(erlow+ertarget*(max(forecast[i]*(1-ertarget)-lrp,0)/(usr-lrp)),ertarget),erlow)
            
      if(forecast[i]*(1-er[i])<=lrp){
        er[i]<-erlow  
      }
      
    }
}



catch<-forecast*er
spawners<-forecast*(1-er)

catchstep<-forecast*erstep
spawnersstep<-forecast*(1-erstep)


erretro<-NULL
erretrostep<-NULL

retrospawner<-seq(0,200,by=2)

for(i in seq_along(retrospawner)){
    if(retrospawner[i]>usr){
      erretro[i]<-ertarget
      erretrostep[i]<-ertarget
      scaler[i]<-1
      scaler2[i]<-1
     
    }else if(retrospawner[i]<=lrp){
      erretro[i]<-erlow
      erretrostep[i]<-erlow 
      scaler[i]<-0 
      scaler2[i]<-0 
    }else{
      
        erretrostep[i]<-ermed
     
      scaler[i]<-max(forecast[i]*(1-ertarget)-lrp,0)/(usr-lrp)
      scaler2[i]<-forecast[i]*(1-ertarget)/(usr) 
      #erretro[i]<-ertarget*retrospawner[i]/(usr)
      erretro[i]<-erlow+(ertarget-erlow)*((retrospawner[i]-lrp)/(usr-lrp))
            

      
    }
}









par(mfrow=c(2,2))

plot(spawners,er,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
lines(spawnersstep,erstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
lines(retrospawner,erretrostep,type="b",lwd=2,col=adjustcolor("darkblue", alpha.f = 0.6))
lines(retrospawner,erretro,type="b",lwd=2,col=adjustcolor("green", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,er,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
lines(forecast,erstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(spawners,catch,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
lines(spawnersstep,catchstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
lines(retrospawner,catchretro,type="b",lwd=2,col=adjustcolor("green", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,catch,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
lines(forecast,catchstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")




par(mfrow=c(2,2))

plot(spawners,er,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
#lines(spawnersstep,erstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
#lines(retrospawner,erretrostep,type="b",lwd=2,col=adjustcolor("darkblue", alpha.f = 0.6))
#lines(retrospawner,erretro,type="b",lwd=2,col=adjustcolor("green", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,er,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
#lines(forecast,erstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(spawners,catch,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
#lines(spawnersstep,catchstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
#lines(retrospawner,catchretro,type="b",lwd=2,col=adjustcolor("green", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")
plot(forecast,catch,type="b",lwd=2,col=adjustcolor("orange", alpha.f = 0.6))
#lines(forecast,catchstep,type="b",lwd=2,col=adjustcolor("darkgreen", alpha.f = 0.6))
abline(v=usr, col="blue")
abline(v=lrp, col="darkred")



ERtarget <- 0.72
ERlow <- 0.1
LRP <- 14000  # S low
URP <- 1000000 # S high 
forecast <- seq(0, 5000000, 5000)
 
ER_function <- function(forecast, ERtarget, ERlow, LRP, URP){
  slope <- (ERtarget-ERlow)/(URP - LRP)
  (ERlow + (forecast-LRP)*slope)/(1+forecast*slope)
}
 
plot(forecast, ER_function(forecast, ERtarget, ERlow, LRP, URP))
 
ER <- ER_function(forecast, ERtarget, ERlow, LRP, URP)
S <- forecast*(1-ER)
plot(S, ER)
 
## Linear in terms of Forecast:
ER_function_lf <- function(forecast, ERtarget, ERlow, LRP, URP){
  lower <- LRP/(1-ERlow)
  upper <- URP/(1-ERtarget)
  slope <- (ERtarget-ERlow)/(upper-lower)
  (ERlow + (forecast-lower)*slope)
}
 
ER2 <- ER_function_lf(forecast, ERtarget, ERlow, LRP, URP)
plot(forecast, ER2, type = "l")
plot(forecast*(1-ER2), ER2, type = "l", xlab = "Spawners")
 