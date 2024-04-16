# Explanatory figures

setwd("C:/Github/Ecotest/")

MMSE = readRDS("C:/temp/Ecotest/Batching/Independent_F/MMSE_1.rda")

dat = MMSE@PPD[[sno]][[fno]][[1]]
nsim<-MMSE@nsim
nyears<-MMSE@nyears
proyears = MMSE@proyears
allyears = nyears+proyears

datlist<-list()

Iobs<-dat@Ind
Cobs<-dat@Cat
CAL = dat@CAL
mids = dat@CAL_mids

Bt = MMSE@SB_SBMSY
I_cur<-I_mu<-I_rel<-I_s5<-I_s10<- C_cur<- C_mu<-C_rel<-C_s5<-C_s10<- 
  ML_cur<-ML_mu<-ML_rel<-ML_s5<-ML_s10<- MV_cur<- MV_mu<-MV_rel<- MV_s5<-MV_s10<-
  FM_cur <- FM_mu <- FM_rel <- FM_s5 <- FM_s10 <- rep(NA,nsim)

matage = MMSE@multiHist[[sno]][[fno]]@AtAge$Maturity[,,1]
lenage = MMSE@multiHist[[sno]][[fno]]@AtAge$Length[,,1]

L50 = sapply(1:nsim,function(X,matage,lenage)approx(x=matage[X,],y=lenage[X,],xout=0.5)$y,matage=matage,lenage=lenage)

yrs = 2013+(-nyears+1):(proyears-1)


exp_trend = function(x, Iyr=95, yrs, nyears,ylab="",
                     curplot=F){
  
  Is<-smooth2(x,plot=F)
  ind = 1:Iyr
  plot(yrs[ind],x[ind],pch=19,ylim=c(0,quantile(x,0.99)),col="#999999",xlab="",ylab=""); grid()
  abline(v=yrs[nyears]+0.5,lty=2,lwd=2)
  
  lines(yrs[ind],Is[ind],col="#ff000099",lwd=2)
  cur<-Is[Iyr]
  if(curplot){
    points(yrs[Iyr],cur,pch=2,cex=2,lwd=2,col="#0000ff99")
    mtext
  }
  if(slpplot){
    abline(v=yrs[Iyr]-c(1,6,11),col="#0000ff99",lwd=2,lty=2)
    
  }
  mu = mean(x[nyears:Iyr])

  rel<-cur / mean(Is[1:Iyr])
  s5<-slp3(Iobs[i,Iyr-(5:1)])
  s10<-slp3(Iobs[i,Iyr-(10:1)])
  
  
  
}



# Catches

x = Cobs








