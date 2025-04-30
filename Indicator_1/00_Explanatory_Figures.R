# Explanatory figures

setwd("C:/Github/Ecotest/")
library(miceadds)
source.all("Source")

MMSE = readRDS("C:/temp/Ecotest/Batching/Independent_F/MMSE_1.rda")

nsim<-MMSE@nsim
nyears<-MMSE@nyears
proyears = MMSE@proyears
allyears = nyears+proyears
yrs = 2013+(-nyears+1):(proyears-1)

# --------- Time series data ---------------------------------------------------------------------------

sno = 1; fno = 1
dat = MMSE@PPD[[sno]][[1]][[1]]

Iobs<-dat@Ind
Cobs<-dat@Cat
CAL = dat@CAL
mids = dat@CAL_mids

matage = MMSE@multiHist[[sno]][[fno]]@AtAge$Maturity[,,1]
lenage = MMSE@multiHist[[sno]][[fno]]@AtAge$Length[,,1]

L50 = sapply(1:nsim,function(X,matage,lenage)approx(x=matage[X,],y=lenage[X,],xout=0.5)$y,matage=matage,lenage=lenage)

# Catches

jpeg("Figures/Presentation 1 April 24/Catch_explanatory.jpg",res=600,width=7,height=6,units="in")
  simno = 4; x = Cobs[simno,]
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.05,0.05),omi=c(0,0.5,0,0))
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",relplot=T )
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",slpplot=T )
  mtext("Catch",2,line=0.3,outer=T)
dev.off()

# CPUE

jpeg("Figures/Presentation 1 April 24/CPUE_explanatory.jpg",res=600,width=7,height=6,units="in")
  simno = 4; x = Iobs[simno,]
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.05,0.05),omi=c(0,0.5,0,0))
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",relplot=T )
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",slpplot=T )
  mtext("Nominal CPUE",2,line=0.3,outer=T)
dev.off()

# Mean length

jpeg("Figures/Presentation 1 April 24/ML_explanatory.jpg",res=600,width=7,height=6,units="in")
  simno = 4; 
  CALi = CAL[simno,,]
  totlen = CALi*t(array(mids,dim(t(CALi))))
  x = apply(totlen,1,sum)/apply(CALi,1,sum)
  
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.05,0.05),omi=c(0,0.5,0,0))
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",relplot=T )
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",slpplot=T )
  mtext("Mean length (cm)",2,line=0.3,outer=T)
dev.off()

# var length

jpeg("Figures/Presentation 1 April 24/varlen_explanatory.jpg",res=600,width=7,height=6,units="in")
  simno = 4; 
  CALi = CAL[simno,,]
  x = sapply(1:dim(CALi)[1],function(X,CALi,mids){
    sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
  },CALi=CALi,mids=mids)
  
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.05,0.05),omi=c(0,0.5,0,0))
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",relplot=T )
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",slpplot=T )
  mtext("Variability length (CV)",2,line=0.3,outer=T)
dev.off()

# Fraction matuire

jpeg("Figures/Presentation 1 April 24/fracmat_explanatory.jpg",res=600,width=7,height=6,units="in")
  simno = 4; 
  CALi = CAL[simno,,]
  L50i = L50[simno]
  x = sapply(1:dim(CALi)[1],function(X,CALi,mids,L50i){
    indivs = rep(mids,CALi[X,])
    mean(indivs > L50i,na.rm=T)
  },CALi=CALi,mids=mids,L50i=L50i)
  
  par(mfrow=c(2,1),mai=c(0.5,0.5,0.05,0.05),omi=c(0,0.5,0,0))
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",relplot=T )
  exp_trend(x, Iyr=95, yrs, nyears,ylab="",slpplot=T )
  mtext("Fraction mature in catch",2,line=0.3,outer=T)
dev.off()


# --------- Catch ratios -----------------------------------------------------------------------------------

sno1 = 5; sno2 = 6; fno = 1; simno = 1; Iyr = 95; ind = 1:Iyr
Cobs1 = MMSE@PPD[[sno1]][[1]][[1]]@Cat[simno,]
Cobs2 = MMSE@PPD[[sno2]][[1]][[1]]@Cat[simno,]
Snames = sapply(MMSE@Snames,function(x)substr(x,1,3))

CS1 = smooth2(Cobs1)
CS2 = smooth2(Cobs2)
rat = CS1/CS2


jpeg("Figures/Presentation 1 April 24/catch_ratio_explanatory.jpg",res=600,width=8,height=6,units="in")

  par(mfcol=c(3,2),mai=c(0.3,0.3,0.05,0.05),omi=c(0,0.5,0,0))
  plot(yrs[ind],Cobs1[ind],ylim=c(0,quantile(Cobs1[ind],0.99)),col = "#99999999",pch=19,xlab="",ylab=""); grid()
  lines(yrs[ind],CS1[ind],col="#ff000099",lwd=2); legend('topleft',Snames[sno1],text.col="#ff000099", bty="n")
  mtext("Catch",2,line=2.5)
  plot(yrs[ind],Cobs2[ind],ylim=c(0,quantile(Cobs2[ind],0.99)),col = "#99999999",pch=19,xlab="",ylab=""); grid()
  lines(yrs[ind],CS2[ind],col="#0000ff99",lwd=2); legend('topleft',Snames[sno2],text.col="#0000ff99", bty="n")
  mtext("Catch",2,line=2.5)
  plot(yrs[ind],rat[ind],ylim=c(0,max(rat[ind])),type="l",col="purple",xlab="",ylab="",lwd=2); grid()
  legend('topleft',legend=paste0("Catch ",Snames[sno1]," / Catch ",Snames[sno2]),text.col="purple",bty="n")
  mtext("Catch ratio",2,line=2.5)
  
  exp_trend(rat, Iyr=95, yrs, nyears,ylab="Catch ratio",relplot=T,pointcol='white',smoothcol="purple")
  exp_trend(rat, Iyr=95, yrs, nyears,ylab="Catch ratio",slpplot=T,pointcol='white',smoothcol="purple")
  exp_trend(rat, Iyr=95, yrs, nyears,ylab="Catch ratio",muplot=T,pointcol='white',smoothcol="purple")

dev.off()






nsim<-MMSE@nsim
nyears<-MMSE@nyears
proyears = MMSE@proyears
allyears = nyears+proyears

datlist<-list()

Iobs<-dat@Ind
Cobs<-dat@Cat
