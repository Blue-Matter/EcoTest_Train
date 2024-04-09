
setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_Indicators.R")

MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/Hist_1.rda")
Catchdat = ICCATtoGEO(read.csv("G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped.csv"))

logit_list = list()

nstocks = MMSE@nstocks
SCodes = sapply(MMSE@Snames,function(x)substr(x,1,3))
SSBrel = list()
for(ss in 1:nstocks) SSBrel[[ss]] = apply(MMSE@multiHist[[ss]]$Longline@TSdata$SBiomass,1:2,sum) / MMSE@RefPoint$ByYear$SSBMSY[,ss,1,1:MMSE@nyears]

lapply(1:nstocks,getlogitmod,SSBrel,Catchdat,SCodes,ploty=T)
 

getlogitmod = function(x,SSBrel,Catchdat,SCodes, Pyrs = seq(1975,2010,by=5), ploty=F, yrs = 1950:2013){
  
  col = match(SCodes[x],names(Catchdat))
  SSBr = SSBrel[[x]][1,]
  Cat = as.numeric(Catchdat[,col])
  Cagg = aggregate(Cat, by = list(Yr = Catchdat$YearC, Lat = Catchdat$Lat, Lon = Catchdat$Lon, Cell = Catchdat$CellID), FUN = sum, na.rm=T)
  CR = aggregate(Catchdat$Eff1, by = list(Yr = Catchdat$YearC, Lat = Catchdat$Lat, Lon = Catchdat$Lon, Cell =Catchdat$CellID), sum, na.rm=T)
  CR$x = Cagg$x / CR$x
  CR = CR[CR$x > 0,]
  
  Nbyyr = aggregate(rep(1,nrow(CR)),by=list(Yr = CR$Yr),sum)
  firstYr = min(Nbyyr$Yr[Nbyyr$x > (0.95 * max(Nbyyr$x))])
  yreval = firstYr:max(yrs)
  yind = yrs%in%yreval
  minCR = mean(CR$x)/10
  meanCR = mean(CR$x)
  if(ploty){
    np = length(Pyrs)+3
    par(mfrow=c(ceiling(np/3),3),mai=c(0.5,0.5,0.1,0.1))
    plot(yrs[yind],SSBr[yind],type="l",lwd=2,ylim=c(0,max(SSBr)),ylab = "SSB / SSBMSY",xlab=""); grid(); abline(v=Pyrs,col="green") 
    legend('top',legend=SCodes[x],bty="n")
    sapply(1:length(Pyrs),plotCR,CR=CR,Pyrs=Pyrs,minCR=minCR,meanCR=meanCR)
  }
    
  CRt = CR[CR$Yr>=firstYr & CR$Yr <=max(yrs),]
  
  cells = unique(CRt$Cell)
  ncell = length(cells)
 
  ny = length(yreval)
  CRarr = CRmin = CRmean = array(0,c(ncell,ny))
  CRarr[as.matrix(cbind(match(CRt$Cell,cells),match(CRt$Yr,yreval)))]=CRt$x
  
  # frac cells with CR over min

  CRmin[CRarr>minCR] = 1
  CRminy = apply(CRmin,2,mean)
  
  # frac cells with CR over mean
  
  CRmean[CRarr>meanCR] = 1
  CRmeany = apply(CRmean,2,mean)
 
  logit = function(p)log(p/(1-p))
  ilogit = function(x)exp(x)/(1+exp(x))
  
  mod = lm(y~x,data=data.frame(x=log(SSBr[yind]),y=logit(CRminy)))
  
  if(ploty){
    plot(yreval,CRminy,col="blue",ylim=c(0,max(CRminy)),type="l",xlab="",ylab="Fraction above 1/10 hist. mean");grid()
    plot(SSBr[yind],CRminy,col="white",pch=19,xlim = c(0,max(SSBr[yind])),ylim=c(0,max(CRminy)),xlab = "SSB / SSBMSY",ylab="Fraction above 1/10 hist. mean");grid()
    text(SSBr[yind],CRminy,yreval,cex=0.8)
    
    SBseq=seq(0.01,5,length.out=1000)
    pred = predict(mod,newdata = data.frame(x=log(SBseq)))
    lines(SBseq,ilogit(pred),col="blue",lwd=2)
    #plot(yreval,CRmeany,col="red",ylim=c(0,max(CRmeany)),type="l",xlab="",ylab="Fraction above hist. mean")
    #plot(SSBr[yind],CRmeany,col="white",pch=19,xlab = "SSB / SSBMSY",ylab="Fraction above 1/10 hist. mean");text(SSBr[yind],CRmeany,yreval,cex=0.8)
  }
  
  mod
}


plotCR = function(x,CR,Pyrs,mincex = 0.2,addcex=1.8,minCR,meanCR){
  
  CRrng = CR[CR$Yr %in% Pyrs,]
  CRs = CR[CR$Yr == Pyrs[x],]
  
  xlim = range(CRrng$Lon)+c(-2.5,2.5)
  ylim = range(CRrng$Lat)+c(-2.5,2.5)
  CRmax = quantile(CRrng$x^0.5,0.99)
  
  plot(xlim,ylim,col="white",xlab="",ylab="")
  ls = (-100:100)*5
  abline(h = ls, v=ls,col="grey")
  
  cex = 0.2 + ((CRs$x^0.5)/CRmax)*addcex
  cols = rep("#0000ff95",nrow(CRs))
  cols[CRs$x < minCR] = "black"
  cols[CRs$x > meanCR] = '#ff000095'
  points(CRs$Lon, CRs$Lat, cex=cex,pch=19,col=cols)
  legend('top',legend=Pyrs[x],bty="n")
  
}

