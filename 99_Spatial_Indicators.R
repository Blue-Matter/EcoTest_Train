# Functions for extracting indices



ICCATtoGEO<-function(Catchdat, SquareTypeCode = "5x5", FleetCode = "JPN"){   # Convert ICCAT format (corner closest to GMT/equator) to South West corner

  ICCATdat = Catchdat[Catchdat$SquareTypeCode == SquareTypeCode,]
  if(FleetCode != "all") ICCATdat = ICCATdat[ICCATdat$FleetCode == FleetCode,]
  
  Latadj<-as.numeric(unlist(strsplit(type,"x"))[1])/2
  Lonadj<-as.numeric(unlist(strsplit(type,"x"))[2])/2
                           
  ICCATdat$Lon[ICCATdat$QuadID==1] <-  as.numeric(ICCATdat$Lon[ICCATdat$QuadID==1]+Lonadj)
  ICCATdat$Lat[ICCATdat$QuadID==1] <- as.numeric(ICCATdat$Lat[ICCATdat$QuadID==1]+Latadj)
  
  ICCATdat$Lon[ICCATdat$QuadID==2] <- as.numeric(ICCATdat$Lon[ICCATdat$QuadID==2]+Lonadj)
  ICCATdat$Lat[ICCATdat$QuadID==2] <-- as.numeric(ICCATdat$Lat[ICCATdat$QuadID==2]-Latadj)
  
  ICCATdat$Lon[ICCATdat$QuadID==3] <-- as.numeric(ICCATdat$Lon[ICCATdat$QuadID==3]-Lonadj)
  ICCATdat$Lat[ICCATdat$QuadID==3] <-- as.numeric(ICCATdat$Lat[ICCATdat$QuadID==3]-Latadj)
  
  ICCATdat$Lon[ICCATdat$QuadID==4] <-- as.numeric(ICCATdat$Lon[ICCATdat$QuadID==4]-Lonadj)
  ICCATdat$Lat[ICCATdat$QuadID==4] <- as.numeric(ICCATdat$Lat[ICCATdat$QuadID==4]+Latadj)
  
  keep = !is.na(ICCATdat$Lat) & !is.na(ICCATdat$Lon); print(mean(keep)*100)
  
  CellID = paste(ICCATdat$Lat,ICCATdat$Lon,sep="_")
  ICCATdat=cbind(ICCATdat,CellID)
  ICCATdat
  
}


getquant = function(SCode, Catchdat, type="CR"){
  
  col = match(SCodes[x],names(Catchdat))
  keep = Catchdat[,col] != "NULL"
  cdat = Catchdat[keep,]
  Cat = as.numeric(Catchdat[keep,col])
  Cagg = aggregate(Cat, by = list(Yr = cdat$YearC, Lat = cdat$Lat, Lon = cdat$Lon, Cell = cdat$CellID), FUN = sum, na.rm=T)
  Eagg = aggregate(cdat$Eff1, by = list(Yr = cdat$YearC, Lat = cdat$Lat, Lon = cdat$Lon, Cell =cdat$CellID), sum, na.rm=T)
  
  if(type == "CR"){
    dat = Eagg
    dat$x = Cagg$x / Eagg$x
  }else if(type == "Catch"){
    dat = Cagg
  }else if(type == "Effort"){
    dat = Eagg
  }
  
  dat=dat[dat$x > 0,]
  dat
  
}


normdat = function(dat,ntype="none"){
  
  if(ntype == "annual fraction"){
    anntot = aggregate(dat$x,by=list(Yr=dat$Yr),sum)
    dat$x = dat$x/anntot$x[match(dat$Yr,anntot$Yr)]
  }
  dat
  
}

getindicator = function(datt, ref_lev = 0.5, itype = "percentage",yreval){
  
  cells = unique(datt$Cell)
  ncell = length(cells)
  ny = length(yreval)
  datarr = indicator = array(0,c(ncell,ny))
  datarr[as.matrix(cbind(match(datt$Cell,cells),match(datt$Yr,yreval)))]=datt$x
  
  if(itype == "percentage"){ # ncells that make up annual percentage
    for(yy in 1:ny){
      vec = datarr[,yy]
      ord = order(vec,decreasing=T)
      cums = cumsum(vec[ord])
      indy = as.numeric(cums<(ref_lev*max(cums)))
      indicator[ord,yy]<-indy
    }
  }else if(itype == "abs_mean"){
    indicator[datarr>(mean(datarr)*ref_lev)] = 1
  }
  
  indicator
}

probust = function(p, tiny = 1E-6){
  p[p<tiny] = tiny
  p[p>(1-tiny)] = (1-tiny)
  p
}
logit = function(p)log(probust(p)/(1-probust(p)))
ilogit = function(x)exp(x)/(1+exp(x))


# Pyrs = seq(1960,2010,by=5); ploty=T; yrs = 1950:2013
# type = "Catch"; ntype = "annual fraction" 
getlogitmod = function(x,SSBrel,Catchdat,SCodes, type = "CR",ntype = "none", itype = "percentage",
                       ref_lev = 0.5, 
                       Pyrs = seq(1960,2010,by=5), ploty="all", yrs = 1950:2013){
  
  dat = getquant(SCode = SCodes[x],Catchdat,type=type)
  dat = normdat(dat,ntype=ntype)
  dat=dat[!is.na(dat$x)&!(dat$x == "Inf"),]
  SSBr = SSBrel[[x]][1,]
  
  firstYr = min(Pyrs) #min(Nbyyr$Yr[Nbyyr$x > (0.95 * max(Nbyyr$x))])
  yreval = firstYr:max(yrs)
  yind = yrs%in%yreval
  minval = mean(dat$x)/10
  meanval = mean(dat$x)

  datt = dat[dat$Yr>=firstYr & dat$Yr <=max(yrs),]
  datt = datt[!is.na(datt$Yr),]
  indicator = getindicator(datt,itype = itype,ref_lev=ref_lev,yreval = yreval)
  itrend = apply(indicator,2,mean)
  
  mod = lm(y~x,data=data.frame(x=log(SSBr[yind]),y=logit(itrend)))
  SBseq=seq(0.01,5,length.out=1000)
  pred = predict(mod,newdata = data.frame(x=log(SBseq)))
 
  
  if(ploty == "all"){
    np = length(Pyrs)+3
    par(mfrow=c(ceiling(np/3),3),mai=c(0.5,0.5,0.1,0.1))
    plot(yrs[yind],SSBr[yind],type="l",lwd=2,ylim=c(0,max(SSBr)),ylab = "SSB / SSBMSY",xlab=""); grid(); abline(v=Pyrs,col="green") 
    legend('top',legend=SCodes[x],bty="n")
    
    plot(yreval,itrend,col="blue",ylim=c(0,max(itrend)),type="l",xlab="",ylab="Indicator (fraction)");grid();abline(v=Pyrs,col="green") 
    plot(SSBr[yind],itrend,col="white",pch=19,xlim = c(0,max(SSBr[yind])),ylim=c(0,max(itrend)),xlab = "SSB / SSBMSY",ylab="Fraction");grid()
    text(SSBr[yind],itrend,yreval,cex=0.8)
    lines(SBseq,ilogit(pred),col="blue",lwd=2)
    
    sapply(1:length(Pyrs),plotdat,dat=dat,Pyrs=Pyrs,minval=minval,meanval=meanval)
 
  }else if(ploty == "rel"){
    plot(SSBr[yind],itrend,col="white",pch=19,xlim = c(0,max(SSBr[yind])),ylim=c(0,max(itrend)),xlab = "SSB / SSBMSY",ylab="Fraction");grid()
    text(SSBr[yind],itrend,yreval,cex=0.8)
    legend('topleft',legend=SCodes[x],bty="n")
    legend('left',legend=c(paste("type:",type),paste("ntype:",ntype),paste("itype:",itype), paste("reflev:", ref_lev)),bty='n')
    lines(SBseq,ilogit(pred),col="blue",lwd=2)
  }
  
  mod
}


plotdat = function(x,dat,Pyrs,mincex = 0.2,addcex=1.8,minval,meanval){
  
  datrng = dat[dat$Yr %in% Pyrs,]
  dats = dat[dat$Yr == Pyrs[x],]
  
  xlim = range(datrng$Lon)+c(-2.5,2.5)
  ylim = range(datrng$Lat)+c(-2.5,2.5)
  datmax = quantile(datrng$x^0.5,0.99)
  
  plot(xlim,ylim,col="white",xlab="",ylab="")
  ls = (-100:100)*5
  abline(h = ls, v=ls,col="grey")
  
  cex = 0.2 + ((dats$x^0.5)/datmax)*addcex
  cols = rep("#0000ff95",nrow(dats))
  cols[dats$x < minval] = "black"
  cols[dats$x > meanval] = '#ff000095'
  points(dats$Lon, dats$Lat, cex=cex,pch=19,col=cols)
  legend('top',legend=Pyrs[x],bty="n")
  
}

