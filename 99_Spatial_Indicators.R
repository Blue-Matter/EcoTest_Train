# Functions for extracting indices



ICCATtoGEO<-function(Catchdat, SquareTypeCode = "5x5", FleetCode = "JPN"){   # Convert ICCAT format (corner closest to GMT/equator) to South West corner

  ICCATdat = Catchdat[Catchdat$SquareTypeCode == SquareTypeCode,]
  if(FleetCode != "all") ICCATdat = ICCATdat[ICCATdat$FleetCode == FleetCode,]
  
  Latadj<-as.numeric(unlist(strsplit(SquareTypeCode,"x"))[1])/2
  Lonadj<-as.numeric(unlist(strsplit(SquareTypeCode,"x"))[2])/2
                           
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
  
  col = match(SCode,names(Catchdat))
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
  datarr = array(NA,c(ncell,ny))
  indicator = array(0,c(ncell,ny))
  datarr[as.matrix(cbind(match(datt$Cell,cells),match(datt$Yr,yreval)))]=datt$x
  Lat = datt$Lat[match(cells,datt$Cell)]
  Lon = datt$Lon[match(cells,datt$Cell)]
  #cbind(cells, Lat, Lon)
  coords = data.frame(Lat = Lat, Lon = Lon)
  
  if(itype == "percentage" | itype == "positive percentage"){ # ncells that make up annual percentage
    for(yy in 1:ny){
      vec = datarr[,yy]
      if(!all(is.na(vec))){
        ord = order(vec,decreasing=T)
        cums = cumsum(vec[ord])
        indy = as.numeric(cums<(ref_lev*max(cums,na.rm=T)))
        if(itype == "percentage")indy[is.na(indy)] = 0 # positive percentage carried NAs and is based on the spatial range of observations in that year
        indicator[ord,yy]<-indy
      }else{
        indicator[,yy] = 0
      }
    }
  }else if(itype == "abs_mean" | itype == "positive abs_mean"){
    
    cond = datarr>(mean(datarr,na.rm=T)*ref_lev)
    indicator[] = as.numeric(cond)
    if(itype == "abs_mean")indicator[is.na(indicator)] = 0 
    for(yy in 1:ny){
      if(all(is.na(datarr[,yy])))indicator[,yy]=0
    }  
  }
  
  list(indicator = indicator, datarr = datarr, coords = coords)
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
  yobs = yrs[yind]
 
  datt = dat[dat$Yr>=firstYr & dat$Yr <=max(yrs),]
  datt = datt[!is.na(datt$Yr),]
  indicators = getindicator(datt, itype = itype, ref_lev = ref_lev, yreval = yreval)
  itrend = apply(indicators$indicator,2,mean,na.rm=T)
  
  obsyrs = range(datt$Yr)
  yfit = yrs>=obsyrs[1] & yrs<=obsyrs[2]
  yfiteval = yreval>=obsyrs[1] & yreval<=obsyrs[2]
  mod = lm(y~x,data=data.frame(x=log(SSBr[yfit]),y=logit(itrend[yfiteval])))
  SBseq=seq(0.01,5,length.out=1000)
  pred = predict(mod,newdata = data.frame(x=log(SBseq)))
 
  
  if(ploty == "all"){
    
    np = length(Pyrs)+3
    par(mfrow=c(ceiling(np/3),3),mai=c(0.5,0.5,0.1,0.1))
    plot(yrs[yind],SSBr[yind],type="l",lwd=2,ylim=c(0,max(SSBr)),ylab = "SSB / SSBMSY",xlab=""); grid(); abline(v=Pyrs,col="green") 
    legend('top',legend=SCodes[x],bty="n")
     
    plot(yreval,itrend,col="blue",ylim=c(0,max(itrend)),type="l",xlab="",ylab="Indicator (fraction)");grid();abline(v=Pyrs,col="green") 
    legend('bottomleft',text.col = "red",legend=c(paste("type:",type),paste("ntype:",ntype),paste("itype:",itype), paste("reflev:", ref_lev)),bty='n')
    
    plot(SSBr[yind],itrend,col="white",pch=19,xlim = c(0,max(SSBr[yfit])),ylim=c(0,max(itrend)),xlab = "SSB / SSBMSY",ylab="Fraction");grid()
    text(SSBr[yfit],itrend[yfiteval],yreval[yfiteval],cex=0.8)
    lines(SBseq,ilogit(pred),col="blue",lwd=2)
    
    sapply(1:length(Pyrs),plotdat,dat=dat,Pyrs=Pyrs,indicators=indicators, yreval=yreval)
 
  }else if(ploty == "rel"){
    
    plot(SSBr[yind],itrend,col="white",pch=19,xlim = c(0,max(SSBr[yfit])),ylim=c(0,max(itrend)),xlab = "SSB / SSBMSY",ylab="Fraction");grid()
    text(SSBr[yfit],itrend[yfiteval],yreval[yfiteval],cex=0.8)
    legend('topleft',legend=SCodes[x],bty="n")
    legend('left',text.col='red',legend=c(paste("type:",type),paste("ntype:",ntype),paste("itype:",itype), paste("reflev:", ref_lev)),bty='n')
    lines(SBseq,ilogit(pred),col="blue",lwd=2)
    
  }
  
  mod
}


plotdat = function(x,dat,Pyrs,mincex = 0.4,addcex=2.6,indicators, yreval){
  
  datarr = indicators$datarr
  indicator = indicators$indicator
  coords = indicators$coords
  
  dats = datarr[,yreval == Pyrs[x]]
  
  inds = indicator[,yreval == Pyrs[x]]
  
  xlim = range(coords$Lon)+c(-2.5,2.5)
  ylim = range(coords$Lat)+c(-2.5,2.5)
  datmax = quantile(dats,0.99,na.rm=T)
  
  plot(xlim,ylim,col="white",xlab="",ylab="")
  ls = (-100:100)*5
  abline(h = ls, v=ls,col="grey")
  
  cex = mincex + ((dats/datmax)^0.5)*addcex
  cols = rep("#0000ff95",length(dats))
  cols[inds==1] = '#ff000095'
  points(coords$Lon, coords$Lat, cex=cex,pch=19,col=cols)
  legend('top',legend=Pyrs[x],bty="n")
  
}

