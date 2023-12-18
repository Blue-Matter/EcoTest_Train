plotB<- function(x=NULL, maxcol = 6, qcol = rgb(0.4, 0.8, 0.95), lcol = "dodgerblue4",
                      quants = c(0.05, 0.25, 0.75, 0.95), curyr = 2018, addline = FALSE, ...) {
  MMSE <- x
  if (!methods::is(MMSE, 'MMSE')) stop('Object must be class `MMSE`')
  if(is.na(maxcol))maxcol=ceiling(length(MMSE@MPs)/0.5) # defaults to portrait 1:2
  MPs<-MMSE@MPs
  MPrefs<-MMSE@MPrefs
  nMPs<-length(MPrefs[,1,1])
  yrs<-curyr+(1:MMSE@proyears)
  ns<-MMSE@nstocks
  nf<-MMSE@nfleets
  
  plots<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxcol))
  
  # --- Biomass projection ---------------------------------------------------
  B_BMSY<-MMSE@SB_SBMSY
  Blims <- c(0,quantile(B_BMSY[B_BMSY!="Inf"],0.95))
  
  for(pp in 1:length(plots)){
    
    toplot<-plots[[pp]]
    nt<-length(toplot)
    par(mfcol=c(ns,nt),mai=c(0.3,0.3,0.3,0.01),omi=c(0.4,0.5,0.05,0.05))
    
    for(MP in toplot){
      
      for(ss in 1:ns){
        plot(range(yrs),Blims,col="white",yaxs="i")
        plotquant(B_BMSY[,ss,MP,],p=quants,yrs,qcol,lcol,ablines=c(0.5,1),addline=addline)
        mtext(paste(paste0("F",1:nf),MPrefs[MP,,ss],collapse=", "),3,line=0.2,font=2,cex=0.7)
        
        if(MP==toplot[1])mtext(MMSE@Snames[ss],2,line=2.5)
      }
    }
  }
  
  mtext("Projection Year",1,font=2,outer=T,line=1.2)
  mtext("Biomass / BMSY",2,font=2,outer=T,line=1.9)
  
}

plotF<- function(x=NULL, maxcol = 6, qcol = rgb(0.4, 0.8, 0.95), lcol = "dodgerblue4",
                 quants = c(0.05, 0.25, 0.75, 0.95), curyr = 2018, addline = FALSE, ...) {
  MMSE <- x
  if (!methods::is(MMSE, 'MMSE')) stop('Object must be class `MMSE`')
  if(is.na(maxcol))maxcol=ceiling(length(MMSE@MPs)/0.5) # defaults to portrait 1:2
  MPs<-MMSE@MPs
  MPrefs<-MMSE@MPrefs
  nMPs<-length(MPrefs[,1,1])
  yrs<-curyr+(1:MMSE@proyears)
  ns<-MMSE@nstocks
  nf<-MMSE@nfleets
  
  plots<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxcol))
 
  # --- F projection -----------------------------------------------------------
  
  F_FMSY<-MMSE@F_FMSY
  F_FMSYsum<-apply(F_FMSY,c(1,2,4,5),sum,na.rm=T)
  Flims<- c(0,quantile(F_FMSYsum,0.95,na.rm=T))
  
  for(pp in 1:length(plots)){
    
    toplot<-plots[[pp]]
    nt<-length(toplot)
    par(mfcol=c(ns,nt),mai=c(0.3,0.3,0.3,0.01),omi=c(0.4,0.5,0.05,0.05))
    
    for(MP in toplot){
      
      for(ss in 1:ns){
        plot(range(yrs),Flims,col="white",yaxs="i")
        plotquant(F_FMSYsum[,ss,MP,],p=quants,yrs,qcol,lcol,ablines=c(0.5,1),addline=addline)
        mtext(paste(paste0("F",1:nf),MPrefs[MP,,ss],collapse=", "),3,line=0.2,font=2,cex=0.7)
        
        if(MP==toplot[1])mtext(MMSE@Snames[ss],2,line=2.5)
      }
    }
  }
  
  mtext("Projection Year",1,font=2,outer=T,line=1.2)
  mtext("Total (all fleets) F / FMSY",2,font=2,outer=T,line=1.9)
  
}

