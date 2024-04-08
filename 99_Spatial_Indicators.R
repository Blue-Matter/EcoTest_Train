# Functions for extracting indices



ICCATtoGEO<-function(Catchdat, type = "5x5"){   # Convert ICCAT format (corner closest to GMT/equator) to South West corner
  
  ICCATdat = Catchdat[Catchdat$SquareTypeCode == type,]
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




fit_logit = function(x,Iind,Catchdat){
  
  
  
  
}




process_spatial_data = function(MMSE, Catchdat, Iind){
  
   dim(MMSE@SB_SBMSY)
 
  
  cat("\n")
  allout
  
}



get_spatial_sim_data = function(ff,filelocs){
  
  set.seed(ff)
  MMSE = readRDS(filelocs[ff])
  
  nsim <- MMSE@nsim
  nyears <- MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
  Iyr <- sample((nyears+1):(allyears-2),nsim,replace=T)
  Byr = Iyr - nyears
  Iind <- cbind(1:nsim,Iyr)
  
  nf=MMSE@nfleets
  ns=MMSE@nstocks
  outs = list()
  for(ss in 1:ns)outs[[ss]] = proc_dat(MMSE,Iind=Iind,sno=ss)
  outs[[ns+1]] = cat_ratios(MMSE, Iind=Iind)
  
  for(i in 1:(ns+1)){
    out = outs[[i]]
    if(i <=ns) names(out) = paste0(names(out),"_",i)
    if(i==1)one_tab=out
    if(i>1)one_tab=cbind(one_tab,out)
  } 
  
  cat(".")
  one_tab
  
}

