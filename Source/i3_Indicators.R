# Functions for extracting indices
#  xx=exp(rlnorm(70,0,0.5)*sin(seq(0,12,length.out=70))); enp.mult=0.2; nint=30; plotname=""; plot=T
interpolate<-function(xx,plot=F,enp.mult=0.2,nint=30,plotname=""){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  allpredout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  
  allpredout[tofill]<-exp(predict(out))
  intx = seq(1,length(xx),length.out=nint)
  interpolated = exp(predict(out,newdata=data.frame(x=intx)))
  out = interpolated/mean(interpolated)
  if(plot){
    plot(xx,type="p",xlab="x",ylab="y",main=plotname)
    lines(allpredout,col="#ff000090",lwd=2)
    points(intx,interpolated,col="blue",pch=19)
  }
  out
  
}

# MMSE = readRDS('C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_801_1000/MMSE_958.rda');Iind=matrix(c(1,89),nrow=1); sno = 1; fno=1; plotsmooth=T; nint = 40

proc_dat_F3<-function(MMSE, Iind=NA, sno = 1, fno=1, plotsmooth=F, nint = 40){
  
  dat = MMSE@PPD[[sno]][[fno]][[1]]
  stock = MMSE@Stocks[[sno]]
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
  
  datlist<-list()
  
  Iobs<-dat@VInd
  Cobs<-dat@Cat
  mids = dat@CAL_mids
  
  Iyr = Iind[,2]
  Byr = Iyr - nyears

  Isd1 <- Isd2 <- Isd3 <- Csd1 <- Csd2 <- Csd3 <- CF <- ML_Linf <- rep(NA, nsim)
  Iint = Cint = MAint = MLint = MVint = FMint= array(NA,c(nsim,nint))
  
  matage = MMSE@multiHist[[sno]][[fno]]@AtAge$Maturity[,,1]
  lenage = MMSE@multiHist[[sno]][[fno]]@AtAge$Length[,,1]
  
  OM = MMSE@OM[[sno]][[fno]]
  
  L50 = OM$L50
  Linf = OM$Linf
  M = OM$M
  K = OM$K
  M_K = M/K
  maxa = -log(0.05)/M # age at 5% cumulative survival
  
  L5 = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$L5_y[,1]
  LFS = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$LFS_y[,1]
  VML = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$Vmaxlen_y[,1]
  L5_L50 = L5/L50
  LFS_L50 = LFS/L50
  
  for(i in 1:nsim){
    ny = Iind[i,2]
    tind = 1:ny
    # Index
    Iint[i,] = interpolate(Iobs[i,tind],plotsmooth,enp.mult=0.2,nint=nint,plotname="Index")
    ind = Iind[i,2] - (20:1)
    Isd1[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
    Isd2[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
    Isd3[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
    
    # Catch
    Cint[i,] = interpolate(Cobs[i,tind],plotsmooth,enp.mult=0.2,nint=nint,plotname="Catch")
    allfleetC = sapply(MMSE@PPD[[sno]],function(x,i)sum(x[[1]]@Cat[i,]),i=i)
    CF[i] = allfleetC[fno]/sum(allfleetC)  
    Csd1[i] = sd(smooth2(Cobs[i,ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
    Csd2[i] = sd(smooth2(Cobs[i,ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
    Csd3[i] = sd(smooth2(Cobs[i,ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
    
    # Age
    CAAi = dat@CAA[i,,]
    totage = CAAi*t(array(1:ncol(CAAi),dim(t(CAAi))))
    muage = apply(totage,1,sum)/apply(CAAi,1,sum)
    if(all(is.na(muage)))muage[]=NA
    MAint[i,] = interpolate(muage[tind],plotsmooth,enp.mult=0.125,nint=nint,plotname="Mean Age")
   
    # Length
    CALi = dat@CAL[i,,]
    totlen = CALi*t(array(mids,dim(t(CALi))))
    mulen = apply(totlen,1,sum)/apply(CALi,1,sum)
    if(all(is.na(mulen)))mulen[]<-NA 
    MLint[i,] = interpolate(mulen[tind],plotsmooth,enp.mult=0.125,nint=nint,plotname="Mean Length")
    Ls = smooth2(mulen[tind],plotsmooth,enp.mult = 0.125) 
    ML_Linf[i]= Ls[ny]/Linf[i]
    
    # standard deviation in length
    CVlen = sapply(1:dim(CALi)[1],function(X,CALi,mids){
      sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
    },CALi=CALi,mids=mids)
    if(all(is.na(CVlen)))CVlen[]=NA
    MVint[i,] =  interpolate(CVlen[tind],plotsmooth,enp.mult=0.125,nint=nint,plotname="SD Length")
   
    # Fraction mature
    L50i = L50[i]
    Fmat = sapply(1:dim(CALi)[1],function(X,CALi,mids,L50i){
      indivs = rep(mids,floor(CALi[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50i,na.rm=T)
    },CALi=CALi,mids=mids,L50i=L50i)
    if(all(is.na(Fmat)))Fmat[]=NA
    FMint[i,] =  interpolate(Fmat[tind],plotsmooth,enp.mult=0.125,nint=nint,plotname="Fraction mature")
   
  }
 
  dfsingle = data.frame(Isd1, Isd2, Isd3, Csd1, Csd2, Csd3, CF, ML_Linf, L5_L50, LFS_L50, VML)
  dfmulti = data.frame(Iint, Cint, MAint, MLint, MVint, FMint)
  names(dfmulti) = paste0(rep(c("I","C","MA","ML","MV","FM"),each=nint),"_",rep(1:nint,6))
  df=cbind(dfsingle, dfmulti)
  names(df) = paste0(names(df),"_s",sno,"_f",fno)
  df
}




get_sim_data_3 = function(ff,filelocs){
  
  set.seed(ff)
  MMSE = readRDS(filelocs[ff])
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
  Iyr<-sample((nyears+1):(allyears-2),nsim,replace=T)
  Byr = Iyr - nyears
  Iind<-cbind(1:nsim,Iyr)
  
  ns=MMSE@nstocks
  nf = MMSE@nfleets
  outs = list()
  
  for(sno in 1:ns)outs[[sno]] = temp = proc_dat_LH(MMSE,Iind=Iind,sno=sno)
  i=sno
  for(sno in 1:ns){
    for(fno in 1:nf){
      i=i+1
      cat(".")
      outs[[i]] = proc_dat_F3(MMSE,Iind=Iind,sno = sno, fno=fno)
    } 
  }
  
  for(j in 1:i){
    out = outs[[j]]
    #if(j <=ns) names(out) = paste0(names(out),"_",i)
    if(j==1) one_tab=out
    if(j>1) one_tab=cbind(one_tab,out)
  }
  ny = Iyr
  one_tab = cbind(one_tab,ny)
  #one_tab = data.table::cbindlist(outs)
  row.names(one_tab) = paste0("sim_",1:nrow(one_tab))
  cat(".")
  one_tab
  
}

process_sim_data_3 = function(MSEdir, parallel=T, cores = NA){
  
  files = list.files(MSEdir)
  keep = grepl("MMSE",files)
  filelocs = list.files(MSEdir,full.names=T)[keep]
  nfile = length(filelocs)
  
  if(parallel){
    library(snowfall)
    library(parallel)
     if(is.na(cores))cores = detectCores()/2
     sfInit(parallel=T,cpus = cores)
     sfExport(list = c("proc_dat_LH","proc_dat_F3","smooth2","interpolate","slp3","smooth3"))
     allout = sfLapply(1:nfile,get_sim_data_3,filelocs=filelocs)
     # test = sfLapply(1, get_sim_data_2, filelocs=filelocs)
  }else{
    allout = lapply(1:nfile, get_sim_data_3, filelocs=filelocs)
  }
  
  cat("\n")
  allout
  
}

debugy = function(){
  for(ff in 12:nfile) test = get_sim_data_3(ff,filelocs)
  
}


