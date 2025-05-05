# Functions for extracting indices

proc_dat_LH<-function(MMSE, Iind=NA, sno = 1, fno=1, plotsmooth=F){
  
  dat = MMSE@PPD[[sno]][[fno]][[1]]
  stock = MMSE@Stocks[[sno]]
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
 
  datlist<-list()
  Iyr = Iind[,2]
  Byr = Iyr - nyears
  
  Bt = MMSE@SB_SBMSY
  OM = MMSE@OM[[sno]][[fno]]
  L50 = OM$L50
  Linf = OM$Linf
  L50_Linf = L50/Linf
  M = OM$M
  K = OM$K
  M_K = M/K
  maxa = -log(0.05)/M # age at 5% cumulative survival
  
  Bind<-cbind(1:nsim,rep(sno,nsim),rep(1,nsim),Byr)
  Brel = MMSE@SB_SBMSY[Bind]
  
  df = data.frame(Brel, K, M_K, maxa, L50_Linf)
  names(df) = paste0(names(df),"_s",sno)
  df
  
}

proc_dat_F<-function(MMSE, Iind=NA, sno = 1, fno=1, plotsmooth=F){
  
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
 
 
  I_cur<-I_mu<-I_rel<-I_g5<-I_g10<-I_g20<- 
    Isd1 <- Isd2 <- Isd3 <-
    C_cur<- C_mu<-C_rel<-C_g5<-C_g10<- C_g20<- 
    Csd <-
    ML_cur<-ML_mu<-ML_rel<-ML_g5<-ML_g10<-ML_g20<- ML_L50 <- ML_Linf <-
    MA_cur<-MA_mu<-MA_rel<-MA_g5<-MA_g10<-MA_g20<-
    MV_cur<- MV_mu<-MV_rel<- MV_g5<-MV_g10<-MV_g20<-
    FM_cur <- FM_mu <- FM_rel <- FM_g5 <- FM_g10 <-FM_g20 <- rep(NA,nsim)
  
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
    # Index
    Is<-smooth2(Iobs[i,],plot=plotsmooth)
    I_cur[i]<-Is[Iind[i,2]]
    I_mu[i] = mean(Iobs[i,nyears:Iind[i,2]])
    I_rel[i]<-I_cur[i] / mean(Is[1:Iind[i,2]])
    I_g5[i]<-slp3(Iobs[i,Iind[i,2]-(5:1)])
    I_g10[i]<-slp3(Iobs[i,Iind[i,2]-(10:1)])
    I_g20[i]<-slp3(Iobs[i,Iind[i,2]-(20:1)])
    
    ind = Iind[i,2] - (20:1)
    Isd1[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
    Isd2[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
    Isd3[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
    
    # Catch
    Cs = smooth2(Cobs[i,],plot=plotsmooth)
    Csd[i] = sd(smooth2(Cobs[i,],plot=plotsmooth,ret='resid'))
    C_cur[i]<-Cs[Iind[i,2]]
    C_mu[i]<-mean(Cobs[i,nyears:Iind[i,2]])
    C_rel[i] <- C_cur[i] / mean(Cs[1:Iind[i,2]],na.rm=T)
    C_g5[i]<-slp3(Cobs[i,Iyr[i]-(5:1)])
    C_g10[i]<-slp3(Cobs[i,Iyr[i]-(10:1)])
    C_g20[i]<-slp3(Cobs[i,Iyr[i]-(20:1)])
    
    # Age
    CAAi = dat@CAA[i,,]
    totage = CAAi*t(array(1:ncol(CAAi),dim(t(CAAi))))
    muage = apply(totage,1,sum)/apply(CAAi,1,sum)
    if(all(is.na(muage))){
      As = rep(NA,length(muage))
    }else{
      As = smooth2(muage,plot=plotsmooth, enp.mult = 0.125)
    }
    
    MA_cur[i]=As[Iind[i,2]]
    MA_mu[i]=mean(muage[nyears:Iind[i,2]])
    MA_rel[i] <- MA_cur[i] / mean(As[1:Iind[i,2]],na.rm=T)
    MA_g5[i]<-slp3(muage[Iyr[i]-(5:1)])
    MA_g10[i]<-slp3(muage[Iyr[i]-(10:1)])
    MA_g20[i]<-slp3(muage[Iyr[i]-(20:1)])
   
    # Length
    CALi = dat@CAL[i,,]
    totlen = CALi*t(array(mids,dim(t(CALi))))
    mulen = apply(totlen,1,sum)/apply(CALi,1,sum)
    if(all(is.na(mulen))){
      Ls = rep(NA,length(mulen))
    }else{
      Ls = smooth2(mulen,plot=plotsmooth, enp.mult = 0.125)
    }
    
    ML_cur[i]=Ls[Iind[i,2]]
    ML_mu[i]=mean(mulen[nyears:Iind[i,2]])
    ML_rel[i] <- ML_cur[i] / mean(Ls[1:Iind[i,2]],na.rm=T)
    ML_g5[i]<-slp3(mulen[Iyr[i]-(5:1)])
    ML_g10[i]<-slp3(mulen[Iyr[i]-(10:1)])
    ML_g20[i]<-slp3(mulen[Iyr[i]-(20:1)])
    ML_Linf[i]= ML_cur[i]/Linf[i]
    ML_L50[i]= ML_cur[i]/L50[i]
    
    # standard deviation in length
    CVlen = sapply(1:dim(CALi)[1],function(X,CALi,mids){
      sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
    },CALi=CALi,mids=mids)
    
    if(all(is.na(CVlen))){
      CVs = rep(NA,length(CVlen))
    }else{
      CVs = smooth2(CVlen, plot=plotsmooth,enp.mult = 0.1)
    }
    
    MV_cur[i]=CVs[Iind[i,2]]
    MV_mu[i] = mean(CVlen[nyears:Iind[i,2]],na.rm=T)
    MV_rel[i] <- MV_cur[i] / mean(CVs[1:Iind[i,2]])
    MV_g5[i]<-slp3(CVlen[Iyr[i]-(5:1)])
    MV_g10[i]<-slp3(CVlen[Iyr[i]-(10:1)])
    MV_g20[i]<-slp3(CVlen[Iyr[i]-(20:1)])
    
    L50i = L50[i]
    Fmat = sapply(1:dim(CALi)[1],function(X,CALi,mids,L50i){
      indivs = rep(mids,floor(CALi[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50i,na.rm=T)
    },CALi=CALi,mids=mids,L50i=L50i)
    
    if(all(is.na(Fmat))){
      FMs = rep(NA,length(Fmat))
    }else{
      FMs = smooth2(Fmat, plot=plotsmooth, enp.mult=0.1)
    }
    
    FM_cur[i]= FMs[Iind[i,2]]
    FM_mu[i] = mean(FMs[nyears:Iind[i,2]])
    FM_rel[i] <- FM_cur[i] / mean(FMs[1:Iind[i,2]])
    FM_g5[i]<-slp3(Fmat[Iyr[i]-(5:1)])
    FM_g10[i]<-slp3(Fmat[Iyr[i]-(10:1)])
    FM_g20[i]<-slp3(Fmat[Iyr[i]-(20:1)])
    
  }
  
  df = data.frame(I_rel, I_g5, I_g10, I_g20, 
             Isd1, Isd2, Isd3,
             C_rel, C_g5, C_g10, C_g20, Csd, 
             ML_cur, ML_rel, ML_g5, ML_g10, ML_g20, ML_L50, ML_Linf,
             MV_cur, MV_rel, MV_g5, MV_g10, MV_g20,
             FM_cur, FM_rel, FM_g5, FM_g10, FM_g20,
             MA_cur, MA_rel, MA_g5, MA_g10, MA_g20, 
             L5_L50, LFS_L50, VML)
  names(df) = paste0(names(df),"_s",sno,"_f",fno)
  df
}

cat_ratios_2 = function(MMSE, Iind, fno=1, plotsmooth=F, refyrs = 20){
  
  ns = MMSE@nstocks
  nsim = MMSE@nsim
  allyears = MMSE@proyears+MMSE@nyears
  nyears = MMSE@nyears
  
  CRtab = NULL
  
  for(i in 1:nsim){ 
    
    catsmth = catsmth2 = array(NA,c(ns,allyears-1))
    
    for(sno in 1:ns){
      catsmth[sno,] = smooth2(MMSE@PPD[[sno]][[fno]][[1]]@Cat[i,],plot=plotsmooth)
      catsmth2[sno,] = smooth2(MMSE@PPD[[sno]][[fno]][[1]]@Cat[i,],enp.mult = 0.075,ret='resid',plot=plotsmooth) # log residuals after broad detrending
    }
    
    ncomp = ((ns-1)*ns)/2
    
    CR_rel =  CR_g5 = CR_g10 = CR_g20 = CR_mu = CR = namy = rep(NA,ncomp)
    
    sampyr = Iind[i,2]
    j = 0
    for(ss in 1:(ns-1)){
      for(s2 in (ss+1):ns){
        
        j=j+1
        namy[j] = paste0("s",ss,"_s",s2)
        rat = catsmth[ss,]/catsmth[s2,]
        CR[j] = rat[sampyr]
        invalidrat = rat<(0.01*mean(rat))&((1:length(rat))<refyrs)
        rat2 = rat; rat2[invalidrat] = NA
        CR_mu[j] = mean(rat2[1:refyrs],na.rm=T)
        CR_rel[j] = rat2[sampyr] / mean(rat2[1:refyrs],na.rm=T)
        CR_g5[j] = slp3(rat2[sampyr-(5:1)])
        CR_g10[j] = slp3(rat2[sampyr-(10:1)])
        CR_g20[j] = slp3(rat2[sampyr-(20:1)])
        
      }
    }  
    
    CRvec = c(CR,CR_mu,CR_rel,CR_g5,CR_g10, CR_g20)
    names(CRvec) = paste(rep(c("CR","CR_mu","CR_rel","CR_g5","CR_g10","CR_g20"),each=ncomp),
                       rep(namy,5),sep="_")
    
    
    CC20 = CC40 = CC60 = rep(NA, ncomp)
    
    j = 0
    ind20 = Iind[i,2] - (20:1)
    ind40 = Iind[i,2] - (40:21)
    ind60 = Iind[i,2] - (60:41)
    
    for(ss in 1:(ns-1)){
      for(s2 in (ss+1):ns){
        
        j=j+1
        CC20[j] = cor(catsmth2[ss,ind20],catsmth2[s2,ind20])
        CC40[j] = cor(catsmth2[ss,ind40],catsmth2[s2,ind40]) 
        CC60[j] = cor(catsmth2[ss,ind60],catsmth2[s2,ind60]) 
        
      }
    }  
    
    CCvec = c(CC20,CC40,CC60)
    names(CCvec) = paste(rep(c("CC20","CC40","CC60"),each=ncomp),
                         rep(namy,3),sep="_")
    
    
                   
    CRtab = rbind(CRtab, c(CRvec,CCvec))  
    #cat(".")
    
  } # simulation
  #cat("\n")
  CRtab = as.data.frame(CRtab)
  names(CRtab) = paste0(names(CRtab),"_f",fno)
  CRtab
  
}






get_sim_data_2 = function(ff,filelocs){
  
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
  for(sno in 1:ns)outs[[sno]] = proc_dat_LH(MMSE,Iind=Iind,sno=sno)
  i=sno
  for(sno in 1:ns){
    for(fno in 1:nf){
      i=i+1
      outs[[i]] = proc_dat_F(MMSE,Iind=Iind,sno = sno, fno=fno)
    } 
  }
  for(fno  in 1:nf){
    i=i+1
    outs[[i]] = cat_ratios_2(MMSE, Iind=Iind, fno=fno)
  }
  
  for(j in 1:i){
    out = outs[[j]]
    #if(j <=ns) names(out) = paste0(names(out),"_",i)
    if(j==1) one_tab=out
    if(j>1) one_tab=cbind(one_tab,out)
  }
  row.names(one_tab) = paste0("sim_",1:nrow(one_tab))
  
  cat(".")
  one_tab
  
}

process_sim_data_2 = function(MSEdir, parallel=T, cores = NA){
  
  files = list.files(MSEdir)
  keep = grepl("MMSE",files)
  filelocs = list.files(MSEdir,full.names=T)[keep]
  nfile = length(filelocs)
  
  if(parallel){
    library(snowfall)
    library(parallel)
     if(is.na(cores))cores = detectCores()/2
     sfInit(parallel=T,cpus = cores)
     sfExport(list = c("proc_dat_LH","proc_dat_F","cat_ratios_2","smooth2","slp3","smooth3"))
     allout = sfLapply(1:nfile,get_sim_data_2,filelocs=filelocs)
     # test = sfLapply(1, get_sim_data_2, filelocs=filelocs)
  }else{
    allout = lapply(1:nfile, get_sim_data_2, filelocs=filelocs)
  }
  
  cat("\n")
  allout
  
}



