

library(keras3)
library(r4ss)
library(data.table)
library(miceadds)
library(readxl)
library(writexl)

fdir = "C:/GitHub/EcoTest_train"
tdir = "C:/GitHub/EcoTest"
ddir = "C:/GitHub/EcoTestData"
source.all(paste0(fdir,"/Source"))
source.all(paste0(tdir,"/R"))

xldir = "C:/GitHub/EcoTest_Train/Indicator_4/Real_data"
xlfiles = list.files(xldir)
dirs = paste0(xldir,"/",xlfiles)

TD[1:3,grepl("_s1",names(TD))][,1:15]

make_i4 = function(xlfile){
  
  # pars sheet
  pars_sheet = as.data.frame(read_xlsx(xlfile,1))
  maxa = as.numeric(pars_sheet$Value[2])
  M = log(0.05)/-maxa
  K = as.numeric(pars_sheet$Value[3])
  L50_Linf = as.numeric(pars_sheet$Value[4])
  M_K = M/K
  
  
  
  
  
  ts_sheet = read_xlsx(xlfile,2)
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
  
}


