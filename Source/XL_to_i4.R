
proc_i4=function(xlfile, CVs = list(K = 0.1, maxa = 0.1, L50_Linf = 0.1, Linf = 0.025, 
                          Index = 0.2, Catch = 0.2, ML = 0.1, MA = 0.1, MV = 0.05, FM = 0.01,
                          L5_L50 = 0.1, LFS_L50 = 0.1, VML = 0.1), 
        plotsmooth = T, nint = 40, nsim = 100, retros = 6, extrapolate = T){
  
  out=list()
  for(i in 1:(retros+1)){
    out[[i]] <- make_i4(xlfile, CVs, plotsmooth, nint, nsim, peel=i-1, extrapolate)
  }
  names(out) = paste("peel", 0:retros, sep="_")
  out
}


make_i4 = function(xlfile, CVs = list(K = 0.1, maxa = 0.1, L50_Linf = 0.1, Linf = 0.025, 
                                      Index = 0.2, Catch = 0.2, ML = 0.1, MA = 0.1, MV = 0.05, FM = 0.01,
                                      L5_L50 = 0.1, LFS_L50 = 0.1, VML = 0.1), 
                   plotsmooth = T, nint = 40, nsim = 100, peel = 0, extrapolate = T){
  
  # Internal functions
  plu = function(ps, nam) as.numeric(ps[match(nam, ps[,1]), 2])  # parameter look up
  
  tsresamp = function(ats, CV, nsim, peel, normalz = F){ # time series resample
    nt = length(ats); 
    atp = ats[1:(nt-peel)]
    if(normalz)atp=atp/mean(atp,na.rm=T)
    np = length(atp)
    matrix(atp * trlnorm(nsim * np, 1, CV), byrow=T, nrow=nsim)
  } 
  
  getsdsmooth = function(x, enp.mult){
    if(sum(!is.na(x))>5){
      return(sd(smooth2(x, ret = 'resid', enp.mult = enp.mult),na.rm=T))
    }else{
      return(NA)
    }
  }
  
  ifunc = function(inp, enp.mult, nint = 40, extrapolate = T){ # function to interpolate (and potentially extrapolate) time series
    if(sum(!is.na(inp[1,]))>5){
      out = t(apply(inp,1,interpolate, enp.mult=enp.mult, nint=nint, extrapolate=extrapolate))
      return(out)
    }else{
      return(NULL)
    }
  }
  
  
  ps = as.data.frame(read_xlsx(xlfile,1))  # Parameter sheet
  tss = as.data.frame(read_xlsx(xlfile,2)) # Time series sheet
  
  # Get all params
  
  K = trlnorm(nsim, plu(ps, "K"),CVs$K)
  maxa = trlnorm(nsim, plu(ps, "max_age"), CVs$maxa)
  M = log(0.05)/-maxa
  M_K = M/K
  L50_Linf = trlnorm(nsim, plu(ps,"L50_Linf"), CVs$L50_Linf)
  Linf = trlnorm(nsim, plu(ps,"Linf"), CVs$Linf)
  
  CF_f1 = plu(ps, "Fleet_1_CF")
  CF_f2 = plu(ps, "Fleet_2_CF")
  CF_f3 = plu(ps, "Fleet_3_CF")
  
  ns = data.frame(K, M_K, maxa, L50_Linf, CF_f1, CF_f2, CF_f3)
  
  # Time series stuff
  
  fnams = c("Total_","Fleet_1_", "Fleet_2_","Fleet_3_")
  flabs = c("T", "f1", "f2","f3")
  
  allts = NULL

  for(ff in 1:length(fnams)){
    
    # Load time series
    Iobs = tsresamp(tss[,paste0(fnams[ff],"Index")], CV = CVs$Index, nsim, peel = peel, normalz = T)    # Mean 1
    Cobs = tsresamp(tss[,paste0(fnams[ff],"catch")], CV = CVs$Catch, nsim, peel = peel, normalz = T)    # Mean 1
    MLobs = tsresamp(tss[,paste0(fnams[ff],"mean_len")], CV = CVs$ML, nsim, peel = peel) / Linf         # Fraction of Linf
    MAobs = tsresamp(tss[,paste0(fnams[ff],"mean_age")], CV = CVs$MA, nsim, peel = peel) / maxa         # Fraction of maxage
    MVobs = tsresamp(tss[,paste0(fnams[ff],"CV_len")], CV = CVs$MV, nsim, peel = peel)                  # Unscaled (CV)
    FMobs = tsresamp(tss[,paste0(fnams[ff],"frac_mat")], CV = CVs$FM, nsim, peel = peel)                # Unscaled frac mature
    
    # Variability properties
    ind = ncol(Iobs) - (19:0)
    Isd1 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.4)
    Isd2 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.2)
    Isd3 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.1)
    Csd1 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.4)
    Csd2 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.2)
    Csd3 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.1)
    
    VMLmu = plu(ps, paste0(fnams[ff],"VML"))
    VMLsd = CVs$VML
    
    df = data.frame(Isd1, Isd2, Isd3, Csd1, Csd2, Csd3, 
                    L5_L50 = trlnorm(nsim, plu(ps, paste0(fnams[ff],"L5_L50")), CVs$L5_L50),
                    LFS_L50 = trlnorm(nsim, plu(ps, paste0(fnams[ff],"LFS_L50")), CVs$LFS_L50),
                    VML = rbeta(nsim, alphaconv(VMLmu, VMLsd),betaconv(VMLmu, VMLsd)))
    colnames(df) = paste0(c("Isd1","Isd2","Isd3","Csd1","Csd2","Csd3","L5_L50","LFS_L50","VML"),"_",flabs[ff])
    df = df[,!is.na(df[1,])] # remove NA columns (has to be all rows)
    
    # Interpolated time series
    if(plotsmooth)  par(mfrow=c(2,3))
   
    Iint = ifunc(Iobs, 0.2, nint)
    Cint = ifunc(Cobs, 0.2, nint)
    MLint = ifunc(MLobs, 0.125, nint)
    MAint = ifunc(MAobs, 0.125, nint)
    MVint = ifunc(MVobs, 0.125, nint)
    FMint = ifunc(FMobs, 0.125, nint)
    dimo = function(x){if(class(x)[1]=="matrix")return(ncol(x)); if(is.null(x))return(0)}
    times=c(dimo(Iint), dimo(Cint), dimo(MAint), dimo(MLint), dimo(MVint), dimo(FMint))
    labs = paste(rep(c("I", "C", "MA", "ML", "MV","FM"), times),1:40, flabs[ff],sep="_")
    tsdat = as.data.frame(cbind(Iint,Cint,MAint, MLint, MVint, FMint))
    names(tsdat) = labs
    if(ff==1) allts = cbind(df, tsdat)
    if(ff>1) allts = cbind(allts, df, tsdat)
    
  }   
  
  out = cbind(ns,allts)
  out
  
}
