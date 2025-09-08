
what_fleet = function(dir){
  replist= SS_output(dir)
  SSplotCatch(replist)
  SSplotIndices(replist)
  SSplotData(replist)
  
}

ss_names=function(io){
  io$inputs$dat$fleetinfo$fleetname
  #unique(io$inputs$dat$CPUE$index)
  #unique(io$inputs$dat$catch$fleet)
}


getmucv = function(comps,lbins, ndraw=20){
  
  nrec = nrow(comps)
  nbin = length(lbins)
  suml = t(t(comps)*lbins)
  totl = apply(suml,1,sum)
  mu = totl / apply(comps,1,sum)
  cv = rep(NA,length(mu))
  mlsamp = array(NA,c(length(ml),ndraw))
  
  for(rr in 1:nrec){
    nobs = sum(comps[rr,])
    samp=rmultinom(ndraw,nobs,prob=comps[rr,]/sum(comps[rr,]))
    sums = samp*lbins
    tots = apply(sums,2,sum)
    cv[rr] = sd(tots/nobs)/mu[rr]
  }
  
  list(mu=mu,cv=cv)
}

#  i= 1; io = ios[[i]]; Fnam = Fnams[[i]]; Inam = Inams[[i]]; nsamp = 10; catch_CV = 0.1; L50_CV = 0.05; Linf_CV = 0.025; K_CV = 0.05; M_CV = 0.05; plotsmooth=T
SS_2_ET = function(io, Fnam = c("F4_JPN_LL","F8_ESP_LL"), Inam = c("S2_JPN_LATE","S4_EU_ESP"),
                   nsamp=10, catch_CV = 0.05, L50_CV = 0.05, Linf_CV = 0.025, K_CV = 0.05, M_CV = 0.05, plotsmooth=T){
  
  control = io$inputs$control
  Fnames = io$inputs$dat$fleetinfo$fleetname 
  Find = match(Fnam, Fnames)
  Iind = match(Inam, Fnames)
 
  # get catches -----------
  dat = io$inputs$dat
  ct = dat$catch
  nc = nrow(ct)
  newind = rep(3,nc) # default to fleet 3
  for(ff in 1:2)newind[ct$fleet == Find[ff]]=ff  
  aggcatch = aggregate(ct$catch,by=list(year = ct$year, fleet = newind),sum)
  
  # get indices ----------
  ind = dat$CPUE
  ni = nrow(ind)
  newind = rep(3,ni)
  for(ff in 1:2)newind[ind$index == Iind[ff]]=ff  
  aggind = aggregate(ind$obs, by=list(year = ind$year, index = newind),mean)
  aggindcv = aggregate(ind$se_log, by=list(year = ind$year, index = newind),mean)
  
  
  lengthswitch = !is.null(dat$lencomp)
  if(lengthswitch){
    lbins = dat$lbin_vector
    lencomp = dat$lencomp[dat$lencomp$year>0,]
    
    labs = names(lencomp)
    bysex = paste0("f",lbins[1])%in%labs
    if(bysex){
      comps = lencomp[,match(paste0("f",lbins),labs)]+lencomp[,match(paste0("m",lbins),labs)]
    }
    
    # get mean length ----------
    mucv = getmucv(comps,lbins)  
    nl = nrow(lencomp)
    newind = rep(3,nl)
    for(ff in 1:2)newind[lencomp$fleet == Find[ff]]=ff  
    aggml = aggregate(mucv$mu, by=list(year = lencomp$year, fleet = newind),mean)
    aggmlcv = aggregate(mucv$cv, by=list(year = lencomp$year, fleet = newind),mean) # this isn't ideal as the aggregate fleet 3 should have a CV to reflect combination of many mean lengths
    #aggmlcv$x[is.na(aggmlcv$x)]=0
    
    # get cv in lengths
    CVlen = sapply(1:nl,function(X,comps,lbins){
      sd(rep(lbins,comps[X,]))/mean(rep(lbins,comps[X,]))
    },comps=comps,lbins=lbins) 
    
    aggsdl = aggregate(CVlen, by=list(year = lencomp$year, fleet = newind),mean)
    aggsdl=aggsdl[!is.na(aggsdl$x),]
     
    # get fraction mature
    Fmat = sapply(1:nl,function(X,comps,lbins,L50){
      indivs = rep(lbins,floor(comps[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50,na.rm=T)
    },comps=comps,lbins=lbins,L50=L50)
    
    aggFM = aggregate(Fmat, by=list(year = lencomp$year, fleet = newind),mean)
    
  }
  
  # selectivities
  selex = io$outputs$sizeselex
  sel = selex[selex$Yr==max(selex$Yr) & selex$Sex==1,]
  fleets = sel$Fleet
  nf = length(fleets)
  for(ff in 1:nf){
    ssel = out$sel[out$sel$Fleet ==fleets[ff] & out$sel$Factor =="Asel" ,]
    cind = match(ages, names(ssel))
    isel = ssel[ssel$Yr == min(ssel$Yr),]
    sel[ff, , ] = as.numeric(unlist(isel[,cind]))
  }
  
  
  L5 = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$L5_y[,1]
  LFS = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$LFS_y[,1]
  VML = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$Vmaxlen_y[,1]
  L5_L50 = L5/L50
  LFS_L50 = LFS/L50
  
  
  
  
  
  
  
  I_cur<-I_mu<-I_rel<-I_g5<-I_g10<-I_g20<- 
    Isd1 <- Isd2 <- Isd3 <-
    C_cur<- C_mu<-C_rel<-C_g5<-C_g10<- C_g20<- 
    Csd <-
    ML_cur<-ML_mu<-ML_rel<-ML_g5<-ML_g10<-ML_g20<- ML_L50 <- ML_Linf <-
    MA_cur<-MA_mu<-MA_rel<-MA_g5<-MA_g10<-MA_g20<-
    MV_cur<- MV_mu<-MV_rel<- MV_g5<-MV_g10<-MV_g20<-
    FM_cur <- FM_mu <- FM_rel <- FM_g5 <- FM_g10 <-FM_g20 <- rep(NA,nsamp)
  
  MGp = control$MG_parms
  L50 = rlnorm(nsamp, MGp[grepl("Mat50",rownames(MGp)),3],L50_CV)
  Linf = rlnorm(nsamp, MGp[(1:nrow(MGp))[grepl("L_at_Amax",rownames(MGp))][1],3],Linf_CV)
  K =  rlnorm(nsamp, MGp[(1:nrow(MGp))[grepl("VonBert_K",rownames(MGp))][1],3], K_CV)
  if("natM"%in%names(control)){
    Ma = as.numeric(control$natM[1,,drop=T])
    Sa = exp(-cumsum(Ma))
    M = rlnorm(nsamp,sum(Sa*Ma)/sum(Sa),M_CV) # survival weighted M
  }else{
    M =  rlnorm(nsamp,MGp[(1:nrow(MGp))[grepl("natM",rownames(MGp))][1],3],M_CV)
  }
  M_K = M/K
  maxa = -log(0.05)/M # age at 5% cumulative survival
  
  
  
  fleetdat = list()
  
  for(ff in 1:3){
    
    # catch
    fcat = aggcatch$x[aggcatch$fleet==ff]
    nycat = length(fcat)
    
    # cpue
    find = aggind$x[aggind$index==ff]
    findcv = aggindcv$x[aggind$index==ff]
    nyind = length(find)
    
    
    if(lengthswitch){
      # mean length
      fml = aggml$x[aggml$fleet ==ff]
      fmlcv = aggmlcv$x[aggmlcv$fleet==ff]
      nyml = length(fml)
      
      #sd length
      fsdl = aggsdl$x[aggsdl$fleet ==ff]
      nysdl = length(fsdl)
      
      # frac mature
      fFM = aggFM$x[aggFM$fleet==ff]
      nyFM = length(fFM)
      
    }  
    
    for(i in 1:nsamp){
      # Indices
      findsim = rlnorm(nyind, log(find),findcv)
      Is = smooth2(findsim,plot=plotsmooth)
      Isd[i] = sd(smooth2(findsim,plot=plotsmooth,ret='resid'))
      I_cur[i]<-Is[ny] 
      #C_mu[i]<-mean(Cobs[i,ny:Iind[i,2]]) # no info
      I_rel[i] <- I_cur[i] / mean(Is[1:nyind],na.rm=T)
      I_g5[i]<-slp3(findsim[nyind-(5:1)])
      if(nyind>=10)I_g10[i]<-slp3(findsim[nyind-(10:1)])
      if(nyind>=20)I_g20[i]<-slp3(findsim[nyind-(20:1)])
      ind = nyind - (min(20,nyind):1) +1
      Isd1[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
      Isd2[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
      Isd3[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
      
      # Catches 
      fcatsim = rlnorm(nycat,log(fcat),catch_CV)
      Cs = smooth2(fcatsim,plot=plotsmooth)
      Csd[i] = sd(smooth2(fcatsim,plot=plotsmooth,ret='resid'))
      C_cur[i]<-Cs[ny] 
      #C_mu[i]<-mean(Cobs[i,ny:Iind[i,2]]) # no info
      C_rel[i] <- C_cur[i] / mean(Cs[1:nycat],na.rm=T)
      C_g5[i]<-slp3(fcatsim[nycat-(5:1)])
      C_g10[i]<-slp3(fcatsim[nycat-(10:1)])
      C_g20[i]<-slp3(fcatsim[nycat-(20:1)])
      
      if(lengthswitch){
        # Mean length
        mlsim = rlnorm(nyml,log(fml),fmlcv)
        Ls = smooth2(mlsim,plot=plotsmooth, enp.mult = 0.125)
        ML_cur[i]=Ls[nyml]
        #ML_mu[i]=mean(mulen[nyears:Iind[i,2]])
        ML_rel[i] <- ML_cur[i] / mean(Ls[1:nyml],na.rm=T)
        ML_g5[i]<-slp3(mlsim[nyml-(5:1)])
        if(nyml>=10) ML_g10[i]<-slp3(mlsim[nyml-(10:1)])
        if(nyml>=20) ML_g20[i]<-slp3(mlsim[nyml-(20:1)])
        ML_Linf[i]= ML_cur[i]/Linf[i]
        ML_L50[i]= ML_cur[i]/L50[i]
        
        # cv length
        CVs = smooth2(fsdl, plot=plotsmooth,enp.mult = 0.1)
        MV_cur[i]=CVs[nysdl]
        #MV_mu[i] = mean(fsdl[nyears:Iind[i,2]],na.rm=T)
        MV_rel[i] <- MV_cur[i] / mean(CVs[1:nysdl])
        MV_g5[i]<-slp3(fsdl[nysld-(5:1)])
        if(nysdl>=10)MV_g10[i]<-slp3(fsdl[nysdl-(10:1)])
        if(nysdl>=20)MV_g20[i]<-slp3(fsdl[nysdl-(20:1)])
        
        # fraction mature
        FMs = smooth2(fFM, plot = plotsmooth, enp.mult = 0.1)
        FM_cur[i]= FMs[nyFM]
        #FM_mu[i] = mean(FMs[nyears:Iind[i,2]])
        FM_rel[i] <- FM_cur[i] / mean(FMs[1:nyFM])
        FM_g5[i]<-slp3(Fmat[nyFM-(5:1)])
        if(nyFM >=10) FM_g10[i]<-slp3(Fmat[nyFM-(10:1)])
        if(nyFM >=20) FM_g20[i]<-slp3(Fmat[nyFM-(20:1)])
        
        
      }
      
    }
    
    find = 
  }
  
  
}



all_ss3 = function(dir){
  inputs = get_ss3_inputs(dir)
  SSout = SS_output(dir)
  outputs = extract_SS(SSout)
  list(inputs=inputs,outputs=outputs)
}

get_ss3_inputs = function(dir){

  dat = SS_readdat(paste0(dir,"/data.ss"))
  control = SS_readctl(paste0(dir,"/control.ss"),datlist=dat)
  starter = SS_readstarter(paste0(dir,"/starter.ss"))
  forecast = SS_readforecast(paste0(dir,"/forecast.ss"))
  list(dat=dat, control=control, starter = starter, forecast = forecast)

}


# fit = SS_output(dir = here('assessment/runs/straw_dog'))    
extract_SS<-function(fit, name="Unnamed"){
  if(class(fit) != "list")  fit = readRDS(fit) 
  
  out<-list()
  out$model_name<-name
  # Estimates
  
  out$dep = fit$current_depletion
  
  # dimensions
  ts = fit$timeseries
  out$startyr = fit$startyr
  out$endyr = fit$endyr
  out$ages = fit$agebins
  out$yrs = ts$Yr
  
  # MLE time series
  out$SSB = ts$SpawnBio
  out$B = ts$Bio_all
  out$SSN = ts$mature_num
  out$ind = fit$cpue[,names(fit$cpue) %in% c("Fleet_name","Yr","Vuln_bio","Obs","Exp","Calc_Q","SE","Dev","Like")]
  out$FM = fit$exploitation[,c(1,4:ncol(fit$exploitation))]
  out$sel = fit$ageselex
  out$M = fit$Natural_Mortality
  out$catch = fit$catch[,names(fit$catch) %in% c("Fleet_name","Yr","Obs","Exp","se","F","Like")]
  NAA = fit$natage[,names(fit$natage)%in%c("Sex","Yr","Beg/Mid",0,out$ages)]
  out$NAA = NAA[NAA[,3]=="B",c(1:2,4:ncol(NAA))]
  
  out$npar = fit$N_estimated_parameters
  
  # Estimates with uncertainty
  dq=fit$derived_quants[,1:3]
  out$FM_est = dq[grepl("F_",dq[,1]),]
  out$SSB_est = dq[grepl("SSB_",dq[,1]),]
  out$Rec_est = dq[grepl("Recr_",dq[,1]),]
  out$Dep_est = dq[grepl("Bratio_",dq[,1]),]
  
  
  # Reference points
  out$Refs=list()
  out$Refs$SSB_MSY = dq[dq[,1]=="SSB_MSY",]
  out$Refs$Brel = dq[dq[,1]=="B_MSY/SSB_unfished",]
  out$Refs$R0 = dq[dq[,1] =="Recr_unfished",]
  out$Refs$SSB0 = dq[dq[,1]=="SSB_unfished",]
  out$Refs$B0 = dq[dq[,1]=="Totbio_unfished",]
  out$Refs$OFL =  dq[grepl("OFLCatch",dq[,1]),]
  out$Refs$MSY = dq[dq[,1]=="Ret_Catch_MSY",]
  
  # Objective function
  out$nll = fit$likelihoods_used
  out$nll_ft = fit$likelihoods_by_fleet
  
  out$AIC  <- 2 * (out$nll$values[1] + out$npar) 
  out$CoVar = fit$CoVar
  # Retro
  #if(!is.null(Fit$peels))out$rho<-mohns_rho(Fit)
  out$RecDevs = fit$recruitpars
  #out$conv = fit$opt$convergence
  out$RunTime = fit$RunTime
  out$Nwarnings = fit$Nwarnings
  out$Warnings = fit$warnings
  out$maxgrad = fit$maximum_gradient_component
  out$sigmaR = fit$sigma_R_info
  
  out$par = fit$parameters
  out$par_h = fit$parameters_with_highest_gradients
  
  out$agedbase = fit$agedbase
  out$sizeselex = fit$sizeselex
  
  out$FleetNames = fit$FleetNames
  # object.size(out)/object.size(fit)
  out
  
}

extract_SS_list<-function(mdirs){
  
  out<-list()
  #mods<-list()
  for(i in 1:length(mdirs)){
    Fit<-readRDS(mdirs[[i]])
    out[[i]]<-extract_SS(Fit)
  }
  names(out)<-names(mdirs)
  
  out
  
}


get_ss3_outputs = function(dir){
  fit = SS_output(dir)
  extract_SS(fit)
}