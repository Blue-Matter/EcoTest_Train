# gets first two fleets selectivity and weighted selectivity of all other fleets (Catch in last year)
getsels_i4 = function(io,Find,dat,ploty){
  # selectivities
  
  selex = io$outputs$sizeselex
  sel = selex[selex$Yr==max(selex$Yr) & selex$Sex==1 & selex$Factor=="Lsel",]
  
  # you were here - need to locate fleet selectivities
  fleets = sel$Fleet
  nf = length(fleets)
  sind = 6:ncol(sel)
  slen =as.numeric(names(sel)[sind])
  sels = array(NA,c(4, length(sind)))
  for(ff in 1:2) sels[ff+1,] = as.numeric(sel[sel$Fleet==Find[ff],sind,drop=T])
  
  frest = (1:nf)[!(sel$Fleet%in%Find)]
  clast = dat$catch[dat$catch$year==max(dat$catch$year),]
  swt = rep(0,nf); swt[match(clast$fleet,fleets)]= clast$catch
  srest = apply(sel[frest,sind]*swt[frest],2,sum)
  sels[4,]=srest # weighted by catch in last year for now
  
  frest = (1:nf)
  swt = rep(0,nf); swt[match(clast$fleet,fleets)]= clast$catch
  srest = apply(sel[frest,sind]*swt[frest],2,sum)
  sels[1,]=srest # weighted by catch in last year for now
  
  
  sels=sels/apply(sels,1,max)
  if(ploty){matplot(slen,t(sels),type="l",lty=1,col=c("black","red","green","blue"),lwd=2); abline(h=0.05)}
  L5s = LFSs = rep(NA,4)
  for(ff in 1:4){
    ascend = 1:which.max(sels[ff,])
    L5s[ff] = approx(sels[ff,ascend],slen[ascend],0.05)$y
    LFSs[ff] = approx(sels[ff,ascend],slen[ascend],0.95)$y
  }
  VMLs = sels[,ncol(sels)]
  
  list(L5s = L5s,LFSs=LFSs, VMLs = VMLs)
  
}


#  dd= 1; io = ios[[dd]]; Fnam = Fnams[[dd]]; Inam = Inams[[dd]]; nsamp = 10; catch_CV = 0.1; L50_CV = 0.05; Linf_CV = 0.025; K_CV = 0.05; M_CV = 0.05; VML_CV = 0.05; peel=0; plotsmooth=T
SS_2_ET_raw = function(io, Fnam = c("F4_JPN_LL","F8_ESP_LL"), Inam = c("S2_JPN_LATE","S4_EU_ESP"),peel=0){
  
  # So this function just gets total (1) then fleet-specific time series data into list format
  # Later this can be converted to and from an excel file
  
  lastyr = io$inputs$dat$endyr - peel # Retro
  control = io$inputs$control
  Fnames = io$inputs$dat$fleetinfo$fleetname
  Find = match(Fnam, Fnames)
  Iind = match(Inam, Fnames)
  
  # === Biological info ======================================
  MGp = control$MG_parms
  Linf_mu = MGp[(1:nrow(MGp))[grepl("L_at_Amax",rownames(MGp))][1],3]
  #Linf = rep(100, nsamp) #rlnorm(nsamp, Linf_mu ,Linf_CV)
  K_mu = MGp[(1:nrow(MGp))[grepl("VonBert_K",rownames(MGp))][1],3]
  K =  rlnorm(nsamp, log(K_mu), K_CV)
  L50_mu = MGp[grepl("Mat50",rownames(MGp)),3]
  if(L50_mu==0 & "Age_Maturity" %in% names(control)){
    Mata = control$Age_Maturity
    Lata = (Linf_mu*(1-exp(-K_mu*(0:100))))[1:length(Mata)]
    L50_mu = approx(Mata,Lata,0.5)$y
  }
  
  #L50 = rlnorm(nsamp, log(L50_mu/Linf_mu*100), L50_CV) # all lengths are % Linf
  L50_Linf = L50_mu/Linf_mu
  
  if("natM"%in%names(control)){
    Ma = as.numeric(control$natM[1,,drop=T])
    Sa = exp(-cumsum(Ma))
    M_mu = sum(Sa*Ma)/sum(Sa) # survival weighted M
  }else{
    M_mu =  MGp[(1:nrow(MGp))[grepl("NatM",rownames(MGp))][1],3]
  }
  #M = rlnorm(nsamp,log(M_mu),M_CV)
  M_K = M_mu/K
  maxa = -log(0.05)/M_mu # age at 5% cumulative survival
  
  
  # === Fleet specific stuff ================================
  
  # get catches -----------
  dat = io$inputs$dat
  ct = dat$catch
  nc = nrow(ct)
  newind = rep(4,nc) # default to fleet 3
  for(ff in 1:2)newind[ct$fleet == Find[ff]]=ff+1
  aggcatch = aggregate(ct$catch,by=list(year = ct$year, fleet = newind),sum)
  totcatch = aggregate(ct$catch,by=list(year = ct$year, fleet = rep(1,nc)),sum)
  aggcatch = rbind(totcatch, aggcatch)
  aggcatch = aggcatch[aggcatch$year > 0 & aggcatch$year<=lastyr,] # Retro
  allfleet = aggregate(ct$catch,by=list(fleet = ct$fleet),sum)
  fleetcatch = aggregate(ct$catch,by=list(fleet = newind),sum)
  fleetcatch$x = fleetcatch$x/sum(fleetcatch$x)
  #print(allfleet)
  #print(fleetcatch)
  
  # get indices ----------
  ind = dat$CPUE
  ni = nrow(ind)
  newind = rep(4,ni)
  for(ff in 1:2)newind[ind$index == Iind[ff]]=ff+1
  aggind = aggregate(ind$obs, by=list(year = ind$year, index = newind),mean)
  aggindcv = aggregate(ind$se_log, by=list(year = ind$year, index = newind),mean)
  SSB = io$outputs$SSB
  yrs = unique(aggcatch$year)
  tind = data.frame(year = yrs,index = rep(1,length(yrs)),x = SSB[1:length(yrs)]); names(tind) = names(aggind)
  aggind = rbind(tind,aggind)
  aggind=aggind[aggind$year>0 & aggind$year<=lastyr,] # Retro
  aggindcv=aggindcv[aggindcv$year>0 & aggindcv$year<=lastyr,] # Retro
  
  lengthswitch = !is.null(dat$lencomp)
  if(lengthswitch){
    lbinlab = dat$lbin_vector
    lbins = dat$lbin_vector  # all lengths are scaled to a % of Linf
    lencomp = dat$lencomp[dat$lencomp$year>0,]
    
    labs = names(lencomp)
    bysex = paste0("f",lbinlab[1])%in%labs
    if(bysex){
      comps = lencomp[,match(paste0("f",lbinlab),labs)]+lencomp[,match(paste0("m",lbinlab),labs)]
    }else{
      comps = lencomp[,match(paste0("l",lbinlab),labs)]
    }
    inflate=1E4/max(apply(comps,1,sum))
    comps =ceiling(comps*inflate)
    
    # get mean length ----------
    mucv = getmucv(comps,lbins)
    nl = nrow(lencomp)
    newind = rep(4,nl)
    for(ff in 1:2)newind[lencomp$fleet == Find[ff]]=ff+1
    aggml = aggregate(mucv$mu, by=list(year = lencomp$year, fleet = newind),mean)
    aggmlcv = aggregate(mucv$cv, by=list(year = lencomp$year, fleet = newind),mean) # this isn't ideal as the aggregate fleet 3 should have a CV to reflect combination of many mean lengths
    totml = aggregate(mucv$mu, by=list(year = lencomp$year, fleet = rep(1,nrow(lencomp))),mean)
    totmlcv = aggregate(mucv$cv, by=list(year = lencomp$year, fleet = rep(1,nrow(lencomp))),mean) # this isn't ideal as the aggregate fleet 3 should have a CV to reflect combination of many mean lengths
    aggml = rbind(totml,aggml)
    aggmlcv = rbind(totmlcv, aggmlcv)
    
    aggml=aggml[aggml$year>0 & aggml$year<=lastyr,] # Retro
    aggmlcv=aggmlcv[aggmlcv$year>0 & aggmlcv$year<=lastyr,] # Retro
    
    # get cv in lengths
    CVlen = sapply(1:nl,function(X,comps,lbins){
      sd(rep(lbins,comps[X,]),na.rm=T)/mean(rep(lbins,comps[X,]),na.rm=T)
    },comps=comps,lbins=lbins)
    
    aggsdl = aggregate(CVlen, by=list(year = lencomp$year, fleet = newind),mean)
    totsdl = aggregate(CVlen, by=list(year = lencomp$year, fleet = rep(1,nrow(lencomp))),mean)
    aggsdl = rbind(totsdl,aggsdl)
    aggsdl=aggsdl[!is.na(aggsdl$x),]
    aggsdl=aggsdl[aggsdl$year>0 & aggsdl$year<=lastyr,] # Retro
    
    # get fraction mature
    Fmat = sapply(1:nl,function(X,comps,lbins,L50_mu){ # all these are phrased as % Linf
     
      indivs = rep(lbins,comps[X,]*1E3) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50_mu,na.rm=T)
    },comps=comps,lbins=lbins,L50_mu=L50_mu)
    
    aggFM = aggregate(Fmat, by=list(year = lencomp$year, fleet = newind),mean)
    totFM = aggregate(Fmat, by=list(year = lencomp$year, fleet = rep(1,nrow(lencomp))),mean)
    aggFM = rbind(totFM, aggFM)
    aggFM=aggFM[aggFM$year>0 & aggFM$year<=lastyr,] # Retro
  }
  
  selpar = getsels_i4(io,Find,dat,plotsmooth) # selectivity by fleet (MLE, not sampled)
  
  par = list(maxage = maxa, K = K_mu, L50_Linf = L50_Linf, selpar = selpar, Linf = Linf_mu, L50=L50_mu)
  ts = list(Catch = aggcatch, Meanlen = aggml, CVlen = aggsdl, Fracmat = aggFM, Index = aggind)
 
  yrs = io$input$dat$styr:io$input$dat$endyr
  SSBlab = paste0("SSB_",yrs)
  SSBest = io$outputs$SSB_est
  Brel = SSBest[match(SSBlab,rownames(SSBest)),2:3] / SSBest[match("SSB_MSY",rownames(SSBest)),2]
  
  list(par = par, ts = ts, Brel = Brel)
  
}


getlastml = function(ft,ts){
  sub = ts$Meanlen[ts$Meanlen$fleet==ft,]
  sub$x[sub$year==max(sub$year)]
}

fillxl = function(sum, xlfile,tofile,spec){
  
  pars_sheet = read_xlsx(xlfile,1)
  ts_sheet = read_xlsx(xlfile,2)
  
  par = sum$par
  ts = sum$ts
  
  # pars sheet
  L50 = par$L50
  Linf = par$Linf
  sels = sapply(par$selpar,function(x)x)/matrix(c(rep(L50,8),rep(1,4)),ncol=3)
  Cfrac = aggregate(ts$Catch$x,by=list(fleet=ts$Catch$fleet),sum)
  ML_Linf = sapply(1:4,getlastml,ts=ts)/Linf
  dats = round(c(par$maxa, par$K, par$L50_Linf, par$Linf, as.vector(t(sels)), Cfrac$x[2:4]/Cfrac$x[1] ),3)
  pars_sheet[,2] = c(spec,dats) 
  
  # time series sheet
  years = unique(ts$Catch$year)
  ny = length(years)
  tsmat = array(NA,c(ny,25))
  tsmat[,1] = years
  j = 1
  for(ff in 1:4){
    j=j+1
    scat = ts$Catch[ts$Catch$fleet == ff,]
    yind = match(scat$year,years)
    tsmat[yind,j] = scat$x
    j=j+1
    
    #ml
    sml = ts$Meanlen[ts$Meanlen$fleet == ff & ts$Meanlen$year>=min(years),]
    yind = match(sml$year,years)
    tsmat[yind,j] = sml$x
    j=j+2
    #ma
    
    #cvlen
    sml = ts$CVlen[ts$CVlen$fleet == ff & ts$CVlen$year>=min(years),]
    yind = match(sml$year,years)
    tsmat[yind,j] = sml$x
    j=j+1
    
    #fracmat
    sml = ts$Fracmat[ts$Fracmat$fleet == ff & ts$Fracmat$year>=min(years),]
    yind = match(sml$year,years)
    tsmat[yind,j] = sml$x
    j=j+1
    
    #index
    sml = ts$Index[ts$Index$index == ff& ts$Index$year>=min(years),]
    yind = match(sml$year,years)
    tsmat[yind,j] = sml$x/mean(sml$x)
    
  }
  tsmat = round(tsmat,3)
  tsmat=as.data.frame(tsmat)
  names(tsmat) = names(ts_sheet)
  writexl::write_xlsx(list(Parameters=pars_sheet, Time_series = tsmat),path=tofile)
}

