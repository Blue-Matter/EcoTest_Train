# Batching functions for generating large simulated datasets


overwritePE = function(temp){ # recasts projected process errors given historical rec dev SD and lag 1 AC

  set.seed(temp@seed)
  proyears = temp@proyears
  nstock = length(temp@cpars)
  nfleet = length(temp@cpars[[1]])
  nsim = temp@nsim
  
  for(ss in 1:nstock){
    
    Perr = temp@cpars[[ss]][[1]]$Perr_y
    nPE = ncol(Perr)
    edind = nPE - ((proyears-1):0)
    preind = (1:nPE)[!((1:nPE)%in%edind)]
    AC = apply(log(Perr[,preind]),1,function(x)acf(x,1,plot=F)[[1]][2,1,1])
    SD = apply(log(Perr[,preind]),1,sd)
    
    for(y in edind){ # overwrite projection process errors
      Perr[,y] <- exp(AC*log(Perr[,y-1]) + ((1-AC^2)^0.5)* rnorm(n=nsim, mean=0, sd=SD) - (1-AC)* (SD^2)/2)
    }
    
    for(ff in 1:nfleet)  temp@cpars[[ss]][[ff]]$Perr_y = Perr
    
  }
  
  temp
  
}


make_sim_dataset = function(MMSE){
  nstock = length(MMSE@PPD)
  outs = list()
  for(ss in 1:nstock)outs[[ss]] = proc_dat(MMSE,sno=ss)
  
}

# x = 1; MOM = readRDS("./MOM/MOM_stitch_100sim.rds"); MPs = "Frand_MMP"; doPE = T; largedir = "C:/temp/Ecotest/batching/Independent_F"

runbatch = function(x, MOM,  MPs, largedir, doPE=T, dostoch = T){ # x is the batch number of 100 simulations
  
  temp = MOM
  temp@seed = x
  set.seed(x)
  Effmat <<-totEffmat[(x-1)*100+(1:100),]  
 
  if(doPE) temp = overwritePE(temp)                           # sapply(temp@cpars,function(x)x[[1]]$Perr_y[1,])
  if(dostoch) temp = add_stochasticity(temp)                  # adds stochasticity in M, K, Linf, and stock depletion
  
  #histfile = paste0(largedir,"/Hist_",x,".rda")
  
  #if(!(file.exists(histfile))){
     Hist = SimulateMOM(temp, parallel = FALSE)
  #  saveRDS(Hist,histfile)
  #}else{
  #  Hist = readRDS(histfile)
  #}
                   # saveRDS(Hist,"C:/temp/Ecotest/dump/Hist.rda")
  MMSE = ProjectMOM(Hist, MPs = MPs[[1]], checkMPs = FALSE)  # saveRDS(MMSE2,"C:/temp/Ecotest/dump/MMSE2.rda")
  MMSE = trim_MMSE(MMSE)
  saveRDS(MMSE, paste0(largedir,"/MMSE_",x,".rda"))
  
  
}

gettodosims=function(largdir, maxsim=500){
  (1:maxsim)[!((1:maxsim) %in% sapply(list.files(largedir),function(x)as.numeric(strsplit(strsplit(x,"_")[[1]][2],".rda")[[1]][1])))]
}


trim_MMSE = function(MMSE){

  ns = length(MMSE@multiHist)
  nf = length(MMSE@multiHist[[1]])
  for(ss in 1:ns){
    for(ff in 1:nf){
      
      # Samppars Stock
      sloty = c("SSB","N","Biomass","mov", "surv",  "Len_age",  "M_ageArray", "Wt_age", "Mat_age", "LatASD","Fec_Age")
      for(snam in sloty)MMSE@multiHist[[ss]][[ff]]@SampPars$Stock[[snam]] = NULL
      
      # Samppars Fleet
      sloty =c("CBret", "CB", "retL_real", "SLarray_real", "retL", "Fdisc_array2", "SLarray", "retA_real" ,"V_real_2", "retA_real_2", "retA", "V_real", "Wt_age_C" , "Fdisc_array1","V") 
      for(snam in sloty)MMSE@multiHist[[ss]][[ff]]@SampPars$Fleet[[snam]] = NULL
      
      # Samppars Obs
      MMSE@multiHist[[ss]][[ff]]@SampPars$Obs$Sample_Area = NULL
      
      # Data Misc
      MMSE@multiHist[[ss]][[ff]]@Data@Misc$StockPars=NULL
      MMSE@multiHist[[ss]][[ff]]@Data@Misc$FleetPars=NULL
      
      # TSdat 
      sloty = c("Number", "Biomass", "VBiomass", "Removals", "Landings", "Discards", "Find", "RecDev", "SPR", "Unfished_Equilibrium")
      for(snam in sloty)MMSE@multiHist[[ss]][[ff]]@TSdata[[snam]] = NULL
      
      # AtAge
      sloty = c("Z.Mortality", "F.Mortality", "Fret.Mortality", "Number", "Biomass", "VBiomass", "SBiomass", "Removals", "Landings", "Discards", "Select", "Retention")  
      for(snam in sloty)MMSE@multiHist[[ss]][[ff]]@AtAge[[snam]] = NULL
      
    }
  }
  
  MMSE
  
}

