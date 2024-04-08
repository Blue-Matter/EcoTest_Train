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

runbatch = function(x, MOM,  MPs, largedir, doPE=T){ # x is the batch number of 100 simulations
  
  temp = MOM
  temp@seed = x
  Effmat <<-totEffmat[(x-1)*100+(1:100),]  
  
  if(doPE) temp = overwritePE(temp)                           # sapply(temp@cpars,function(x)x[[1]]$Perr_y[1,])
  
  histfile = paste0(largedir,"/Hist_",x,".rda")
  
  if(!(file.exists(histfile))){
    Hist = SimulateMOM(temp, parallel = FALSE)
    saveRDS(Hist,histfile)
  }else{
    Hist = readRDS(histfile)
  }
                   # saveRDS(Hist,"C:/temp/Ecotest/dump/Hist.rda")
  MMSE = ProjectMOM(Hist, MPs = MPs[[1]], checkMPs = FALSE)  # saveRDS(MMSE2,"C:/temp/Ecotest/dump/MMSE2.rda")
  saveRDS(MMSE, paste0(largedir,"/MMSE_",x,".rda"))
  
  
}


# [100,] 0.7573089

