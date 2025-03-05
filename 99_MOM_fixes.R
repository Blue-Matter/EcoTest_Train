



stoch_SLarray= function(Linf, L5, LFS, Vmaxlen, lens){
  
  srs <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5)
  sls <- (LFS - L5) /((-log(0.05,2))^0.5)
  t(sapply(1:length(Linf), MSEtool::getsel, lens, LFS, sls=sls, srs=srs))
  
}  


# MOM=MOM3; Mcv = 0.15; Kcv = 0.15; Linfcv = 0.025; MKcor = 0.8; qcv = 0.1; hcv=0.1; L50cv = 0.05; L5cv = 0.1; LFScv = 0.1; Vmaxlencv = 0.1; ploty = T


add_stochasticity = function(MOM, Mcv = 0.15, Kcv = 0.15, Linfcv = 0.025, MKcor = 0.8, qcv = 0.1, hcv=0.1, 
                             L50cv = 0.05, L5cv = 0.1, LFScv = 0.1, Vmaxlencv = 0.1, ploty=F){
  
  ns = length(MOM@Stocks)
  nf = length(MOM@Fleets[[1]])
  nsim = MOM@nsim
  MKcov = MKcor * sqrt(Mcv^2 * Kcv^2)
  allyears = dim(MOM@cpars[[1]][[1]]$M_ageArray)[3]
  nsim = dim(MOM@cpars[[1]][[1]]$M_ageArray)[1]
  na = dim(MOM@cpars[[1]][[1]]$M_ageArray)[2]
  calbins = MOM@cpars[[1]][[1]]$CAL_binsmid
  nl = length(calbins)
  
  #cor = cov/sqrt(varx, *vary)
  #cov = cor * sqrt(varx, *vary)
  
  bounds = log(c(0.66,1.5))
  
  for(ss in 1:ns){
    
    MKerr = rmvnorm(nsim*10,mean = c(0,0),sigma = matrix(c(Mcv^2,MKcov,MKcov,Kcv^2),nrow=2)); if(ploty){plot(MKerr[,1],MKerr[,2])}
    MKerr = MKerr[MKerr[,1]>bounds[1] & MKerr[,1]<bounds[2]  & MKerr[,2]> bounds[1] & MKerr[,2]< bounds[2],][1:nsim,]; if(ploty){ points(MKerr[,1],MKerr[,2],pch = 19, col="#0000ff60"); abline(h=bounds,v=bounds,col="#0000ff60")}
    MOM@cpars[[ss]][[1]]$M_ageArray = MOM@cpars[[ss]][[1]]$M_ageArray * exp(MKerr[,1]); 
    if(ploty)hist(MOM@cpars[[ss]][[1]]$M_ageArray[,1,1])
    
    Ks =  MOM@Stocks[[ss]]@K[1] * exp(MKerr[,2])
    Linferr = rnorm(nsim*10,0,Linfcv)
    Linftrial = MOM@Stocks[[ss]]@Linf[1] * exp(Linferr[Linferr>bounds[1] & Linferr < bounds[2]])
    Linf = Linftrial[Linftrial < (max(MOM@cpars[[ss]][[1]]$CAL_bins)-1)][1:nsim]
    t0 = MOM@Stocks[[ss]]@t0[1]
    a = MOM@Stocks[[ss]]@a
    b = MOM@Stocks[[ss]]@b
    
    dims = dim(MOM@cpars[[ss]][[1]]$Len_age)
    na = dims[2]; ny = dims[3]
    ind = as.matrix(expand.grid(1:nsim,1:na, 1:ny))
    Len_age = array(NA, dims)
    Len_age[ind] = Linf[ind[,1]]*(1-exp(-Ks[ind[,1]]*((ind[,2]-0.5)-t0)))
    if(ploty){matplot(t(Len_age[1:10,,1]),type="l"); lines(MOM@cpars[[ss]][[1]]$Len_age[1,,1],lwd=2)}
    MOM@cpars[[ss]][[1]]$Len_age = Len_age
    MOM@cpars[[ss]][[1]]$K = Ks
    MOM@cpars[[ss]][[1]]$Linf = Linf
    
    Wt_age = a*Len_age^b
    if(ploty){matplot(t(Wt_age[1:10,,1]),type="l"); lines(MOM@cpars[[ss]][[1]]$Wt_age[1,,1],lwd=2)}
    
    MOM@cpars[[ss]][[1]]$Wt_age = Wt_age
    MOM@cpars[[ss]][[1]]$qs = rlnorm(nsim,0,qcv)
    #MOM@cpars[[ss]][[1]]$D =  MOM@Stocks[[ss]]@D[1]*rlnorm(nsim,0,Dcv)
    
    # Steepness 
    hmu = MOM@Stocks[[ss]]@h[1]
    betatest = rbeta(nsim*100, MSEtool::alphaconv(hmu,hcv), MSEtool::betaconv(hmu,hcv))
    MOM@cpars[[ss]][[1]]$h = betatest[betatest > 0.22 & betatest < 0.99][1:nsim]
    
    # Selectivity uncertainty 
    simno = 1 # previously everything was deterministic so we take the first simulation
    La = MOM@cpars[[ss]][[1]]$Len_age
    
    for(ff in 1:nf){
      
      Va = MOM@cpars[[ss]][[ff]]$V
      SLarray = array(NA,c(nsim,nl,allyears))
      tofill = apply(Va[1,,],2,function(x)all(x==1))
      if(any(tofill)){
        first = max((1:dim(Va)[2])[tofill])+1
        tofill = (1:dim(Va)[2])[tofill]
        for(jj in tofill)Va[,,jj] = Va[,,first]
      }
      L5str = -1
      
      for(y in 1:allyears){
        
        amax = which.max(Va[simno,,y])
        L5mu = approx(Va[simno,1:amax,y],La[simno,1:amax,y],0.05,ties = "ordered")$y #; plot(Va[simno,,y],La[simno,,y], type="l");abline(v=0.05,h=L5mu,col="red")
        LFSmu = La[simno,,y][match(max(Va[simno,,y]),Va[simno,,y])]
        if(is.na(L5mu))L5mu = 0.5 * LFSmu
        Vmmu = Va[simno,na,y]
        Vmmu[Vmmu >0.98] = 0.98
        Vmmu[Vmmu < 0.02] = 0.02
        
        if(L5mu != L5str){
          
          L5 = rlnorm(nsim, log(L5mu), L5cv)
          LFS = rlnorm(nsim, log(LFSmu), LFScv)
          L5[L5 > (0.95 * LFS)] = LFS[L5 > (0.95 * LFS)] * 0.95 # can't be larger than LFS
          LFS[LFS > (0.95*Linf)] = Linf[LFS > (0.95*Linf)]*0.95 # can't be larger than Linf
          Vmstoch = rbeta(nsim, MSEtool::alphaconv(Vmmu,Vmaxlencv), MSEtool::betaconv(Vmmu,Vmaxlencv))
          SLa = stoch_SLarray(Linf, L5, LFS, Vmstoch, lens = calbins)
          L5str = L5mu
          #print(paste(y,L5mu))
        }
        
        SLarray[,,y] = SLa
        
      }  
     
      MOM@cpars[[ss]][[ff]]$SLarray = SLarray
      
    }
  
  }
  
  #MOM = add_SL_array(MOM) # lapply(MOM2@cpars,function(x)x[[1]]$SLarray[1,1:25,1:5])
  MOM
}



Len_from_age = function(Vav,Lav,CAL_mu,lenCV = 0.15){
  nl = length(CAL_mu)
  na = length(Vav)
  iALK<-array(NA,c(na,nl))
  ind=as.matrix(expand.grid(1:na,1:nl))
  iALK[ind]<-dlnorm(CAL_mu[ind[,2]],log(Lav[ind[,1]]),lenCV) # SD is determined by linear model of Allioud et al. 2017
  #contour(x=1:na,y=CAL_mu,iALK,nlevels=10)
  #plot()
  SelL = apply(Vav*iALK,2,sum) #; plot(SelL)
  SelL/max(SelL)
}

makeSLarray = function(Va, La, CAL_mu){
  nCALb = length(CAL_mu)
  nsim=dim(La)[1]
  na = dim(La)[2]
  allyears = dim(La)[3]
  SLarray = array(NA,c(nsim,nCALb,allyears))
  for(i in 1:nsim){
    for(y in 1:allyears){
      SLarray[i,,y]= Len_from_age(Va[i,,y],La[i,,y],CAL_mu)
    }
  }
  SLarray
}

makeSLarray_LI = function(Va, La, CAL_mu){
  nCALb = length(CAL_mu)
  nsim=dim(La)[1]
  na = dim(La)[2]
  allyears = dim(La)[3]
  SLarray = array(NA,c(nsim,nCALb,allyears))
  for(i in 1:nsim){
    for(y in 1:allyears){
      Stemp = approx(La[i,,y],Va[i,,y],CAL_mu)$y
      withval = Stemp[!(is.na(Stemp))]
      Stemp[is.na(Stemp)] = withval[length(withval)]
      SLarray[i,,y]= Stemp
    }
  }
  SLarray
}

add_SL_array = function(MOM){
  
  ns = length(MOM@Stocks)
  nf = length(MOM@Fleets[[1]])
  
  for(ss in 1:ns){
    CAL_mu = MOM@cpars[[ss]][[1]]$CAL_binsmid
    La = MOM@cpars[[ss]][[1]]$Len_age
    for(ff in 1:nf){
      Va = MOM@cpars[[ss]][[ff]]$V
      SLarray = makeSLarray_LI(Va,La, CAL_mu)
      MOM@cpars[[ss]][[ff]]$SLarray = SLarray
      cat(".")
      # plot(La[1,,10],Va[1,,10],type="l"); lines(CAL_mu,SLarray[1,,10],col="red")
    }
  }
  cat("\n")
  MOM
}


fix_selectivity_1 = function(MOM){
  
  ns = length(MOM@Stocks)
  nf = length(MOM@Fleets[[1]])
  nsim = MOM@nsim
  nyears = MOM@Fleets[[1]][[1]]@nyears
  proyears = MOM@proyears
 
  CAL_bins_all = lapply(MOM@cpars,function(x)x[[1]]$CAL_bins)
  stock_with_longest = which.max(sapply(CAL_bins_all,max))
  CAL_bins = unlist(CAL_bins_all[stock_with_longest])
  CAL_mids_all = lapply(MOM@cpars,function(x)x[[1]]$CAL_binsmid)
  CAL_binsmid = unlist(CAL_mids_all[stock_with_longest])
  newnl = length(CAL_binsmid)
  
  for(ss in 1:ns){
    for(ff in 1:nf){
      Va = MOM@cpars[[ss]][[ff]]$V
      areall1s = function(x)all(x==1)
      the1s = apply(Va[1,,],2,areall1s)
      if(any(the1s)){
        indfill = (1:dim(Va)[3])[the1s]
        indfrom = max(indfill)+1
        Va[,,indfill] = Va[,,indfrom]
      }
      MOM@cpars[[ss]][[ff]]$V = Va
      MOM@cpars[[ss]][[ff]]$CAL_bins = CAL_bins
      MOM@cpars[[ss]][[ff]]$CAL_binsmid = CAL_binsmid
      #MOM@cpars[[ss]][[ff]]$retL = NULL
      #MOM@cpars[[ss]][[ff]]$
     
      if("retL" %in% names(MOM@cpars[[ss]][[ff]])){
        
        Fdisc_array2 = retL = array(NA, c(nsim,newnl,nyears+proyears))  
        refmids = c(0,unlist(CAL_mids_all[ss]),1E10)
        
        for(sim in 1:nsim){
          for(y in 1:(proyears+nyears)){
            old_retL = MOM@cpars[[ss]][[ff]]$retL[sim,,y]
            old_retL = c(0,old_retL,old_retL[length(old_retL)]) 
            retL[sim,,y] = approx(refmids,old_retL,CAL_binsmid)$y
            
            old_Fd2 = MOM@cpars[[ss]][[ff]]$Fdisc_array2[sim,,y]
            old_Fd2 = c(0,old_Fd2,old_Fd2[length(old_Fd2)]) 
            Fdisc_array2[sim,,y] = approx(refmids,old_Fd2,CAL_binsmid)$y
            
          }
        }
        
        MOM@cpars[[ss]][[ff]]$Fdisc_array2 = Fdisc_array2
        MOM@cpars[[ss]][[ff]]$retL = retL
        
      }
          
    }
  }
  
  MOM
  
}

