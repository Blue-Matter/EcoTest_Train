
add_stochasticity = function(MOM, Mcv = 0.15, Kcv = 0.15, Linfcv = 0.025, MKcor = 0.8, Dcv = 0.2, ploty=F){
  
  ns = length(MOM@Stocks)
  nsim = MOM@nsim
  MKcov = MKcor * sqrt(Mcv^2 * Kcv^2)
  
  #cor = cov/sqrt(varx, *vary)
  #cov = cor * sqrt(varx, *vary)
  
  bounds = log(c(0.66,1.5))
  
  for(ss in 1:ns){
    
    MKerr = rmvnorm(nsim*10,mean = c(0,0),sigma = matrix(c(Mcv^2,MKcov,MKcov,Kcv^2),nrow=2)); if(ploty){plot(MKerr[,1],MKerr[,2])}
    MKerr = MKerr[MKerr[,1]>bounds[1] & MKerr[,1]<bounds[2]  & MKerr[,2]> bounds[1] & MKerr[,2]< bounds[2],][1:nsim,]; if(ploty){ points(MKerr[,1],MKerr[,2],pch = 19, col="#0000ff60"); abline(h=bounds,v=bounds,col="#0000ff60")}
    MOM@cpars[[ss]][[1]]$M_ageArray = MOM@cpars[[ss]][[1]]$M_ageArray * exp(MKerr[,1]); if(ploty)hist(MOM@cpars[[ss]][[1]]$M_ageArray[,1,1])
    
    Ks =  MOM@Stocks[[ss]]@K[1] * exp(MKerr[,2])
    Linferr = rnorm(nsim*10,0,Linfcv)
    Linf = MOM@Stocks[[ss]]@Linf[1] * exp(Linferr[Linferr>bounds[1] & Linferr < bounds[2]][1:100])
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
    
    Wt_age = a*Len_age^b
    if(ploty){matplot(t(Wt_age[1:10,,1]),type="l"); lines(MOM@cpars[[ss]][[1]]$Wt_age[1,,1],lwd=2)}
    
    MOM@cpars[[ss]][[1]]$Wt_age = Wt_age
    
    MOM@cpars[[ss]][[1]]$qs = NULL
    MOM@cpars[[ss]][[1]]$D =  MOM@Stocks[[ss]]@D[1]*rlnorm(nsim,0,Dcv)
    
  }
  
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
  
  for(ss in 1:ns){
    for(ff in 1:nf){
      Va = MOM@cpars[[ss]][[ff]]$V
      areall1s = function(x)all(x==1)
      the1s = apply(Va[1,,],2,areall1s)
      indfill = (1:dim(Va)[3])[the1s]
      indfrom = max(indfill)+1
      Va[,,indfill] = Va[,,indfrom]
      MOM@cpars[[ss]][[ff]]$V = Va
    }
  }
  
  MOM
  
}

