

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


fix_maturity=function(MOM2){
  
  cat("Blue marlin 50% mature at 256cm http://www.iccat.int/Documents/SCRS/Manual/CH2/2_1_6_BUM_ENG.pdf \n")
  
  
  
  
}