
get_CurB = function(MOM){
  
  nstocks = length(MOM@Stocks)
  R0 = MOM@Stocks[[1]]@R0
  a = MOM@Stocks[[1]]@a
  b = MOM@Stocks[[1]]@b
  nsim = length(MOM@cpars[[1]][[1]]$L50)
  maxage = MOM@Stocks[[1]]@maxage
  age_mat = matrix(rep(0:maxage,each=nsim),nrow=nsim)
  CurB = array(NA,c(nstocks,nsim))
  
  for(ss in 1:nstocks){
    
    cpars = MOM@cpars[[ss]][[1]]
    len_age = cpars$Linf * (1-exp(-cpars$K * age_mat))
    wt_age = a*len_age^b
    surv = exp(-age_mat*cpars$M)
    BpR = apply(surv*wt_age,1,sum) 
    CurB[ss,] = R0 * BpR * cpars$D
    
  }
  
  CurB
  
}

get_CurF = function(MOM){
  nstocks = length(MOM@Stocks)
  nfleets = length(MOM@Fleets[[1]])
  nsim = length(MOM@cpars[[1]][[1]]$L50)
  nyears = ncol(MOM@cpars[[1]][[1]]$Find)
  CurF = array(NA, c(nstocks, nfleets, nsim))
  for(ss in 1:nstocks){
    for(ff in 1:nfleets){
      CurF[ss,ff,] = MOM@cpars[[ss]][[ff]]$Find[,nyears]
   }
  }
  CurF
}

calc_Catch_Frac = function(nstocks, nfleets, nsim, CurF, Fmult){ #}, CurB, Ffrac){
  
  Catch_Frac =  list()
  for(ss in 1:nstocks){
    CF = array(NA, c(nsim, nfleets))
    for(ff in 1:nfleets){
      CF[,ff] =  CurF[ss,ff,] * Fmult[ff,] # * CurB[ss, ]Ffrac[ss,] *
    }
    CF = CF / apply(CF,1,sum)
    Catch_Frac[[ss]] = CF
  }
  Catch_Frac
  
}


make_Catch_Frac = function(MOM, Frange = c(0.2, 1)){
 
  nstocks = length(MOM@Stocks)
  nfleets = length(MOM@Fleets[[1]])
  nsim = length(MOM@cpars[[1]][[1]]$L50)
  make_Fmult = function(i, Frange, nfleets)c(1, runif(nfleets-1,Frange[1],Frange[2]))
  Fmult = sapply(1:nsim, make_Fmult, Frange = Frange, nfleets = nfleets)
  #CurB = get_CurB(MOM)
  CurF = get_CurF(MOM)
  calc_Catch_Frac(nstocks, nfleets, nsim, CurF, Fmult)
  
}


cat("Multi-stock, multi-fleet catch fraction calculator for Indicator 2 loaded \n")
