



make_Catch_Frac3 = function(MOM, Frange = c(0.05, 1)){
 
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
