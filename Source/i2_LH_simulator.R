



add_F_dynamics=function(cpars0, Effs, ss, L5_p = c(0.5, 0.1), LFS_p = c(0.1, 0.7), plot=F){
  # cpars0 = cpars[[ss]][[ff]]; L5_p = c(0.5, 0.1); LFS_p = c(0.1, 0.7)
  nsim = length(cpars0$L50)
  L50 = cpars0$L50
  Linf = cpars0$Linf
  LB = L5_p[1] * L50
  UB = L50 + L5_p[2] * (Linf - L50)
  L5 = runif(nsim,LB,UB)
  # cbind(LB,L5,UB)
  LB_LFS = (L5*1.05) + LFS_p[1] * (Linf - (L5*1.05))
  UB_LFS = (L5*1.05) + LFS_p[2] * (Linf - (L5*1.05))
  LFS = runif(nsim, LB_LFS, UB_LFS)
  
  if(plot){
    tab = as.data.frame(cbind(round(LB,2),round(L5,2), round(UB,2),rep("-",nsim), round(L50,2),rep("-",nsim), round(LB_LFS,2), round(LFS,2), round(UB_LFS,2), rep("-",nsim), Linf))
    names(tab) = c("LB_L5","L5","UB_L5","bla","L50","bla2","LB_LFS","LFS","UB_LFS","bla3","Linf")
    print(tab)
  }
  
  Vmaxlen = runif(nsim,0,1)
  Find = t(sapply(1:nsim,function(x,Effs,ss)Effs[[x]][,ss],Effs=Effs, ss=ss))
  
  cpars0$L5 = L5
  cpars0$LFS = LFS
  cpars0$Vmaxlen = Vmaxlen
  cpars0$Find = Find
  cpars0
}

EcoTest_sample = function(vec){
  nv = length(vec)
  pnts = seq(0,1,length.out=nv+1)
  LB = pnts[1:nv]
  UB = pnts[2:(nv+1)]
  rnd= runif(1)
  print(rnd)
  pos = (1:nv)[rnd >=LB & rnd <UB]
  vec[pos]
}


getrelF = function(cpars_ss, ny=15){
  
  Fsum = cpars_ss[[1]]$Find
  nyears = ncol(Fsum)
  for(ff in 2:length(cpars_ss))Fsum = Fsum + cpars_ss[[1]]$Find
  Frec = apply(Fsum[,nyears-(0:(ny-1))],1,mean)
  Fmu = apply(Fsum,1,mean)
  Frec / Fmu
  
}

add_dep = function(cpars_ss, dep_rng = c(0.025,0.5), plot=F){
  nsim = length(cpars_ss[[1]]$L50)
  Frel = getrelF(cpars_ss)
  Dsamp = runif(nsim,dep_rng[1],dep_rng[2])
  Dord = Dsamp[order(Dsamp)]
  Dnew = rep(NA,nsim)
  ord = order(Frel,decreasing=T)
  Dnew[ord] = Dord
  if(plot) {cbind(Frel,ord, Dnew); plot(Frel, Dnew,pch=19)}
  Dnew
}


get_seed = function()as.numeric(format(Sys.time(), "%OS3"))*1000

add_M_K_L50_Linf = function(cpars0, nsim, M_max = 0.6, M_min = 0.1, Linf_max = 1000, K_min = 0.15, plot=F, seed=1){
  OM= new('OM')
  OM@nsim = nsim * 10 # make too many so you can filter out very high Linf and M samples
  OM@Species = ""
  OM@seed = get_seed() # LH2OM uses the OM random seed. if this is kept the same you get nstocks identical stocks!
  Class = "predictive"; Order = "predictive"; genus = "predictive"; species = "predictive"
  set.seed(get_seed()) # LH2OM sets the seed, so if you want sampling of family you need to reset it here to something new and changeable
  Family = sample(c("Scombridae","Istiophoridae","Carcharhinidae"),1); print(Family)
  cpars = MSEtool::LH2OM(OM, Class = "predictive", Order = "predictive", Family = Family, plot=plot)@cpars
  
  keep = ((1:OM@nsim)[cpars$M < M_max &
                      cpars$M > M_min &
                      cpars$Linf < Linf_max &
                      cpars$L50/cpars$Linf < 0.9 &
                      cpars$K > K_min])[1:nsim]
  
  cpars2 = cpars; cpars2$M = cpars$M[keep]; cpars2$K = cpars$K[keep]; cpars2$L50 = cpars$L50[keep]; cpars2$Linf = cpars$Linf[keep]; 
  cpars2$L50 = cpars2$L50 /cpars2$Linf *100; cpars2$Linf[] = 100
  if(plot) {par(mfrow = c(1,3)); plot(cpars2$M, cpars2$K); plot(cpars2$M, cpars2$L50); plot(cpars2$K, cpars2$L50)}# plot(cpars2$L50, cpars2$Linf); plot(cpars2$M,cpars2$Linf); plot(cpars2$K, cpars2$Linf)}
  cpars0$Linf = cpars2$Linf
  cpars0$L50 = cpars2$L50
  cpars0$K = cpars2$K
  cpars0$M = cpars2$M
  cpars0
}


dummy_effort = function(Fleet){
  tiny = 1E-6
  Fleet@EffYears = c(1,Fleet@nyears); Fleet@EffLower = c(1,1); Fleet@EffUpper = c(1+tiny, 1+tiny);  Fleet@Esd = c(0.1, 0.1)
  Fleet
}

dummy_selectivity = function(Fleet){
  Fleet@L5 = c(10,30);  Fleet@LFS = c(40, 60);  Fleet@Vmaxlen = c(0,1)
  Fleet@LR5 = c(0,0);  Fleet@LFR = c(0,0);  Fleet@Rmaxlen = c(1,1)
  Fleet
}


make_Imp = function(i){MSEtool::Perfect_Imp}

make_Imps = function(nstocks, nfleets){
  Imps = list()
  for(ss in 1:nstocks)Imps[[ss]]=lapply(1:nfleets,make_Imp)
  Imps
}

make_Ob = function(i){
  Oby = MSEtool::Generic_Obs
  Oby@CAL_nsamp = c(200,300)
  Oby@CAL_ESS = c(100,200)
  Oby
}

make_Obs = function(nstocks, nfleets){
  Obs = list()
  for(ss in 1:nstocks)Obs[[ss]]=lapply(1:nfleets,make_Ob)
  Obs
}

make_Fleets = function(nstocks, nfleets, nyears, CurrentYr){
  Fleets= list()
  for(ss in 1:nstocks)Fleets[[ss]]=lapply(1:nfleets,make_Fleet,nyears=nyears, CurrentYr = CurrentYr)
  Fleets
}

make_Fleet = function(i, nyears, CurrentYr){
  tiny = 1E-6
  Fleet = new('Fleet')
  Fleet@Name = "Generic Simulated Fleet"
  Fleet@nyears = nyears
  Fleet@CurrentYr = CurrentYr
  Fleet = dummy_effort(Fleet)
  Fleet@qinc = c(0,tiny)
  Fleet@qcv = c(0,tiny)
  Fleet = dummy_selectivity(Fleet)
  Fleet@isRel = FALSE
  Fleet@DR = c(0.1, 0.9)
  Fleet@LR5 = c(0,0)
  Fleet@LFR = c(0,0)
  Fleet@Rmaxlen = c(1, 1)
  Fleet@Spat_targ = c(1,1)
  Fleet
}

make_Stock = function(i, maxage){
  tiny = 1E-6
  Stock = new('Stock')
  Stock@Name = "Generic Simulated Fish"
  Stock@maxage = maxage # assumed to be asymptotic dynamics after 
  Stock@R0 = 1E6
  Stock@M = c(0.2, 0.2)
  Stock@Msd = c(0,tiny)
  Stock@h = c(0.6,0.95)
  Stock@SRrel = 1
  Stock@Perr = c(0.35, 0.7)
  Stock@AC = c(0, 0.6)
  Stock@Linf = c(100, 100)
  Stock@Linfsd = c(tiny, tiny)
  Stock@K = c(0.2, 0.2)
  Stock@Ksd = c(0, tiny)
  Stock@t0 = c(0, 0)
  Stock@LenCV = c(0.1, 0.15)
  Stock@L50 = c(60,70)
  Stock@L50_95 = c(5, 10)
  Stock@D = c(0.025, 0.6)
  Stock@a = 1E-5
  Stock@b = 3
  Stock@Size_area_1 = c(0.5, 0.5)
  Stock@Frac_area_1 = c(0.5, 0.5)
  Stock@Prob_staying = c(0.5, 0.5)
  Stock@Fdisc = c(0.1, 0.9)
  Stock@Source = "EcoTest: simulated generic fish"
  Stock
}


cat("Multi-stock, multi-fleet simulator code for Indicator 2 loaded \n")

