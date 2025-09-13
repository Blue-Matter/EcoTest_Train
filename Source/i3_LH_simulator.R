



add_F_dynamics3=function(cpars0, Effs, ss, L5_p = c(0.25, 0.1), LFS_p = c(0.1, 0.9), plot=F){
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


add_dep3 = function(cpars_ss, dep_rng = c(0.025,0.65), plot=F){
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


add_M_K_L50_Linf3 = function(cpars0, nsim, M_max = 0.6, M_min = 0.05, Linf_max = 1000, K_min = 0.075, plot=F, seed=1){
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



cat("Multi-stock, multi-fleet simulator code for Indicator 3loaded \n")

