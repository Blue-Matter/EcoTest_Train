make_MOM = function(Name = "Simulated for AI training", nsim = 100, 
                    proyears = 50, interval = 1, seed = 1, 
                    nstocks = 3, nfleets = 3){
  
  # Name = "Simulated for AI training"; nsim = 100; proyears = 50; interval = 1; seed = 1; nstocks = 3; nfleets = 3
  
  MOM = new('MOM')
  MOM@Name = Name
  MOM@Agency = "Blue Matter"
  MOM@Region = "Global"
  MOM@Sponsor = "ICCAT"
  MOM@Latitude = 0
  MOM@Longitude = 0
  MOM@nsim = nsim
  MOM@proyears = proyears
  MOM@interval = interval
  MOM@pstar = 0.5
  MOM@maxF = 5
  MOM@reps = 1
  MOM@cpars = list()
  MOM@seed = seed
  MOM@Source = "EcoTest Homepage (coming soon)"
  
  MOM@Stocks = lapply(1:nstocks,make_Stock,nsim=nsim)
  MOM@Fleets = make_Fleets(nfleets)
  
  
}

make_Stock = function(i,nsim){
  tiny = 1E-6
  Stock = new('Stock')
  Stock@Name = "Generic Simulated Fish"
  Stock@maxage = 30 # assumed to be asymptotic dynamics after 30
  Stock@R0 = 1000
  Stock@M = c(0.2, 0.2)
  Stock@Msd = c(0,tiny)
  Stock@h = 0.75
  Stock@SRrel = 1
  Stock@Perr = c(0.3, 0.7)
  Stock@AC = c(0.3 ,0.6)
  Stock@Linf = c(1,1)
  Stock@Linfsd = c(0.01,0.025)
  Stock@K = c(0.2, 0.2)
  Stock@Ksd = c(0, tiny)
  Stock@t0 = c(-0.5,-0.5)
  Stock@LenCV = c(0.1, 0.15)
  Stock@L50 = c(0.6,0.7)
  Stock@L50_95 = c(0.05, 0.1)
  Stock@D = c(0.05, 0.5)
  Stock@a = 1
  Stock@b = 3
  Stock@Size_area_1 = c(0.5, 0.5)
  Stock@Frac_area_1 = c(0.5, 0.5)
  Stock@Prob_staying = c(0.5, 0.5)
  Stock@Fdisc = c(0.1, 0.9)
  Stock@Source = "EcoTest: simulated generic fish"
  Stock
}

make_cpars = function(nsim, M_max = 1, Linf_max = 1000, plot=F){
  OM= new('OM')
  OM@nsim = nsim * 2 # make too many so you can filter out very high Linf and M samples
  OM@Species = ""
  cpars = LH2OM(OM, Class = "Actinopterygii", Order = "Scombriformes", Family = "Scombridae")@cpars
  keep = ((1:OM@nsim)[cpars$M < M_max & cpars$Linf < Linf_max])[1:nsim]
  cpars2 = cpars; cpars2$M = cpars$M[keep]; cpars2$K = cpars$K[keep]; cpars2$L50 = cpars$L50[keep]; cpars2$Linf = cpars$Linf[keep]; 
  cpars2$L50 = cpars2$L50 /cpars2$Linf; cpars2$Linf[] = 1
  if(plot) {par(mfrow = c(1,3)); plot(cpars2$M, cpars2$K); plot(cpars2$M, cpars2$L50); plot(cpars2$K, cpars2$L50)}# plot(cpars2$L50, cpars2$Linf); plot(cpars2$M,cpars2$Linf); plot(cpars2$K, cpars2$Linf)}
  cpars2
}
