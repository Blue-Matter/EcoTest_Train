make_MOM3 = function(Name = "Simulated for AI training", nsim = 100, 
                    proyears = 40, interval = 1, seed = 1, 
                    nstocks = 3, nfleets = 3, nyears = 70,
                    maxage = 30,
                    CurrentYr = 2024, plot = F){
  
  # Name = "Simulated for AI training"; nsim = 10; proyears = 50; interval = 1; seed = 1; nstocks = 3; nfleets = 3; nyears = 75; CurrentYr = 2024; plot=F
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
  MOM@Stocks = lapply(1:nstocks,make_Stock, maxage = maxage)
  MOM@Fleets = make_Fleets(nstocks, nfleets, nyears, CurrentYr)
  MOM@Obs = make_Obs(nstocks, nfleets)
  MOM@Imps = make_Imps(nstocks, nfleets)
  MOM@cpars = make_cpars3(nstocks, nfleets, nsim, nyears, plot = plot, seed = seed)
  MOM@CatchFrac = make_Catch_Frac3(MOM)
  MOM
}


make_cpars3 = function(nstocks, nfleets, nsim, nyears, M_max = 0.6, M_min = 0.05, Linf_max = 1000, 
                      plot=F, Ecor = 0.6, seed){
  cpars = make_blank_cpars(nstocks, nfleets)
  set.seed(seed)
  
  # have to make sure its the same stock dynamics across fleets
  cpars_by_stock = list()
  
  for(ss in 1:nstocks){
    cpars_by_stock[[ss]] = list()
    cpars_by_stock[[ss]] = add_M_K_L50_Linf3(cpars_by_stock[[ss]], nsim, M_max, M_min, Linf_max, plot=plot)
  }
  
  for(ff in 1:nfleets){
    Effs = list()
    for(i in 1:nsim) Effs[[i]]= Stoch_effort3(ny=nyears, plot=F, nstocks = nstocks, Ecor = Ecor)
    for(ss in 1:nstocks){
      cpars[[ss]][[ff]] = cpars_by_stock[[ss]]
      cpars[[ss]][[ff]] = add_F_dynamics3(cpars[[ss]][[ff]],Effs,ss)
    }
  }
  
  ff=1
  for(ss in 1:nstocks)  cpars[[ss]][[1]]$D = add_dep3(cpars_ss = cpars[[ss]])
  
  
  cpars
}


cat("Multi-stock, multi-fleet simulator code for Indicator 3 loaded \n")

