

Ind2Export3=function(){
  sfLibrary(MSEtool)
  sfLibrary(mvtnorm)
  sfExport(list = c("ConstE_MMP","trim_MMSE","make_MOM3","make_Stock","make_Fleet",
                    "make_Fleets","dummy_effort","dummy_selectivity","make_Obs",
                    "make_Ob","make_Imps","make_Imp","make_cpars3","make_blank_cpars",
                    "add_M_K_L50_Linf3","get_seed","Stoch_effort3","int_effort3",
                    "add_F_dynamics3","add_dep3","getrelF", "make_Catch_Frac3",
                    "get_CurF","calc_Catch_Frac","largedir"))
}


runbatch3 = function(x, nsim=5){
  MOM = make_MOM3(nsim = nsim, seed=x)
  hist = SimulateMOM(MOM) 
  # saveRDS(hist,paste0(largedir,"/hist/hist_",i,".rds"))
  MMSE = ProjectMOM(hist, "ConstE_MMP", checkMPs = FALSE)
  MMSE = trim_MMSE(MMSE)
  saveRDS(MMSE, paste0(largedir,"/MMSE/MMSE_",x,".rda"))
}
