


dolog_3=function(dat){
  dolog(dat, types = c("I_rel","C_rel","FM_cur","FM_rel","ML_cur",
                             "ML_rel","MV_cur","MV_rel","ML_Linf", "MA_cur", "MA_rel",
                             "ML_L50","CR_rel","maxa","K","M_K","L50_Linf","L5_L50", 
                             "LFS_L50","Csd","CF","Isd"))
}

fixnames3 = function(allout3){
  weirdnames = paste0(c("Brel","K","M_K","maxa","L50_Linf")
  
}
                          
                         
makerawdata_3 = function(allout3, sno=1, isBrel = F, clean = T, 
                       inc_Irel = T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T,
                       stock_in = NA, fleet_in = NA, Bmin = 0.1){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; inc_CAL = T; inc_CAA = T
  cdat = as.data.frame(rbindlist(allout3))
  if("Brel_s1_12"%in%names(cdat))allout3 = fixnames3(allout3)
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  stocks = sapply(dnames[Brelcols],function(x)strsplit(x,"_")[[1]][2])
  if(is.na(stock_in[1]))stock_in = stocks[stocks!=sno]
  if(is.na(fleet_in[1]))fleet_in = 1:3
  
  cdat = stock_subsetter(cdat, sno, stock_in) # keep only those listed in sno and stock_in
  cdat = fleet_subsetter(cdat, fleet_in)                # keep only those fleets in nfleet
  cdat = data_subsetter(cdat,  inc_Irel, inc_I, inc_CR, inc_CAL, inc_CAA) # only include those data that are specified to be available
  dat = resp_subsetter(cdat, sno, isBrel)                     # gets the right response variable according to sno
  dat = dolog_2(dat)                                  # log imperfect fractions 
  dat = dologit(dat,types = "VML")                    # logit proportions (but rescaled 0.05 - 0.95 prior to logit)
  if(clean) dat=cleandat_2(dat)                       # clean NAs and Infs 
  
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat = dat[dat$Res > log(Bmin),]
  dat
  
}


