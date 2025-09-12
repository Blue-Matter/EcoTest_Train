


dotrans=function(dat){
  # log things
  nolog = c("Res","VML")
  boolarray = sapply(nolog,function(x)grepl(x,names(dat)))
  dolog = !apply(boolarray,1,any) 
  dat[,dolog] = log(dat[,dolog])
  # logit things
  types = "VML"; bnds = c(0.025,0.975)
  boolarray = sapply(types,function(x)grepl(x,names(dat)))
  tologit = apply(boolarray,1,any) 
  rng = bnds[2] - bnds[1]
  dat[,tologit] = bnds[1]+dat[,tologit]*rng
  dat[,tologit] = MSEtool:::logit(dat[,tologit])
  dat
}
                      
                         
makerawdata_3 = function(allout3, Brange = c(0.025,4)){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; inc_CAL = T; inc_CAA = T
  cdat = as.data.frame(rbindlist(allout3))
  dnames = names(cdat)
  cdat = cdat[,!(names(cdat)%in%c("Brel_s2","Brel_s3"))]
  dat = resp_subsetter(cdat, sno=1, isBrel=F)                     # gets the right response variable according to sno
  dat = dotrans(dat)                                  # log imperfect fractions 
  dat = cleandat_2(dat,resprange = log(Brange))                       # clean NAs and Infs 
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat
  
}


