



makerawdata_4 = function(cdat, Brange = c(0.025,4)){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; inc_CAL = T; inc_CAA = T
  #cdat = as.data.frame(rbindlist(allout3))
  cnams = names(cdat)
  s1dat = cdat[,grepl("_s1",cnams)]
  s2dat = cdat[,grepl("_s2",cnams)]
  s3dat = cdat[,grepl("_s3",cnams)]
  colnms = sapply(names(s1dat),function(x)stringr::str_remove(x,"_s1"))
  names(s1dat) = names(s2dat) = names(s3dat) = colnms
  dat = rbind(s1dat,s2dat,s3dat)
  dat = dat[,!(names(dat)%in%c("CF_T"))]
  
  dat = log(dat)
  dat = cleandat_2(dat,resprange = log(Brange))                       # clean NAs and Infs
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat
  
}

