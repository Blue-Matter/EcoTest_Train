


cleandat = function(dat,resprange=c(0.05,2)){
  goodrow = function(x)(sum(is.na(x)| x==0 | x==1 | x == 'Inf' | is.null(x)| x=="NULL")==0)&x[1]>resprange[1] & x[1]<resprange[2]
  dat[apply(dat,1,goodrow),]
}

reformat = function(x){
  
  renam = grepl("-",x)
  charmatch("-",x[162])
  
}

makerawdata = function(allout, sno=1, isBrel = F, clean = T, 
                       inc_spat = T, inc_Irel = T, inc_L50_Linf = F){
  
  cdat = as.data.frame(rbindlist(allout))
  
  
  #for(i in 1:ncol(cdat)) names(cdat)[i] = gsub("-","_",names(cdat)[i])
  dnames = names(cdat)
 
  Brelcols = grepl("Brel",dnames)
  ns = sum(Brelcols)
  
  keepcol = grep(paste0("Brel_",sno),dnames)
  
  Res = cdat[,keepcol]
  
  # level relative to mean level of other stocks
  if(isBrel){Res = ((ns-1) * Res) /(apply(cdat[,Brelcols],1,sum)-Res)}
  
  dat = cbind(Res,cdat[,!Brelcols])
  
  if(!inc_spat) dat = dat[,!grepl('spat',names(dat))]
  if(!inc_Irel) dat = dat[,!grepl('I_rel',names(dat))]
  if(!inc_L50_Linf){
    cond = grepl("L50",names(dat)) & !grepl("ML_L50",names(dat))
    dat= dat[,!cond]
    cond = grepl("Linf",names(dat)) & !grepl("ML_Linf",names(dat))
    dat= dat[,!cond]
  }
  
  # log ratios
  cond = grepl("C_rel", names(dat)) | grepl("CR_mu",names(dat)) | grepl("FM_rel", names(dat)) |
    grepl("ML_cur", names(dat)) | grepl("ML_rel" , names(dat)) | grepl("MV_cur",names(dat)) |
    grepl("MV_rel",names(dat))| grepl("ML_Linf" ,names(dat)) |   grepl("ML_L50" , names(dat)) |
    grepl("CR_rel",names(dat))
  
  dat[,cond] = log(dat[,cond])
  
  if(clean)dat=cleandat(dat)
  dat
  
}