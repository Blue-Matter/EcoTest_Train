
plot_at_age = function(hist, ss, ff, slot= "Weight"){
  dat = hist[[ss]][[ff]]@AtAge[[slot]]
  matplot(t(dat[,,1]),type="l",ylab=slot,xlab="Age")
}  

doSageplot = function(hist, slot){
  nstocks = length(hist);par(mfrow=c(1,nstocks))
  for(ss in 1:nstocks) plot_at_age(hist, ss, 1, slot)
}

plot_SP = function(hist, ss, ff, slot= "V"){
  dat = hist[[ss]][[ff]]@SampPars$Fleet[[slot]]
  matplot(t(dat[,,1]),type="l",ylab=slot,xlab="Age")
}  

doFageplot = function(hist, slot, at_age=T){
  nstocks = length(hist); nfleets = length(hist[[1]]); par(mfrow=c(nfleets,nstocks))
  for(ff in 1:nfleets){
    for(ss in 1:nstocks){
      if(at_age)plot_at_age(hist, ss, ff, slot)
      if(!at_age)plot_SP(hist, ss, ff, slot)
    } 
  }
}

testy = function(){
  SB = apply(hist[[ss]][[ff]]@TSdata$SBiomass,1:2,sum)
  matplot(t(SB),type="l")
  MOMtest@cpars[[ss]][[ff]]$D
  SB[,ncol(SB)]/SB[,1]
}

quick_stock_plots = function(hist){
  doSageplot(hist, "Weight")
  doSageplot(hist, "Length")
  doSageplot(hist, "Maturity")
}


quick_fleet_plots = function(hist){
  doFageplot(hist, "Select")
  doFageplot(hist, "V",at_age=F)
  
  
  
}
