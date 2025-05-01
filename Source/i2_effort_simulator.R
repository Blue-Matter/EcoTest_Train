# Indicator 2 simulation set up



int_effort = function(ny){  
  inc = runif(1,0,0.1)
  dec = runif(1,-0.05,0)
  pow = runif(1,0.7,1.5)
  itc = runif(1,0,0.1)
  amp = runif(1,0.1, 2)
  
  slp = runif(1, -1/75, 1/75)
  per = runif(1, 0.02, 0.13)
  yind = 0:(ny-1)
  yindp = (yind*per)^pow
  cyc = cyc2 = exp(sin(yindp)) * (1:ny)*inc
  
  decpos = sample(seq(ceiling(ny/3),ceiling(2*ny/3)),1)
  cind = decpos:ny
  deci = (1:length(cind))*dec
  cyc2[cind] = cyc[cind] +deci
  rnd = cyc2
  rnd2 =rnd
  if(which.min(rnd)>10){
    grad = (rnd[1]-rnd[which.min(rnd)])/which.min(rnd)
    rnd2 = rnd + (1:ny)*grad
  }
  rndi = rnd2 - min(rnd2)+itc
  rndi=rndi/max(rndi)
 
}  


Stoch_effort = function(ny=75, plot=F, nstocks = 3, Ecor = 0.6){
  # ny=75; plot=F; nstocks = 3; cor = 0.5
  int1 = int_effort(ny)
  int2 = int_effort(ny)
  rat = runif(1)
  eff = (int1*rat)+(int2*(1-rat))
  eff = eff/max(eff)
  
  CV = runif(1, 0.1, 0.3) # this is the normal CV before creating lognornal dist
 
  strt = sample(1:40,1)
  smult = ((1:strt)^2.5)/strt^2.5
  effi = eff
  effi[1:strt]= eff[1:strt]*smult
  
  sigma = array(Ecor,c(nstocks, nstocks))
  diag(sigma) = 1
  errs = rmvnorm(ny, mean=rep(0, nstocks), sigma = sigma )
  errs = errs * CV
  errs = exp(errs)
  effsi = effi * errs
  
  if(plot){
    plot(effsi[,1],pch=19,col="grey")
    points(effsi[,2],pch=1,col="darkgrey")
    grid(col="red")
    lines(eff, col="black",lwd=2)
    lines(int1,col="green",lwd=2,lty=2)
    lines(int2,col="blue",lwd=2)
    legend('topleft',legend=c(paste0("strt = ",strt),
                              paste0("rat = ",round(rat,2)),
                              paste0("Ecor = ",round(Ecor,2)),
                              paste0("CV =", round(CV,2))))
  }
  
  effsi
}

cat("Effort simulator code for Indicator 2 loaded \n")

