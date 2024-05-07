# Tools for making correlated F trends
rnorm_T95<-function(n=1,mean=0,sd=1){
  ntrial=n*ceiling(48*(1/n)^0.5+2)  #prevents over sampling at high n; 0.05^50 = 1E-65 chance of not sampling in the interval given n=1
  trial<-rnorm(ntrial,mean,sd)
  bound<-sd*1.959964
  if(sd[1]==0){
    return(rep(mean,n)) # exception for sd=0 (Deterministic == T)
  }else{
    return(trial[trial>(mean-bound) & trial < (mean+bound)][1:n])
  }
}



smooth4<-function(xx,plot=F,enp.mult=0.1,plotname=""){
  tofill<-!is.na(xx)
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  predout[tofill]<-exp(predict(out))
  if(plot){
    plot(xx,type="p",xlab="x",ylab="y",main=plotname)
    lines(predout,col="#ff000090",lwd=2)
  }
  predout
}

ACfunc = function(arr,autocor=0.98, inds = 3:6, enp.mult=0.2,ploty = F){
  
  dims = dim(arr)
  arrAC = arr
  if(length(autocor)==1)autocor=rep(autocor,dims[2])
  
  repi = ceiling(seq(1E-5,1,length.out=20)*dims[1])
  
  for(j in inds){
    for(i in 1:dims[1]){
      var = arr[i,j,]
      smt = smooth4(var, plot=ploty, enp.mult=enp.mult)
      res = resAC = log(var/smt)
      rescv = sd(res)
      for(y in 2:dims[3])resAC[y] <- autocor[j]*resAC[y-1] + ((1-autocor[j]^2)^0.5)* rnorm_T95(1,0,rescv) - (1-autocor[j])* (rescv^2)/2
      arrAC[i,j,] = smt * exp(resAC)
      if(i %in% repi)cat("*")
    }
    
  }
  cat("\n")
  
  par(mfrow=c(2,2)) 
  sims = sample(1:dims[1],2);stock = sample(3:dims[2], 1)
  matplot(t(arr[sims,stock,]),type="l"); matplot(t(arrAC[sims,stock,]),type="l")
  plot(arr[sims[1],3,],arrAC[sims[1],3,]); plot(arr[sims[1],4,],arrAC[sims[1],4,])
  arrAC
}

Simfunc = function(object,newdat, slopeadj = 1,RSD_adj = 5){
  
  rel_sum <- summary(object)
  SE = sapply(rel_sum,function(x)x$coefficients[,2])*slopeadj
  f0 <- predict(object,newdat = newdat)
  simdat = array(NA,dim(f0))
  nc = ncol(f0)
  nm <- dimnames(f0)
  n <- nrow(f0)
  nin = ncol(newdat)
  RSD0 <- sapply(rel_sum, getElement, "sigma") # Residual standard deviation (normal space)
  RSD_LN <- sdconv(1, RSD0) * RSD_adj # Residual standard deviation (lognormal space)
  
  for(s in 1:nc){
    interr = rep(SE[,s],each=n)
    adj = -0.5 * interr^2
    newerdat = newdat * rlnorm(nin*n,adj,interr)
    ftemp = predict(object,newdata=newerdat)
    simdat[,s] = ftemp[,s]* rlnorm(n, -0.5 * RSD_LN^2, RSD_LN)
  }
  
  simdat
}  



plot_an_F_sim = function(muE_AC,simno=1,snames){
  par(mfcol=c(4,2))
  for(i in 1:2){
    for(j in 3:6){
      plot(muE_AC[simno,i,],muE_AC[simno,j,],type="l")
      if(i == 1)mtext(snames[j],2,line=2.6)
      if(j == 3)mtext(snames[i],3,line=0.6)
    }
  }
  
}



