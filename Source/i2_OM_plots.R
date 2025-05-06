
dohist = function(x, xlab, panellab){
  hist(x, freq=FALSE, breaks=20, xlab=xlab, main="", ylab = "", col="slategrey", border="white")
  mtext(panellab,adj=0.01,line=-0.45,cex=0.65)
} 

plot_OM_LH = function(MOM){
  
  par(mfrow=c(2,3),mai = c(0.7,0.4,0.05,0.05),omi=c(0.025,0.4,0.025,0.025))
  
  Ms = unlist(lapply(MOM@cpars,function(x)x[[1]]$M))
  Ks = unlist(lapply(MOM@cpars,function(x)x[[1]]$K))
  Ds = unlist(lapply(MOM@cpars,function(x)x[[1]]$D))
  L50s = unlist(lapply(MOM@cpars,function(x)x[[1]]$L50))/100
  hs = runif(length(Ms),MOM@Stocks[[1]]@h[1],MOM@Stocks[[1]]@h[2])
  
  dohist(Ms, "Natural Mortality Rate (M)(per year)","(a)")
  dohist(Ks, "Somatic Growth (K)(per year)","(b)")
  dohist(Ms/Ks, "M/K ratio","(c)")
  dohist(L50s, "Length at 50% Maturity / Asymptotic Length","(d)")
  dohist(Ds, "SSB(2024) / SSB(unfished)","(e)")
  dohist(hs, "Steepness of Stock Recruitment Function","(f)")
  
  mtext("Relative Frequency",2,outer=T)
  
}


plot_OM_Find = function(MOM){
  
  par(mfrow=c(2,2),mai = c(0.7,0.4,0.05,0.05),omi=c(0.025,0.4,0.025,0.025))
  cols = c("black","darkgrey","blue","red","green","orange", "purple" )
  Fs = rbindlist(lapply(MOM@cpars[[1]], function(x)as.data.frame(x$Find)))
  Fs = rbindlist(lapply(MOM@cpars[[2]], function(x)as.data.frame(x$Find)))
  #Fs3 = rbindlist(lapply(MOM@cpars[[3]], function(x)as.data.frame(x$Find)))
  #allFs = rbind(Fs1, Fs2, Fs3)
  ny = dim(Fs)[2]
  yrlab = 2025-ny:1
  nc = length(cols)
  for(i in 1:4){
    ind = nc*(i-1)+1:nc
    allFs = Fs[ind,]
    matplot(yrlab,t(allFs),col="white",lty=1,type="l",lwd=2,xlab="",ylab="")
    grid()
    matplot(yrlab,t(allFs),col=cols,lty=1,type="l",lwd=2,add=T)
    mtext(paste0("(",letters[i],")"),adj=0.02,line=-1.1,cex=0.85)
  }
 
  mtext("Relative Apical Fishing Mortality Rate",2,outer=T)
  
}

plot_OM_Vuln = function(MOM){
  
  L5s = unlist(lapply(MOM@cpars,function(x)x[[1]]$L5))
  LFSs = unlist(lapply(MOM@cpars,function(x)x[[1]]$LFS))
  Vmaxlens = unlist(lapply(MOM@cpars,function(x)x[[1]]$Vmaxlen))
  L50s = unlist(lapply(MOM@cpars,function(x)x[[1]]$L50))
  
  par(mfrow=c(2,2),mai = c(0.7,0.4,0.05,0.05),omi=c(0.025,0.4,0.025,0.025))
  
  dohist(L5s / L50s, "Length at 5% Vuln. / Length at 50% Maturity","(a)")
  dohist(LFSs / L50s, "Length at Full Vuln. / Length at 50% Maturity","(b)")
  dohist(L5s/LFSs, "Length at 5% Vuln. / Length at Full Vuln.","(c)")
  dohist(Vmaxlens, "Vuln. at Asymptotic Length","(d)")
  
  mtext("Relative Frequency",2,outer=T)
  
}

plot_OM_F_cor = function(MOM){
  Fs1 = MOM@cpars[[1]][[1]]$Find
  Fs2 = MOM@cpars[[2]][[1]]$Find
  ny = dim(Fs)[2]
  yrlab = 2025-ny:1
  
  inds = c(5,6,7)
  inds = 8:10
  par(mfrow=c(3,2),mai = c(0.5,0.5,0.05,0.05),omi=c(0.025,0.3,0.025,0.025))
  for(i in inds){
    plot(yrlab,Fs1[i,],pch=19,col="#ff000099",xlab="",ylab="");grid()
    mtext("Fishing Effort",2,line=2.6,cex=0.9)
    mtext("Year",1,line=2.6,cex=0.9)
    if(i == inds[1])legend('topleft',legend=c("Stock 1","Stock 2"),text.col = c("red","blue"),bty='n')
    points(yrlab,Fs2[i,],pch=19,col="#0000ff99")  
    plot(Fs1[i,],Fs2[i,],pch=19,col="#99999999",xlab="",ylab=""); grid()
    mtext("Effort on Stock 1",1,line=2.6,cex=0.9)
    mtext("Effort on Stock 2",2,line=2.6,cex=0.9)
  }  
  
}


