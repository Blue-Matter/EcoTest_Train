# curplot=T; slpplot=T; x=Cobs[1,]; anncol = "#0000ff99"; anncol2 ="green"
exp_trend = function(x, Iyr=95, yrs, nyears,ylab="",
                     curplot=F, slpplot=F, muplot = F, relplot=F, 
                     anncol = "#0000ff99",anncol2 ="green",pos="left",
                     pointcol="#999999", smoothcol ="#ff000099" ){
  
  Is<-smooth2(x,plot=F)
  ind = 1:Iyr
  plot(yrs[ind],x[ind],pch=19,ylim=c(0,quantile(x,0.99)),col=pointcol,xlab="",ylab=""); grid()
  abline(v=yrs[nyears]+0.5,lty=2,lwd=2)

  lines(yrs[ind],Is[ind],col=smoothcol,lwd=2)
  cur<-Is[Iyr]
  if(curplot){
    points(yrs[Iyr],cur,pch=2,cex=2,lwd=2,col=anncol)
    legend(pos,legend="Current smoothed value (c)",bty='n',cex=1,text.col=anncol)
  }
  if(slpplot){
    s5<-slp3(Is[Iyr-(5:1)])
    s10<-slp3(Is[Iyr-(10:1)])
    abline(v=yrs[Iyr]-c(6,11),col=c(anncol,anncol2),lwd=2,lty=2)
    xs= c(1,Iyr)
    
    liney = exp(s5*xs) *mean(Is[Iyr-(5:1)])
    mumean = mean(Is[Iyr-(5:1)])
    yatx = approx(c(1,Iyr),liney,Iyr-3.5)$y
    lines(yrs[xs],liney*mumean/yatx,col=anncol)
    
    liney = exp(s10*xs) *mean(Is[Iyr-(10:1)])
    mumean = mean(Is[Iyr-(10:1)])
    yatx = approx(c(1,Iyr),liney,Iyr-6)$y
    lines(yrs[xs],liney*mumean/yatx,col=anncol2)
    
    legend(pos,legend=c("Smoothed slope last 5 years (5)", "Smoothed slope last 10 years (10)"),
           text.col=c(anncol, anncol2),bty="n")
  }
  if(muplot){
    mu = mean(x[nyears:Iyr])
    lines(yrs[c(nyears,Iyr)],rep(mu,2),col=anncol,lwd=2)
    legend(pos,legend=c("Mean value in projection (m)"),bty="n",text.col=anncol)
  }
  if(relplot){
    allmu = mean(x[1:Iyr])
    abline(h=allmu,col=anncol,lwd=2)
    points(yrs[Iyr],cur,pch=2,cex=2,lwd=2,col=anncol)
    legend(pos,legend=c("Current value relative to series mean (r)"),bty="n",text.col=anncol)
  }
  
}

