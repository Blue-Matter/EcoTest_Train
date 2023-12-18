# Functions for extracting indices

slp3<-function(y0){
  y=log(y0)
  x1<-1:length(y0)
  x1[is.na(y)]=NA
  mux<-mean(x1,na.rm=T)
  muy<-mean(y,na.rm=T)
  SS<-sum((x1-mux)^2,na.rm=T)
  (1/SS)*sum((x1-mux)*(y-muy),na.rm=T)
  
}


smooth2<-function(xx,plot=F,enp.mult=0.3,plotname=""){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  #xx[tofill] = log(xx[tofill])
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  #out<-loess(y~x,dat=dat,span=0.5)
  predout[tofill]<-exp(predict(out))
  if(plot){
    plot(xx,type="p",xlab="x",ylab="y",main=plotname)
    lines(predout,col="#ff000090",lwd=2)
  }
  predout
}

proc_dat<-function(MMSE,sno = 1, fno=1){
  
  dat = MMSE@PPD[[sno]][[fno]][[1]]
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
 
  datlist<-list()
  
  Iobs<-dat@Ind
  Cobs<-dat@Cat
  CAL = dat@CAL
  mids = dat@CAL_mids
 
  #prob<-seq(1,0.5,length.out=length(allyears))
  Iyr<-sample((nyears+1):(allyears-2),nsim,replace=T)#,prob=prob)
  Byr = Iyr - nyears
  
  Iind<-cbind(1:nsim,Iyr)
 
  Bt = MMSE@SB_SBMSY
  I_cur<-I_mu<-I_rel<-I_s5<-I_s10<- C_cur<- C_mu<-C_rel<-C_s5<-C_s10<- 
    ML_cur<-ML_mu<-ML_rel<-ML_s5<-ML_s10<- MV_cur<- MV_mu<-MV_rel<- MV_s5<-MV_s10<-
    FM_cur <- FM_mu <- FM_rel <- FM_s5 <- FM_s10 <- rep(NA,nsim)
  
  plotsmooth=F
  
  L50 = MMSE@OM[[sno]][[fno]]$L50
  
  for(i in 1:nsim){
    # Index
    Is<-smooth2(Iobs[i,],plot=plotsmooth)
    I_cur[i]<-Is[Iind[i,2]]
    I_mu[i] = mean(Iobs[i,nyears:Iind[i,2]])
    I_rel[i]<-I_cur[i] / mean(Is[1:Iind[i,2]])
    I_s5[i]<-slp3(Iobs[i,Iind[i,2]-(5:1)])
    I_s10[i]<-slp3(Iobs[i,Iind[i,2]-(10:1)])
    
    # Catch
    Cs = smooth2(Cobs[i,],plot=plotsmooth)
    C_cur[i]<-Cs[Iind[i,2]]
    C_mu[i]<-mean(Cobs[i,nyears:Iind[i,2]])
    C_rel[i] <- C_cur[i] / mean(Cs[1:Iind[i,2]])
    C_s5[i]<-slp3(Cobs[i,Iyr[i]-(5:1)])
    C_s10[i]<-slp3(Cobs[i,Iyr[i]-(10:1)])
    
    # Length
    CALi = CAL[i,,]
    totlen = CALi*array(mids,dim(CALi))
    mulen = apply(totlen,1,sum)/apply(CALi,1,sum)
    Ls = smooth2(mulen,plot=plotsmooth)
    ML_cur[i]=Ls[Iind[i,2]]
    ML_mu[i]=mean(mulen[nyears:Iind[i,2]])
    ML_rel[i] <- ML_cur[i] / mean(Ls[1:Iind[i,2]])
    ML_s5[i]<-slp3(mulen[Iyr[i]-(5:1)])
    ML_s10[i]<-slp3(mulen[Iyr[i]-(10:1)])
    
    # standard deviation in length
    CVlen = sapply(1:dim(CALi)[1],function(X,CALi,mids){
      sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
     },CALi=CALi,mids=mids)
    CVs = smooth2(CVlen, plot=plotsmooth)
    MV_cur[i]=CVs[Iind[i,2]]
    MV_mu[i] = mean(CVlen[nyears:Iind[i,2]])
    MV_rel[i] <- MV_cur[i] / mean(CVs[1:Iind[i,2]])
    MV_s5[i]<-slp3(CVlen[Iyr[i]-(5:1)])
    MV_s10[i]<-slp3(CVlen[Iyr[i]-(10:1)])
    
    L50i = L50[i]
    Fmat = sapply(1:dim(CALi)[1],function(X,CALi,mids,L50i){
      indivs = rep(mids,CALi[X,])
      mean(indivs > L50i,na.rm=T)
    },CALi=CALi,mids=mids,L50i=L50i)
    FMs = smooth2(Fmat, plot=plotsmooth)
    FM_cur[i]= FMs[Iind[i,2]]
    FM_mu[i] = mean(FMs[nyears:Iind[i,2]])
    FM_rel[i] <- FM_cur[i] / mean(FMs[1:Iind[i,2]])
    FM_s5[i]<-slp3(Fmat[Iyr[i]-(5:1)])
    FM_s10[i]<-slp3(Fmat[Iyr[i]-(10:1)])

  }
  
  Bind<-cbind(1:nsim,rep(sno,nsim),rep(1,nsim),Byr)
  Brel = MMSE@SB_SBMSY[Bind]
 
  data.frame(Brel, I_rel, I_s5, I_s10, C_rel, C_s5, C_s10, 
             ML_rel, ML_s5, ML_s10, MV_rel, MV_s5, MV_s10,
             FM_rel, FM_s5, FM_s10)
  
}

filt=function(x){
  x[x%in%c("Inf","-Inf")] = NA
  x
} 


corplot = function(out,quanty=0.5,ptcex=0.8,maxn=20, lasinv=F,lab=""){
  
  nams = names(out)
  ndat = dim(out)[2]
  nobs = dim(out[1])
  ni <- min(ndat,maxn)
  
  ncol = 200
  cols <- rainbow(ncol,start=0.1,end=0.3,alpha=0.7)
  cols <- c(cols,rep(cols[ncol],1000-ncol))
  Brel = out$Brel
  Bcol = cols[ceiling(Brel*100)]; plot(Brel[order(Brel)],col=Bcol[order(Brel)],pch=19)
  par(mfrow=c(ni-1,ni-1),mai=rep(0,4),omi=c(0.55,0.75,0.05,0.05))
  cutoff= c(quanty/100,(100-quanty)/100)
  # labys = c(letters,paste0(letters,letters),paste0(letters,letters,letters),paste0(letters,letters,letters,letters))
  labys = paste0(rep(letters,9),rep(1:9,each=length(letters)))
  labys = paste0("(",labys,")")
  k = 0
  
  for(i in 2:ni){ # i<-i+1
    
    for(j in 1:(ni-1)){
      
      if(j==i|j>i){
        
        plot(1,1,col='white',axes=F)
        
      }else{
        xs = filt(out[,j])
        ys = filt(out[,i])
        xlim<-range(quantile(xs,cutoff,na.rm=T))
        ylim<-range(quantile(ys,cutoff,na.rm=T))
        plot(xs,ys,pch=19,xlim=xlim,ylim=ylim,cex=ptcex,col=Bcol,axes=F)
        k=k+1; mtext(labys[k],3,adj=0.02,line=-1,cex=0.6)
        axis(1,c(-1000,1000))
        axis(2,c(-1000,1000))
        axis(3,c(-1000,1000))
        axis(4,c(-1000,1000))
        #points(obs[j],obs[i],pch=19,cex=1.2,col=cols[2])
        if(nams[j]=="Brel")abline(v=c(0.5,1),lty=c(2,1),col="#99999999")
        
      }
      if(i==2&j==(ni-1)){
        legend('center',legend=lab,text.col="blue",bty='n',text.font=2,cex=1.2)
      }
      
      if(!lasinv){
        if(j==1)mtext(nams[i],2,line=2,cex=0.6,las=2)
        if(i==ni)mtext(nams[j],1,line=1,cex=0.6,las=2)
      }else{
        if(j==1)mtext(nams[i],2,line=2,cex=0.6,srt=90)
        if(i==ni)mtext(nams[j],1,line=1,cex=0.6,las=1)
      }
      #if(j==1)mtext(i,2,line=2,cex=0.5,las=2)
      #if(i==nplotted)mtext(j,1,line=1,cex=0.5,las=2)
      
    }
    
  }
  
#}

  
  
  
  
}



