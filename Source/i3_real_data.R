# MMSE = readRDS("C:/Users/tcar_/Dropbox/temp/Ecotest/Ind3/MMSE/MMSE_1001.rda")
# multiHist = MMSE@multiHist
#  L50 = 29; Linf = 51.47; M = 0.75; K = 0.32; Mat_slope = 0.8; Len_age = 51.47 * (1-exp(-0.32 * (c(0.5,1.5,2.5,3.5)+0.83)))
#  Cat_f = list(dat_FRI$Catch); CAL_f = list(dat_FRI$CAL); I_f = NA; L5_f = list(32); LFS_f = list(42); Vmaxlen_f = list(0.99); CAL_mids = dat_FRI$Lbins
#  ff=1
make_MMSE_I3 = function(L50, Linf, M, K, Mat_slope, Len_age,
                        Cat_f, CAL_f, I_f, L5_f, LFS_f, Vmaxlen_f,
                        CAL_mids, plot=T){

  nyears = length(Cat_f[[1]])
  nf = length(Cat_f)
  nages = length(Len_age)
  MMSE = new('MMSE')
  MMSE@nsim = 1
  MMSE@nyears = nyears-2
  MMSE@proyears = 2

  PPD = OM = multiHist = list() # hierarchical lists
  PPD[[1]] =  multiHist[[1]] = OM[[1]] = list() # one stock
  Mat_age = domat(Len_age, L50, Mat_slope=Mat_slope, plot=F)
  Sel_len = DoubleNormal(CAL_mids,L5_f[[1]],LFS_f[[1]],Vmaxlen_f[[1]])
  Sel_age = approx(CAL_mids,Sel_len,Len_age)$y
  #if(plot){
   # plot(CAL_mids,Sel_len); abline(v=Len_age); abline(h=Sel_age)
  #}
  Sel_age = Sel_age/max(Sel_age)


  # Fleet specific
  for(ff in 1:nf){

    # Posterior predicted data
    Data = new("Data")
    Data@Cat  = matrix(Cat_f[[ff]],nrow=1)  # nsim, nyear
    Data@VInd = matrix(I_f[[ff]],  nrow=1)  # nsim, nyear
    Data@CAL = CAL_f[[ff]]  # dat@CAL[i,,] # nsim, nyear, n CAL
    Data@CAL_mids = CAL_mids

    PPD[[1]][[ff]] = list()
    PPD[[1]][[ff]][[1]] = Data  # only one MP

    # Sampled parameters
    multiHist[[1]][[ff]]=new("Hist")
    multiHist[[1]][[ff]]@SampPars = list()
    multiHist[[1]][[ff]]@SampPars$Fleet$L5_y = matrix(L5_f[[ff]], nrow=1)
    multiHist[[1]][[ff]]@SampPars$Fleet$LFS_y = matrix(LFS_f[[ff]], nrow=1)
    multiHist[[1]][[ff]]@SampPars$Fleet$Vmaxlen_y = matrix(Vmaxlen_f[[ff]], nrow=1)
    multiHist[[1]][[ff]]@AtAge$Maturity = array(Mat_age, c(1, nages, 1))
    multiHist[[1]][[ff]]@AtAge$Length = array(Len_age, c(1, nages, 1))

    multiHist[[1]][[ff]]@AtAge$Selectivity = array(Sel_age, c(1, nages, 1))

  }
  class(multiHist) = "multiHist"

  OM[[1]][[1]] = list()
  OM[[1]][[1]]$L50 = L50
  OM[[1]][[1]]$Linf = Linf
  OM[[1]][[1]]$M = M
  OM[[1]][[1]]$K = K

  MMSE@multiHist = multiHist
  MMSE@OM = OM
  MMSE@PPD = PPD

  MMSE
}


prep_NN_input=function(L50, Linf, M, K, Mat_slope, Len_age ,Cat_f, CAL_f, I_f, L5_f, LFS_f, Vmaxlen_f, CAL_mids, yrs, plot=T){
  MMSE = make_MMSE_I3(L50, Linf, M, K, Mat_slope, Len_age, Cat_f, CAL_f, I_f, L5_f, LFS_f, Vmaxlen_f,CAL_mids,plot=plot)
  make_real_NN_input(MMSE, yrs=yrs, plotsmooth=plot)
}

bootcal = function(CAL_f){
  newCAL = CAL_f
  ny = nrow(CAL_f[[1]][1,,])
  for(i in 1:ny){
    ns = CAL_f[[1]][1,i,]
    if(!all(ns==0)){
      prob = ns/sum(ns)
      newCAL[[1]][1,i,]=  rmultinom(1,sum(ns),prob)[,1]
    }
  }
  newCAL
}

bootcat = function(Cat_f, catf_CV){
  newCat = Cat_f
  ny = ncol(Cat_f[[1]])
  newCat[[1]][1,] = MSEtool::trlnorm(ny,1,catf_CV) *  Cat_f[[1]][1,]
  newCat
}

bootint = function(x,L50, L50_CV,M,K,Mat_slope, Len_age, MK_CV,M_CV,Linf,Linf_CV,LFS_f,LFS_CV,L5_f,CAL_f, I_f,LFS_L5_CV,Vmaxlen_f=Vmaxlen_f,Cat_f,catf_CV,CAL_mids,yrs){
  set.seed(x)
  L50val =trl(2,L50,L50_CV)[1]
  MKrat = trl(2,M/K,MK_CV)[1]
  Mval = trl(2,M,M_CV)[1]
  Kval = Mval / MKrat
  Linfval = trl(2,Linf,Linf_CV)[1]
  LFSval = list(trl(2,LFS_f[[1]],LFS_CV)[1])
  LFS_L5 = LFS_f[[1]]/L5_f[[1]]
  L5val = list(LFSval[[1]] / trl(2,LFS_L5,LFS_L5_CV)[1])
  MMSE = make_MMSE_I3(L50=L50val,Linf= Linfval , M=Mval, K=Kval, Mat_slope, Len_age,
                      Cat_f = bootcat(Cat_f,catf_CV), CAL_f = bootcal(CAL_f),
                      I_f, L5_f=L5val, LFS_f=LFSval, Vmaxlen_f,CAL_mids,plot=F)
  make_real_NN_input(MMSE, yrs=yrs, plotsmooth=F)
}


bootstrap_NN_input = function(L50, Linf, M, K, Mat_slope, Len_age ,Cat_f, CAL_f, I_f, L5_f, LFS_f, Vmaxlen_f, CAL_mids, yrs, plot=T,
                              nboot = 200, L50_CV = 0.1, Linf_CV = 0.025, MK_CV = 0.1, M_CV=0.15,
                              catf_CV = 0.075, LFS_L5_CV = 0.05, LFS_CV = 0.05){
  out = list()
  set.seed(1)
  trl = MSEtool::trlnorm

  setup(cpus=8)
  sfExport(list=c("trl","make_MMSE_I3","make_real_NN_input","bootcal","bootcat","domat","proc_dat_LH_3",
                  "proc_allFdat_F3","calcTsel","intselpars","interpolate","smooth2"))
  out = sfSapply(1:nboot,bootint, L50=L50, L50_CV=L50_CV,M=M,K=K,Mat_slope=Mat_slope, Len_age=Len_age, MK_CV=MK_CV,M_CV=M_CV,Linf=Linf,Linf_CV=Linf_CV,
                 LFS_f=LFS_f,LFS_CV=LFS_CV,L5_f=L5_f, CAL_f= CAL_f, I_f=I_f, LFS_L5_CV = LFS_L5_CV,Vmaxlen_f=Vmaxlen_f,Cat_f=Cat_f,catf_CV = catf_CV, CAL_mids = CAL_mids, yrs=yrs)
  sfStop()
  temp = as.data.frame(t(out))
  for(cc in 1:ncol(temp))temp[,cc] = as.numeric(temp[,cc])
  temp
}



make_real_NN_input = function(MMSE, yrs, plotsmooth=T){
  outs = list()
  Cat = MMSE@PPD[[1]][[1]][[1]]@Cat
  ny = ncol(Cat)
  Iind = matrix(c(1,ny),nrow=1)
  sno = 1
  outs[[1]] = proc_dat_LH_3(MMSE,Iind=Iind,sno=sno,plotsmooth=plotsmooth)
  outs[[2]] = proc_allFdat_F3(MMSE,Iind=Iind,sno=sno,yrs=yrs,plotsmooth=plotsmooth)
  one_tab = cbind(outs[[1]],outs[[2]],ny)
  keep = (1:ncol(one_tab))[!is.na(one_tab[1,])&!grepl("CF",names(one_tab))]
  one_tab = log(one_tab)
  one_tab[,keep]
}

domat = function(Len_age,L50,Mat_slope=0.8,plot=T){
  nage = length(Len_age)
  aa = (1:nage)-0.5
  Len_mat =
    a50 = approx(Len_age,aa,L50)$y
  mat_age = 1/(1+exp(-Mat_slope*(Len_age-L50)))
  if(plot)plot(aa,mat_age,ylab="Spawn. Frac.",xlab="Age")
  mat_age
}


concat = function(MMSE, Ind = 2, sno = 1, plotsmooth = T, nint = 40){
  outs = list()
  out[[1]] = proc_allFdat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  out[[2]] = proc_dat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  out[[3]] = proc_dat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  one_tab = cbind(out[[1]], out[[2]], out[[3]])
  one_tab
}


statplot = function(pred_boot, pred1, yr="2025", main="",xlim=c(0,2),doxlab=T,doylab=T,bty="n",pos='topleft'){
  #par(mai=c(0.7,0.7,0.4,0.05))
  plot(density(pred_boot),col="blue",main="",xlab="",ylab="",xlim=xlim); grid()
  abline(v=c(0.5,1),lty=2,lwd=2,col=c("red","orange"))
  abline(v=pred1,lwd=3,col="black")

  qs = quantile(pred_boot,c(0.05,0.5,0.95)); qsr = round(qs,2)
  abline(v=qs,col="blue",lty=2)
  POFed = mean(pred_boot<1)
  PCS = mean(pred_boot<0.5)
  lines(density(pred_boot),col="blue",lwd=2)
  legend(pos,legend=c(paste0("P(SSB",yr," < SSBMSY) = ",round(POFed*100,2),"%"),
                            paste0("P(SSB",yr," < 0.5 SSBMSY) = ",round(PCS*100,2),"%"),
                            "Deterministic estimate",
                            paste0("Stochastic estimate (",qsr[2]," [",qsr[1],",",qsr[3],"])")),
         text.col=c("orange","red","black","blue"),text.font=c(1,1,2,1),bty=bty)
  if(doxlab)mtext(paste0("SSB/SSBMSY (",yr,")"),1,line=2.4)
  if(doylab)mtext("Relative Freq.",2,line=2.4)
  mtext(main,3,line=0.4,font=2)
}


do_sense = function(boot_ind, pars = c("K_s1","M_K_s1","maxa_s1","L50_Linf_s1"),
                              parnam = c("K","M","Max. Age","L50"),
                              nsense = 9, range = c(0.25,0.25,0.25,0.25)){
  tboot = boot_ind
  tdata = boot_ind$sdata
  np = length(pars)
  res = list(); 
  
  for(pp in 1:np){
    
    loglim = log(1-range[pp])
    sseq = seq(loglim,-loglim,length.out=nsense)
    mults = exp(sseq)
    res[[pp]] = list()
    pind = match(pars[pp],colnames(tdata))
    
    for(ss in 1:nsense){
      
      tempdata = tdata
      tempdata[,pind] = log(exp(tempdata[,pind])*mults[ss])
      tboot$sdata = tempdata
      res[[pp]][[ss]] = pred_ind(tboot)
      
    }
    
    names(res[[pp]])=round(mults,2)
    
  }
  
  names(res) = parnam
  res
  
}

sensepanel = function(subsense, lab){
  levs = sapply(subsense,quantile,p=c(0.05,0.25,0.5,0.75,0.95))
  nlev = ncol(levs)
  ylim = c(0,max(levs)*1.025)
  plot(c(0.5,nlev+0.5),ylim, col="white",xlab="",ylab="",axes=F);grid()
  for(i in 1:nlev){
    lines(rep(i,2),levs[c(1,5),i],lwd=2,col="blue")
    lines(rep(i,2),levs[c(2,4),i],lwd=5,col="blue")
  }
  points(levs[3,],pch=19,cex=1.1,col="blue")
  lines(levs[3,],col="blue")
  axis(2,c(-1000,1000))
  reflev = ceiling(nlev/2)
  abline(h=levs[3,reflev],lty=2)
  abline(v=reflev,lty=2)
  axis(2); axis(3,c(-1000,1000)); axis(4,c(-1000,1000))
  axis(1,at = 0:(nlev+1),c("",round(1+(as.numeric(colnames(levs))-1)/2,2),""))
  mtext(lab,font=2,line=0.3)
}

plot_sense=function(senseobj,speclab){
  np = length(senseobj)
  parnams = names(senseobj)
  ncol = ceiling(np^0.5)
  nrow = ceiling(np/ncol)
  par(mfrow=c(nrow,ncol),mai = c(0.5,0.3,0.3,0.05), omi=c(0.5,0.5,0.4,0.05))
  for(pp in 1:np) sensepanel(senseobj[[pp]],lab=parnams[pp])
  mtext("Multiplier on Parameter",1,outer=T,line=0.3)
  mtext("SSB/SSBMSY (2025)",2,outer=T,line=0.3)
  mtext(speclab,3,line=0.3,outer=T,font=2,cex=1.1)
}



