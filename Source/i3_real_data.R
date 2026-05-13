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
  Mat_age = domat(Len_age, L50, Mat_slope=0.8, plot=plot)
  Sel_len =

  # Fleet specific
  for(ff in 1:nf){

    # Posterior predicted data
    Data = new("Data")
    Data@Cat  = matrix(Cat_f[[ff]],nrow=1)  # nsim, nyear
    Data@VInd = matrix(I_f[[ff]],  nrow=1)  # nsim, nyear
    Data@CAL = matrix(CAL_f[[ff]],nrow=1)  # dat@CAL[i,,] # nsim, nyear, n CAL
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

make_real_NN_input = function(MMSE){
  outs = list()
  Cat = MMSE@PPD[[1]][[1]][[1]]@Cat
  ny = ncol(Cat)
  Iind = matrix(c(1,ny),nrow=1)
  sno = 1
  outs[[1]] = proc_dat_LH_3(MMSE,Iind=Iind,sno=sno)
  outs[[2]] = proc_allFdat_F3(MMSE,Iind=Iind,sno=sno)


    for(fno in 1:nf){
      i=i+1
      if(!quiet)cat(".")
      outs[[i]] = proc_dat_F3(MMSE,Iind=Iind,sno = sno, fno=fno)
    }
  }

  for(j in 1:i){
    out = outs[[j]]
    #if(j <=ns) names(out) = paste0(names(out),"_",i)
    if(j==1) one_tab=out
    if(j>1) one_tab=cbind(one_tab,out)
  }
  ny = Iyr
  one_tab = cbind(one_tab,ny)


}

domat = function(Len_age,L50,Mat_slope=0.8,plot=T){
  nage = length(Len_age)
  aa = (1:nage)-0.5
  Len_mat =
    a50 = approx(Len_age,aa,L50)$y
  mat_age = 1/(1+exp(-Mat_slope*(Len_age-L50)))
  plot(aa,mat_age,ylab="Spawn. Frac.",xlab="Age")
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
