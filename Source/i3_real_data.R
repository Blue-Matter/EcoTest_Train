

make_MMSE_I3 = function(L50, Linf, M, K, Mat_age, Len_age,
                        Cat_f, CAL_f, I_f, L5_f, LFS_f, Vmaxlen_f,
                        CAL_mids){

  nyears = length(Cat_f[[1]])
  nages = Mat_age
  nf = length(Cat_f)
  MMSE = new('MMSE')
  MMSE@nsim = 1
  MMSE@nyears = nyears-2
  MMSE@proyears = 2

  PPD = multiHist = OM = list() # hierarchical lists
  PPD[[1]] =  multiHist[[1]] = OM[[1]] = list() # one stock

  # Fleet specific
  for(ff in 1:nf){

    # Posterior predicted data
    Data = new("Data")
    Data@Cat  = matrix(Cat_f[[ff]],nrow=1)
    Data@VInd = matrix(I_f[[ff]],  nrow=1)
    Data@CAL = matrix(CAL_f[[ff]],nrow=1)  # dat@CAL[i,,] # nsim, nyear, n CAL
    Data@CAL_mids = CAL_mids
    MMSE@PPD[[sno]][[fno]] = list()
    MMSE@PPD[[sno]][[fno]][[1]] = Data  # only one MP

    # Sampled parameters
    multiHist[[1]][[ff]]=new("Hist")
    multiHist[[1]][[ff]]@SampPars = list()
    multiHist[[1]][[fno]]@SampPars$Fleet$L5_y = matrix(L5_f[[fno]], nrow=1)
    multiHist[[1]][[fno]]@SampPars$Fleet$LFS_y = matrix(LFS_f[[fno]], nrow=1)
    multiHist[[1]][[fno]]@SampPars$Fleet$Vmaxlen_y = matrix(Vmaxlen_f[[fno]], nrow=1)
    multiHist[[1]][[fno]]@AtAge$Maturity = array(Mat_age, c(1, nages, 1))
    multiHist[[1]][[fno]]@AtAge$Length = array(Len_age, c(1, nages, 1))

  }

  OM[[1]][[1]] = list()
  OM[[1]][[1]]$L50 = L50
  OM[[1]][[1]]$Linf = Linf
  OM[[1]][[1]]$M = M
  OM[[1]][[1]]$K = K

  MMSE@multiHist = multiHist
  MMSE@OM = OM

  MMSE

}


concat = function(MMSE, Ind = 2, sno = 1, plotsmooth = T, nint = 40){
  outs = list()
  out[[1]] = proc_allFdat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  out[[2]] = proc_dat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  out[[3]] = proc_dat_F3(MMSE, Iind = Iind, sno = sno, plotsmooth = plotsmooth, nint = nint)
  one_tab = cbind(out[[1]], out[[2]], out[[3]])
  one_tab
}
