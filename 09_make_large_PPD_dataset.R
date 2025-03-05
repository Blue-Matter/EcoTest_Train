
library(openMSE)
library(MSEtool)
library(dplyr)
library(mvtnorm)
library(parallel)


setwd("C:/GitHub/Ecotest")

source("99_Indicators.R")
source("99_make_MMP.R")
source("99_batching.R")
source("99_MOM_fixes.R")

MOM = readRDS("./MOM/MOM_latest.rds")

largedir = "C:/temp/Ecotest/batching/Dependent_F" # "C:/temp/Ecotest/batching/Independent_F"
totEffmat <<- readRDS("./Batch/totEffmat_cor.rda") 
# totEffmat <<- readRDS("./Batch/totEffmat.rda")



sfInit(cpus=parallel::detectCores(),parallel=T) 
sfLibrary(MSEtool); sfLibrary(mvtnorm)
sfExport("overwritePE"); sfExport("totEffmat"); sfExport("Frand_MMP"); sfExport("Ftv_MMP");
sfExport("add_stochasticity"); sfExport("trim_MMSE"); sfExport("stoch_SLarray")
sfExport("makeCAL"); sfExport("remakeCAL")

todosims = gettodosims(largedir,maxsim=1000)
sfSapply(todosims, runbatch, MOM=MOM, MPs = "Ftv_MMP", largedir)
#sfSapply(todosims, runbatch, MOM=MOM, MPs = "Frand_MMP", largedir)

#sfSapply(1:20, runbatch, MOM=MOM, MPs = "Ftv_MMP", largedir)


# === End of script ==============================================================




# vvvvvvvv Code for debugging vvvvvvvvVVVVVV

runbatch(1,MOM=MOM,MPs = "Frand_MMP", largedir=largedir)

dim(MOM@cpars[[3]][[2]]$Fdisc_array2)

lapply(MOM@cpars,function(x)dim(x[[1]]$SLarray))

lapply(MOM@cpars,function(x)length(x[[2]]$CAL_binsmid))
lapply(MOM@cpars,function(x)length(x[[2]]$CAL_bins))

test = SampleStockPars(MOM@Stocks[[3]],nsim=100,nyears=MOM@Fleets[[1]][[1]]@nyears,cpars=MOM@cpars[[1]][[1]])
test = SampleFleetPars(MOM@Fleets[[1]][[1]],Stock=MOM@Stocks[[1]],nsim=100,nyears=MOM@Fleets[[1]][[1]]@nyears,cpars=MOM@cpars[[1]][[1]])



# vvvvvvvv more disused code vvvvvvvvvvvvvv

nyears = MOM@Fleets[[1]]$Longline@nyears
proyears = dim(PPD[[1]][[1]][[1]]@Cat)[2] - nyears
CurrentYr = MOM@Fleets[[1]]$Longline@CurrentYr

nf=MMSE@nfleets
ns=MMSE@nstocks

Fname <- sapply(PPD,function(x,sloty)MSEtool:::SIL(x,sloty), sloty="Name") %>% matrix(nf) # nf x np matrix
Sname <- substr(Fname[1, ], 1, 3)

outs = list()
for(ss in 1:ns)outs[[ss]] = proc_dat(MMSE,sno=ss)


sind = match(unique(Sname),Sname) # first is the female for which maturity indicator works

for(ss in sind){
  jpeg(paste0("Figures/Indicators/Cor_Plots/Cor_Plot_",Sname[ss],".jpg"),height=11,width=11,units='in',res=600)
    corplot(outs[[ss]],maxn=20,lab = Sname[ss])
  dev.off()
}



