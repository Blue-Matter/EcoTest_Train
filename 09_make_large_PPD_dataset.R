
library(MSEtool)
library(dplyr)
library(mvtnorm)

setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_Indicators.R")
source("99_make_MMP.R")
source("99_batching.R")
source("99_MOM_fixes.R")

MOM0 = readRDS("./MOM/MOM_stitch_100sim_simplified.rds")
MOM1 = fix_selectivity_1(MOM0) #  MOM1@cpars[[1]][[1]]$V[1,1:4,1:10] # first 10 years of sim 1
MOM2 = add_stochasticity(MOM1) # adds stochasticity in M, K, Linf, and stock depletion
MOM = add_SL_array(MOM2) # lapply(MOM2@cpars,function(x)x[[1]]$SLarray[1,1:25,1:5])
saveRDS(MOM,"./MOM/MOM_stoch.rds")

MOM = readRDS("./MOM/MOM_stoch.rds")

largedir = "C:/temp/Ecotest/batching/Independent_F"
totEffmat <<- readRDS("./Batch/totEffmat.rda")

sfInit(cpus=parallel::detectCores()/2,parallel=T)
sfLibrary(MSEtool)
sfExport("overwritePE");sfExport("totEffmat"); sfExport("Frand_MMP")

todosims = gettodosims(largedir)
sfSapply(todosims, runbatch, MOM=MOM, MPs = "Frand_MMP", largedir)







# vvvvvvvv disused code vvvvvvvvvvvvvv

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