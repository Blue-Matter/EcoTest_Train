
library(MSEtool)

setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_Indicators.R")
source("99_make_MMP.R")
source("99_batching.R")

MOM = readRDS("./MOM/MOM_stitch_100sim.rds")
proyears = MOM@proyears

nbatch = 10

largedir = "C:/temp/Ecotest/batching/Independent_F"

totEffmat <<- readRDS("./Batch/totEffmat.rda")

sapply(1:nbatch, runbatch,vMPs = "Frand_MMP", largedir)



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