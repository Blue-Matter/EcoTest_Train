
library(MSEtool)

setwd("C:/GitHub/Ecotest")
source("99_Indicators.R")
MMSE = readRDS("./MOM/MMSE_100sim.rds")

jpeg("C:/temp/Ecotest/dump/Fproj_base.jpg",res=600,height=12,width=8,units="in"); plotF(MMSE); dev.off()
PPD = readRDS("./MOM/PPD_100sim.rds")
MOM = readRDS("./MOM/MOM_stitch_100sim.rds")

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