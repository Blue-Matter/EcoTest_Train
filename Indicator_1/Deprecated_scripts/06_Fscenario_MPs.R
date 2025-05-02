

library(dplyr)
setwd("C:/GitHub/Ecotest")
source("99_make_MMP.R")
source("99_plotting.R")

# Make a DataList for PM development

MSEforPPD = readRDS("MOM/MSEforPPD.rds")
DataList0 = MSEforPPD@PPD
ns = length(DataList0)
nf = length(DataList0[[1]])

DataList=list()
for(s in 1:ns){
  DataList[[s]]=list()
  for(f in 1:nf){
    DataList[[s]][[f]] = DataList0[[s]][[f]][[1]]
  }
}


# Correlation of % cover with historical biomass
# Correlation of CPUE with historical biomass
# 
MOM <- readRDS("./MOM/MOM_stitch_100sim.rds")
multiHist <- SimulateMOM(MOM, parallel = FALSE)
saveRDS(multiHist,"./MOM/multiHist_100sim.rds")
MMSE <- ProjectMOM(multiHist, MPs = "Frand_MMP", checkMPs = FALSE)
saveRDS(MMSE,"./MOM/MMSE_100sim.rds")
PPD = MMSE@PPD
saveRDS(PPD,"./MOM/PPD_100sim.rds")

jpeg("Figures/Projections/Fscenario_B.jpg",height=15,width=10,units='in',res=400)
  plotB(MMSE)
dev.off()
