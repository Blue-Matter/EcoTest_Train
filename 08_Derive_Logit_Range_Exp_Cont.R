
setwd("C:/GitHub/Ecotest"); MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_1.rda")

setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest"); MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_21.rda")

source("99_Indicators.R")
source("99_Spatial_Indicators.R")
Catchdat = ICCATtoGEO(read.csv("G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped_fleet.csv"), FleetCode = "JPN")

nstocks = MMSE@nstocks
SCodes = sapply(MMSE@Snames,function(x)substr(x,1,3))
SSBrel = list()
for(ss in 1:nstocks) SSBrel[[ss]] = apply(MMSE@multiHist[[ss]]$Longline@TSdata$SBiomass,1:2,sum) / MMSE@RefPoint$ByYear$SSBMSY[,ss,1,1:MMSE@nyears]




par(mfrow=c(3,2),mai=c(0.5,0.5,0.1,0.1))
xs = c(1,2,4,5,7,9)

type = "Catch"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.5
type = "CR"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.2
type = "CR"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.8

for(x in xs) {getlogitmod(x,SSBrel,Catchdat,SCodes, type = type,
            ntype = ntype, itype = itype, ref_lev = ref_lev, Pyrs = seq(1970,2010,by=5), 
            ploty="rel", yrs = 1950:2013)}



getlogitmod(x,SSBrel,Catchdat,SCodes, type = type,
            ntype = ntype, itype = itype, ref_lev = ref_lev, Pyrs = seq(1970,2010,by=5), 
            ploty="all", yrs = 1950:2013)
  
type = "CR"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.5



lapply(1:nstocks,getlogitmod,SSBrel,Catchdat,SCodes,ploty=T)


