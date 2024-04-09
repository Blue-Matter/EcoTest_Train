
setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_Indicators.R")
source("99_Spatial_Indicators.R")

MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE.rda")
Catchdat = ICCATtoGEO(read.csv("G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped.csv"))
#Catchdat = ICCATtoGEO(read.csv("G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped.csv"),type="1x1")

nstocks = MMSE@nstocks
SCodes = sapply(MMSE@Snames,function(x)substr(x,1,3))
SSBrel = list()
for(ss in 1:nstocks) SSBrel[[ss]] = apply(MMSE@multiHist[[ss]]$Longline@TSdata$SBiomass,1:2,sum) / MMSE@RefPoint$ByYear$SSBMSY[,ss,1,1:MMSE@nyears]


type = "Catch"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.5


type = "CR"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.2

type = "CR"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.8


mod = getlogitmod(x=x,SSBrel,Catchdat,SCodes, type = type,
            ntype = ntype, itype = itype, ref_lev = ref_lev, Pyrs = seq(1970,2010,by=5), 
            ploty=T, yrs = 1950:2013)

  
type = "CR"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.5



lapply(1:nstocks,getlogitmod,SSBrel,Catchdat,SCodes,ploty=T)


