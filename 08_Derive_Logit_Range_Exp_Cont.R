
# ---- Prerequisites --------------------------------------------------------------------------------------

setwd("C:/GitHub/Ecotest"); MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_1.rda")

setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest"); MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_1.rda")

source("99_Indicators.R")
source("99_Spatial_Indicators.R")
Catchdat = ICCATtoGEO(read.csv("G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped_fleet.csv"), FleetCode = "JPN")

nstocks = MMSE@nstocks
SCodes = sapply(MMSE@Snames,function(x)substr(x,1,3))
SSBrel = list()
for(ss in 1:nstocks) SSBrel[[ss]] = apply(MMSE@multiHist[[ss]]$Longline@TSdata$SBiomass,1:2,sum) / MMSE@RefPoint$ByYear$SSBMSY[,ss,1,1:MMSE@nyears]



# ----------- Exploration ---------------------------------------------------------------------------------

par(mfrow=c(3,2),mai=c(0.5,0.5,0.1,0.1))

# catch data    normalized to ann frac     percentage cells that make up   top 50%
type = "Catch"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.5  # all cells
type = "Catch"; ntype = "annual fraction"; itype = "positive percentage"; ref_lev = 0.5 # only those positive in each year

# catch data no normalization  relative to abs mean    values above 50% abs mean
type = "Catch"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.5  # all cells
type = "Catch"; ntype = "none"; itype = "positive abs_mean"; ref_lev = 0.5 # only those positive in each year

# CPUE data    normalized to ann frac     percentage cells that make up   top 50%
type = "CR"; ntype = "annual fraction"; itype = "percentage"; ref_lev = 0.5  # all cells
type = "CR"; ntype = "annual fraction"; itype = "positive percentage"; ref_lev = 0.5 # only those positive in each year

# CPUE data no normalization  relative to abs mean    values above 50% abs mean
type = "CR"; ntype = "none"; itype = "abs_mean"; ref_lev = 0.5  # all cells
type = "CR"; ntype = "none"; itype = "positive abs_mean"; ref_lev = 0.5 # only those positive in each year

for(x in 1:6) {getlogitmod(x,SSBrel,Catchdat,SCodes, type = type,
            ntype = ntype, itype = itype, ref_lev = ref_lev, Pyrs = seq(1970,2010,by=5), 
            ploty="rel", yrs = 1950:2013)}

getlogitmod(x=6,SSBrel,Catchdat,SCodes, type = type,
            ntype = ntype, itype = itype, ref_lev = ref_lev, Pyrs = seq(1970,2010,by=5), 
            ploty="all", yrs = 1950:2013)


# ---------- Correlations / Hypothesized relationships ---------------------------------------------------  

configs = lapply(1:6,function(x)new('list'))

# BET
configs[[1]][[1]] = data.frame(type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)

# SWO
configs[[2]][[1]] = data.frame(type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)
configs[[2]][[2]] = data.frame(type = "CR", ntype = "none", itype = "positive abs_mean", ref_lev = 0.5) # only those positive in each year
configs[[2]][[3]] = data.frame(type = "CR", ntype = "none", itype = "abs_mean", ref_lev = 0.5) 

# BSH
configs[[3]][[1]] = data.frame(type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)

# SMA
configs[[4]][[1]] = data.frame(type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)

# WHM
configs[[5]][[1]] = data.frame(type = "CR", ntype = "none", itype = "abs_mean", ref_lev = 0.5) # only those positive in each year
configs[[5]][[2]] = data.frame(type = "Catch", ntype = "none", itype = "abs_mean", ref_lev = 0.5)  # all cells

# BUM
configs[[6]][[1]] = data.frame(type = "CR", ntype = "none", itype = "abs_mean", ref_lev = 0.5) 
configs[[6]][[2]] = data.frame(type = "CR", ntype = "none", itype = "positive abs_mean", ref_lev = 0.5) # only those positive in each year

# plots the first for each
par(mfrow=c(3,2),mai=c(0.5,0.5,0.1,0.1))
for(ss in 1:6)  getlogitmod(x=ss,SSBrel,Catchdat,SCodes, type = configs[[ss]][[1]]$type,ntype = configs[[ss]][[1]]$ntype, itype = configs[[ss]][[1]]$itype, ref_lev = configs[[ss]][[1]]$ref_lev, Pyrs = seq(1970,2010,by=5), ploty="rel", yrs = 1950:2013)

# plots all 
np = sum(sapply(configs,function(x)length(x)))
par(mfrow=c(ceiling(np/2),2),mai=c(0.5,0.5,0.1,0.1))
for(ss in 1:6)for(mm in 1:length(configs[[ss]]))  getlogitmod(x=ss,SSBrel,Catchdat,SCodes, type = configs[[ss]][[mm]]$type,ntype = configs[[ss]][[mm]]$ntype, itype = configs[[ss]][[mm]]$itype, ref_lev = configs[[ss]][[mm]]$ref_lev, Pyrs = seq(1970,2010,by=5), ploty="rel", yrs = 1950:2013)


#  -------- Models -------------------------------------------

spat_mods = configs = new('list')

# reusing the 'no relationship' fit of SMA for BSH

# BET
configs[[1]] = data.frame(ss = 1, type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)
configs[[2]] = data.frame(ss = 2, type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)
configs[[3]] = data.frame(ss = 4, type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)
configs[[4]] = data.frame(ss = 4, type = "Catch", ntype = "annual fraction", itype = "positive percentage", ref_lev = 0.5)
configs[[5]] = data.frame(ss = 5, type = "CR", ntype = "none", itype = "abs_mean", ref_lev = 0.5) 
configs[[6]] = data.frame(ss = 6, type = "CR", ntype = "none", itype = "abs_mean", ref_lev = 0.5) 

par(mfrow=c(3,2),mai=c(0.5,0.5,0.1,0.1))
for(ss in 1:6)spat_mods[[ss]] = getlogitmod(x=configs[[ss]]$ss,SSBrel,Catchdat,SCodes, type = configs[[ss]]$type,ntype = configs[[ss]]$ntype, itype = configs[[ss]]$itype, ref_lev = configs[[ss]]$ref_lev, Pyrs = seq(1970,2010,by=5), ploty="rel", yrs = 1950:2013)

saveRDS(spat_mods,file="Spatial Models/spat_mods.rds")




jpeg("Figures/Presentation 1 April 24/spatial_model_explanatory_BET.jpg",res=600,width=8,height=7,units="in")
  ss=1
  getlogitmod(x=configs[[ss]]$ss,SSBrel,Catchdat,SCodes, type = configs[[ss]]$type,ntype = configs[[ss]]$ntype, itype = configs[[ss]]$itype, ref_lev = configs[[ss]]$ref_lev, Pyrs = seq(1970,2010,by=5), ploty="all", yrs = 1950:2013)
dev.off()


jpeg("Figures/Presentation 1 April 24/spatial_model_explanatory_WHM.jpg",res=600,width=8,height=7,units="in")
  ss=6
  getlogitmod(x=configs[[ss]]$ss,SSBrel,Catchdat,SCodes, type = configs[[ss]]$type,ntype = configs[[ss]]$ntype, itype = configs[[ss]]$itype, ref_lev = configs[[ss]]$ref_lev, Pyrs = seq(1970,2010,by=5), ploty="all", yrs = 1950:2013,pos="topright")
dev.off()


