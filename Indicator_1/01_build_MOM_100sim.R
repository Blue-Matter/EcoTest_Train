
library(openMSE)

#### Import SS3 assessments into multi-fleet, 2-sex operating models (MOM)
#### multiHist contains the historical reconstruction from the MOM

nsim <- 100

### BUM
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2018 BUM SS3\\SS_BASE_2018_v3"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BUM_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_BUM_100sim.rds")

### WHM
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2019 WHM SS3\\StockSynthesis\\Model_6"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_WHM_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_WHM_100sim.rds")


### SMA
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2019 SMA SS3\\run_1_try_09_0"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_SMA_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_SMA_100sim.rds")


### BSH
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\BSH\\Preliminary_Run_6_input"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
for(p in 1:length(MOM@Stocks)) {
  for(f in 1:length(MOM@Fleets[[1]])) MOM@cpars[[p]][[f]]$Fec_age <- 39 * MOM@cpars[[p]][[f]]$Mat_age
}

multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BSH_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_BSH_100sim.rds")

# SWO
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\ICCAT_SWO_Assessment\\NSWO_MSE_SS3_Base_v2"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_SWO_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_SWO_100sim.rds")

# BET
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2021 BET\\M20_h0.8_sigmaR0.4"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = nsim)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BET_100sim.rds")
saveRDS(multiHist, file = "MOM/multiHist_BET_100sim.rds")

#MOM@cpars[[1]][[1]]$V[1,1:4,1:10]
#MOM@cpars[[1]][[1]]$SLarray[1,,1:10]
