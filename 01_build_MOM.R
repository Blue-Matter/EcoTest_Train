
library(MSEtool)

#### Import SS3 assessments into multi-fleet, 2-sex operating models (MOM)
#### multiHist contains the historical reconstruction from the MOM


### BUM
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2018 BUM SS3\\SS_BASE_2018_v3"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BUM.rds")
saveRDS(multiHist, file = "MOM/multiHist_BUM.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_BUM", dir = getwd())

### WHM
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2019 WHM SS3\\StockSynthesis\\Model_6"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_WHM.rds")
saveRDS(multiHist, file = "MOM/multiHist_WHM.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_WHM", dir = getwd())


### SMA
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2019 SMA SS3\\run_1_try_09_0"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_SMA.rds")
saveRDS(multiHist, file = "MOM/multiHist_SMA.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_SMA", dir = getwd())



### BSH
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\BSH\\Preliminary_Run_6_input"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
for(p in 1:length(MOM@Stocks)) {
  for(f in 1:length(MOM@Fleets[[1]])) MOM@cpars[[p]][[f]]$Fec_age <- 39 * MOM@cpars[[p]][[f]]$Mat_age
}

multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BSH.rds")
saveRDS(multiHist, file = "MOM/multiHist_BSH.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_BSH", dir = getwd())

# SWO
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\ICCAT_SWO_Assessment\\NSWO_MSE_SS3_Base_v2"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_SWO.rds")
saveRDS(multiHist, file = "MOM/multiHist_SWO.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_SWO", dir = getwd())

# BET
ssdir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Assessments\\2021 BET\\M20_h0.8_sigmaR0.4"
replist <- r4ss::SS_output(ssdir)
MOM <- SS2MOM(replist, nsim = 2)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(MOM, file = "MOM/MOM_BET.rds")
saveRDS(multiHist, file = "MOM/multiHist_BET.rds")

plot_SS2MOM(multiHist, replist, filename = "SS2MOM_BET", dir = getwd())


