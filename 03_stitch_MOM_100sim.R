

library(MSEtool)
library(dplyr)

source("MOM_fn.R")

byc <- c("BSH", "SMA", "WHM", "BUM")
targ <- c("BET", "SWO")

multiHist_byc <- lapply(paste0('MOM/multiHist_', byc, '_100sim.rds'), readRDS)
multiHist_targ <- lapply(paste0('MOM/multiHist_', targ, '_100sim.rds'), readRDS)

MOM_byc <- lapply(paste0('MOM/MOM_', byc, '_100sim.rds'), readRDS)
MOM_targ <- lapply(paste0('MOM/MOM_', targ, '_100sim.rds'), readRDS)

LL_byc <- list(1:9, c(1:12)[-c(8, 10, 12)], 2, 2) %>% structure(names = byc)
LL_targ <- list(10:18, 1:10) %>% structure(names = targ)

M2_byc <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_byc, y = LL_byc)
M2_targ <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_targ, y = LL_targ)

MOM <- local({
  args <- c(M2_targ, M2_byc) %>% structure(names = c(targ, byc))
  do.call("MOM_stitch", args)
})


lapply(args,function(x)x@cpars[[1]][[1]]$Mat_age[1,,1])
lapply(MOM@cpars,function(x)print(x[[1]]$Mat_age[1,,1]))
lapply(args,function(x)x@cpars[[1]][[1]]$Len_age[1,,1])
lapply(MOM@cpars,function(x)print(x[[1]]$Len_age[1,,1]))

saveRDS(MOM, file = "MOM/MOM_stitch_100sim.rds")

MOM= readRDS("MOM/MOM_stitch_100sim.rds")
MOM2 = MOM_simplify(MOM)

saveRDS(MOM2, file = "MOM/MOM_stitch_100sim_simplified.rds")



#MOM <- readRDS(file = "MOM/MOM_stitch_100sim.rds")

setup(8)
multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = TRUE)
sfStop()
saveRDS(multiHist, file = "MOM/multiHist_stitch_100sim.rds")
#multiHist <- readRDS(file = "MOM/multiHist_stitch_100sim.rds")

