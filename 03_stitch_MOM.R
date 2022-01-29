

library(MSEtool)
library(dplyr)

source("MOM_fn.R")

byc <- c("BSH", "SMA", "WHM", "BUM")
targ <- c("BET", "SWO")

multiHist_byc <- lapply(paste0('multiHist_', byc, '.rds'), readRDS)
multiHist_targ <- lapply(paste0('multiHist_', targ, '.rds'), readRDS)

MOM_byc <- lapply(paste0('MOM_', byc, '.rds'), readRDS)
MOM_targ <- lapply(paste0('MOM_', targ, '.rds'), readRDS)

LL_byc <- list(1:9, c(1:12)[-c(8, 10, 12)], 2, 2) %>% structure(names = byc)
LL_targ <- list(10:18, 1:11) %>% structure(names = targ)

M2_byc <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_byc, y = LL_byc)
M2_targ <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_targ, y = LL_targ)


MOM <- local({
  args <- c(M2_targ, M2_byc) %>% structure(names = c(targ, byc))
  do.call("MOM_stitch", args)
})
saveRDS(MOM, file = "MOM_stitch.rds")

multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)

