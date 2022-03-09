

library(MSEtool)
library(dplyr)

source("MOM_fn.R")

byc <- c("BSH", "SMA", "WHM", "BUM")
targ <- c("BET", "SWO")

multiHist_byc <- lapply(paste0('MOM/multiHist_', byc, '.rds'), readRDS)
multiHist_targ <- lapply(paste0('MOM/multiHist_', targ, '.rds'), readRDS)

MOM_byc <- lapply(paste0('MOM/MOM_', byc, '.rds'), readRDS)
MOM_targ <- lapply(paste0('MOM/MOM_', targ, '.rds'), readRDS)

LL_byc <- list(1:9, c(1:12)[-c(8, 10, 12)], 2, 2) %>% structure(names = byc)
LL_targ <- list(10:18, 1:11) %>% structure(names = targ)

M2_byc <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_byc, y = LL_byc)
M2_targ <- Map(function(x, y) aggregate_fleet(x = x, LL = y), x = MOM_targ, y = LL_targ)

MOM <- local({
  args <- c(M2_targ, M2_byc) %>% structure(names = c(targ, byc))
  do.call("MOM_stitch", args)
})
saveRDS(MOM, file = "MOM/MOM_stitch.rds")
#MOM <- readRDS(file = "MOM/MOM_stitch.rds")

multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
saveRDS(multiHist, file = "MOM/multiHist_stitch.rds")
#multiHist <- readRDS(file = "MOM/multiHist_stitch.rds")


# Compare SSB between multiHist_stitch and individual multiHist
stitch_SSB <- local({
  SSBs <- sapply(c(1, 2, 4, 5, 7, 9), function(p) multiHist[[p]][[1]]@TSdata$SBiomass[1, , ] %>% rowSums()) %>%
    structure(dimnames = list(Year = MOM@Fleets[[1]][[1]]@CurrentYr - 64:1 + 1, 
                              Stock = c(targ, byc))) %>%
    reshape2::melt() %>% mutate(Type = "Stitch")
  
  SSBi <- lapply(1:6, function(i) {
    m <- c(multiHist_targ, multiHist_byc)[[i]]
    SSB <- m[[1]][[1]]@TSdata$SBiomass[1, , ] %>% rowSums()
    data.frame(Year = m[[1]][[1]]@Data@OM$CurrentYr[1] - length(SSB):1 + 1, 
               Stock = c(targ, byc)[i], value = SSB, Type = "Indiv.")
  }) %>% bind_rows()
  
  ggplot(rbind(SSBs, SSBi), aes(Year, value, colour = Type)) + 
    facet_wrap(~Stock, scales = "free_y") + geom_line() + expand_limits(y = 0)
})

