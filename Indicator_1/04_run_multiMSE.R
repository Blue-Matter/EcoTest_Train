
library(MSEtool)

source("99_make_MMP.R")

MOM <- readRDS("MOM/MOM_stitch_100sim.rds")
MOM@interval <- 1

#multiHist <- multiMSE(MOM, Hist = TRUE, checkMPs = FALSE, parallel = FALSE)
#multiHist <- readRDS("MOM/multiHist_stitch_100sim.rds")

# Add FMSY to Data@Misc
np <- length(MOM@Stocks)
nf <- length(MOM@Fleets[[1]])

#FMSY <- lapply(1:np, function(p) {
#  lapply(1:nf, function(f) multiHist[[p]][[f]]@Ref$ByYear$FMSY)
#})

#FinF <- lapply(1:np, function(p) {
#  lapply(1:nf, function(f) multiHist[[p]][[f]]@AtAge$F.Mortality[, , MOM@Fleets[[p]][[f]]@nyears, 1])
#})
#rm(multiHist)

Frel <- readRDS("Frel/F_lm.rds")
FMSY_rel_linear <- make_MMP(Frel, do_sim_byc = FALSE, Ftarg_CV = 0.05, rel_log_log = FALSE)

Frel_power <- readRDS("Frel/F_lm_power.rds")
FMSY_rel_power <- make_MMP(Frel_power, do_sim_byc = FALSE, Ftarg_CV = 0.05, rel_log_log = TRUE)

# Run MSE in batches to reduce memory
batch <- list(1:25, 26:50, 51:75, 76:100)

for(i in 2:length(batch)) {
  MOM_batch <- SubCpars(MOM, batch[[i]])
  message("Starting batch ", i)
  multiHist <- SimulateMOM(MOM_batch, parallel = FALSE)
  
  for(p in 1:np) {
    for(f in 1:nf) {
      multiHist[[p]][[f]]@Data@Misc <- lapply(1:MOM_batch@nsim, function(x) {
        out <- list()
        if(f == 1) {
          out$FMSY <- multiHist[[p]][[f]]@Ref$ByYear$FMSY[x, ]
        }
        out$FinF <- multiHist[[p]][[f]]@AtAge$F.Mortality[x, , MOM@Fleets[[p]][[f]]@nyears, 1]
        return(out)
      })
    }
  }
  
  MMSE <- ProjectMOM(multiHist, MPs = c("FMSY_rel_linear", "FMSY_rel_power", "NoFishing"), checkMPs = FALSE)
  saveRDS(MMSE, file = paste0("MSE/MMSE_100sim_batch_", i, ".rds"))
}

