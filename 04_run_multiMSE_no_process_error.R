
library(MSEtool)

source("99_make_MMP.R")


MOM <- readRDS("MOM/MOM_stitch.rds")

# Add FMSY to Data@Misc
np <- length(MOM@Stocks)
nf <- length(MOM@Fleets[[1]])

Histnpe <- local({
  maxage <- MOM@Stocks[[1]]@maxage
  nyears <- MOM@Fleets[[1]][[1]]@nyears
  proyears <- MOM@proyears
  
  for(p in 1:np) {
    for(f in 1:nf) {
      MOM@cpars[[p]][[f]]$Perr_y[, maxage + nyears + 1:proyears] <- 1
    }
  }
  SimulateMOM(MOM, parallel = FALSE, silent = TRUE)
})
Histnpe[[1]][[1]]@Misc$MOM@interval <- 1
saveRDS(Histnpe, file = "MOM/multiHist_npe.rds")

for(p in 1:np) {
  for(f in 1:nf) {
    
    Histnpe[[p]][[f]]@Data@Misc <- lapply(1:MOM@nsim, function(x) {
      out <- list()
      if(f == 1) {
        out$FMSY <- Histnpe[[p]][[f]]@Ref$ByYear$FMSY[x, ]
      }
      out$FinF <- Histnpe[[p]][[f]]@AtAge$F.Mortality[x, , MOM@Fleets[[p]][[f]]@nyears, 1]
      return(out)
    })
    
  }
}

Frel <- readRDS("Frel/F_lm.rds")
FMSY_rel_linear <- make_MMP(Frel, do_sim_byc = FALSE, Ftarg_CV = 0, rel_log_log = FALSE)

Frel_power <- readRDS("Frel/F_lm_power.rds")
FMSY_rel_power <- make_MMP(Frel_power, do_sim_byc = FALSE, Ftarg_CV = 0, rel_log_log = TRUE)

# No process error
npe <- ProjectMOM(Histnpe, MPs = c("FMSY_rel_linear", "FMSY_rel_power", "NoFishing"), checkMPs = FALSE)
saveRDS(npe, file = "MSE/MMSE_npe.rds")

