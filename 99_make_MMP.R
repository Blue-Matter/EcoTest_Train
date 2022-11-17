
library(MSEtool)

make_MMP <- function(Frel, do_sim_byc = TRUE, Ftarg_CV = 0, rel_log_log = FALSE, frac = 1) {
  
  rel_sum <- summary(Frel)
  RSD <- sapply(rel_sum, getElement, "sigma") # Residual standard deviation (normal space)
  RSD_LN <- sdconv(1, RSD) # Residual standard deviation (lognormal space)
  
  RSD_LN_targ <- sdconv(1, Ftarg_CV)
  force(Frel)
  force(rel_sum)
  force(do_sim_byc)
  force(rel_log_log)
  force(frac)
  
  MMP <- function(x, DataList, reps = 1, ...) {
    ny <- ncol(DataList[[1]][[1]]@Cat)
    nsim <- nrow(DataList[[1]][[1]]@Cat)
    
    np <- length(DataList)
    nf <- length(DataList[[1]])
    
    Fname <- MSEtool:::SIL(DataList, "Name") %>% matrix(nf) # nf x np matrix
    Sname <- substr(Fname[1, ], 1, 3)
    
    unisex <- table(Sname) == 1
    Sunique <- names(unisex)
    
    # Terminal F at age by fleet
    FM <- sapply(1:np, function(p) { # age x fleet x stock
      sapply(1:nf, function(f) DataList[[p]][[f]]@OM$qs[x] * DataList[[p]][[f]]@Misc[[x]]$FinF)
    }, simplify = "array")
    
    #HistE <- sapply(1:np, function(p) { 
    #  sapply(1:nf, function(f) DataList[[p]][[f]]@OM$qs[x] * DataList[[p]][[f]]@OM$FinF[x])
    #})
    
    # Fleet apical terminal F
    FMfleet <- apply(FM, c(2, 3), max)
    
    # Aggregate apical terminal F
    FMaggmax <- apply(FM, c(1, 3), sum) %>% apply(2, max) %>% structure(names = Sname)
    
    # Ratio of apical F between fleets and aggregate - identical apical age among fleets implies colSums = 1 
    relF <- t(t(FMfleet)/FMaggmax)
    
    targ <- c("SWO", "BET")
    byc <- c("BSH", "SMA", "WHM", "BUM")
    
    # Aggregate FMSY for target sp.
    FMSY_targ <- sapply(targ, function(sp) {
      ind <- which(sp == Sname)[1] # Assume females have higher catchability than males
      val <- DataList[[ind]][[1]]@Misc[[x]]$FMSY[ny + 1]
      return(val)
    })
    
    # F by fleet for target sp. (row 1 = longline, row 2 = other)
    Fout_targ <- lapply(targ, function(sp) {
      ind <- which(sp == Sname)
      val <- relF[, ind, drop = FALSE] * frac * FMSY_targ[sp] 
      val * rlnorm(length(val), -0.5 * RSD_LN_targ^2, RSD_LN_targ)
    }) %>% structure(names = targ)
    
    # F other target sp.
    # Calculate bycatch F longline from LL (use linear model)
    FLL_byc <- local({
      # Longline F corresponding to aggregate FMSY
      FLL_targ <- sapply(Fout_targ, function(xx) xx[1, 1])
      sapply(targ, function(xx) relF[1, which(xx == Sname)[1]])
      newdata <- FLL_targ %>% t() %>% as.data.frame()
      if(rel_log_log) {
        Fpred <- predict(Frel, newdata = log(newdata))[1, ] %>% exp()
      } else {
        Fpred <- predict(Frel, newdata = newdata)[1, ]
      }
      
      if(do_sim_byc) {
        dev <- rlnorm(length(byc), -0.5 * RSD_LN^2, RSD_LN)
        Fpred <- Fpred * dev
      }
      Fpred
    })
    
    RecList <- lapply(1:np, function(p) replicate(nf, new("Rec")))
    for(sp in Sunique) {
      pind <- which(sp == Sname) # Assign population to species
      
      # Target stocks fished at F = FMSY for both longline and other fleet at same ratio as in terminal year
      if(any(sp == targ)) {
        # Effort advice is the ratio of new F to terminal year F (FMSY_targ/FMaggmax)
        # Assuming effort = F (Data@OM$qs = 1)
        # Should account for catchability differences between sexes
        for(p in pind) { 
          for(f in 1:nf) { # Specify relative F 
            RecList[[p]][[f]]@Effort <- rep(Fout_targ[[sp]][f, p == pind]/FMfleet[f, p], reps)
          }
        }
      } else {
        
        # Bycatch longline F based on linear model, status quo for other
        for(p in pind) { 
          RecList[[p]][[1]]@Effort <- rep(FLL_byc[sp]/FMfleet[1, p], reps)
          RecList[[p]][[2]]@Effort <- rep(1, reps)
        }
      }
      
    }
    
    for(p in 1:np) { 
      for(f in 1:nf) { # Specify relative F 
        RecList[[p]][[f]]@Misc <- DataList[[p]][[f]]@Misc[[x]]
      }
    }
    return(RecList)
  }
  structure(MMP, class = "MMP")
}


NoFishing <- function(x, DataList, reps) {
  RecList <- lapply(1:10, function(p) {
    lapply(1:2, function(f) {
      NFref(x, DataList[[p]][[f]])
    })
  })
  return(RecList)
}
class(NoFishing) <- "MMP"


