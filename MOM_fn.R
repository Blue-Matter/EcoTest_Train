
#' Need identical fleet structure for each MOM
MOM_stitch <- function(..., silent = FALSE) {
  
  dots <- list(...) 
  
  # Check maxage within MOM
  if (!silent) message("\nChecking maxage within individual MOMs...")

  dots <- lapply(dots, MOM_std_maxage, silent = TRUE)
  
  # Check maxage among MOM
  maxage <- vapply(dots, function(MOM) {
    vapply(MOM@Stocks, "slot", numeric(1), "maxage") %>% unique()
  }, numeric(1))
  
  if (!silent) {
    message("\nChecking maxage among MOMs...")
    message("Maximum age will be: ", max(maxage))
  }
  
  MOM_list <- lapply(dots, MOM_std_maxage, silent = silent, maxage_new = max(maxage))
  
  # Align years among MOM
  if (!silent) {
    message("Checking historical and projection years among MOMs...")
  }
  
  nyears <- vapply(MOM_list, function(MOM) MOM@Fleets[[1]][[1]]@nyears, numeric(1))
  CurrentYr <- vapply(MOM_list, function(MOM) MOM@Fleets[[1]][[1]]@CurrentYr, numeric(1))
  StartYr <- CurrentYr - nyears + 1
  proyears <- vapply(MOM_list, slot, numeric(1), "proyears")
  
  StartYr_new <- min(StartYr)
  CurrentYr_new <- min(CurrentYr)
  nyears_new <- length(StartYr_new:CurrentYr_new)
  proyears_new <- min(proyears)
  
  #if (StartYr != StartYr_new || CurrentYr != CurrentYr_new || any(proyears_new != proyears)) {
    
    if (!silent) {
      message("Historical period of MOMs will be adjusted to years: ", StartYr_new, " - ", CurrentYr_new)
      message("Number of projection years: ", proyears_new)
    }
    
    MOM_list <- lapply(MOM_list, MOM_std_years, silent = silent, CurrentYr_new = CurrentYr_new,
                       StartYr_new = StartYr_new, proyears_new = proyears_new)
  #}
  
  # Checking simulations among MOM
  if (!silent) {
    message("Checking number of simulations among MOMs...")
  }
  
  nsim <- vapply(MOM_list, slot, numeric(1), "nsim")
  
  if (length(unique(nsim)) > 1) {
    stop("Number of simulations in MOMs do not match")
  }
  
  # Checking fleets among MOM
  if (!silent) {
    message("Checking number of fleets among MOMs...")
  }
  
  nf <- sapply(MOM_list, function(x) length(x@Fleets[[1]]))
  
  if (length(unique(nf)) > 1) {
    stop("Number of fleets in MOMs do not match")
  }
  
  # Stitch MOM
  if (!silent) {
    message("Stitching together MOMs...")
  }
  
  MOM <- suppressMessages(new("MOM")) 
  MOM@nsim <- unique(nsim) 
  MOM@proyears <- proyears_new 
  MOM@reps <- 1
  MOM@maxF <- 5
  MOM@pstar <- 0.5
  MOM@interval <- 1
  MOM@Name <-  "MOM created by MOM_stitch"
 
  if(is.null(names(MOM_list))) names(MOM_list) <- paste0("MOM_", 1:length(MOM_list))
  Sname <- lapply(1:length(MOM_list), function(m) {
    if(length(MOM_list[[m]]@Stocks) > 1) {
      paste0(names(MOM_list)[m], ".", names(MOM_list[[m]]@Stocks))
    } else {
      names(MOM_list)[m]
    }
  })

  # SexPars
  MOM@Stocks <- do.call(c, lapply(MOM_list, slot, "Stocks")) %>% structure(names = do.call(c, Sname))
  for(p in 1:length(MOM@Stocks)) MOM@Stocks[[p]]@Name <- names(MOM@Stocks)[p]
  MOM@Fleets <- do.call(c, lapply(MOM_list, slot, "Fleets")) %>% structure(names = do.call(c, Sname))
  MOM@Obs <- do.call(c, lapply(MOM_list, slot, "Obs"))
  MOM@Imps <- do.call(c, lapply(MOM_list, slot, "Imps"))
  MOM@CatchFrac <- do.call(c, lapply(MOM_list, slot, "CatchFrac")) %>% structure(names = do.call(c, Sname))
  MOM@cpars <- do.call(c, lapply(MOM_list, slot, "cpars")) %>% structure(names = do.call(c, Sname))
  
  SSBfrom_exist <- sapply(MOM_list, function(x) length(x@SexPars$SSBfrom))
  
  if(any(SSBfrom_exist)) {
    
    MOM@SexPars$SSBfrom <- local({
      np <- sapply(MOM_list, function(x) length(x@Stocks))
      SSBfrom <- matrix(0, sum(np), sum(np))
      
      for(i in 1:length(np)) {
        if(i == 1) {
          ind <- 1:np[i]
        } else {
          ind <- (cumsum(np)[i-1]+1):cumsum(np)[i]
        }
        if(SSBfrom_exist[i]) {
          SSBfrom[ind, ind] <- MOM_list[[i]]@SexPars$SSBfrom
        } else {
          SSBfrom[ind, ind] <- diag(np[i])
        }
        
      }
      structure(SSBfrom, dimnames = list(names(MOM@Stocks), names(MOM@Stocks)))
    })
    
  }
  
  
  if(!silent) message("Complete.")
  
  return(MOM)
}

# For stocks with multiple maxages
#' @importFrom abind abind
MOM_std_maxage <- function(MOM, silent = FALSE, maxage_new) {
  
  maxage <- vapply(MOM@Stocks, "slot", numeric(1), "maxage")
  if (missing(maxage_new)) {
    maxage_new <- max(maxage)
  } else if(any(maxage_new < maxage)) {
    stop("maxage_new must be greater than: ", max(maxage))
  }
  
  if (all(maxage == maxage_new)) {
    return(MOM)
  }
  
  MOM2 <- MOM
  
  n_age_new <- maxage_new + 1
  
  np <- length(MOM2@Stocks)
  nf <- length(MOM2@Fleets[[1]])
  
  Sname <- MOM2@Stocks %>% names()
  if (is.null(Sname)) Sname <- 1:np
  
  Fname <- MOM2@Fleets[[1]] %>% names()
  if (is.null(Fname)) Fname <- 1:nf
  
  for (p in 1:np) {
    MOM2@Stocks[[p]]@maxage <- maxage_new
    for (f in 1:nf) {
      if (length(MOM@cpars[[p]][[f]])) {
        if (!silent) message("Updating cpars for Stock: ", Sname[p], ", Fleet: ", Fname[f])
        cpars <- MOM@cpars[[p]][[f]]
        if(!is.null(cpars$plusgroup) && !cpars$plusgroup && is.null(cpars$M_ageArray)) {
          stop("cpars$plusgroup is FALSE. Provide cpars$M_ageArray")
        }
        MOM2@cpars[[p]][[f]] <- lapply(names(cpars), cpars_extend_age, cpars = cpars, n_age_new = n_age_new,
                                       n_age = maxage[p] + 1) %>% structure(names = names(cpars))
      }
    }
  }
  
  return(MOM2)
}




MOM_std_years <- function(MOM, silent = FALSE, CurrentYr_new, StartYr_new, proyears_new) {
  
  nyears <- MOM@Fleets[[1]][[1]]@nyears
  CurrentYr <- MOM@Fleets[[1]][[1]]@CurrentYr
  StartYr <- CurrentYr - nyears + 1
  proyears <- MOM@proyears
  HistYr <- seq(StartYr, CurrentYr)
  
  HistYr_new <- seq(StartYr_new, CurrentYr_new)
  nyears_new <- length(HistYr_new)
  
  if (StartYr == StartYr_new && CurrentYr == CurrentYr_new && proyears == proyears_new) {
    return(MOM)
  } else if(CurrentYr_new > CurrentYr) {
    stop("CurrentYr_new must be less than :", CurrentYr)
  } else if (length(intersect(HistYr, HistYr_new)) < 2) {
    stop("No years in common between proposed MOM historical years.")
  } else if(proyears_new > proyears) {
    stop("Proposed proyears needs to be less than or equal to MOM proyears.")
  }
  
  MOM2 <- MOM
  
  np <- length(MOM2@Stocks)
  nf <- length(MOM2@Fleets[[1]])
  
  Sname <- MOM2@Stocks %>% names()
  if (is.null(Sname)) Sname <- 1:np
  
  Fname <- MOM2@Fleets[[1]] %>% names()
  if (is.null(Fname)) Fname <- 1:nf
  
  MOM2@proyears <- proyears_new
  for (p in 1:np) {
    for (f in 1:nf) {
      MOM2@Fleets[[p]][[f]]@nyears <- nyears_new
      MOM2@Fleets[[p]][[f]]@CurrentYr <- CurrentYr_new
      
      if (length(MOM@cpars)) { # Assume Find and qs is specified in cpars
        if (!silent) message("Updating cpars for Stock: ", Sname[p], ", Fleet: ", Fname[f])
        cpars <- MOM@cpars[[p]][[f]]
        MOM2@cpars[[p]][[f]] <- lapply(names(cpars), cpars_subset_years, cpars_list = MOM@cpars, HistYr = HistYr,
                                       HistYr_new = HistYr_new, proyears_new = proyears_new, p = p, f = f,
                                       maxage = MOM2@Stocks[[p]]@maxage) %>% structure(names = names(cpars))
      }
    }
  }
  
  return(MOM2)
}


cpars_extend_age <- function(xx, cpars, n_age_new, n_age) {
  cpars_with_age <- c("Len_age", "LatASD", "Wt_age", "Mat_age", "M_ageArray", "Wt_age_C", "Fec_age", "V", "retA", "Fdisc_array1",
                      "mov", "Perr_y")
  
  x <- cpars[[xx]]
  if (length(dim(x)) && any(xx == cpars_with_age)) {
    
    if (xx == "mov") {
      
      if (length(dim(x)) == 4) {
        par_maxage <- x[, n_age, , ]
        x_new <- array(par_maxage, c(dim(par_maxage), n_age_new - n_age)) %>% aperm(c(1, 4, 2, 3))
      }
      if (length(dim(x)) == 5) {
        par_maxage <- x[, n_age, , , ]
        x_new <- array(par_maxage, c(dim(par_maxage), n_age_new - n_age)) %>% aperm(c(1, 5, 2:4))
      }
      x_out <- abind::abind(x, x_new, along = 2)
      
    } else if (xx == "Perr_y") {
      
      par_maxage <- x[, 1]
      x_new <- matrix(par_maxage, length(par_maxage), n_age_new - n_age)
      x_out <- cbind(x_new, x)
      
    } else if (length(dim(x)) == 3 && dim(x)[2] == n_age) {
      
      par_maxage <- x[, n_age, ]
      if(!is.null(cpars$plusgroup) && !cpars$plusgroup && xx == "M_ageArray") { # Set very high M if all animals die after n_age
        par_maxage[] <- 1e8
      }
      x_new <- array(par_maxage, c(dim(par_maxage), n_age_new - n_age)) %>% aperm(c(1, 3, 2))
      x_out <- abind::abind(x, x_new, along = 2)
      
    } else {
      stop("Don't know what to do with cpars: ", xx)
    }
    
  } else {
    x_out <- x
    if (xx == "Data") {
      x_out@MaxAge <- NA
      if (length(x_out@CAA) > 1 && dim(x_out@CAA)[3] == n_age) {
        x_out@CAA <- abind::abind(x_out@CAA, array(0, c(dim(x_out@CAA)[1:2], n_age_new - n_age)), along = 3)
      }
      
      if (length(x_out@Vuln_CAA) > 1 && dim(x_out@Vuln_CAA)[2] == n_age) {
        x_out@Vuln_CAA <- cbind(x_out@Vuln_CAA, 
                                matrix(x_out@Vuln_CAA[, n_age], nrow(x_out@Vuln_CAA), n_age_new - n_age))
      }
      if (length(x_out@AddIndV) > 1 && dim(x_out@AddIndV)[3] == n_age) {
        x_out@AddIndV <- abind::abind(x_out@AddIndV, 
                                      array(x_out@AddIndV[, , n_age], c(dim(x_out@AddIndV)[1:2], n_age_new - n_age)), 
                                      along = 3)
      }
    }
  }
  return(x_out)
}

cpars_subset_years <- function(xx, cpars_list, HistYr, HistYr_new, proyears_new, p, f, maxage) {
  
  if (max(HistYr_new) > max(HistYr)) {
    stop("The new last historical year can not be later than the previous last historical year.")
  }
  
  cpars_array_fullyear <- c("Len_age", "LatASD", "Wt_age", "Mat_age", "M_ageArray", "ageM", "age95", "Wt_age_C", "Fec_age",
                            "V", "SLarray", "retA", "retL", "Fdisc_array1", "Fdisc_array2", "AddIerr")
  cpars_matrix_fullyear <- c("DR_y", "Cerr_y", "Cobs_y", "Ierr_y", "SpIerr_y", "VIerr_y", "Derr_y", "Aerr_y", "Recerr_y",
                             "Eerr_y", "Eobs_y", "Mrand", "Krand", "Linfrand", "Linfarray", "Karray", "t0array", "ageMarray",
                             "age95array", "Marray")
  cpars_matrix_proyears <- c("qvar", "TAC_y", "E_y", "SizeLim_y")
  
  nyears <- length(HistYr)
  nyears_new <- length(HistYr_new)
  
  HistYr_ind <- HistYr_NA <- match(HistYr_new, HistYr)
  extend_init <- min(HistYr_new) < min(HistYr)
  if(extend_init && any(is.na(HistYr_ind))) {
    HistYr_ind[is.na(HistYr_ind)] <- 1
  }
  
  cpars <- cpars_list[[p]][[f]]
  x <- cpars[[xx]]
  if (xx == "mov" && length(dim(x)) == 5) {
    
    x_hist <- x[, , , , HistYr_ind]
    x_pro <- x[, , , , nyears + 1:proyears_new]
    x_out <- abind::abind(x_hist, x_pro, along = 5)
    
  } else if (xx == "Perr_y") {
    x_out <- matrix(1, nrow(x), maxage + nyears_new + proyears_new)
    
    if (extend_init) { # extend init period, assume F = 0 in missing years but abundance should still match
      nadd <- nyears_new - nyears
      
      x_out[, nadd + 1:(maxage + nyears)] <- x[, 1:(maxage + nyears)]
      
    } else if (min(HistYr_new) > min(HistYr)) {  # cut model short - assume constant M-at-age during missing years
      
      # First year abundance
      # Need to deal with plusgroup esp. if F in missing years is high
      x_out[, maxage + 1:nyears_new] <- x[, maxage + HistYr_ind] # Historical recruitment during new historical period
      ncut <- min(HistYr_ind) - 1
      if(ncut <= maxage + 1) {
        FM <- sapply(cpars_list[[p]], function(f) {
          array(f$qs * f$Find[, 1:ncut], c(length(f$qs), ncut, maxage + 1)) * 
            aperm(f$V[, , 2:min(HistYr_ind) - 1], c(1, 3, 2))
        }, simplify = "array") %>% apply(1:3, sum)
        
        Perr_init <- x[, ncut + 1:maxage] 
        
        # Additional mortality from cumulative F
        for(a in 1:maxage - 1) { # Loop over a = 0, ... , maxage - 1
          for(y in 1:ncut) {
            if (a - y + 1 > 0) {
              aind <- a - y + 1
              yind <- ncut - y + 1
              Perr_init[, maxage - a] <- Perr_init[, maxage - a] * exp(-FM[, yind, aind])
            }
          }
        }
        x_out[, 1:maxage] <- Perr_init
        
      } else {
        stop("Number of years removed is greater than maxage + 1")
      }
      
    }
    
    x_out[, maxage + nyears_new + 1:proyears_new] <- x[, maxage + nyears + 1:proyears_new]
    
  } else if (xx == "Find") { 
    
    if (extend_init) { # extend init period, assume F = 0 in missing years
      x_out <- x[, HistYr_NA]
      x_out[is.na(x_out)] <- 0
    } else {
      x_out <- x[, HistYr_ind]
    }
    
  } else if (xx == "MPA") {
    
    x_hist <- x[HistYr_ind, ]
    x_pro <- x[nyears + 1:proyears_new, ]
    x_out <- rbind(x_hist, x_pro)
    
  } else if (is.matrix(x)) {
    
    if (any(xx == cpars_matrix_fullyear)) {
      x_hist <- x[, HistYr_ind]
      x_pro <- x[, nyears + 1:proyears_new]
      x_out <- cbind(x_hist, x_pro)
    } else if (any(xx == cpars_matrix_proyears)) {
      x_out <- x[, 1:proyears_new]
    }
    
  } else if (is.array(x) && any(xx == cpars_array_fullyear)) {
    
    x_hist <- x[, , HistYr_ind]
    x_pro <- x[, , nyears + 1:proyears_new]
    x_out <- abind::abind(x_hist, x_pro, along = 3)
    
  } else if (xx == "Data") {
    x_out <- x
    
    if (length(x_out@Year) > 1) x_out@Year <- HistYr_new
    if (!is.na(x@LHYear)) x_out@LHYear <- max(HistYr_new)
    
    Data_matrix <- c("Cat", "CV_Cat", "Effort", "CV_Effort", "Ind", "CV_Ind", "SpInd", "CV_SpInd", "VInd", "CV_VInd",
                     "ML", "Lc", "Lbar")
    Data_array <- c("CAA", "CAL", "AddInd", "CV_AddInd")
    
    for(i in Data_matrix) {
      if (length(slot(x_out, i)) > 1 && dim(slot(x_out, i))[2] == nyears) {
        slot(x_out, i) <- slot(x_out, i)[, HistYr_NA, drop = FALSE]
      }
    }
    for(i in Data_array) {
      if (length(slot(x_out, i)) > 1 && dim(slot(x_out, i))[3] == nyears) {
        slot(x_out, i) <- slot(x_out, i)[, , HistYr_NA, drop = FALSE]
      }
    }
  } else {
    x_out <- x
  }
  
  return(x_out)

}

calculate_F_at_age <- function(x) {
  V <- sapply(x, getElement, "V", simplify = "array")
  retA <- sapply(x, getElement, "retA", simplify = "array")
  Find <- sapply(x, function(xx) xx$qs * xx$Find, simplify = "array")
  Fdisc1 <- sapply(x, getElement, "Fdisc_array1", simplify = "array")
  
  nsim <- dim(Find)[1]
  nyears <- dim(Find)[2]
  proyears <- dim(V)[3] - nyears
  
  F_at_age <- MSEtool:::single_fleet_F_at_age(Find, V, retA, Fdisc1, nsim, nyears) %>% pmax(1e-8)
  Find_new <- apply(F_at_age, c(1, 3), max)
  V_new <- apply(F_at_age, c(1, 3), function(x) x/max(x)) %>% aperm(c(2, 1, 3))
  V_pro <- replicate(proyears, V_new[, , nyears])
  
  list(Find = Find_new, V = abind::abind(V_new, V_pro, along = 3))
}



aggregate_fleet <- function(x, LL, f_name = c("Longline", "Other"), silent = FALSE, model_discards = FALSE) {
  np <- length(x@Stocks)
  nf <- length(x@Fleets[[1]])
  
  if (all(LL <= nf)) {
    if (!silent) message("Fleets to aggregate as Longline: ", paste0(LL, collapse = ", "))
  } else {
    stop("Can't find Fleets: ", paste0(LL[LL > nf], collapse = ", "))
  }
  
  nf_new <- 2
  MOM <- x
  
  MOM@CatchFrac <- local({
    Cat <- sapply(1:np, function(p) {
      sapply(1:nf, function(f) {
        xx <- x@cpars[[p]][[f]]$Data@Cat[1, ]
        xx[length(xx)]
      })
    }) # nf x np matrix
    
    tot_catch <- colSums(Cat)
    LL <- colSums(Cat[LL, , drop = FALSE])
    Cat2 <- cbind(LL, tot_catch - LL) %>% structure(dimnames = list(NULL, f_name)) # np x 2 matrix
    
    lapply(1:np, function(p) {
      matrix(Cat2[p, ]/sum(Cat2[p, ]), MOM@nsim, 2, byrow = TRUE) %>% 
        structure(dimnames = list(NULL, colnames(Cat2)))
    }) %>% structure(names = names(x@Stocks))
  })
  
  fLL <- list(LL, setdiff(1:nf, LL)) # length(fLL[[2]]) = 0 if there is no other fleet
  LL_only <- !length(fLL[[2]])
  
  if(LL_only) {
    MOM@Fleets <- lapply(1:np, function(p) x@Fleets[[p]][c(1, 1)] %>% structure(names = f_name))
    MOM@Obs <- lapply(1:np, function(p) x@Obs[[p]][c(1, 1)])
    MOM@Imps <- lapply(1:np, function(p) x@Imps[[p]][c(1, 1)])
  } else {
    MOM@Fleets <- lapply(1:np, function(p) x@Fleets[[p]][sapply(fLL, getElement, 1)] %>% structure(names = f_name))
    MOM@Obs <- lapply(1:np, function(p) x@Obs[[p]][sapply(fLL, getElement, 1)])
    MOM@Imps <- lapply(1:np, function(p) x@Imps[[p]][sapply(fLL, getElement, 1)])
  }
  
  MOM@cpars <- lapply(1:np, function(p) {
    vars <- c("Fdisc", "Fdisc_array1", "Fdisc_array2", "retL", "Fdisc", "SLarray")
    
    # Longline
    cpars_LL <- x@cpars[[p]][getElement(fLL, 1)]
    if(model_discards) {
      LL_update <- MSEtool:::calculate_single_fleet_dynamics(cpars_LL)
      cpars_LL[[1]][names(LL_update)] <- LL_update
    } else {
      LL_update <- calculate_F_at_age(cpars_LL)
      cpars_LL[[1]][vars] <- NULL
      cpars_LL[[1]]$Find <- LL_update$Find
      cpars_LL[[1]]$V <- LL_update$V
      cpars_LL[[1]]$retA[] <- 1
    }
    
    if(LL_only) { # Other
      cpars_oth <- x@cpars[[p]][1]
      cpars_oth[[1]]$Find[] <- 1e-8
    } else {
      cpars_oth <- x@cpars[[p]][getElement(fLL, 2)]
      
      if(model_discards) {
        oth_update <- MSEtool:::calculate_single_fleet_dynamics(cpars_oth)
        cpars_oth[[1]][names(oth_update)] <- oth_update
      } else {
        oth_update <- calculate_F_at_age(cpars_oth)
        cpars_oth[[1]][vars] <- NULL
        cpars_oth[[1]]$Find <- oth_update$Find
        cpars_oth[[1]]$V <- oth_update$V
        cpars_oth[[1]]$retA[] <- 1
      }
    }
    cpars_LL[[1]]$Data <- cpars_oth[[1]]$Data <- NULL
    list(cpars_LL[[1]], cpars_oth[[1]]) %>% structure(names = f_name)
  })
  
  if (length(x@Efactor)) {
    if(LL_only) {
      MOM@Efactor <- lapply(1:np, function(p) {
        matrix(c(1, 1e-8), x@nsim, ncol = 2, byrow = TRUE) %>% structure(dimnames = list(NULL, f_name))
      }) %>% structure(names = names(x@Efactor))
    } else {
      MOM@Efactor <- lapply(1:np, function(p) {
        cbind(apply(x@Efactor[[p]][, fLL[[1]], drop = FALSE], 1, mean), apply(x@Efactor[[p]][, fLL[[2]], drop = FALSE], 1, mean)) %>% 
          structure(dimnames = list(NULL, f_name))
      }) %>% structure(names = names(x@Efactor))
    }
  }
  
  if (length(x@CatchFrac)) {
    MOM@CatchFrac <- lapply(1:np, function(p) {
      cbind(rowSums(x@CatchFrac[[p]][, fLL[[1]], drop = FALSE]), rowSums(x@CatchFrac[[p]][, fLL[[2]], drop = FALSE])) %>% 
        structure(dimnames = list(NULL, f_name))
    }) %>% structure(names = names(x@CatchFrac))
  }
  if (length(x@Allocation)) {
    MOM@Allocation <- lapply(1:np, function(p) {
      cbind(rowSums(x@Allocation[[p]][, fLL[[1]], drop = FALSE]), rowSums(x@Allocation[[p]][, fLL[[2]], drop = FALSE])) %>% 
        structure(dimnames = list(NULL, f_name))
    }) %>% structure(names = names(x@Efactor))
  }
  
  #for(p in 1:np) {
  #  for(f in 1:nf) {
  #    MOM@Fleets[[p]][[f]]@Name <- names(MOM@Fleets[[p]])[f]
  #  }
  #}
  
  # Ignore SexPars, they are independent of fleet structure
  # Ignore Complexes and Rel for now
  return(MOM)
}
