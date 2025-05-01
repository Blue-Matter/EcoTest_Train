
EcoTest_gettaxa = function (Class = "predictive", Order = "predictive", Family = "predictive", 
          Genus = "predictive", Species = "predictive", ParentChild_gz = MSEtool::LHdatabase$ParentChild_gz, 
          msg = TRUE) 
{
  Taxa_Table <- MSEtool::Taxa_Table
  Species2 <- strsplit(Taxa_Table$Species, " ")
  Match = 1:nrow(Taxa_Table)
  if (Class != "predictive") 
    Match = Match[which(tolower(Taxa_Table$Class[Match]) == 
                          tolower(Class))]
  if (Order != "predictive") 
    Match = Match[which(tolower(Taxa_Table$Order[Match]) == 
                          tolower(Order))]
  if (Family != "predictive") 
    Match = Match[which(tolower(Taxa_Table$Family[Match]) == 
                          tolower(Family))]
  if (Genus != "predictive") 
    Match = Match[which(tolower(Taxa_Table$Genus[Match]) == 
                          tolower(Genus))]
  if (Species != "predictive") {
    Species2 <- paste(Genus, Species)
    Match = Match[which(tolower(Taxa_Table$Species[Match]) == 
                          tolower(Species2))]
  }
  full_taxonomy <- c(Class, Order, Family, Genus, Species)
  spIn <- trimws(paste(gsub("predictive", "", full_taxonomy), 
                       collapse = " "))
  if (length(Match) == 0) {
    if (msg) 
      message(spIn, " not found in FishBase database")
    Class <- Order <- Family <- Genus <- Species <- "predictive"
  }
  full_taxonomy <- c(Class, Order, Family, Genus, Species)
  if (!all(Species == "predictive")) {
    if (length(unique(Taxa_Table$Species[Match])) != 1) 
      stop("inputs are not unique")
    if (length(unique(Taxa_Table$Species[Match])) == 1) {
      tmp = unique(Taxa_Table$Species[Match])[1]
      full_taxonomy[5] <- strsplit(tmp, " ")[[1]][2]
    }
  }
  if (!all(c(Species, Genus) == "predictive")) {
    if (length(unique(Taxa_Table$Genus[Match])) != 1) 
      stop("inputs are not unique")
    if (length(unique(Taxa_Table$Genus[Match])) == 1) 
      full_taxonomy[4] = unique(Taxa_Table$Genus[Match])[1]
  }
  if (!all(c(Species, Genus, Family) == "predictive")) {
    if (length(unique(Taxa_Table$Family[Match])) != 1) 
      stop("inputs are not unique")
    if (length(unique(Taxa_Table$Family[Match])) == 1) 
      full_taxonomy[3] = unique(Taxa_Table$Family[Match])[1]
  }
  if (!all(c(Species, Genus, Family, Order) == "predictive")) {
    if (length(unique(Taxa_Table$Order[Match])) != 1) 
      stop("inputs are not unique")
    if (length(unique(Taxa_Table$Order[Match])) == 1) 
      full_taxonomy[2] = unique(Taxa_Table$Order[Match])[1]
  }
  if (!all(c(Species, Genus, Family, Order, Class) == "predictive")) {
    if (length(unique(Taxa_Table$Class[Match])) != 1) 
      stop("inputs are not unique")
    if (length(unique(Taxa_Table$Class[Match])) == 1) 
      full_taxonomy[1] = unique(Taxa_Table$Class[Match])[1]
  }
  match_taxonomy = full_taxonomy
  fam_gen_sp <- tolower(paste(match_taxonomy[3:5], collapse = "_"))
  nm_ind <- which(grepl(fam_gen_sp, tolower(ParentChild_gz$ChildName)))
  if (length(nm_ind) == 0) {
    temp <- strsplit(fam_gen_sp, "_")[[1]]
    temp[3] <- "predictive"
    fam_gen_sp <- paste0(temp, collapse = "_")
    nm_ind <- which(grepl(fam_gen_sp, tolower(ParentChild_gz$ChildName)))
  }
  if (length(nm_ind) == 0) {
    temp <- strsplit(fam_gen_sp, "_")[[1]]
    temp[2] <- "predictive"
    fam_gen_sp <- paste0(temp, collapse = "_")
    nm_ind <- which(grepl(fam_gen_sp, tolower(ParentChild_gz$ChildName)))
  }
  fullname <- gsub("_", " ", ParentChild_gz$ChildName[nm_ind])
  if (length(fullname) > 1) 
    fullname <- fullname[length(fullname)]
  ind <- !grepl("predictive", strsplit(fullname, " ")[[1]])
  if (all(!ind)) {
    if (msg) 
      message_info("Predicting from all species in FishBase")
  }
  else if (any(!ind)) {
    if (msg) 
      message_info("Closest match: ", fullname)
  }
  else {
    if (msg) 
      message_info("Species match: ", fullname)
  }
  match_taxonomy = unique(as.character(Add_predictive(ParentChild_gz$ChildName[nm_ind])))
  if (length(match_taxonomy) > 1) 
    match_taxonomy <- match_taxonomy[length(match_taxonomy)]
  match_taxonomy
}



EcoTest_predictLH = function (inpars = list(), Genus = "predictive", Species = "predictive", 
          nsamp = 100, db = MSEtool::LHdatabase, dist = c("unif", "norm"), 
          filterMK = TRUE, plot = TRUE, Class = "predictive", Order = "predictive", 
          Family = "predictive", msg = TRUE) 
{
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  dist <- match.arg(dist)
  inpars_1 <- inpars
  names <- names(inpars)
  lens <- lapply(inpars, length) == 2
  valnames <- c("Linf", "L50", "K", "M")
  if (!prod(names %in% valnames)) 
    stop("invalid names in inpars. Valid names are: ", paste0(valnames, 
                                                              " "))
  notnames <- !valnames %in% names
  chkNA <- as.logical(unlist(lapply(lapply(inpars_1, is.na), 
                                    prod)))
  notnames[chkNA] <- TRUE
  inpars_1[valnames[!valnames %in% names]] <- NA
  if (any(notnames)) {
    for (x in seq_along(valnames)) {
      if (notnames[x]) {
        nm <- valnames[x]
        inpars[[nm]] <- NA
        if (msg) 
          message("Predicting ", nm)
      }
    }
  }
  if (prod(valnames[1:2] %in% names)) {
    if (lens["Linf"]) {
      if (prod(valnames[1:2] %in% names)) {
        if (lens["Linf"]) {
          if (msg) 
            message("Predicting L50 from Linf")
          lens["L50"] <- FALSE
          inpars$L50 <- NA
        }
        if (lens["L50"]) {
          if (msg) 
            message("Predicting Linf from L50")
          lens["Linf"] <- FALSE
          inpars$Linf <- NA
        }
      }
      lens["L50"] <- FALSE
      inpars$L50 <- NA
    }
    if (lens["L50"]) {
      if (msg) 
        message("Predicting Linf from L50")
      lens["Linf"] <- FALSE
      inpars$Linf <- NA
    }
  }
  if (prod(valnames[3:4] %in% names)) {
    if (lens["M"]) {
      if (msg) 
        message("Predicting K from M")
      lens["K"] <- FALSE
      inpars$K <- NA
    }
    if (lens["K"]) {
      if (msg) 
        message("Predicting M from K")
      lens["M"] <- FALSE
      inpars$M <- NA
    }
  }
  multi <- 100
  if (is.logical(filterMK)) {
    filter <- "OM"
    filterM <- filterK <- FALSE
    if (prod(c("K", "M") %in% names) & filterMK & !(all(is.na(inpars_1$K)) || 
                                                    all(is.na(inpars_1$M)))) {
      if (all(is.na(inpars$M))) {
        filterM <- TRUE
        if (msg) 
          message_info("Filtering predicted M within bounds:", 
                       paste0(inpars_1$M, collapse = "-"))
      }
      if (all(is.na(inpars$K))) {
        filterK <- TRUE
        if (msg) 
          message_info("Filtering predicted K within bounds:", 
                       paste0(inpars_1$K, collapse = "-"))
      }
      multi <- 500
    }
  }
  if (is.numeric(filterMK)) {
    filter <- "perc"
    if (length(filterMK) != 2) 
      stop("filterMK must be numeric values of length 2 (lower and upper percentiles) OR logical")
    if (all(is.na(inpars$M))) {
      filterM <- TRUE
      if (msg) 
        message_info("Filtering predicted M within percentiles:", 
                     paste0(filterMK, collapse = "-"))
    }
    if (all(is.na(inpars$K))) {
      filterK <- TRUE
      if (msg) 
        message_info("Filtering predicted K within percentiles:", 
                     paste0(filterMK, collapse = "-"))
    }
    multi <- 500
  }
  taxa <- MSEtool:::gettaxa(Class, Order, Family, Genus, Species, msg = msg)
  if (is.null(taxa)) 
    return(NULL)
  if (!methods::is(db, "list")) 
    stop("db must be database list from FishLife", call. = FALSE)
  Which <- grep(taxa, db$ParentChild_gz[, "ChildName"])
  mu <- db$ParHat$beta_gj[Which, ]
  covar <- db$Cov_gjj[Which, , ]
  names(mu) <- gsub("Loo", "Linf", names(mu))
  names(mu) <- gsub("Lm", "L50", names(mu))
  sampvals <- exp(as.data.frame(MASS::mvrnorm(nsamp * multi, 
                                              mu, covar))) %>% dplyr::select(Linf, K, M, L50)
  sampvals$relLm <- sampvals$L50/sampvals$Linf
  sampvals$MK <- sampvals$M/sampvals$K
  outpars <- inpars
  if (dist == "unif") {
    for (x in names(outpars)) {
      if (!all(is.na(outpars[[x]]))) 
        outpars[[x]] <- myrunif(nsamp * multi, min(outpars[[x]]), 
                                max(outpars[[x]]))
    }
  }
  else {
    for (x in names(outpars)) {
      if (!all(is.na(outpars[[x]]))) {
        varsd <- (max(outpars[[x]]) - mean(outpars[[x]]))/2
        varmean <- mean(outpars[[x]])
        outpars[[x]] <- rnorm(nsamp * multi, varmean, 
                              varsd)
      }
    }
  }
  missing <- lapply(lapply(outpars, is.na), prod) == 1
  missnm <- names(outpars)[missing]
  for (x in missnm) {
    if (x == "L50") {
      outpars$L50 <- outpars$Linf * sampvals$relLm
    }
    if (x == "Linf") {
      outpars$Linf <- outpars$L50/sampvals$relLm
    }
    if (x == "M") {
      outpars$M <- outpars$K * sampvals$MK
    }
    if (x == "K") {
      outpars$K <- outpars$M/sampvals$MK
    }
  }
  missing <- lapply(lapply(outpars, is.na), prod) == 1
  missnm <- names(outpars)[missing]
  for (x in missnm) {
    outpars[[x]] <- sampvals[[x]]
  }
  Out <- as.data.frame(do.call("cbind", outpars))
  ind <- Out$L50 > 0.95 * Out$Linf
  Out <- Out[!ind, ]
  if (filter == "OM") {
    if (filterK) {
      ind <- Out$K > min(inpars_1$K) & Out$K < max(inpars_1$K)
      if (sum(ind) < 2) {
        warning("No samples of K within bounds: ", paste(as.character(inpars_1$K), 
                                                         collapse = " "), "\nIgnoring bounds on K")
      }
      else {
        Out <- Out[ind, ]
      }
    }
    if (filterM) {
      ind <- Out$M > min(inpars_1$M) & Out$M < max(inpars_1$M)
      if (sum(ind) < 2) {
        warning("No samples of M within bounds: ", paste(as.character(inpars_1$M), 
                                                         collapse = " "), "\nIgnoring bounds on M")
      }
      else {
        Out <- Out[ind, ]
      }
    }
  }
  if (filter == "perc") {
    if (filterK) {
      ind <- Out$K > quantile(Out$K, filterMK[1]) & Out$K < 
        quantile(Out$K, filterMK[2])
      Out <- Out[ind, ]
    }
    if (filterM) {
      ind <- Out$M > quantile(Out$M, filterMK[1]) & Out$M < 
        quantile(Out$M, filterMK[2])
      Out <- Out[ind, ]
    }
  }
  if (nrow(Out) < nsamp) {
    warning("Could not generate ", nsamp, " samples within specified bounds. Sampling with replacement")
    rows <- sample(1:nrow(Out), nsamp, replace = TRUE)
    Out <- Out[rows, ]
  }
  else {
    Out <- Out[1:nsamp, ]
  }
  Out <- Out %>% dplyr::select(valnames)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(4, 4), mai = c(0.3, 0.3, 0.4, 0.05), omi = c(0.02, 
                                                               0.02, 0.3, 0.02))
    colline = MSEtool::makeTransparent("blue", 60)
    lwdline = 4
    histcol = "black"
    if (length(inpars_1) > 1) {
      dobounds <- TRUE
      bounds <- matrix(NA, nrow = 2, ncol = 4)
      count <- 0
      for (nm in valnames) {
        count <- count + 1
        tempval <- inpars_1[[nm]]
        if (length(tempval) > 0) 
          bounds[, count] <- range(tempval)
      }
    }
    else {
      dobounds <- FALSE
    }
    for (i in 1:4) {
      if (length(inpars_1[[valnames[i]]]) > 1) {
        rng <- range(inpars_1[[valnames[i]]])
      }
      else {
        rng <- c(NA, NA)
      }
      rng2 <- range(Out[, i])
      if (!all(is.na(rng))) {
        rng2[1] <- min(c(min(rng, na.rm = TRUE), min(rng2, 
                                                     na.rm = TRUE)), na.rm = TRUE)
        rng2[2] <- max(c(max(rng, na.rm = TRUE), max(rng2, 
                                                     na.rm = TRUE)), na.rm = TRUE)
      }
      for (j in 1:4) {
        if (i == j) {
          if (i == 1) {
            hist(Out[, 1], main = "Asymptotic length (Linf)", 
                 col = histcol, border = "white", xlab = "", 
                 axes = F, xlim = rng2, ylab = "")
            axis(1)
            abline(v = rng, col = colline, lwd = lwdline)
          }
          else if (i == 2) {
            hist(Out[, 2], main = "Length at 50% maturity (L50)", 
                 col = histcol, border = "white", xlab = "", 
                 axes = F, xlim = rng2, ylab = "")
            axis(1)
            abline(v = rng, col = colline, lwd = lwdline)
          }
          else if (i == 3) {
            hist(Out[, 3], main = "Growth rate (K)", 
                 col = histcol, border = "white", xlab = "", 
                 axes = F, xlim = rng2, ylab = "")
            axis(1)
            abline(v = rng, col = colline, lwd = lwdline)
          }
          else {
            hist(Out[, 4], main = "Natural mortality rate (M)", 
                 col = histcol, border = "white", xlab = "", 
                 axes = F, xlim = rng2, ylab = "")
            axis(1)
            abline(v = rng, col = colline, lwd = lwdline)
          }
        }
        else {
          plot(Out[, j], Out[, i], axes = F, col = "white", 
               xlab = "", ylab = "")
          if (dobounds) 
            polygon(bounds[c(1, 1, 2, 2), j], bounds[c(1, 
                                                       2, 2, 1), i], col = colline, border = "white")
          points(Out[, j], Out[, i], pch = 19)
          axis(1)
          axis(2)
        }
      }
    }
  }
  Out
}



EcoTest_LH2OM = function (OM, dist = c("unif", "norm"), filterMK = FALSE, plot = TRUE, 
          Class = "predictive", Order = "predictive", Family = "predictive", 
          msg = TRUE, db = MSEtool::LHdatabase) 
{

 # dist = c("unif", "norm"); filterMK = FALSE; plot = TRUEClass = "Actinopterygii"; Order = "Scombriformes"; Family = "Scombridae"; msg = TRUE; db = MSEtool::LHdatabase
  
    if (!methods::is(OM, "OM")) 
    stop("OM must be class 'OM'")
  dist <- match.arg(dist)
  set.seed(OM@seed)
  if (length(OM@nsim) < 1) 
    OM@nsim <- 48
  if (length(OM@cpars) > 0) {
    cnames <- names(OM@cpars)
    if (any(c("Linf", "L50", "M", "K") %in% cnames)) {
      message("Life-history parameters already in OM@cpars.\nReturning original OM")
      return(OM)
    }
  }
  sls <- c("Linf", "L50", "K", "M")
  for (sl in sls) {
    slval <- slot(OM, sl)
    if (any(is.na(slval)) | all(slval == 0) | length(slval) < 
        1) {
      assign(sl, NA)
    }
    else {
      assign(sl, slval)
    }
  }
  if (length(OM@M) > 2) {
    message("Age-dependant M has been set in OM. (length(OM@M) > 2)\nReturning original OM")
    return(OM)
  }
  Genus <- unlist(strsplit(OM@Species, " "))[1]
  Species <- unlist(strsplit(OM@Species, " "))[2]
  if (is.na(Genus) || nchar(Genus) < 1) 
    Genus <- "predictive"
  if (is.na(Species) || nchar(Species) < 1) 
    Species <- "predictive"
  Out <- EcoTest_predictLH(inpars = list(Linf = Linf, L50 = L50, K = K, 
                                 M = M), Genus, Species, nsamp = OM@nsim, db = db, dist = dist, 
                   filterMK = filterMK, plot = plot, Class = Class, Order = Order, 
                   Family = Family, msg = msg)
  if (is.null(Out)) {
    message("Could not complete prediction. Returning original OM")
    return(OM)
  }
  OM@Linf <- c(0, 0)
  OM@L50 <- c(0, 0)
  OM@M <- c(0, 0)
  OM@K <- c(0, 0)
  if (any(Out$Linf <= 0)) 
    warning("Some Linf values are <= 0")
  if (any(Out$M <= 0)) 
    warning("Some M values are <= 0")
  if (any(Out$K <= 0)) 
    warning("Some K values are <= 0")
  if (any(Out$L50 <= 0)) 
    warning("Some L50 values are <= 0")
  OM@cpars$Linf <- Out$Linf
  OM@cpars$M <- Out$M
  OM@cpars$K <- Out$K
  OM@cpars$L50 <- Out$L50
  OM
}


cat("Instructive (depricated) EcoTest LH2OM code loaded \n")

