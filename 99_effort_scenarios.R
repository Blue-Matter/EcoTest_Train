
# Make a set of effort simulations for projections that are roughly around FMSY levels for each stock

library(MSEtool)
library(dplyr)
library(mvtnorm)
library(PerformanceAnalytics)
library(ecostats)

setwd("C:/GitHub/Ecotest")

source("99_batching.R")
source("99_plotting.R")
source("99_MOM_fixes.R")
source("99_make_effort.R")


# ---- Correlated F and F variability ---------------------------------------------------------

# === 1 === run a basic current effort projection to get approximate ratio of F/FMSY =======================          1

MOM =  readRDS("./MOM/MOM_latest.rds") %>% overwritePE %>% add_stochasticity
Hist = SimulateMOM(MOM)
saveRDS(Hist,"C:/temp/Ecotest/Hist_Base.rds")

E1 <- function(x, DataList, reps = 1, ...) {
  np <- length(DataList)
  nf <- length(DataList[[1]])
  RecList <- lapply(1:np, function(p) replicate(nf, new("Rec")))
  for(p in 1:np) { 
    for(f in 1:nf) { 
      RecList[[p]][[f]]@Effort <- 1
    }
  }
  RecList
}
class(E1) = "MMP"

MMSE = ProjectMOM(Hist, MPs = "E1", checkMPs = FALSE) 
MMSE = trim_MMSE(MMSE)
saveRDS(MMSE,"C:/temp/Ecotest/CalibrationFs/MMSE_Base.rds")



# === 2 === look at F/FMSY and remake effort ================================================================        2

MMSE = readRDS("C:/temp/Ecotest/MMSE_Base.rds")

curFs = apply(MMSE@FM[,,1,1,2:5],2,mean) # Fleet 1, MP 1, proyear 1 (time invariant)
out = plotF(MMSE)
Fadj = 1/out[[1]]

# create the matrices effmod per sim and year

Frel = readRDS("Frel/F_lm.rds") # multivariate linear model of Fs (real space)

nsim = 50000
proyears = 50
ns = 6; nt = 2

# - Errors -------------------------------------------------

targCVs = c(0.1,0.1)
SWOslopeCV = 0.05

# Fleet 1 - BET effort trajectory -------------------------

plotind = 1:10
trad = runif(nsim,-0.02,0.02) # BET % change
tarr = array(NA,c(nsim,ns,proyears))

yrind = rep(1:proyears,each=nsim)
slp = array((1+rep(trad,proyears))^yrind,c(nsim, proyears)); matplot(t(slp[plotind,]),type="l")
soffset = array(runif(nsim,0,proyears),c(nsim,proyears))
sdiv = array(runif(nsim,5,10),c(nsim,proyears))
smag = array(runif(nsim,0,0.3),c(nsim,proyears))
sina = array(exp(sin(((yrind-soffset) / sdiv))*smag),c(nsim,proyears)); matplot(t(sina[plotind,]),type="l")

muE = array(NA,c(nsim,ns,proyears))
muE[,1,] = curFs[1]* sina * slp * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears))
matplot(t(muE[plotind,1,]),type="l",ylim=c(0,max(muE[plotind,1,])))

# Fleet 2 Swordfish with error            error in slope by sim                          error in slope by sim                    
muE[,2,] = curFs[2] * sina * slp * array(trlnorm(nsim,1,SWOslopeCV),c(nsim,proyears)) * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears)) 
matplot(t(muE[plotind,2,]),type="l",ylim=c(0,max(muE[plotind,2,]))); plot(muE[plotind,1,],muE[plotind,2,])

ord = c(3,6,4,5) # position in the muE array
ind = as.matrix(expand.grid(1:nsim,ord,1))

for(y in 1:proyears){
  newdat = data.frame(BET = muE[,1,y], SWO = muE[,2,y])
  ind[,3] = y
  muE[ind] = Simfunc(Frel,newdat)
  cat("*")
}
cat("\n")
plot_an_F_sim(muE,1,MMSE@Snames)

muE_AC = ACfunc(muE,0.95,enp.mult=0.2,ploty=F)
plot_an_F_sim(muE_AC,7,MMSE@Snames)

#muEtot = array(aperm(muE_AC,c(1,3,2)),c(nsim*proyears,ns)); colnames(muEtot) = MMSE@Snames
#chart.Correlation(muEtot)

relEff = muE_AC * array(rep(Fadj,each=nsim),dim(muE_AC))
saveRDS(relEff,"./Batch/totEffmat_cor.rda")




# === 3 iteratively rerun analyses with new F adjustment ==========================================================================         3

ratey = 1 # the extend of changes to the F adjustment multiplier

for(i in 1:5){
  totEffmat <<- readRDS("./Batch/totEffmat_cor.rda") 
  
  runbatch(1,MOM=MOM,MPs = "Ftv_MMP", largedir=largedir)
  
  MMSE = readRDS("C:/temp/Ecotest/batching/Dependent_F/MMSE_1.rda")
  out = plotF(MMSE)
  Fadj = 1/out[[1]]
  
  Fadj2 = exp(log(Fadj)*ratey)
  
  relEff =relEff * array(rep(Fadj2,each=nsim),dim(muE_AC))
  saveRDS(relEff,"./Batch/totEffmat_cor.rda")

}

# save reference sims

saveRDS(relEff,"./Batch/totEffmat_calib.rda")


# make big set and rescale


Frel = readRDS("Frel/F_lm.rds") # multivariate linear model of Fs (real space)

nsim = 60000
proyears = 50
ns = 6; nt = 2

# - Errors -------------------------------------------------

targCVs = c(0.1,0.1)
SWOslopeCV = 0.05

# Fleet 1 - BET effort trajectory -------------------------

plotind = 1:10
trad = runif(nsim,-0.02,0.02) # BET % change
tarr = array(NA,c(nsim,ns,proyears))

yrind = rep(1:proyears,each=nsim)
slp = array((1+rep(trad,proyears))^yrind,c(nsim, proyears)); matplot(t(slp[plotind,]),type="l")
soffset = array(runif(nsim,0,proyears),c(nsim,proyears))
sdiv = array(runif(nsim,5,10),c(nsim,proyears))
smag = array(runif(nsim,0,0.3),c(nsim,proyears))
sina = array(exp(sin(((yrind-soffset) / sdiv))*smag),c(nsim,proyears)); matplot(t(sina[plotind,]),type="l")

muE = array(NA,c(nsim,ns,proyears))
muE[,1,] = curFs[1]* sina * slp * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears))
matplot(t(muE[plotind,1,]),type="l",ylim=c(0,max(muE[plotind,1,])))

# Fleet 2 Swordfish with error            error in slope by sim                          error in slope by sim                    
muE[,2,] = curFs[2] * sina * slp * array(trlnorm(nsim,1,SWOslopeCV),c(nsim,proyears)) * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears)) 
matplot(t(muE[plotind,2,]),type="l",ylim=c(0,max(muE[plotind,2,]))); plot(muE[plotind,1,],muE[plotind,2,])

ord = c(3,6,4,5) # position in the muE array
ind = as.matrix(expand.grid(1:nsim,ord,1))

for(y in 1:proyears){
  newdat = data.frame(BET = muE[,1,y], SWO = muE[,2,y])
  ind[,3] = y
  muE[ind] = Simfunc(Frel,newdat)
  cat("*")
}
cat("\n")
plot_an_F_sim(muE,1,MMSE@Snames)

muE_AC = ACfunc(muE,0.95,enp.mult=0.2,ploty=F)
plot_an_F_sim(muE_AC,7,MMSE@Snames)


muE_AC_calib = readRDS("./Batch/totEffmat_calib.rda")

rat = apply(muE_AC[1:100,,2:20],2,quantile,p=0.5) / apply(muE_AC_calib[1:100,,2:20],2,quantile,p=0.5)

calibE = muE_AC / array(rep(rat,each=nsim),dim(muE_AC))

saveRDS(calibE,"./Batch/totEffmat_cor.rda")















