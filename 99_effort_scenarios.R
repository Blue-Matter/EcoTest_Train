
library(MSEtool)
library(dplyr)
library(mvtnorm)
library(PerformanceAnalytics)
library(ecostats)

setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_batching.R")
source("99_plotting.R")
source("99_MOM_fixes.R")
source("99_make_effort.R")

# ---- constant F policy (a preliminary analysis and a hack and used totEffmat) --------------

# totEffmat <- readRDS("./Batch/totEffmat.rda")



# ---- Correlated F and F variability ---------------------------------------------------------

# run a basic current effort projection to get approximate ratio of F/FMSY

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
saveRDS(MMSE,"C:/temp/Ecotest/MMSE_Base.rds")



# look at F/FMSY

MMSE = readRDS("C:/temp/Ecotest/MMSE_Base.rds")
curFs = apply(MMSE@FM[,,1,1,1],2,mean) # Fleet 1, MP 1, proyear 1 (time invariant)
out = plotF(MMSE)
Fadj = 1/out[[1]]

# create the matrices effmod per sim and year

Frel = readRDS("Frel/F_lm.rds") # multivariate linear model of Fs (real space)

nsim = 100
proyears = 50
ns = 6; nt = 2

# - Errors -------------------------------------------------

targCVs = c(0.1,0.1)
SWOslopeCV = 0.05


# Fleet 1 - BET effort trajectory -------------------------

trad = runif(nsim,-0.02,0.02) # BET % change
tarr = array(NA,c(nsim,ns,proyears))

yrind = rep(1:proyears,each=nsim)
slp = array((1+rep(trad,proyears))^yrind,c(nsim, proyears)); matplot(t(slp),type="l")
soffset = array(runif(nsim,0,proyears),c(nsim,proyears))
sdiv = array(runif(nsim,5,10),c(nsim,proyears))
smag = array(runif(nsim,0,0.3),c(nsim,proyears))
sina = array(exp(sin(((yrind-soffset) / sdiv))*smag),c(nsim,proyears)); matplot(t(sina),type="l")

muE = array(NA,c(nsim,ns,proyears))
muE[,1,] = curFs[1]* sina * slp * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears))
matplot(t(muE[,1,]),type="l",ylim=c(0,max(muE[,1,])))

# Fleet 2 Swordfish with error            error in slope by sim                          error in slope by sim                    
muE[,2,] = curFs[2] * sina * slp * array(trlnorm(nsim,1,SWOslopeCV),c(nsim,proyears)) * array(trlnorm(nsim*proyears,1,targCVs[1]),c(nsim,proyears)) 
matplot(t(muE[,2,]),type="l",ylim=c(0,max(muE[,2,]))); plot(muE[,1,],muE[,2,])

ord = c(3,6,4,5) # position in the muE array
ind = as.matrix(expand.grid(1:nsim,ord,1))

for(y in 1:proyears){
  newdat = data.frame(BET = muE[,1,y], SWO = muE[,2,y])
  ind[,3] = y
  muE[ind] = Simfunc(Frel,newdat)
}
plot_an_F_sim(muE,1,MMSE@Snames)

muE_AC = ACfunc(muE,0.95,enp.mult=0.2,ploty=F)
plot_an_F_sim(muE_AC,7,MMSE@Snames)

muEtot = array(aperm(muE_AC,c(1,3,2)),c(nsim*proyears,ns)); colnames(muEtot) = MMSE@Snames
chart.Correlation(muEtot)

relEff = muE_AC / array(rep(curFs,each=nsim),dim(muE_AC)) * array(rep(Fadj,each=nsim),dim(muE_AC))

saveRDS(relEff,"./Batch/totEffmat_cor.rda")






























totEffmat <- readRDS("./Batch/totEffmat.rda")

source("99_Indicators.R")
source("99_make_MMP.R")
source("99_batching.R")
source("99_MOM_fixes.R")
source("99_plotting.R")

MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_2.rda")
out = plotF(MMSE)
Fadj = rep(NA,6)
Fadj[1:6] = 1/out[[1]][1:6]
# round(Fadj,2)

totsims = 250000
cv = 0.5
totEffmat = array(NA,c(totsims,10))
for(i in 1:10)totEffmat[,i] = rlnorm(totsims,log(Fadj[i]),cv)

varcov = matrix(0,10,10)
diag(varcov) = 0.50
totEffmat = exp(mvtnorm::rmvnorm(totsims,mean=log(Fadj),sigma = varcov))
Frelmat = totEffmat / array(rep(Fadj,each=totsims),dim(totEffmat))
isout = function(x)sum(!(x>0.35 & x<3))==0
keep = apply(Frelmat,1,isout); round(sum(keep)/length(keep)*100,3)
totEffmat = totEffmat[keep,]; print(nrow(totEffmat))
saveRDS(totEffmat[1:50000,],"./Batch/totEffmat.rda")

apply(totEffmat,2,quantile,p=0.5)
apply(totEffmat,2,mean)
apply(totEffmat,2,function(x)sd(x)/mean(x))
hist(totEffmat[,3],30)



