
library(MSEtool)
library(dplyr)

setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

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