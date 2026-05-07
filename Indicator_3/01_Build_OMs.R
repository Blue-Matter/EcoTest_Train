# Indicator 3

# Jan 2026

library(openMSE)
library(miceadds)
library(SimDesign)

setwd("C:/GitHub/EcoTest")
source.all("Source")

largedir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind3"
setup(cpus=48)

for(x in 1001:1950) runbatch3(x, nsim=48)
