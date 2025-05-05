# Indicator 2

# November 2024

library(openMSE)
library(miceadds)
library(SimDesign)

setwd("C:/GitHub/EcoTest")
source.all("Source")

largedir = "C:/Users/Admin/Dropbox/temp/Ecotest/Ind2"
setup(cpus=48)
#Ind2Export()

for(x in 1:100) runbatch2(x, nsim=48)


for(x in 101:500) runbatch2(x, nsim=48)

for(x in 501:800) runbatch2(x, nsim=48)

#system.time({  sfSapply(1:20, runbatch2, nsim = 4)  })


# can probably do this with depletion x VB0 x Find[lastyr] or there abouts (with some error maybe)
