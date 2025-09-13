# Indicator 2

# November 2024

library(openMSE)
library(miceadds)
library(SimDesign)

setwd("C:/GitHub/EcoTest")
source.all("Source")

largedir = "C:/Users/Admin/Dropbox/temp/Ecotest/Ind3"
setup(cpus=48)
#Ind2Export()

for(x in 1:500) runbatch3(x, nsim=48)





for(x in 101:500) runbatch3(x, nsim=48)

for(x in 501:800) runbatch3(x, nsim=48)

for(x in 601:800) runbatch3(x, nsim=48)

for(x in 801:1000) runbatch3(x, nsim=48)

for(x in 1001:1200) runbatch3(x, nsim=48)


#system.time({  sfSapply(1:20, runbatch2, nsim = 4)  })


# can probably do this with depletion x VB0 x Find[lastyr] or there abouts (with some error maybe)
