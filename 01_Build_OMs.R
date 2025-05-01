# Indicator 2

# November 2024

library(openMSE)
library(miceadds)
setwd("C:/GitHub/EcoTest")
source.all("Source")

MOMtest = make_MOM()
histtest = SimulateMOM(MOMtest) # need catch frac

# can probably do this with depletion x VB0 x Find[lastyr] or there abouts (with some error maybe)
