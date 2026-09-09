
library(EcoTest)
library(EcoTestData)

# Indicator 3
i3 = readRDS("C:/GitHub/EcoTest_Train/2026 Analyses/i3_indicators_test.rds")
nams = sapply(i3,function(x)x$name)

# Indicator 4
fs = list.files("C:/GitHub/EcoTest_Train/2026 Analyses")
fsf = list.files("C:/GitHub/EcoTest_Train/2026 Analyses",full.names = T)
keep = grepl('i4',fs)
i4f = fsf[keep]
i4 = list()
for(ii in 1:length(if4))  i4[[ii]] = readRDS(iff[ii])

i3 = readRDS("C:/GitHub/EcoTest_Train/2026 Analyses/i3_indicators_test.rds" )

# Indicator 3 predictions
p3 = lapply(i3, pred_ind)


# Do retro plots
par(mfrow=c(3,3),par(mai=c(0.8,0.8,0.2,0.05)))
for(dd in 1:length(preds)){
  retro_ind(preds[[dd]])
  mtext(dats[dd],3)
}





