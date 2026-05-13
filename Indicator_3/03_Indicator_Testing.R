

library(EcoTest)
library(EcoTestData)
packageVersion('EcoTest')
packageVersion('EcoTestData')

TD = prepare_TD()

setwd("G:/Shared drives/BM shared/1. Projects/EcoTest/SCRS 2026/Figures")

#NN = train_NN(TD, nodes=c(100,10), nepoch = 100)

keep = grepl("_s1", names(TD))
keep[1] = TRUE
TDs1 = TD[,keep]

NN_s1 = train_NN(TDs1, nodes=c(30,10), nepoch = 30)


keep = grepl("_T", names(TDs1))
keep[1] = TRUE
TDs1T = TDs1[,keep]


NN_T = train_NN(TDs1T, nodes=c(30,10), nepoch = 30)


keep = !grepl("I_",names(TDs1T))
keep[1]=TRUE
NI = TDs1T[,keep]


NN_NI = train_NN(NI, nodes=c(20,5), nepoch = 50)

Cind = (1:ncol(NI))[grepl("C_",names(NI))]
CO = NI[,1:Cind[2]]

NN_CO= train_NN(CO, nodes=c(10,5), nepoch = 30)


CLind = c(1:8,(1:ncol(NI))[grepl("C_",names(NI))|grepl("ML_",names(NI))|grepl("MV_",names(NI))|grepl("FM_",names(NI))])
CL = NI[,CLind]
NN_CL= train_NN(CL, nodes=c(50,10), nepoch = 30)


png("PP.png",units='in', width=10,height=12,res=400)
 par(mfrow=c(3,2),mai=c(0.5,0.5, 0.4, 0.2),omi=c(0.3,0.3,0.05,0.05))
 pred_plot(NN_s1,newplot=F, lab = "(A) All single species data")
 pred_plot(NN_T,newplot=F, lab = "(B) As A but with fleet data aggregated")
 pred_plot(NN_NI,newplot=F, lab = "(C) As B but without relative abundance data")
 pred_plot(NN_CL,newplot=F, lab = "(D) As C but without age data (only catch and length)")
 pred_plot(NN_CO,newplot=F, lab = "(E) Only life-history and catch data")
dev.off()
