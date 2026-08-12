
library(openMSE)
library(r4ss)
library(data.table)
library(miceadds)

tdir = "C:/GitHub/EcoTest"
fdir = "C:/GitHub/EcoTest_train"
source.all(paste0(fdir,"/Source"))
source.all(paste0(tdir,"/R"))

# Indicator 2 data generation


ssfold = c(rep("G:/Shared drives/BM shared/1. Projects/TOF Advisory/Blue Shark MSE Workshop/Assessment_files_provided/",5),
           rep("G:/Shared drives/BM shared/1. Projects/EcoTest/Assessments/",4))

ssdirs = c("IO_Joel_Rice/SS3","NAtl_Nathan_Taylor/SS3","SAtl_Nathan_Taylor/SS3","SWPac_Philipp_Neubauer/SS3","NP_Nicholas_D_Barth/SS3",
           "SALB","NSWO","2019 WHM SS3/StockSynthesis/Model_6","2021 BET/M20_h0.8_sigmaR0.4")

ssnams = c("BSH_IO", "BSH_NA", "BSH_SA", "BSH_S{","BSH_NP","ALB_SA","SWO_N","WHM_NA","BET_N")


nss = length(ssnams)

ios = list()
for(dd in 1:nss)ios[[dd]] = all_ss3(dir=paste0(ssfold[dd],ssdirs[dd]))
names(ios) = ssnams



oms = list()
for(dd in 1:nss)oms[[dd]] = SS2MOM(paste0(ssfold[dd],ssdirs[dd]))
names(oms) = ssnams

OMv2 = oms
#for(dd in 1:nss)OMv2[[dd]] = Convert(oms[[dd]])
#names(OMv2) = ssnams

hists = list()
for(dd in 1:nss)hists[[dd]] = Simulate(oms[[dd]])

OM = oms[[1]]
saveRDS(OM,"C:/temp/OM.rds")

temp = Convert(OM)







