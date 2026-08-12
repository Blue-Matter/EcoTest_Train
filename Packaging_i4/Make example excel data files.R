

library(keras3)
library(r4ss)
library(data.table)
library(miceadds)

fdir = "C:/GitHub/EcoTest_train"
tdir = "C:/GitHub/EcoTest"
ddir = "C:/GitHub/EcoTestData"
source.all(paste0(fdir,"/Source"))
source.all(paste0(tdir,"/R"))



# Make an example EcoTest data input file from an SS assessment


ssfold = c(rep("G:/Shared drives/BM shared/1. Projects/TOF Advisory/Blue Shark MSE Workshop/Assessment_files_provided/",5),
           rep("G:/Shared drives/BM shared/1. Projects/EcoTest/Assessments/",4))

ssdirs = c("IO_Joel_Rice/SS3","NAtl_Nathan_Taylor/SS3","SAtl_Nathan_Taylor/SS3","SWPac_Philipp_Neubauer/SS3","NP_Nicholas_D_Barth/SS3",
           "SALB","NSWO","2019 WHM SS3/StockSynthesis/Model_6","2021 BET/M20_h0.8_sigmaR0.4")
ssnams = c("Shark_1", "Shark_2", "Shark_3", "Shark_4","Shark_5","Tuna_1","Billfish_1","Billfish_2","Tuna_2")
nss = length(ssnams)
mss = 1

ios = list()

for(dd in mss:nss)ios[[dd]] = get_ss3_inputs(dir=paste0(ssfold[dd],ssdirs[dd]))
names(ios) = ssnams









lapply(ios,ss_fleet_helper) # what are primary fleets / surveys

Fnams = list(c("F4_JPN_LL","F8_ESP_LL"),   c("F1_EU-ESP","F2_JPN"),
             c("FS2_BRA","FS4_JPN"),             c("F_Tar_NZ","F_Tar_EU"),
             c("F1_MEX","F19_TAIW_LG"),
             c("Fleet_01","Fleet_02"),  c("US_2","PORT_6"),
             c("LongLine_2","Gill_Net_1"),c("11_Japan_LL_TRO","4_PS_FAD_9119"))
Inams = list(c("S2_JPN_LATE","S4_EU_ESP"), c("S1_ESP-LL-N","S2_JP-LL-N"),
             c("BRA_index_TB2","JPN_index_TB2"), c("S_Tar_NZ","S_Tar_EU"),
             c("S10_MEX","S3_TAIW_LG"),
             c("CTP-LL", "JPN-LL3"), c("US_Survey_12","PORT_Survey_13"),
             c("US_LL","Ven_GN"), c("11_Japan_LL_TRO","4_PS_FAD_9119"))

for(dd in mss:nss){
  temp = SS_2_ET(io=ios[[dd]], Fnam = Fnams[[dd]], Inam = Inams[[dd]])
  objname = paste0(ssnams[dd],"_data")
  assign(objname,temp)
  to_file=paste0(tdir,"/data/",ssnams[dd],"_data.rda")
  do.call(save, list(objname,file=to_file))
}

for(dd in mss:nss){
  temp = SS_2_ET_Retro(io=ios[[dd]], Fnam = Fnams[[dd]], Inam = Inams[[dd]], npeels=8)
  objname = paste0(ssnams[dd],"_retro")
  assign(objname,temp)
  to_file=paste0(tdir,"/data/",ssnams[dd],"_retro.rda")
  do.call(save, list(objname,file=to_file))
}
