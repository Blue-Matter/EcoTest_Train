

library(keras3)
library(r4ss)
library(data.table)
library(miceadds)
library(readxl)
library(writexl)

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
anams = c("Blue_Shark_IO", "Blue_Shark_NAtl", "Blue_Shark_SAtl", "Blue_Shark_SWPac","Blue_Shark_NPac",
          "Albacore_SAtl", "Swordfish_NAtl", "White_Marlin_Atl","Bigeye_Tuna_Atl")
nss = length(ssnams)
mss = 1

ios = list()

for(dd in mss:nss)ios[[dd]] = all_ss3(dir=paste0(ssfold[dd],ssdirs[dd]))
names(ios) = ssnams

# ---------------------------------------------------------------------
saveRDS(ios,"C:/GitHub/EcoTest_Train/Packaging_i4/all_ss3_io.rds")

# --------------------------------------------------------------------
ios = readRDS("C:/GitHub/EcoTest_Train/Packaging_i4/all_ss3_io.rds")


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

pre_format = list()
for(dd in mss:nss){
  cat(".")
  pre_format[[dd]] = SS_2_ET_raw(io=ios[[dd]], Fnam = Fnams[[dd]], Inam = Inams[[dd]])
}
names(pre_format) = anams

# ---------------------------------------------------------------------
saveRDS(pre_format,"C:/GitHub/EcoTest_Train/Packaging_i4/pre_format.rds")

# --------------------------------------------------------------------
pre_format = readRDS("C:/GitHub/EcoTest_Train/Packaging_i4/pre_format.rds")


# Make Excel Input files

for(dd in 1:nss){
  
  sum=pre_format[[dd]]
  xlfile = "C:/GitHub/EcoTest_Train/Indicator_4/Real_data/EcoTest_Input.xlsx"
  tofile = paste0("C:/GitHub/EcoTest_Train/Indicator_4/Real_data/",names(pre_format)[dd],".xlsx")
  spec = anams[dd]
  fillxl(sum, xlfile,tofile,spec)

}


