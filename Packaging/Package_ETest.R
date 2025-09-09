

library(keras3)
library(r4ss)
library(data.table)

fdir = "C:/GitHub/EcoTest"
tdir = "C:/GitHub/ETest"


# Run numbers
model_index = readRDS(paste0(fdir,"/Indicator_2/Fitted_Models/runs.rds"))
save(model_index,file = paste0(tdir,"/data/model_index.rda"))

# Fitted models
for(x in 1:nrow(model_index)){
  model = load_model(paste0(fdir,"/Indicator_2/Fitted_Models/Model_",x,".keras"))
  save_model(model,paste0(tdir,"/data/Model_",x,".keras"))
}


# Validation data

allout_1 = readRDS(paste0(fdir,"/Indicator_2/allout_1_200.rds"))
allout_2 = readRDS(paste0(fdir,"/Indicator_2/allout_201_400.rds"))
allout_3 = readRDS(paste0(fdir,"/Indicator_2/allout_401_600.rds"))
TD = rbindlist(c(allout_1, allout_2, allout_3))
save(TD,file=paste0(tdir,"/data/TD.rda"))


# -code

files = c("i2_neural_net.R", "i2_SS_funcs","99_neural_net","99_indicators")
for(ff in 1:length(files)){
  file.copy(paste0(fdir,"/Source/",files[ff]),paste0(tdir,"/R/",files[ff]),overwrite = T)
}

# - Data objects

ssfold = "G:/Shared drives/BM shared/1. Projects/TOF Advisory/Blue_Shark_MSE/Assessment_files_provided/"
ssdirs = c("IO_Joel_Rice/SS3","NAtl_Nathan_Taylor/SS3","SAtl_Nathan_Taylor/SS3","SWPac_Philipp_Neubauer/SS3")
ssnams = c("Shark1", "Shark2", "Shark3", "Shark4")
nss = length(ssnams)

ios = list()
for(dd in 1:nss)ios[[dd]] = all_ss3(dir=paste0(ssfold,ssdirs[dd]))
names(ios) = ssnams
lapply(ios,ss_names) # what are primary fleets / surveys
       
Fnams = list(c("F4_JPN_LL","F8_ESP_LL"),   c("F1_EU-ESP","F2_JPN"),       c("FS2_BRA","FS4_JPN"),             c("F_Tar_NZ","F_Tar_EU"))
Inams = list(c("S2_JPN_LATE","S4_EU_ESP"), c("S1_ESP-LL-N","S2_JP-LL-N"), c("BRA_index_TB2","JPN_index_TB2"), c("S_Tar_NZ","S_Tar_EU"))

for(dd in 1:nss){
  temp = SS_2_ET(ios[[dd]], Fnam =Fnams[[dd]], Inam = Inams[[dd]])
  save(temp,file=paste0(tdir,"/data/",ssnams[dd],".rda"))
}




