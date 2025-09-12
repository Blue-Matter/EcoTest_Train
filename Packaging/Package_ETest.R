

library(keras3)
library(r4ss)
library(data.table)
library(miceadds)

fdir = "C:/GitHub/EcoTest"
tdir = "C:/GitHub/ETest"
source.all(paste0(fdir,"/Source"))


# Run numbers
#model_index = readRDS(paste0(fdir,"/Indicator_2/Fitted_Models/runs.rds"))
#save(model_index,file = paste0(tdir,"/data/model_index.rda"))

# Fitted models
#for(x in 1:nrow(model_index)){
#  model = load_model(paste0(fdir,"/Indicator_2/Fitted_Models/Model_",x,".keras"))
#  save_model(model,paste0(tdir,"/data/Model_",x,".keras"))
#}


# Testing data

allout_1 = readRDS(paste0(fdir,"/Indicator_2/allout_1_200.rds"))
allout_2 = readRDS(paste0(fdir,"/Indicator_2/allout_201_400.rds"))
allout_3 = readRDS(paste0(fdir,"/Indicator_2/allout_401_600.rds"))
allout = c(allout_1, allout_2, allout_3)

# Logged data
TD = makerawdata_2(allout, sno=1, isBrel=F,  inc_Irel = T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T,  stock_in = 1:3, fleet_in = 1:3, Bmin = 0.05)

save(TD,file=paste0(tdir,"/data/TD.rda"))


# indicator 3
cr = 1.5
allout3 = list()
for(i in 1:8){
  temp = readRDS(paste0(fdir,"/Indicator_3/allout_github_",i,".rds"))
  allout3 = c(allout3,temp)
}
TD = makerawdata_3(allout3)


# - code

files = c("i2_neural_net.R", "i2_SS_funcs.R","99_neural_net.R","99_Indicators.R","i2_data_processing.R")
for(ff in 1:length(files)){
  file.copy(paste0(fdir,"/Source/",files[ff]),paste0(tdir,"/R/",files[ff]),overwrite = T)
}

# - Data objects

ssfold = c(rep("G:/Shared drives/BM shared/1. Projects/TOF Advisory/Blue_Shark_MSE/Assessment_files_provided/",5),
           rep("G:/Shared drives/BM shared/1. Projects/EcoTest/Assessments/",2))
           
ssdirs = c("IO_Joel_Rice/SS3","NAtl_Nathan_Taylor/SS3","SAtl_Nathan_Taylor/SS3","SWPac_Philipp_Neubauer/SS3","NP_Nicholas_D_Barth/SS3","SALB","NSWO")
ssnams = c("Shark_1", "Shark_2", "Shark_3", "Shark_4","Shark_5","Tuna_1","Billfish_1")
nss = length(ssnams)

ios = list()
for(dd in 1:nss)ios[[dd]] = all_ss3(dir=paste0(ssfold[dd],ssdirs[dd]))
names(ios) = ssnams

# you were here

for(dd in 1:nss){
  temp = ios[[dd]]
  objname = paste0(ssnams[dd],"_io")
  assign(objname,temp)
  to_file=paste0(tdir,"/data/",ssnams[dd],"_io.rda")
  do.call(save, list(objname,file=to_file))
}


lapply(ios,ss_names) # what are primary fleets / surveys
       
Fnams = list(c("F4_JPN_LL","F8_ESP_LL"),   c("F1_EU-ESP","F2_JPN"),       
             c("FS2_BRA","FS4_JPN"),             c("F_Tar_NZ","F_Tar_EU"),
             c("F1_MEX","F19_TAIW_LG"),
             c("Fleet_01","Fleet_02"),  c("US_2","PORT_6"))
Inams = list(c("S2_JPN_LATE","S4_EU_ESP"), c("S1_ESP-LL-N","S2_JP-LL-N"), 
             c("BRA_index_TB2","JPN_index_TB2"), c("S_Tar_NZ","S_Tar_EU"),
             c("S10_MEX","S3_TAIW_LG"),
             c("CTP-LL", "JPN-LL3"), c("US_Survey_12","PORT_Survey_13"))

for(dd in 1:nss){
  temp = SS_2_ET(io=ios[[dd]], Fnam = Fnams[[dd]], Inam = Inams[[dd]])
  objname = paste0(ssnams[dd],"_data")
  assign(objname,temp)
  to_file=paste0(tdir,"/data/",ssnams[dd],"_data.rda")
  do.call(save, list(objname,file=to_file))
}

for(dd in 1:nss){
  temp = SS_2_ET_Retro(io=ios[[dd]], Fnam = Fnams[[dd]], Inam = Inams[[dd]], npeels=8)
  objname = paste0(ssnams[dd],"_retro")
  assign(objname,temp)
  to_file=paste0(tdir,"/data/",ssnams[dd],"_retro.rda")
  do.call(save, list(objname,file=to_file))
}

# Example data object Some_data

MMSE = readRDS("C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_1_200/MMSE_1.rda")
PPD = MMSE@PPD; PD = PPD[[1]][[1]][[1]]; simno = 1; ay = 70
Catch = getdat(PPD, "Cat", simno)
Index = getdat(PPD, "VInd", simno)
ML = getdat(PPD, "ML",simno)
K = getsdat(PPD, "vbK",simno)
Linf = getsdat(PPD, "vbLinf",simno)
L50 = getsdat(PPD, "L50",simno)
L5 = getfdat(PPD,"LFC",simno)
LFS = getfdat(PPD,"LFS",simno)
VML = getfdat(PPD,"Vmaxlen",simno)
M = getsdat(PPD,"Mort",simno)
maxa = -log(0.05)/M

Some_data  = list(Catch = Catch, Index = Index, Mean_length = ML, K = K, 
                     Linf = Linf, L50 = L50, L5 = L5, LFS = LFS, VML = VML, 
                     M = M, maxa = maxa)
save(Some_data,file = paste0(tdir,"/data/Some_data.rda"))

  
  
  


