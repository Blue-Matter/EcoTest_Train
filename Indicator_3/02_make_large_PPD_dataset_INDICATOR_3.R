
library(openMSE)
library(dplyr)
library(mvtnorm)
library(parallel)
library(miceadds)
library(SimDesign)
library(data.table)

setwd("C:/GitHub/EcoTest_train")
source.all("Source")


# --- Processes simulated data -------------------------------------------------------------

# prediction is 60 minutes for workstation (48 x 1.0)

# 45600 sims

MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind3/MMSE"
system.time({allout = process_sim_data_4(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_4/Allout.rds")


setwd("C:/GitHub/EcoTest_train")

#files = paste0("Indicator_3/",c("allout_1_200.rds","allout_201_400.rds","allout_401_600.rds","allout_601_800.rds","allout_801_1000.rds"))

datlist = readRDS(paste0("Indicator_3/New_allout_2.rds"))
TD = rbindlist(datlist)

sizelim = 100
crat = 1 # approximate compression ratio
fileno = 0

sz = as.numeric(object.size(TD)/crat/1E6)
npack = ceiling(sz/sizelim)
nl = nrow(TD)
chunks <- split(1:nl, cut(seq_along(1:nl), npack, labels = FALSE))
nc = length(chunks)
for(cc in 1:nc){
  fileno=fileno+1
  temp = TD[chunks[[cc]],]
  saveRDS(temp,paste0("Indicator_3/TD_",fileno,".rds"))
}



# END of Train files






# Second 200 runs (9600 sims)
MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_201_400"
system.time({allout = process_sim_data_3(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_3/allout_201_400.rds")

# third 200 runs (9600 sims)
MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_401_600"
system.time({allout = process_sim_data_3(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_3/allout_401_600.rds")

# fourth 200 runs (9600 sims)
MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_601_800"
system.time({allout = process_sim_data_3(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_3/allout_601_800.rds")

# fifth 200 runs (9600 sims)
MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_801_1000"
system.time({allout = process_sim_data_3(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_3/allout_801_1000.rds")

# === End of script =============================================================

# === Testing ==================================================================


files = list.files(MSEdir)
keep = grepl("MMSE",files)
filelocs = list.files(MSEdir,full.names=T)[keep]
MMSE_1 = readRDS(filelocs[1])
MMSE_2 = readRDS(filelocs[2])

head(MMSE_1@OM[[1]][[1]])[,1:8]
head(MMSE_2@OM[[1]][[1]])[,1:8]
