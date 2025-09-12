
library(openMSE)
library(dplyr)
library(mvtnorm)
library(parallel)
library(miceadds)
library(SimDesign)
library(data.table)

setwd("C:/GitHub/EcoTest")
source.all("Source")


# --- Processes simulated data -------------------------------------------------------------

# prediction is 60 minutes for workstation (48 x 1.0)  

# First 200 runs (9600 sims)
MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE_1_200"
system.time({allout = process_sim_data_3(MSEdir, parallel=T, cores = parallel::detectCores())})
saveRDS(allout,"Indicator_3/allout_1_200.rds")

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
