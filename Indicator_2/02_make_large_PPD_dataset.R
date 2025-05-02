
library(openMSE)
library(dplyr)
library(mvtnorm)
library(parallel)
library(miceadds)
library(SimDesign)
library(data.table)

setwd("C:/GitHub/EcoTest")
source.all("Source")

MSEdir = "C:/Users/tcar_/Dropbox/temp/Ecotest/Ind2/MMSE"

files = list.files(MSEdir)
keep = grepl("MMSE",files)
filelocs = list.files(MSEdir,full.names=T)[keep]
MMSE_1 = readRDS(filelocs[1])
MMSE_2 = readRDS(filelocs[2])

head(MMSE_1@OM[[1]][[1]])[,1:8]

head(MMSE_2@OM[[1]][[1]])[,1:8]



# --- Processes simulated data -------------------------------------------------------------

# characterizes by simulation range expansion / contraction and simulates for projections
system.time({allout = process_sim_data(MSEdir = "C:/temp/Ecotest/batching/Dependent_F",
                                       parallel=T, cores = parallel::detectCores())})

saveRDS(allout,"Indicator/Processed_data_stoch_DependentF.rds")



