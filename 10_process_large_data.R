


# --- Prerequisites ------------------------------------------------------------------------

library(data.table)
setwd("C:/GitHub/Ecotest") #setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")
source("99_Indicators.R")


# --- Processes simulated data -------------------------------------------------------------

spat_mods = readRDS("Spatial Models/spat_mods.rds")

# takes all the MSE simulated data
# characterizes by simulation range expansion / contraction and simulates for projections
allout2 = process_sim_data(MSEdir = "C:/temp/Ecotest/batching/Independent_F",
                          spat_mods, parallel=F, cores = 5)

saveRDS(allout2,"Indicator/Processed_data_2.rds")







