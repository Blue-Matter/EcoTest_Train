


# --- Prerequisites ------------------------------------------------------------------------

library(data.table)
setwd("C:/GitHub/Ecotest") #setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")
source("99_Indicators.R")


# --- Processes simulated data -------------------------------------------------------------

spat_mods = readRDS("Spatial Models/spat_mods.rds")

# takes all the MSE simulated data
# characterizes by simulation range expansion / contraction and simulates for projections
allout6 = process_sim_data(MSEdir = "C:/temp/Ecotest/batching/Independent_F",
                          spat_mods, parallel=F, cores = 10)

saveRDS(allout6,"Indicator/Processed_data_stoch.rds")







