


# --- Prerequisites ------------------------------------------------------------------------

library(data.table)
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest") #setwd("C:/GitHub/Ecotest") 
source("99_Indicators.R")


# --- Processes simulated data -------------------------------------------------------------

spat_mods = readRDS("Spatial Models/spat_mods.rds")

# takes all the MSE simulated data
# characterizes by simulation range expansion / contraction and simulates for projections
allout8 = process_sim_data(MSEdir = "C:/temp/Ecotest/batching/Dependent_F",
                          spat_mods, parallel=T, cores = 10)

saveRDS(allout8,"Indicator/Processed_data_stoch_DependentF.rds")







