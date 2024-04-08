
library(data.table)

setwd("C:/GitHub/Ecotest")
setwd("C:/Users/tcarruth/Documents/GitHub/Ecotest")

source("99_Indicators.R")

# takes all the MSE simulated data
# characterizes by simulation range expansion / contraction and simulates for projections
allout = process_sim_data(MSEdir = "C:/temp/Ecotest/batching/Independent_F",
                          Spatial_file = "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases/All_catch_cropped.csv")



