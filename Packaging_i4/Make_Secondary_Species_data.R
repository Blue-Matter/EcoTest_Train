

# ==== Process BON and F tuna data ===========================================================================

library(MSEtool)
library(dplyr)
library(mvtnorm)
library(parallel)
library(miceadds)
library(data.table)
library(EcoTest)
library(EcoTestData)
library(readxl)
library(writexl)

TD = prepare_TD()

setwd("C:/GitHub/EcoTest_train")
source.all("Source")
setwd("G:/Shared drives/BM shared/1. Projects/EcoTest/SCRS 2026/Small_tunas")

t1 = read.csv("t1.csv",header=F)         # total catch
t2 = read.csv("t2.csv",header=T)         # size frequency
t2ce = read.csv("t2ce_LL_curated.csv",header=T)  # catch and effort


# ================== FRIGATE TUNA ============================================================================

dat_FRI = procICCATdat("FRI",t1, t2, t2ce, tF = "EU.ESP-ES-ETRO", tG = "PS" )

png("G:/Shared drives/BM shared/1. Projects/EcoTest/Workshop 2026/Figures/FRI_inputs.png",res=400,width=8,height=6.5,units="in")
  makeXL(xlfile = "C:/GitHub/EcoTest_Train/Indicator_4/Blank_data/EcoTest_Input.xlsx",
         tofile = "C:/GitHub/EcoTest_Train/Indicator_4/Secondary_species/FRI.xlsx",
         spec = "Frigate Tuna",
         L50 = 29.0,
         Linf = 51.47,   # Grudtsev and Korolevich, 1986
         M = 0.48,
         K = 0.32, # Grudtsev and Korolevich, 1986
         C_T = dat_FRI$Catch,
         CAL_T = dat_FRI$CAL,
         CAL_mids = dat_FRI$Lbins,
         I_T = NA,
         L5_T = 32,
         LFS_T = 42,
         VML = 0.975,
         yrs = dat_FRI$yrs,
         plot=T)
dev.off()


# === LITTLE TUNNY ==============================================================================================

dat_LTA = procICCATdat(SpecCode = "LTA",t1, t2, t2ce, tF = "SEN-SN-Art", tG = "HAND", model = 'log(CPUE)~Y + Q + F')

png("G:/Shared drives/BM shared/1. Projects/EcoTest/Workshop 2026/Figures/LTA_inputs.png",res=400,width=8,height=6.5,units="in")
  makeXL(xlfile = "C:/GitHub/EcoTest_Train/Indicator_4/Blank_data/EcoTest_Input.xlsx",
         tofile = "C:/GitHub/EcoTest_Train/Indicator_4/Secondary_species/LTA.xlsx",
         spec = "Little Tunny",
         L50 = 34.4,     # Cruz-Cástan et al., 2019, #
         Linf = 76.45,   # da Silva et al. 2024
         M = 0.392,      # El-Haweet et al. (2013)
         K = 0.553,      # da Silva et al. 2024
         C_T = dat_LTA$Catch,
         CAL_T = dat_LTA$CAL,
         CAL_mids = dat_LTA$Lbins,
         I_T = dat_LTA$Index,
         L5_T = 30,
         LFS_T = 42,
         VML = 0.975,
         yrs = dat_LTA$yrs,
         plot=T)
dev.off()


# === BLF BLACKFIN TUNA  ==============================================================================================

dat_BLF= procICCATdat("BLF",t1, t2, t2ce, tF = "VEN", tG = "PS")

png("G:/Shared drives/BM shared/1. Projects/EcoTest/Workshop 2026/Figures/BLF_inputs.png",res=400,width=8,height=6.5,units="in")
  makeXL(xlfile = "C:/GitHub/EcoTest_Train/Indicator_4/Blank_data/EcoTest_Input.xlsx",
         tofile = "C:/GitHub/EcoTest_Train/Indicator_4/Secondary_species/BLF.xlsx",
         spec = "Blackfin Tuna",
         L50 = 39,       # Valle-Gomez, 1992
         Linf = 82.4,    # Gutierrez et al. 2023
         M = 0.467,      # Gutierrez et al. 2023
         K = 0.365,      # Gutierrez et al. 2023
         C_T = dat_BLF$Catch,
         CAL_T = dat_BLF$CAL,
         CAL_mids = dat_BLF$Lbins,
         I_T = dat_BLF$Index,
         L5_T = 45,
         LFS_T = 52,
         VML = 0.975,
         yrs = dat_BLF$yrs,
         plot=T)
dev.off()




# ============================================================================================================================


# Add data objects to package


XLfiles = list.files("C:/GitHub/EcoTest_Train/Indicator_4/Secondary_species")




