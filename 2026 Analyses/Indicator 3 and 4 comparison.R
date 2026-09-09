
# ===========================================================================================

library(EcoTest)
library(EcoTestData)
TD = prepare_TD()

alldats = findy("ETData")
isassessed = sapply(alldats,function(x){y=get(x);("Brel" %in% names(y))})

# ===== INDICATOR 3 =========================================================================

# === Assessed Species ======================================================================

dats = alldats[isassessed]

# ---- Data in range? -----------------------

for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  data = get(dats[dd])
  chk = check_train_data(data, TD)
  print(chk$fail)
}

# ---- Indicator 3 training ----------------

mod_save_files = paste0("C:/GitHub/EcoTest_Train/2026 Analyses/i3_",dats,".keras")

inds3 = list()
for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  data = get(dats[dd])
  inds3[[dd]] = train_ind_3(data, TD, c(80,40), nepoch = 35, lr = 0.003,
                            model_savefile = mod_save_files[dd], name = dats[dd])
  cat(paste0("MAE: ",round(inds3[[dd]]$MAE,3)," \n"))
}

saveRDS(inds3, "C:/GitHub/EcoTest_Train/2026 Analyses/i3_indicators_assessed.rds" )



# === Non-Assessed Species =============================================================

dats = alldats[!isassessed]

# Indicator 3 training
inds = list()
for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  data = get(dats[dd])
  inds[[dd]] = train_ind_3(data, TD, c(20,10), nepoch = 15)
  cat(paste0("MAE: ",round(inds[[dd]]$MAE,3)," \n"))
}



# ===== INDICATOR 4 =========================================================================

# === Assessed Species ======================================================================

dats = alldats[isassessed]
mod_save_files = paste0("C:/GitHub/EcoTest_Train/2026 Analyses/i4_",dats,".keras")

# training

for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  input = get(dats[dd])
  ind = train_ind_4(input, TD, nodes = c(50,20,60), nepoch = 35, nbatch=8, lr = 0.004,
                    model_savefile = mod_save_files[dd], name = dats[dd])

  saveRDS(ind, paste0("C:/GitHub/EcoTest_Train/2026 Analyses/i4_",dats[dd],".rds" ))

  cat(paste0("MAE: ",round(ind$MAE,3)," \n"))
}



# === Non-Assessed Species =============================================================

dats = alldats[!isassessed]

# Data in range?
for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  data = get(dats[dd])
  chk = check_train_data(data, TD)
  print(chk$fail)
}

# Indicator 3 training
inds = list()
for(dd in 1:length(dats)){
  cat("---------------- \n")
  cat(paste0(dats[dd]," \n"))
  data = get(dats[dd])
  inds[[dd]] = train_ind_3(data, TD, c(20,10), nepoch = 15)
  cat(paste0("MAE: ",round(inds[[dd]]$MAE,3)," \n"))
}

# Indicator 3 predictions
preds = list()
for(dd in 1:length(dats)){
  preds[[dd]] = pred_ind(Ind = inds[[dd]])
}

# Do pred plots
par(mfrow=c(2,2),par(mai=c(0.8,0.8,0.2,0.05)))
for(dd in 1:length(preds)){
  plot_ind(preds[[dd]], F, 2024)
  mtext(dats[dd],3)
}


