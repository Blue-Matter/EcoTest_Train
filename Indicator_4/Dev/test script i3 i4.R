
# ===========================================================================================

library(EcoTest)
library(EcoTestData)
TD = prepare_TD()

alldats = findy("ETData")
isassessed = sapply(alldats,function(x){y=get(x);("Brel" %in% names(y))})



# === Assessed Species ======================================================================

dats = alldats[isassessed]

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
  inds[[dd]] = train_ind_3(data, TD, c(80,40), nepoch = 25)
  cat(paste0("MAE: ",round(inds[[dd]]$MAE,3)," \n"))
}

# Indicator 3 predictions
preds = list()
for(dd in 1:length(dats)){
  preds[[dd]] = pred_ind(Ind = inds[[dd]])
}

# Do retro plots
par(mfrow=c(3,3),par(mai=c(0.8,0.8,0.2,0.05)))
for(dd in 1:length(preds)){
  retro_ind(preds[[dd]])
  mtext(dats[dd],3)
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


