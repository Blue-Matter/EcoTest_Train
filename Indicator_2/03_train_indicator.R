

# --- Prerequisites ----------------------------------------------------------------------

library("Rcpp")
library(processx)
library(keras3)
library(tensorflow)
library(reticulate)
library(dplyr)
library(tfdatasets)
library(ggplot2)
library(progress)
library(data.table)
library(miceadds)

# NOTES ------------------------

# Networks will NOT fit if there is a constant value in an input column (same across simulations)


# --- source code -------------------------------------------------------------------------

setwd("C:/GitHub/Ecotest")
source.all("Source")


# --- makes datasets for AI training ------------------------------------------------------


allout_1 = readRDS("Indicator_2/allout_1_200.rds")
allout_2 = readRDS("Indicator_2/allout_201_400.rds")

allout = c(allout_1, allout_2)

# --- individual test --------------------------------------------------------------------

TD = makerawdata_2(allout, sno=1, isBrel=F, 
                   inc_Irel = T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T, 
                   stock_in = 1:3, fleet_in = 1:3)

out = train_NN(TD, c(3, 2), 50, T, NULL)
 

# --- multiple tests ---------------------------------------------------------------------

runs = expand.grid(list(inc_I = c(F, T), inc_CAL = c(F, T), inc_CAA = c(F, T), nstocks = 1:3, nfleets = 1:3))
saveRDS(runs,"Indicator_2/Fitted_Models/runs.rds")

for(i in 1:nrow(runs)){
  TD = makerawdata_2(allout, sno=1, isBrel=F, inc_Irel = F, inc_CR = T,
                     inc_I = runs$inc_I[i], inc_CAL = runs$inc_CAL[i], inc_CAA = runs$inc_CAA[i], 
                     stock_in = 1:runs$nstocks[i], fleet_in = 1:runs$nfleets[i])
  out = train_NN(TD, c(3, 2), 50, T, NULL, model_savefile = paste0("Indicator_2/Fitted_Models/Model_",i,".keras"))
  saveRDS(out,paste0("Indicator_2/Fitted_Models/Fit_",i,".rds"))
  print(i)
}

pred_plot(out)




# model =  load_model(paste0("Indicator_2/Fitted_Models/Model_",i,".keras"))
# model.save('C:/savepath/savename.hdf5')
# summarize results

  
  
# === Disused code =============================================================
  
# --- individual species / fleet run (tripling data size) ----------------------

TD1 = makerawdata_2(allout, sno=1, isBrel=F, inc_Irel = T, inc_I = F, inc_CR = T, inc_CAL = T, inc_CAA = F, stock_in = 1, fleet_in = 1)
TD2 = makerawdata_2(allout, sno=2, isBrel=F, inc_Irel = T, inc_I = F, inc_CR = T, inc_CAL = T, inc_CAA = F, stock_in = 2, fleet_in = 2)
TD3 = makerawdata_2(allout, sno=3, isBrel=F, inc_Irel = T, inc_I = F, inc_CR = T, inc_CAL = T, inc_CAA = F, stock_in = 3, fleet_in = 3)
names(TD2) = names(TD3) = names(TD1)
TD = rbind(TD1,TD2,TD3)

# --- load.save ----------------------------------------------------------------

save_model_weights_tf(model, './checkpoints/my_checkpoint')

# Create a new model instance
model <- create_model()

# Restore the weights
load_model_weights_tf(model, './checkpoints/my_checkpoint')




  #save_model_hdf5(AIE, "AIE.hdf5", overwrite=T,include_optimizer = FALSE)
  #testAIE <-load_model_tf("AIE.hdf5" )
  
  
  save_model_weights_hdf5(AIE, "AIE_wts.h5")
  
}else{
  print("2/4 Using saved trained model weights")
  AIE %>% load_model_weights_hdf5(filepath = "AIE_wts.h5")
  test_predictions <- AIE %>% predict(test_df %>% dplyr::select(-label))
  
  plot(exp(test_df$label),exp(test_predictions[ , 1]),xlab="East obs Bt",ylab="East Pred Bt",col="blue"); lines(c(0,1E10),c(0,1E10),col='#ff000050',lwd=2)
  
  
}