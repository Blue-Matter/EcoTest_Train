

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

# NOTES ------------------------

# Networks will NOT fit if there is a constant value in an input column


# --- source code -------------------------------------------------------------------------

setwd("C:/GitHub/Ecotest")
source('99_neural_net.R')


# --- makes datasets for AI training ------------------------------------------------------

#  1   2   3   4   5   6
# BET SWO BSH SMA WHM BUM  BSBSWB

allout = readRDS("Indicator/Processed_data_stoch_2.rds")
TD = makerawdata(allout, sno=6, isBrel=F, inc_spat =T, inc_Irel = T)

nr<-nrow(TD)
nc<-ncol(TD)
ind<-(1:nr)%in%sample(1:nr,floor(nr*0.05),replace=FALSE)

train<-TD[!ind,2:nc]
mu = apply(train,2,mean)
sd = apply(train,2,sd)

train = scale(train,center=mu,scale=sd)
train_target<-TD[!ind,1]

testy<-TD[ind,2:nc]
testy=scale(testy,center=mu,scale=sd)

testy_target<-TD[ind,1]

model = keras_model_sequential()

model %>%
  layer_dense(units = 8, activation = "relu", input_shape = c(nc-1)) %>%
  layer_dense(units = 4, activation = "relu") %>%
  layer_dense(units = 1)
  
  
model %>%
  compile(
    loss = 'MSE', #mean_squared_error',#"mse",#mae",#"mse",
    optimizer =  optimizer_adam(), #'adam',#optimizer_rmsprop(),'rmsprop',
    metrics = "MAE"
  )
 
nepoch = 20 

history <- model %>% fit(train, train_target,
  epochs = nepoch,
  validation_split = 0.2,
  verbose = 2
)

pred = exp((model %>% predict(testy))[,1])
sim = exp(testy_target)
NN_fit(sim, pred, history,lev=c(0.5,1), addpow=T)


# some figures

Figloc = "C:/GitHub/EcoTest/Figures/Presentation 1 April 24"

jpeg(paste0(Figloc,"/Importance_BSH.jpg"),res=400,units="in",width=9,height=6)

calc_importance(model,testy,adj=0.25,barno=30)

dev.off()

pred = exp((model %>% predict(testy))[,1])
sim = exp(testy_target)

jpeg(paste0(Figloc,"/NN_training_BSH.jpg"),res=400,units="in",width=9,height=5)
  plot(history)
dev.off()

jpeg(paste0(Figloc,"/NN_fit_WHM_noCPUE_nospat.jpg"),res=400,units="in",width=6,height=5)
  NN_fit(sim, pred, history,lev=c(0.5,1), addpow=T)
  mtext("White Marlin (excluding CPUE and spatial model)",line=0.8)
dev.off()

jpeg(paste0(Figloc,"/NN_fit_BSH.jpg"),res=400,units="in",width=6,height=5)
  NN_fit(sim, pred, history,lev=c(0.5,1), addpow=T)
  mtext("Blue shark (including CPUE and spatial model)",line=0.8)
dev.off()



plot(history)  
saveRDS(history,'Ehistory_26_26.rda')
  
  
  
  
  
  
  
  
  
# vvvvvvvv Disused code vvvvvvvvvvvvvvvvvvv
  
  #save_model_hdf5(AIE, "AIE.hdf5", overwrite=T,include_optimizer = FALSE)
  #testAIE <-load_model_tf("AIE.hdf5" )
  
  
  save_model_weights_hdf5(AIE, "AIE_wts.h5")
  
}else{
  print("2/4 Using saved trained model weights")
  AIE %>% load_model_weights_hdf5(filepath = "AIE_wts.h5")
  test_predictions <- AIE %>% predict(test_df %>% dplyr::select(-label))
  
  plot(exp(test_df$label),exp(test_predictions[ , 1]),xlab="East obs Bt",ylab="East Pred Bt",col="blue"); lines(c(0,1E10),c(0,1E10),col='#ff000050',lwd=2)
  
  
}