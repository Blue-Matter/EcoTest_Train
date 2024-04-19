

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


# --- source code -------------------------------------------------------------------------

setwd("C:/GitHub/Ecotest")
source('99_neural_net.R')


# --- makes datasets for AI training ------------------------------------------------------

allout = readRDS("Indicator/Processed_data_3.rds")
TD = makerawdata(allout, sno=2, isBrel=F, inc_spat = F, inc_Irel = F)
TD[,1] = log(TD[,1])

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
  layer_dense(units = 4, activation = "relu", input_shape = c(nc-1)) %>%
  layer_dense(units = 2, activation = "relu") %>%
  layer_dense(units = 1)
  
  
model %>%
  compile(
    loss = 'MSE', #mean_squared_error',#"mse",#mae",#"mse",
    optimizer =  optimizer_adam(), #'adam',#optimizer_rmsprop(),'rmsprop',
    metrics = "MAE"
  )
 

history <- model %>% fit(train, train_target,
  epochs = 30,
  validation_split = 0.2,
  verbose = 2#,
  #callbacks = list(print_dot_callback)
)

#plot(history)

pred <- model %>% predict(testy)

plot(exp(testy_target),exp(pred[ , 1]),xlab="SSB/SSBMSY (simulated)",ylab="SSB/SSBMSY (pred)"); lines(c(0,1E10),c(0,1E10),col='#ff000050',lwd=2)


power_tab(sim = testy_target, pred = pred)

power_tab(sim, pred, lev = 0.5){
  
  tab = array(NA,c(2,2))
  row.names(tab) = 
  
  
}

  
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