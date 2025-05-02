



# Prerequisites
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

allout = readRDS("Indicator/Processed_data.rds")
TD = makerawdata(allout, sno=6, isBrel=F, inc_spat = T, inc_Irel = T)
TD[,1] = log(TD[,1])

nr<-nrow(TD)
nc<-ncol(TD)
ind<-(1:nr)%in%sample(1:nr,floor(nr*0.05),replace=FALSE)

train_data<-TD[!ind,2:nc]
train_labels<-TD[!ind,1]
test_data<-TD[ind,2:nc]
test_labels<-TD[ind,1]

# Check for NAs # for(i in 1:ncol(train_data))print(sum(is.na(train_data)))

column_names <- colnames(TD)[2:nc]
column_names = paste0("p",1:(nc-1))

train_df <- train_data %>%
  as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = train_labels)

test_df <- test_data %>%
  as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = test_labels)

spec <- feature_spec(train_df, label ~ . ) %>%
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>%
  fit()

layer <- layer_dense_features(
  feature_columns = dense_features(spec),
  dtype = tf$float32
)

layer(train_df)


# model building function for reuse
build_model <- function(spec, train_df) {
  
  input <- layer_input_from_dataset(train_df %>% dplyr::select(-label))
  
  output <- input %>%
    layer_dense_features(dense_features(spec)) %>%
    layer_dense(units = 8, activation = "relu") %>%
    layer_dense(units = 4, activation = "relu") %>%
    layer_dense(units = 1)
  
  mody <- keras_model(input, output)
  
  mody %>%
    compile(
      loss = 'MSE', #mean_squared_error',#"mse",#mae",#"mse",
      optimizer =  optimizer_adam(), #'adam',#optimizer_rmsprop(),'rmsprop',
      metrics = "MAE"
    )
  
  mody
}

# Train model

print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)

early_stop <- callback_early_stopping(monitor = "val_loss", patience = 20)


mody <- build_model(spec,train_df)


history <- mody %>% fit(
  x = train_df %>% dplyr::select(-label),
  y = train_df$label,
  epochs = 30,
  validation_split = 0.2,
  verbose = 2#,
  #callbacks = list(print_dot_callback)
)

plot(history)

test_predictions <- mody %>% predict(test_df %>% dplyr::select(-label))

plot(test_df$label,test_predictions[ , 1],xlab="SSB/SSBMSY (obs)",ylab="SSB/SSBMSY (pred)"); lines(c(0,1E10),c(0,1E10),col='#ff000050',lwd=2)



  
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