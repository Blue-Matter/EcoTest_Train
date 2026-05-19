
library(EcoTest)
library(EcoTestData)


TD = prepare_TD()


# Just time series data Catch and index

Res = TD[,1]
Cs = TD[,grepl("C_",names(TD))&grepl("s1_T",names(TD))]
Is = TD[,grepl("I_",names(TD))&grepl("s1_T",names(TD))]
MLs = TD[,!grepl("Linf",names(TD))&!grepl("VML_",names(TD))&grepl("ML_",names(TD))&grepl("s1_T",names(TD))]

window_size = ncol(Cs)
num_samples <- nrow(Cs) 
n_features <- 3

# Initialize arrays
X_train <- array(0, dim = c(num_samples, window_size, n_features))
y_train <- Res
X_train[,,1] = as.matrix(Cs)
X_train[,,2] = as.matrix(Is)
X_train[,,3] = as.matrix(MLs)


# Split Data into Training and Validation Sets

# Split data
split_index <- floor(0.8 * num_samples)
X_train_split <- X_train[1:split_index, , ]
y_train_split <- y_train[1:split_index]
X_val_split <- X_train[(split_index + 1):num_samples, , ]
y_val_split <- y_train[(split_index + 1):num_samples]

#5. Build and Train the Model
#Define and compile the LSTM model for forecasting.


# Define the LSTM model
model <- keras_model_sequential() %>%
  layer_lstm(units = 20, input_shape = c(window_size, n_features)) %>%
  layer_dense(units = 16) %>%
  layer_dense(units = 1)

# Compile the model
model %>% compile(
  loss = 'mean_squared_error',
  metrics = list("mae"),
  optimizer = optimizer_adam()
)

# Train the model
history <- model %>% fit(
  x = X_train_split,
  y = y_train_split,
  epochs = 200,
  batch_size = 4100,
  validation_data = list(X_val_split, y_val_split)
)

pred = exp((model %>% predict(X_val_split))[,1])
plot(exp(y_val_split), pred,xlab ="Simulated SSB/SSBMSY", ylab = "Predicted SSB/SSBMSY");grid(col="blue")
lines(c(-10,10),c(-10,10),col='red')




# ========= Multi layer approach  =============================================================================================================

# Response variable

Res = TD[,1]

# Time series features --------------------------------------------

Cs = TD[,grepl("C_",names(TD))&grepl("s1_T",names(TD))]
Is = TD[,grepl("I_",names(TD))&grepl("s1_T",names(TD))]
MLs = TD[,!grepl("Linf",names(TD))&!grepl("VML_",names(TD))&grepl("ML_",names(TD))&grepl("s1_T",names(TD))]
MV = TD[,grepl("MV_",names(TD))&grepl("s1_T",names(TD))]
FM = TD[,grepl("FM_",names(TD))&grepl("s1_T",names(TD))]
MA = TD[,grepl("MA_",names(TD))&grepl("s1_T",names(TD))]

window_size = ncol(Cs)
num_samples <- nrow(Cs) 
n_ts_feat <- 6

tsf <- array(0, dim = c(num_samples, window_size, n_ts_feat))
y_train <- Res
tsf[,,1] = as.matrix(Cs)
tsf[,,2] = as.matrix(Is)
tsf[,,3] = as.matrix(MLs)
tsf[,,4] = as.matrix(MV)
tsf[,,5] = as.matrix(FM)
tsf[,,6] = as.matrix(MA)

# Static inputs --------------------------------------------------------


sfind = match(c("K_s1","M_K_s1","maxa_s1","L50_Linf_s1","ML_Linf_s1_T","L5_L50_s1_T", "LFS_L50_s1_T","VML_s1_T"),names(TD))

sf = TD[,sfind]
n_s_feat = ncol(sf)

# Split data 

tind <- 1:(floor(0.9 * num_samples))
vind <- (floor(0.9 * num_samples) + 1):num_samples

tsf_t = tsf[tind,,]
tsf_v = tsf[vind,,]

sf_t = as.matrix(sf[tind,])
sf_v = as.matrix(sf[vind,])

r_t = Res[tind]
r_v = Res[vind]




# Build model --------------------------------------------------------

# Time series features input layer

ts_input <- layer_input(shape = c(window_size, n_ts_feat),name="ts_input")
ts_out <- ts_input %>% layer_lstm(units = 100) 


# Static features input layer

static_input <- layer_input(shape = c(n_s_feat), name = "static_input")
static_out <- static_input %>%  layer_dense(units = 6, activation = "relu")

# Combined layer

combined <- layer_concatenate(list(ts_out, static_out))

# Network hidden layers

output <- combined %>%
  layer_dense(units = 60, activation = "relu") %>%
  layer_dense(units = 1, activation = "linear") 

# Specify model 

model <- keras_model(inputs = c(ts_input, static_input), outputs = output)

# Compile Model

model %>% compile(
  optimizer = "adam",
  loss = "mse",
  metrics = list("mae"),
  run_eagerly = T
)



# Pass both as a list to the 'x' argument
model %>% fit(
  x = list(tsf_t, sf_t),
  y = r_t,
  epochs = 200,
  validation_split = 0.1,
  batch_size = ceiling(length(r_t)/8)
)


pred = exp((model %>% predict(list(tsf_v, sf_v)))[,1])
plot(exp(r_v), pred,xlab ="Simulated SSB/SSBMSY", ylab = "Predicted SSB/SSBMSY",col="#0000ff90",pch=19)
grid(col="red")
lines(c(-10,10),c(-10,10),col='red')






















# Define the LSTM model
model <- keras_model_sequential() %>%
  layer_lstm(units = 20, input_shape = c(window_size, n_features)) %>%
  layer_dense(units = 5) %>%
  layer_dense(units = 1)

# Compile the model
model %>% compile(
  loss = 'mean_squared_error',
  optimizer = optimizer_adam()
)

# Train the model
history <- model %>% fit(
  x = X_train_split,
  y = y_train_split,
  epochs = 20,
  validation_data = list(X_val_split, y_val_split)
)

pred = exp((model %>% predict(X_val_split))[,1])
plot(exp(y_val_split), pred,xlab ="Simulated SSB/SSBMSY", ylab = "Predicted SSB/SSBMSY");grid(col="blue")
lines(c(-10,10),c(-10,10),col='red')




#  ===================== other code

library(keras)
library(tensorflow)

# 1. Define input shapes
batch_size <- 32
timesteps <- 10        # e.g., past 10 days of sequence data
ts_features <- 3       # Number of time-series variables
static_features <- 4   # e.g., categorical embeddings or static variables

# 2. Time Series Branch
ts_input <- layer_input(shape = c(timesteps, ts_features), name = "ts_input")

ts_out <- ts_input %>%
  layer_lstm(units = 64) 

# 3. Static/Other Inputs Branch
static_input <- layer_input(shape = c(static_features), name = "static_input")

static_out <- static_input %>%
  layer_dense(units = 32, activation = "relu")

# 4. Combine and finalize
combined <- layer_concatenate(list(ts_out, static_out))

output <- combined %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "linear") # Predicting a continuous value

# 5. Build and compile the model
model <- keras_model(inputs = c(ts_input, static_input), outputs = output)

model %>% compile(
  optimizer = "adam",
  loss = "mse",
  metrics = list("mae")
)

# --- HOW TO TRAIN ---
# Create dummy data matching shapes
ts_data <- array(runif(100 * timesteps * ts_features), dim = c(100, timesteps, ts_features))
static_data <- array(runif(100 * static_features), dim = c(100, static_features))
target_data <- runif(100)

# Pass both as a list to the 'x' argument
model %>% fit(
  x = list(ts_data, static_data),
  y = target_data,
  epochs = 10,
  batch_size = batch_size
)


