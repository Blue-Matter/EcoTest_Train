



makerawdata_4 = function(cdat, Brange = c(0.025,4)){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; inc_CAL = T; inc_CAA = T
  #cdat = as.data.frame(rbindlist(allout3))
  cnams = names(cdat)
  s1dat = cdat[,grepl("_s1",cnams)]
  s2dat = cdat[,grepl("_s2",cnams)]
  s3dat = cdat[,grepl("_s3",cnams)]
  colnms = sapply(names(s1dat),function(x)stringr::str_remove(x,"_s1"))
  names(s1dat) = names(s2dat) = names(s3dat) = colnms
  dat = rbind(s1dat,s2dat,s3dat)
  dat = dat[,!(names(dat)%in%c("CF_T"))]

  dat = log(dat)
  dat = cleandat_2(dat,resprange = log(Brange))                       # clean NAs and Infs
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat

}




# gets Time series data from TD for total fleet
getTS = function(TD, code = "T"){
  namTD = names(TD)
  codey = paste0("_",code)
  Cs = TD[,grepl("C_",namTD)&grepl(codey,namTD)]
  Is = TD[,grepl("I_",names(TD))&grepl(codey,names(TD))]
  MLs = TD[,!grepl("Linf",names(TD))&!grepl("VML_",names(TD))&grepl("ML_",names(TD))&grepl(codey,names(TD))]
  MAs = TD[,grepl("MA_",names(TD))&grepl(codey,names(TD))]
  FMs = TD[,grepl("FM_",names(TD))&grepl(codey,names(TD))]
  cbind(Cs, Is, MLs, MAs, FMs)
}



train_NN_TS = function(TD, nodes = c(20, 4, 10), nepoch = 20, plot = T, model=NULL, model_savefile=NA,
                    validation_split = 0.1, test_split = 0.1, lev = c(0.5, 1), nbatch = 4){


  Res = TD[,1]

  # ---- Time series inputs ----------

  TS_T = getTS(TD,"T")
  Cs = TS_T[,grepl("C",names(TS_T))]
  TS_f1 = getTS(TD,"f1")
  TS_f2 = getTS(TD,"f2")
  TS_f3 = getTS(TD,"f3")

  allTS = cbind(TS_T, TS_f1, TS_f2, TS_f3)

  window_size = ncol(Cs)
  num_samples <- nrow(Cs)
  n_ts_feat <- ncol(allTS)/window_size

  # Initialize arrays
  tsf <- array(0, dim = c(num_samples, window_size, n_ts_feat))
  for(ff in 1:n_ts_feat){
    ind = ((ff-1) * window_size) + (1:window_size)
    tsf[,,ff] = as.matrix(allTS[,ind])
  }


  # ---- Static Inputs -------------

  sfind = match(c("K","M_K","maxa","L50_Linf","CF_f1","CF_f2","CF_f3",paste(c("L5_L50", "LFS_L50","VML"),rep(c("T","f1","f2","f3"),each=3),sep="_")),names(TD))
  sf = TD[,sfind]
  n_s_feat = ncol(sf)


  # ---- Split data ----------------

  tind <- 1:(floor((1 - test_split) * num_samples))
  vind <- (floor((1 - test_split) * num_samples) + 1):num_samples

  tsf_t = tsf[tind,,]
  tsf_v = tsf[vind,,]

  sf_t = as.matrix(sf[tind,])

  mu = apply(sf_t, 2, mean)
  sd = apply(sf_t, 2, sd)
  sf_t = scale(sf_t, center=mu, scale=sd)
  sf_v = as.matrix(sf[vind,])

  r_t = Res[tind]
  r_v = Res[vind]


  # Build model --------------------------------------------------------

  # Time series features input layer

  ts_input <- layer_input(shape = c(window_size, n_ts_feat),name="ts_input")
  ts_out <- ts_input %>% layer_lstm(units = nodes[1])

  # Static features input layer

  static_input <- layer_input(shape = c(n_s_feat), name = "static_input")
  static_out <- static_input %>%  layer_dense(units = nodes[2], activation = "relu")

  # Combined layer

  combined <- layer_concatenate(list(ts_out, static_out))

  # Network hidden layers

  output <- combined %>%
    layer_dense(units = nodes[3], activation = "relu") %>%
    layer_dense(units = 1)

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
    epochs = nepoch,
    validation_split = validation_split,
    batch_size = ceiling(length(r_t)/nbatch)
  )

  pred = exp((model %>% predict(list(tsf_v, sf_v)))[,1])
  sim = exp(r_v)
  tab = NN_fit(sim, pred, history, lev=lev, addpow=T,nepoch=nepoch,plot=F)

  MAE = history$metrics$val_mae[nepoch]

  if(!is.na(model_savefile)) save_model(model,filepath = model_savefile, overwrite=T)

  outlist=list(MAE = MAE, sim=sim, pred=pred, tab = tab, model=model, train=list(sf = sf, tsf = tsf),
               train_target=r_t, nepoch = nepoch, r2 = cor(sim,pred)^2, history = history,
               lev=lev, mu = mu, sd = sd)

  if(plot)(pred_plot(outlist))

  outlist

}

