
  install.packages("tensorflow")
  install.packages("keras")
  install.packages("Rcpp")
  library("Rcpp")
  install.packages("devtools")
  install.packages("reticulate") #devtools::install_github("rstudio/reticulate", force=TRUE)
  devtools::install_github("r-lib/processx")
  library(processx)
  remotes::install_github("rstudio/tensorflow")
  remotes::install_github("rstudio/keras")
  #devtools::install_github("rstudio/keras")
 
  #install_keras(method = c("auto", "virtualenv", "conda"), conda = "auto",  tensorflow = "gpu", extra_packages = NULL)
  
  #install git
  # add git to path
  
  reticulate::install_python()
  reticulate::install_miniconda()
  
  # install python 64 bit
  
  
  library(tensorflow)
  install_tensorflow(version="cpu")
  
  install.packages("keras3")
  library(keras3)
  install_keras(tensorflow = "cpu")
  
  # might need to get tfdatasets as a tarball
  install.packages('tfdatasets')
  install.packages('progress')