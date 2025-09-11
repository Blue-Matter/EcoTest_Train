
  install.packages("tensorflow")
  install.packages("keras")
  install.packages("Rcpp")
  
  remotes::install_github("rstudio/reticulate") 
  remotes::install_github("rstudio/tensorflow")
  remotes::install_github("rstudio/keras")
  
  reticulate::install_python()
  reticulate::install_miniconda()
  
  # install python 64 bit
  
  tensorflow::install_tensorflow(version="cpu")
  
  keras3::install_keras(tensorflow = "cpu")
  
  # might need to get tfdatasets as a tarball
  install.packages('tfdatasets')
  install.packages('progress')
  
  
  
  remove.packages("tensorflow")
  remove.packages("reticulate")
  remove.packages("keras")
  remove.packages("keras3")
  remove.packages('tfdatasets')
  