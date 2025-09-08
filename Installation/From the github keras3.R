Keras github installation 


install.packages("keras3")
#remotes::install_github("rstudio/keras3")

tensorflow::install_tensorflow(version="2.20-cpu")

keras3::install_keras(backend = "tensorflow", gpu=F) # did not work
library(tensorflow)
py_require("tensorflow")
tf$constant("Hello TensorFlow!")

install_r_tensorflow(python_path = "C:/Users/User_name/AppData/Local/Programs/Python/Python313/python.exe", env_name = "r-tensorflow")

