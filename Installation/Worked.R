
install.packages("keras3")
library(keras3)
install_keras()

# check tensorflow is working
tensorflow::tf$constant("Hello World")


# if you want to take data from fitted with SS3 assessments

install.packages("xfun") # do not install the package that required compilation
install.packages("remotes")
remotes::install_github("r4ss/r4ss")  # v1.52.1


reticulate::py_config()
tensorflow::tf_config()
tensorflow::tf$constant("Hello World")


fdir = "C:/GitHub/EcoTest"
model = load_model(paste0(fdir,"/Indicator_2/Fitted_Models/Model_",1,".keras"))

