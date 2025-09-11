
# ETest Dependencies

# EcoTest uses neural networks conditioned on simulated data. In R this is 
# achieved using various software: Python, Tensorflow and Keras


# Necessary cuda version: https://developer.nvidia.com/cuda-11.3.0-download-archive


# These libraries can be a total PITA to install so I apologize in advance
install.packages('remotes')
library(remotes)

# 1) Python (language for compiling Tensorflow neural networks)
install.packages('reticulate') # install_github("rstudio/reticulate") 
reticulate::install_python()
reticulate::install_miniconda(force=T)


# 2) Tensorflow (libraries for conducting neural network training)
install.packages('tensorflow') # install_github("rstudio/tensorflow")
# https://storage.googleapis.com/tensorflow/versions/2.20.0/tensorflow_cpu-2.20.0-cp313-cp313-manylinux_2_17_x86_64.manylinux2014_x86_64.whl

# 3) Keras (high level API for managing network training)
install.packages('keras3') #install_github("rstudio/keras3")
# < restart R >
keras3::install_keras(tensorflow = "cpu")



# If you are having difficulties with the above try installing the latest version
# of Microsoft Visual C++  (this is needed to install tensorflow correctly below)

"https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170"



# redundant code
reticulate::install_miniconda()

install.packages('tfdatasets')   #  Has code for installing tensorflow
tfdatasets::install_tensorflow() #  Libraries for training neural networks

install.packages('tensorflow')    # libraries for training neural networks
tensorflow::install_tensorflow(version="cpu") # cpu edition install

# Nvidia cuDNN (libraries for deep learning)
https://developer.nvidia.com/cudnn-downloads?target_os=Windows&target_arch=x86_64&target_version=11&target_type=exe_local

# Anaconda
https://www.anaconda.com/download/success


# anaconda powershell 
pip uninstall tensorflow
pip cache purge
pip install tensorflow

# Other code:
library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-tensorflow", python = path_to_python)

library(tensorflow)
install_tensorflow(envname = "r-tensorflow")
use_virtualenv("r-tensorflow")
tf$constant("Hello Tensorflow!")


