
https://www.youtube.com/watch?v=UbpuXWKQJXs&ab_channel=hamannlab

Anaconda 2022.10
Python 3.9.13
cuDNN 8.1.0
CUDAtoolkit 11.2
Tensorflow 2.10.1
R 4.3
Keras 2.13.0

# GPU
python.exe -pip install --upgrad pip
pip install "tensorflow<2.11"
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0

#CPU
python.exe -m pip install--upgrade pip
pip install tensorflow

install.packages("keras")
library(keras)
use_condaenv("base",required=T)

library(tensorflow)
tf$python$client$device$lib$list_local_devices()

ts_config()
system("nvidia-smi")