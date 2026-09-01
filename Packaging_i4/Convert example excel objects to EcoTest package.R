library(miceadds)
library(readxl)

fdir = "C:/GitHub/EcoTest_train"
tdir = "C:/GitHub/EcoTest"
ddir = "C:/GitHub/EcoTestData"
source.all(paste0(fdir,"/Source"))
source.all(paste0(tdir,"/R"))

xldir = "C:/GitHub/EcoTest_Train/Indicator_4/Real_data"
xlfiles = list.files(xldir)
dirs = paste0(xldir,"/",xlfiles)
nams = sapply(xlfiles,function(x)strsplit(x,".xlsx")[[1]][1])

nsim = 32
retros = 6

for(dd in 1:length(dirs)){
  cat(".")
  out = proc_i4(xlfile = dirs[dd], nsim = nsim, retros = retros)
  assign(nams[dd],out)
  do.call(save, list(nams[dd], file = paste0("C:/GitHub/EcoTest/data/",nams[dd],".rda")))
}
xlfile = dirs[1]

