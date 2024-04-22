

# --- Prerequisites ----------------------------------------------------------------------

library(openMSE)

setwd("C:/Users/tcarruth/DocumentsGitHub/Ecotest")

MMSE = readRDS("C:/temp/Ecotest/batching/Independent_F/MMSE_1.rda")


szs4=function(obj){
  slots = slotNames(obj)
  nslot = length(slots)
  sz = rep(NA,nslot)
  for(i in 1:nslot)   sz[i] = object.size(slot(obj,slots[i]))/1E6
  names(sz)=slots
  sz[order(sz,decreasing=T)]
}


szlist = function(obj){
  sz = sapply(obj,object.size)/1E6
  sz[order(sz,decreasing = T)]
}


szs4(MMSE)

szlist(MMSE@multiHist)
szlist(MMSE@multiHist[[1]])
szs4(MMSE@multiHist[[1]][[1]])
szlist(MMSE@multiHist[[1]][[1]]@SampPars)
szlist(MMSE@multiHist[[1]][[1]]@SampPars$Stock)

szlist(MMSE@multiHist[[1]][[1]]@SampPars$Fleet)
szlist(MMSE@multiHist[[1]][[1]]@SampPars$Obs)

szs4(MMSE@multiHist[[1]][[1]]@Data)
szlist(MMSE@multiHist[[1]][[1]]@Data@Misc)

szlist(MMSE@multiHist[[1]][[1]]@AtAge)

temp= MMSE

MMSE = trim_MMSE(temp)

object.size(MMSE) / object.size(temp) *100 


# overwrite massive files

largedir = "C:/temp/Ecotest/batching/Independent_F"

files= list.files(largedir,full.names = T, include.dirs = T)


sfInit(cpus=parallel::detectCores()/2,parallel=T)
sfLibrary(MSEtool) 
sfExport("trim_MMSE");sfExport("files")

overfunc = function(x,files){
  MMSE = readRDS(files[x])
  saveRDS(trim_MMSE(MMSE),file=files[x])
}

sfSapply(1:length(files),overfunc,files=files)

