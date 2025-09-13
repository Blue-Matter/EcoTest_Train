
setwd("C:/GitHub/EcoTest")

files = paste0("Indicator_3/",c("allout_1_200.rds","allout_201_400.rds","allout_401_600.rds","allout_601_800.rds","allout_801_1000.rds"))
sizelim = 100
crat = 220/170 # approximate compression ratio
fileno = 0
for(ff in 1:length(files)){
  datlist = readRDS(files[[ff]])
  sz = as.numeric(object.size(datlist)/crat/1E6)
  npack = ceiling(sz/sizelim)
  nl = length(datlist)
  chunks <- split(1:nl, cut(seq_along(1:nl), npack, labels = FALSE))
  nc = length(chunks)
  for(cc in 1:nc){
    fileno=fileno+1
    temp = datlist[chunks[[cc]]]
    saveRDS(temp,paste0("Indicator_3/allout_github_",fileno,".rds"))
  }
}