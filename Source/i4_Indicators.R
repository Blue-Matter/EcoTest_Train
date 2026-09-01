
process_sim_data_4 = function(MSEdir, parallel=T, cores = NA){

  files = list.files(MSEdir)
  keep = grepl("MMSE",files)
  filelocs = list.files(MSEdir,full.names=T)[keep]
  nfile = length(filelocs)

  if(parallel){
    library(snowfall)
    library(parallel)
    if(is.na(cores))cores = detectCores()/2
    sfInit(parallel=T,cpus = cores)
    sfExport(list = c("proc_dat_LH","proc_dat_F3","proc_allFdat_F3","smooth2","interpolate","slp3","smooth3","calcTsel","intselpars"))
    allout = sfLapply(1:nfile,get_sim_data_3,filelocs=filelocs)
    # test = sfLapply(1, get_sim_data_3, filelocs=filelocs)
  }else{
    allout = lapply(1:nfile, get_sim_data_3, filelocs=filelocs)
    # test = lapply(1, get_sim_data_3, filelocs=filelocs, quiet=F)
  }

  cat("\n")
  allout

}



