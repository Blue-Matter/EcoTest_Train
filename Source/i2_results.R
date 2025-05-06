
getlab = function(i,runs){
  vec = runs[i,]
  txt = paste(c("no index","index")[as.integer(vec[1])+1],
              c("no lens","lens")[as.integer(vec[2])+1],
              c("no ages","ages")[as.integer(vec[3])+1],
              paste0("nstk = ",vec[4]),
              paste0("nflt = ",vec[5]),sep=", ")
  paste0(i,": ",txt)
}