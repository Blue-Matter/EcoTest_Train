
ConstE_MMP <- function(x, DataList, reps = 1, Emult = 1,...) {
  
  np <- length(DataList)
  nf <- length(DataList[[1]])
  RecList <- lapply(1:np, function(p) replicate(nf, new("Rec")))
  
  for(p in 1:np) { 
    for(f in 1:nf) { # Specify relative F 
      RecList[[p]][[f]]@Effort <- Emult
    }
  }
  
  return(RecList)
}
class(ConstE_MMP) = "MMP"
