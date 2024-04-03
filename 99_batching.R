# Batching functions for generating large simulated datasets



solveforR1<-function(v2,y2){
  y1 = (-1+(1 + 8 * max(y2, -1/8) )^0.5)/2
  v1 = 2*v2/(1+y1)
  return(c(v1=v1,y1=y1))
}



runbatch = function(x,MOM,doPE=T){
  
  temp = MOM
  temp@seed = x
  proyears = MOM@proyears
  nstock = length(MOM@cpars)
  nfleet = length(MOM@cpars[[1]])
 
  
  if(doPE){
    for(ss in 1:nstock){
      
      Perr = MOM@cpars[[ss]][[1]]$Perr_y
      nPE = ncol(Perr)
      edind = nPE - ((proyears-1):0)
      preind = (1:nPE)[!((1:nPE)%in%edind)]
      #AC = apply(Perr[,edind],1,function(x)acf(x,1,plot=F)[[1]][2,1,1])
      
      
      MOM@cpars[[ss]][[1]]$
      
    
  }
  
}


