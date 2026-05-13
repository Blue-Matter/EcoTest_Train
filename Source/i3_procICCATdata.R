
procICCATdat = function(SpecCode = "FRI", t1, t2, tF = NA, tG = NA, CALint=T ){
  
  t1s = t1[t1[,1] == SpecCode,c(3,4,5,6)]
  names(t1s) = c("Year","Flag","GearCode","Catch")
  aggF = aggregate(t1s$Catch,by=list(Flag = t1s$Flag, Gear = t1s$GearCode),sum)
  aggF = aggF[order(aggF$x,decreasing = T),]
 
  print(cbind(aggF,cumsum(aggF$x)/sum(aggF$x)*100)[1:10,])
  catch = aggregate(t1s$Catch,by=list(Year=t1s$Year),sum)
  plot(catch$Year, catch$x,type="l",xlab="Year",ylab="Catch (t)")
  
  t2s = t2[t2$SpeciesCode == SpecCode,]
  
  aggF2 = aggregate(t2s$Nr,by=list(Flag = t2s$FleetCode, Gear =t2s$GearCode), sum)
  aggF2 = aggF2[order(aggF2$x,decreasing = T),]
  print(cbind(aggF2,cumsum(aggF2$x)/sum(aggF2$x)*100)[1:10,])
  
  if(is.na(tF)) tF = aggF2$Flag[1]
  if(is.na(tG)) tG = aggF2$Gear[1]
  
  t2s = t2[t2$SpeciesCode == SpecCode & t2$FleetCode == tF & t2$GearCode == tG, ]
  CAL = aggregate(t2s$Nr,by=list(Year = t2s$YearC, LC = ceiling(t2s$ClassFrq/2)*2-1), sum) # cbind(t2s$ClassFrq,ceiling(t2s$ClassFrq/2)*2-1)
  nCAL = nrow(CAL)
  
  cat(paste0("Returning length compositions for ",tF," - ",tG," \n"))
  
  yrs = min(catch$Year):max(catch$Year)
  Lbins = unique(ceiling((1:max(t2s$ClassFrq))/2)*2-1)
  ny = length(yrs)
  nl = length(Lbins)
  CALmat = array(0,c(1,ny,nl))
  ind = as.matrix(cbind(rep(1,nCAL),match(CAL$Year,yrs),match(CAL$LC,Lbins)))
  CALmat[ind] = CAL$x
  CALsumy = apply(CALmat[1,,],1,sum)
  fsy = min((1:ny)[CALsumy!=0])
  lsy = min((1:ny)[CALsumy!=0])
  if(CALint)for(yy in 1:ny)if(CALsumy[yy]==0 & yy<fsy)CALmat[1,yy,] = CALmat[1,fsy,]
  if(CALint)for(yy in 1:ny)if(CALsumy[yy]==0 & yy>lsy)CALmat[1,yy,] = CALmat[1,lsy,]
  
  list(Lbins =Lbins, CAL = CALmat, Catch = matrix(catch$x,nrow=1), yrs = yrs, ny =ny, nl = nl)
  
}