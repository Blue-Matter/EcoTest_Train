


cleandat = function(dat,resprange=c(-3.689,0.693)){
  goodrow = function(x)(sum(is.na(x)| x==0 | x==1 | x == 'Inf' | is.null(x)| x=="NULL")==0)&x[1]>resprange[1] & x[1]<resprange[2]
  dat[apply(dat,1,goodrow),]
}

reformat = function(x){
  
  renam = grepl("-",x)
  charmatch("-",x[162])
  
}

makerawdata = function(allout, sno=1, isBrel = F, clean = T, 
                       inc_spat = T, inc_Irel = T){
  
  cdat = as.data.frame(rbindlist(allout))
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  ns = sum(Brelcols) # os = (1:ns)[!((1:ns)%in%sno)]
  keepcol = grep(paste0("Brel_",sno),dnames)
  Res = log(cdat[,keepcol])
 
  # level relative to mean level of other stocks
  if(isBrel){Res = ((ns-1) * Res) /(apply(cdat[,Brelcols],1,sum)-Res)}
  
  dat = cbind(Res,cdat[,!Brelcols])
  
  if(!inc_spat) dat = dat[,!grepl('spat',names(dat))]
  if(!inc_Irel) dat = dat[,!grepl('I_rel',names(dat))]
  
  # log ratios
  tolog = grepl("C_rel", names(dat)) | grepl("CR_mu",names(dat)) | grepl("FM_cur", names(dat)) | grepl("FM_rel", names(dat)) |
    grepl("ML_cur", names(dat)) | grepl("ML_rel" , names(dat)) | grepl("MV_cur",names(dat)) |
    grepl("MV_rel",names(dat))| grepl("ML_Linf" ,names(dat)) |   grepl("ML_L50" , names(dat)) |
    grepl("CR_rel",names(dat))| grepl("maxa", names(dat)) | grepl("M_K",names(dat))
  
  dat[,tolog] = log(dat[,tolog])
  if(clean)dat=cleandat(dat)
  isconst = apply(dat,2,sd)<1E-10
  if(sum(isconst)>0){
    cat(paste(paste(names(dat)[isconst],collapse=", "), "dropped for sd < 1E-10 \n"))
    dat = dat[,!isconst]
  }
  
  dat
  
}




power_tab = function(sim, pred, lev = c(0.5,1),asprob = T){
  
  nlev = length(lev)
  ncat = nlev+1
  labs=rep(NA,ncat)
  labs[1] = paste("Brel <",lev[1])
  labs[ncat] = paste(lev[nlev], "< Brel")
  if(nlev>1){ for(i in 2:nlev){
    labs[i] = paste(lev[i-1],"< Brel <",lev[i])
  }}
  
  tab = array(NA,c(ncat,ncat))
  row.names(tab) = paste("(Pred)",labs)
  colnames(tab) = paste("(Sim)",labs)
  
  alllev = c(-Inf,lev,Inf)
  for(i in 1:ncat){ #pred
    for(j in 1:ncat){ #sim
      cond = sim > alllev[j] & sim < alllev[j+1] & pred > alllev[i] & pred < alllev[i+1]
      tab[i,j] = sum(cond)
    }  
  }
  
  if(asprob)  tab = tab / rep(colSums(tab),each=ncat)
  tab
}


NN_fit = function(sim, pred, history,lev, addpow=T){
  nlev = length(lev)
  alllev = c(-Inf,lev,Inf)
  ncat = nlev+1
  colt = rainbow(ncat,start=0,end=0.35)
  cols = rep(colt[1],length(sim))
  for(i in 2:ncat){
    cols[pred > alllev[i] & pred < alllev[i+1]] = colt[i]
  }
  plot(sim, pred, xlab="SSB/SSBMSY (simulated)",ylab="SSB/SSBMSY (pred)",pch=19,cex=1.2,col=cols)
  lines(c(0,1E10),c(0,1E10),col='black',lwd=1,lty=2)
  abline(h=levs,v=lev,lty=2)
  
  if(addpow){
    tab = power_tab(sim, pred, lev)
    
    mins_pred = c(min(pred),levs)
    maxs_pred = c(levs,(max(pred)))
    difs_pred = maxs_pred-mins_pred
    muy = mins_pred + difs_pred/4
    
    mins_sim = c(min(sim),levs)
    maxs_sim = c(levs,(max(sim)))
    difs_sim = maxs_sim-mins_sim
    mux = maxs_sim-difs_sim/4
    
    grid = expand.grid(muy,mux)
    text(grid[,2],grid[,1],round(as.vector(tab)*100,1))
  }
  
  legend('topleft',legend = paste("MAE_val:",round(history$metrics$val_MAE[nepoch],4)),bty="n")
  tab
}


calc_importance = function(model, testy, adj = 0.25, barno = 30){
  nin = ncol(testy)
  impdf = rbind(testy,testy,testy,testy,testy)[1:(nin+1),]
  impdf[]=0.0
  ind = as.matrix(cbind(2:(nin+1),1:nin))
  impdf[ind] = adj
  sense = exp((model %>% predict(impdf))[,1])
  dif = abs(sense[2:nin] - sense[1])
  ord = order(dif,decreasing=T)
  out = data.frame(input=colnames(testy)[ord],output_diff = dif[ord])
  y = out[1:barno,2]
  x=barplot(y,xaxt='n');grid(); barplot(y,xaxt='n',add=T)
  yadj = max(y)/8
  text(x=x[,1]-0.25, y=rep(-yadj,barno), labels=out[1:barno,1], xpd=TRUE, srt=90,cex=0.85)
  mtext("Relative weight of input",2,line=2.3)
  out
}

