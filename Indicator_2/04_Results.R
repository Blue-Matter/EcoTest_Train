


library(openMSE)
library(miceadds)
library(SimDesign)
library(data.table)
library(xlsx)

setwd("C:/GitHub/EcoTest")
source.all("Source")


# --- Methods ----------------------------------------

MOM = make_MOM(nsim = 5000, seed=1)

jpeg("Indicator_2/Figures/MOM_LH.jpg",res=400, width=7, height=5, units="in")
  plot_OM_LH(MOM)
dev.off()


jpeg("Indicator_2/Figures/MOM_Find.jpg",res=400, width=9, height=6.5, units="in")
  plot_OM_Find(MOM)
dev.off()


jpeg("Indicator_2/Figures/MOM_Vuln.jpg",res=400, width=7, height=5, units="in")
  plot_OM_Vuln(MOM)
dev.off()


jpeg("Indicator_2/Figures/MOM_F_cor.jpg",res=400, width=6.5, height=8, units="in")
  plot_OM_F_cor(MOM)
dev.off()


# --- Results ------------------------------------------------------------------

runs = readRDS("Indicator_2/Fitted_Models/runs.rds")

vals = sapply(1:nrow(runs),function(x){
  out=readRDS(paste0("Indicator_2/Fitted_Models/Fit_",x,".rds"))
  model = load_model(paste0("Indicator_2/Fitted_Models/Model_",x,".keras"))
  np =model$count_params()
  round(c(npar = np,R_sq = out$r2, MAE = out$MAE, TP_LPR = out$tab[1,1], TP_R = out$tab[2,2]),3)})

res = cbind(runs,t(vals))

saveRDS(res,"Indicator_2/Results/All_models.rds")

 by_fleet_stock(res, quant, 1, 1)
 
 by_fleet_stock=function(res, quant="MAE", fleet=1, stock=1){
    subset(res,res$nstocks==stock & res$nfleets==fleet)
 }

write.xlsx(res,"Indicator_2/Results/Results.xlsx","All")


# --- Demonstration figures ---------------------------------------------------- 
 

runi = c(1,8,11,12,36,62)

jpeg("Indicator_2/Figures/Results_cors.jpg",res=400, width=9.5, height=7, units="in")
  par(mfrow=c(2,3),mai=c(0.6,0.6,0.3,0.05))
  for(j in 1:length(runi)){
    i = runi[j]
    fit =  readRDS(paste0("Indicator_2/Fitted_Models/Fit_",i,".rds"))
    pred_plot(fit, newplot=F, lab = getlab(i,runs))               
  }
dev.off()

