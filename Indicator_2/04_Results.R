


library(openMSE)
library(miceadds)
library(SimDesign)
library(data.table)

setwd("C:/GitHub/EcoTest")
source.all("Source")

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

