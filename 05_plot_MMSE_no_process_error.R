

MMSE <- readRDS("MSE/MMSE_npe.rds")



# F - by fleet
FM <- MMSE@FM %>% structure(
  dimnames = list(
    Sim = 1:MMSE@nsim,
    Stock = MMSE@Snames,
    Fleet = c("Longline", "Other"),
    MP = MMSE@MPs[[1]],
    Year = 2013 + 1:MMSE@proyears
  )
) %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))

g <- ggplot(FM %>% filter(Sim == 1, Sex == "Female"), aes(Year, value, linetype = Fleet)) + 
  facet_grid(paste(Label, Species) ~ MP) + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.6)) + 
  labs(x = "Projection Year", y = "Apical fishing mortality") +
  expand_limits(y = 0)
ggsave("Figures/MMSE/F_npe.png", g, height = 8, width = 8)

# SSB
SSB <- MMSE@SSB %>% structure(
  dimnames = list(
    Sim = 1:MMSE@nsim,
    Stock = MMSE@Snames,
    MP = MMSE@MPs[[1]],
    Year = 2013 + 1:MMSE@proyears
  )
) %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))
g <- ggplot(SSB %>% filter(Sim == 1, Sex == "Female"), aes(Year, value, linetype = MP)) + 
  facet_wrap(~ paste(Label, Species), scales = "free_y") + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) + 
  labs(y = "Spawning biomass", linetype = "Fishing scenario") +
  expand_limits(y = 0)
ggsave("Figures/MMSE/SSB_npe.png", g, height = 4, width = 8)

# Catch
Catch <- MMSE@Catch %>% structure(
  dimnames = list(
    Sim = 1:MMSE@nsim,
    Stock = MMSE@Snames,
    Fleet = c("Longline", "Other"),
    MP = MMSE@MPs[[1]],
    Year = 2013 + 1:MMSE@proyears
  )
) %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3)) %>%
  group_by(Sim, Fleet, Year, Species, MP) %>%
  summarise(value = sum(value)) %>%
  mutate(Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))
g <- ggplot(Catch %>% filter(Sim == 1), aes(Year, value, colour = Fleet, linetype = MP)) + 
  facet_grid(paste(Label, Species) ~ MP, scales = "free_y") + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) + 
  labs(y = "Catch") +
  expand_limits(y = 0)
ggsave("Figures/MMSE/Catch_npe.png", g, height = 8, width = 8)









# B/BMSY 
multiHist <- readRDS("MOM/multiHist_npe.rds")

nsim <- MMSE@nsim
np <- MMSE@nstocks
nf <- MMSE@nfleets
F_at_age <- sapply(1:np, function(p) {
  sapply(1:nf, function(f) {
    sapply(1:nsim, function(x) {
      Fout <- array(NA_real_, c(MMSE@proyears, MMSE@Stocks[[1]]@maxage + 1, MMSE@nMPs))
      for(mm in 1:MMSE@nMPs) {
        Fout[, , mm] <- MMSE@FM[x, p, f, mm, ] * t(multiHist[[p]][[f]]@SampPars$Fleet$V_real[x, , MMSE@nyears + 1:MMSE@proyears])
      }
      Fout
    }, simplify = "array")
  }, simplify = "array")
}, simplify = "array")
F_agg <- F_at_age %>% apply(c(1:4, 6), sum)
Fapic <- apply(F_agg, c(1, 3:5), max) 

%>% aperm(c(2, 1, 3))
V_agg <- apply(F_agg, c(1, 3:5), function(x) x/max(x)) # y, a, mm, x, p

MSYRefs <- sapply(1:np, function(p) {
  StockPars <- multiHist[[p]][[1]]@SampPars$Stock
  sapply(1:MMSE@proyears, function(y) {
    sapply(1:MMSE@nMPs, function(mm) {
      sapply(1:nsim, optMSY_eq,
             M_ageArray=StockPars$M_ageArray[, , -c(1:MMSE@nyears)],
             Wt_age=StockPars$Wt_age[, , -c(1:MMSE@nyears)],
             Mat_age=StockPars$Mat_age[, , -c(1:MMSE@nyears)],
             Fec_age=StockPars$Fec_Age[, , -c(1:MMSE@nyears)],
             V=V_agg[, , mm, , p] %>% aperm(c(3, 1, 2)),
             maxage=StockPars$maxage,
             R0=StockPars$R0,
             SRrel=StockPars$SRrel,
             hs=StockPars$hs,
             SSBpR=StockPars$SSBpR,
             yr.ind=y,
             plusgroup = StockPars$plusgroup,
             simplify = "array")
    }, simplify = "array")
  }, simplify = "array")
}, simplify = "array")



SBMSY <- MSYRefs[3, , , , ] %>% aperm(c(1, 4, 2, 3))
B_BMSY <- structure(MMSE@SSB/SBMSY, dimnames = list(
  Sim = 1:MMSE@nsim,
  Stock = MMSE@Snames,
  MP = MMSE@MPs[[1]],
  Year = 2013 + 1:MMSE@proyears
))  %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))
g <- ggplot(B_BMSY %>% filter(Sim == 1, Sex == "Female"), aes(Year, value, linetype = MP)) + 
  facet_wrap(~ paste(Label, Species)) + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.6)) + 
  #geom_hline(yintercept = 1, linetype = 3) + 
  expand_limits(y = 0) + 
  labs(x = "Projection Year", y = expression(SSB/SSB[MSY]), linetype = "Fishing scenario")
ggsave("Figures/MMSE/BMSY_npe.png", g, height = 4, width = 8)



s########## PPD Longline

# Vulnerable index
VInd <- sapply(1:np, function(p) {
  MMSE@PPD[[p]][[1]][[1]]@VInd[, -c(1:MMSE@nyears)]
}, simplify = "array") %>%
  structure(
    dimnames = list(
      Sim = 1:MMSE@nsim,
      Year = 2013 + 1:(MMSE@proyears-1),
      Stock = MMSE@Snames
    )
  ) %>% 
  reshape2::melt() %>%
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))

VInd_out <- left_join(VInd %>% rename(VInd = value), 
                      B_BMSY %>% rename(B_BMSY = value),
                      by = c("Sim", "Year", "Stock", "Species", "Sex", "Label")) %>%
  filter(Label == "Secondary:", Sim == 1, Sex == "Female")

for(i in unique(VInd_out$MP)) {
  g <- filter(VInd_out, MP == i) %>%
    ggplot(aes(B_BMSY, VInd, colour = Year)) + 
    facet_wrap(~ paste(Label, Species), scales = "free") + 
    geom_point() + 
    geom_path() + 
    theme_bw() +
    labs(x = expression(SSB/SSB[MSY]), y = "Simulated longline CPUE") +
    scale_colour_viridis_c() +
    ggtitle(i)
  ggsave(paste0("Figures/MMSE/VInd_", i, "_npe.png"), g, height = 4, width = 6)
}





# Mean length
ML <- sapply(1:np, function(p) {
  MMSE@PPD[[p]][[1]][[1]]@ML[, -c(1:MMSE@nyears)]
}, simplify = "array") %>%
  structure(
    dimnames = list(
      Sim = 1:MMSE@nsim,
      Year = 2013 + 1:(MMSE@proyears-1),
      Stock = MMSE@Snames
    )
  ) %>% 
  reshape2::melt() %>%
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Primary:", "Secondary:"))

ML_out <- left_join(ML %>% rename(ML = value), 
                    B_BMSY %>% rename(B_BMSY = value),
                    by = c("Sim", "Year", "Stock", "Species", "Sex", "Label")) %>%
  filter(Label == "Secondary:", Sim == 1, Sex == "Female")

for(i in unique(ML_out$MP)) {
  g <- filter(ML_out, MP == i) %>%
    ggplot(aes(B_BMSY, ML, colour = Year)) + 
    facet_wrap(~ paste(Label, Species), scales = "free") + 
    geom_point() + 
    geom_path() + 
    theme_bw() +
    labs(x = expression(SSB/SSB[MSY]), y = "Longline mean length (cm)") +
    scale_colour_viridis_c() +
    ggtitle(i)
  ggsave(paste0("Figures/MMSE/ML_", i, "_npe.png"), g, height = 4, width = 6)
}

