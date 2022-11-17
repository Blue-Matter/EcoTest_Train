

MMSE <- readRDS("MSE/MMSE_npe.rds")



# F
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

g <- ggplot(FM %>% filter(Sim == 1, Sex == "Female"), aes(Year, value, colour = Fleet, linetype = MP)) + 
  facet_wrap(~ paste(Label, Species), scale = "free_y") + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  labs(y = "Apical fishing mortality") +
  expand_limits(y = 0)
ggsave("Figures/MMSE/F_npe.png", g, height = 4, width = 8)

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
  labs(y = "Spawning biomass") +
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
  facet_wrap(~ paste(Label, Species), scales = "free_y") + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  labs(y = "Catch") +
  expand_limits(y = 0)
ggsave("Figures/MMSE/Catch_npe.png", g, height = 4, width = 8)









# F/FMSY - calculate aggregate F then divide by FMSY
nsim <- 2
np <- 10
nf <- 2
F_at_age <- sapply(1:np, function(p) {
  sapply(1:nf, function(f) {
    sapply(1:nsim, function(x) {
      MMSE@FM[x, p, f, 1, ] * t(multiHist[[p]][[f]]@SampPars$Fleet$V_real[x, , MMSE@nyears + 1:MMSE@proyears])
    }, simplify = "array")
  }, simplify = "array")
}, simplify = "array")

F_agg <- F_at_age %>% apply(c(1, 2, 3, 5), sum)
Fapic <- apply(F_agg, c(1, 3, 4), max) %>% aperm(c(2, 1, 3))
V_agg <- apply(F_agg, c(1, 3, 4), function(x) x/max(x))

MSYRefs <- sapply(1:np, function(p) {
  StockPars <- multiHist[[p]][[1]]@SampPars$Stock
  sapply(1:MMSE@proyears, function(y) {
    MSYrefsYr <- sapply(1:nsim, optMSY_eq,
                        M_ageArray=StockPars$M_ageArray[, , -c(1:MMSE@nyears)],
                        Wt_age=StockPars$Wt_age[, , -c(1:MMSE@nyears)],
                        Mat_age=StockPars$Mat_age[, , -c(1:MMSE@nyears)],
                        Fec_age=StockPars$Fec_Age[, , -c(1:MMSE@nyears)],
                        V=V_agg[, , , p] %>% aperm(c(3, 1, 2)),
                        maxage=StockPars$maxage,
                        R0=StockPars$R0,
                        SRrel=StockPars$SRrel,
                        hs=StockPars$hs,
                        SSBpR=StockPars$SSBpR,
                        yr.ind=y,
                        plusgroup=StockPars$plusgroup)
  }, simplify = "array")
}, simplify = "array")

FMSY <- MSYRefs[2, , , ]
F_FMSY <- structure(Fapic/FMSY, dimnames = list(
  Sim = 1:MMSE@nsim,
  Year = 2013 + 1:MMSE@proyears,
  Stock = MMSE@Snames
))  %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Target:", "Bycatch:"))
g <- ggplot(F_FMSY %>% filter(Sim == 1, Sex == "Female"), aes(Year, value)) + 
  facet_wrap(~ paste(Label, Species)) + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 3) + expand_limits(y = 0) + 
  labs(y = expression(F/F[MSY]))
ggsave("Figures/MMSE/FMSY_FMSYtarget.png", height = 4, width = 8)


SBMSY <- MSYRefs[3, , , ]
B_BMSY <- structure(aperm(MMSE@SSB[,,1,], c(1, 3, 2))/SBMSY, dimnames = list(
  Sim = 1:MMSE@nsim,
  Year = 2013 + 1:MMSE@proyears,
  Stock = MMSE@Snames
))  %>% reshape2::melt() %>% 
  mutate(Species = substr(Stock, 1, 3),
         Sex = Stock %>% as.character() %>% 
           sapply(function(x) if(nchar(x) > 3) substr(x, 5, nchar(x)) else "Female"),
         Label = ifelse(Species %in% c("BET", "SWO"), "Target:", "Bycatch:"))
g <- ggplot(B_BMSY %>% filter(Sim == 1, Sex == "Female"), aes(Year, value)) + 
  facet_wrap(~ paste(Label, Species)) + 
  #geom_point() + 
  geom_line() + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 3) + expand_limits(y = 0) + 
  labs(y = expression(SSB/SSB[MSY]))
ggsave("Figures/MMSE/BMSY_FMSYtarget.png", g, height = 4, width = 8)



########## PPD Longline

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
         Label = ifelse(Species %in% c("BET", "SWO"), "Target:", "Bycatch:"))

VInd_out <- left_join(VInd %>% rename(VInd = value), 
                      B_BMSY %>% rename(B_BMSY = value),
                      by = c("Sim", "Year", "Stock", "Species", "Sex", "Label"))

g <- ggplot(VInd_out %>% filter(Sim == 1, Sex == "Female"), aes(B_BMSY, VInd, colour = Year)) + 
  facet_wrap(~ paste(Label, Species), scales = "free") + 
  #geom_vline(xintercept = 0, linetype = 3) + 
  geom_path() + 
  theme_bw() +
  labs(x = expression(SSB/SSB[MSY]), y = "Longline index") +
  scale_colour_viridis_c()
ggsave("Figures/MMSE/VInd_FMSYtarget.png", g, height = 4, width = 8)





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
         Label = ifelse(Species %in% c("BET", "SWO"), "Target:", "Bycatch:"))

ML_out <- left_join(ML %>% rename(ML = value), 
                    B_BMSY %>% rename(B_BMSY = value),
                    by = c("Sim", "Year", "Stock", "Species", "Sex", "Label"))

g <- ggplot(ML_out %>% filter(Sim == 1, Sex == "Female"), aes(B_BMSY, ML, colour = Year)) + 
  facet_wrap(~ paste(Label, Species), scales = "free") + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = expression(SSB/SSB[MSY]), y = "Longline mean length (cm)") +
  scale_colour_viridis_c()
ggsave("Figures/MMSE/ML_FMSYtarget.png", g, height = 4, width = 8)
