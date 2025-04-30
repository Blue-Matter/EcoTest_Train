

library(trelliscopejs)
library(tidyverse)

dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases"
LLCE <- read.csv(file.path(dir, "t2ce_LL.csv"))

# Calculate longline CPUE as catch/Eff1 - aggregated at trip level?
LLN <- LLCE %>%
  filter(CatchUnit == "nr", Eff1Type == "NO.HOOKS") %>%
  mutate(Lat = ifelse(QuadID == 1 | QuadID == 4, Lat, -1 * Lat),
         Lon = ifelse(QuadID == 1 | QuadID == 2, Lon, -1 * Lon),
         BET = BET/Eff1, SWO = SWO/Eff1, WHM = WHM/Eff1, BUM = BUM/Eff1, SMA = SMA/Eff1, BSH = BSH/Eff1) %>%
  select(YearC, Lat, Lon, BET, SWO, BUM, WHM, BSH, SMA, FlagName)

LLW <- LLCE %>%
  filter(CatchUnit == "kg", Eff1Type == "NO.HOOKS") %>%
  mutate(Lat = ifelse(QuadID == 1 | QuadID == 4, Lat, -1 * Lat),
         Lon = ifelse(QuadID == 1 | QuadID == 2, Lon, -1 * Lon),
         BET = BET/Eff1, SWO = SWO/Eff1, WHM = WHM/Eff1, BUM = BUM/Eff1, SMA = SMA/Eff1, BSH = BSH/Eff1) %>%
  select(YearC, Lat, Lon, BET, SWO, BUM, WHM, BSH, SMA, FlagName)

# Species combinations
sp <- c("BET", "SWO", "BSH", "SMA", "WHM", "BUM")

# Calculate joint catches
get_jc <- function(LL, fname) {
  
  tab <- lapply(1:length(sp), function(x) expand.grid(sp[x], sp[x:length(sp)], stringsAsFactors = FALSE)) %>% bind_rows()
  
  lapply(1:nrow(tab), function(x) {
    s1 <- as.symbol(tab$Var1[x])
    s2 <- as.symbol(tab$Var2[x])
    
    ##### Aggregate by year and spatial area - perhaps too fine for calculating correlation among CPUE
    if(s1 == s2) {
      n_records <- LL %>% filter(FlagName == fname) %>% group_by(YearC) %>% summarise(n = n())
      pos <- LL %>% filter(!!s1 > 0, FlagName == fname) %>%
        group_by(YearC) %>%
        summarise(npos = n())
    } else {
      n_records <- LL %>% filter(!!s1 > 0 || !!s2 > 0, FlagName == fname) %>% group_by(YearC) %>% summarise(n = n())
      pos <- LL %>% filter(!!s1 > 0, !!s2 > 0, FlagName == fname) %>%
        group_by(YearC) %>%
        summarise(npos = n())
    }
    
    left_join(pos, n_records, by = "YearC") %>% mutate(S1 = tab$Var1[x], S2 = tab$Var2[x], p_pos = npos/n, FlagName = fname)
      
  }) %>% 
    bind_rows() %>% 
    mutate(S1 = factor(S1, levels = sp), S2 = factor(S2, levels = sp))
  
}

# Calculate spatial correlation in CPUE from non zero catch
get_corr <- function(LL, fname) {
  tab <- lapply(1:length(sp), function(x) expand.grid(sp[x], sp[x:length(sp)], stringsAsFactors = FALSE)) %>% bind_rows()
  
  lapply(1:nrow(tab), function(x) {
    s1 <- as.symbol(tab$Var1[x])
    s2 <- as.symbol(tab$Var2[x])
    
    ##### Aggregate by year and spatial area - perhaps too fine for calculating correlation among CPUE
    LL %>% filter(!!s1 > 0, !!s2 > 0, FlagName == fname) %>%
      group_by(YearC) %>%
      summarise(corr = cor(!!s1, !!s2)) %>%
      mutate(S1 = tab$Var1[x], S2 = tab$Var2[x], FlagName = fname)
    
    
  }) %>% 
    bind_rows() %>% 
    mutate(S1 = factor(S1, levels = sp), S2 = factor(S2, levels = sp))
  
}

# Japan
JPN <- get_jc(LLN, "Japan")
JPN_corr <- get_corr(LLN, "Japan")

# TPE
TPE <- get_jc(LLW, "Chinese Taipei")
TPE_corr <- get_corr(LLW, "Chinese Taipei")

# USA
USA <- get_jc(LLN, "U.S.A.")
USA_corr <- get_corr(LLN, "U.S.A.")

ESP <- get_jc(LLN, "EU.España")
ESP_corr <- get_corr(LLN, "EU.España")

BRA <- get_jc(LLW, "Brazil")
BRA_corr <- get_corr(LLW, "Brazil")

pos <- rbind(JPN, TPE, USA, ESP, BRA)
  
g <- ggplot(pos, aes(YearC, p_pos, colour = FlagName)) + 
  facet_grid(vars(S1), vars(S2)) + 
  geom_line() + 
  #geom_point() + 
  theme_bw() +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  labs(x = "Year", y = "Proportion positive catches of both species", colour = "Flag")
ggsave("Figures/LLCE/LL_spatial_pos.png", g, height = 6, width = 10)

f_corr <- rbind(JPN_corr, 
                TPE_corr, 
                USA_corr, 
                ESP_corr, 
                BRA_corr)

g <- ggplot(f_corr, aes(YearC, corr, colour = FlagName)) + 
  facet_grid(vars(S1), vars(S2)) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_line() + 
  #geom_point() + 
  theme_bw() +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  labs(x = "Year", y = "CPUE correlation", colour = "Flag")
ggsave("Figures/LLCE/LL_spatial_corrs.png", g, height = 6, width = 10)



### Scroll through annual longline CPUE correlations
ggplot(get_corr, aes(Lon, Lat, colour = corr)) + facet_grid(vars(S1), vars(S2)) + theme_bw() + 
  geom_jitter(alpha = 0.2) +
  facet_trelliscope(~ YearC) + 
  theme(panel.spacing = unit(0, "in")) + 
  scale_colour_gradient2(high = "red", low = "blue", mid = "grey90") +
  labs(x = "Longitude", y = "Latitude", colour = "CPUE\ncorrelation")


### Save to disk year-specific map with coastline
year_to_plot <- 2009
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

g <- get_corr %>% filter(YearC == year_to_plot) %>%
  ggplot(aes(Lon, Lat, colour = corr)) + facet_grid(vars(S1), vars(S2)) + theme_bw() + 
  geom_sf(data = coast, inherit.aes = FALSE) +
  coord_sf(xlim = range(get_corr$Lon), ylim = range(get_corr$Lat)) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 3) + 
  geom_jitter(alpha = 0.6) +
  #facet_trelliscope(~ YearC) + 
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  scale_colour_gradient2(high = "red", low = "blue", mid = "grey90") +
  ggtitle(year_to_plot) +
  labs(x = "Longitude", y = "Latitude", colour = "CPUE\ncorrelation")
ggsave(paste0("Figures/LLCE/LL_spatial_", year_to_plot, ".png"), g, height = 10, width = 10)


# Annual correlations across all spatial areas and fleets
get_corr_year <- local({
  tab <- lapply(1:length(sp), function(x) expand.grid(sp[x], sp[x:length(sp)], stringsAsFactors = FALSE)) %>% bind_rows()
  lapply(1:nrow(tab), function(x) {
    s1 <- as.symbol(tab$Var1[x])
    s2 <- as.symbol(tab$Var2[x])
    corr <- group_by(LL, YearC) %>% summarise(corr = cor(!!s1, !!s2, use = "complete.obs")) %>%
      filter(!is.na(corr)) %>%
      mutate(S1 = tab$Var1[x], S2 = tab$Var2[x])
    corr
  }) %>% 
    bind_rows() %>% 
    mutate(S1 = factor(S1, levels = sp), S2 = factor(S2, levels = sp)) %>%
    filter(S1 != S2)
})

g <- ggplot(get_corr_year, aes(YearC, corr)) + facet_grid(vars(S1), vars(S2)) + theme_bw() + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_line() + 
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Year", y = "CPUE correlation") +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.4)) 
ggsave("Figures/LLCE/LLCE_annual_correlation.png", g, height = 4, width = 6)


# Annual correlation in CPUE by country
get_corr_flag <- local({
  tab <- lapply(1:length(sp), function(x) expand.grid(sp[x], sp[-c(1:x)], stringsAsFactors = FALSE)) %>% bind_rows()
  
  lapply(1:nrow(tab), function(x) {
    s1 <- as.symbol(tab$Var1[x])
    s2 <- as.symbol(tab$Var2[x])
    
    ##### Aggregate by year and FlagName
    corr <- group_by(LL, YearC, FlagName) %>% 
      summarise(corr = cor(!!s1, !!s2)) %>%
      filter(!is.na(corr)) %>%
      mutate(S1 = tab$Var1[x], S2 = tab$Var2[x])
    
    corr
    
  }) %>% 
    bind_rows() %>% 
    mutate(S1 = factor(S1, levels = sp), S2 = factor(S2, levels = sp))
  
})


ggplot(get_corr_flag, aes(YearC, corr)) + facet_grid(vars(S1), vars(S2)) + theme_bw() + 
  #geom_hline(yintercept = 0, linetype = 3) + 
  geom_line() + 
  coord_cartesian(ylim = c(-1, 1)) +
  facet_trelliscope(~ FlagName) + 
  labs(x = "Year", y = "CPUE correlation") +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.4)) 
