
library(tidyverse)
theme_set(theme_bw())

LLCE <- local({
  dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases"
  read.csv(file.path(dir, "t2ce_LL.csv")) %>% 
    mutate(Latitude = ifelse(QuadID %in% 2:3, -Lat, Lat), 
           Longitude = ifelse(QuadID %in% 3:4, -Lon, Lon))
})

## Effort by country
top10 <- LLCE %>% 
  filter(#FlagCode %in% c("JPN", "USA"),
    YearC >= 2010,
    TimePeriodID <= 12, 
    Eff1Type == "NO.HOOKS", 
    SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
    ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")
  ) %>%
  group_by(FlagCode) %>% 
  summarise(Eff = sum(Eff1)) %>%
  mutate(Eff_order = rank(Eff)) %>%
  filter(Eff_order >= 22)

g <- LLCE %>% 
  filter(FlagCode %in% top10$FlagCode,
         TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")
  ) %>%
  group_by(YearC, QuadID, FlagCode) %>%
  summarise(Effort = sum(Eff1)) %>%
  group_by(YearC, FlagCode) %>%
  mutate(p = Effort/sum(Effort)) %>%
  ggplot(aes(YearC, p)) + 
  geom_col(aes(colour = factor(QuadID), fill = factor(QuadID))) +
  labs(x = "Year", y = "Proportion of hooks", colour = "Quadrant", fill = "Quadrant") +
  facet_wrap(vars(FlagCode)) +
  theme_bw() + 
  theme(legend.position = "bottom")
ggsave("Figures/LLCE/spatial-effort-country.png", g, height = 4, width = 6)


# Summary table - number of records (exclude < 10 percent of all records)
LLsumry <- LLCE %>% 
  filter(TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         #ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")
         ) %>%
  group_by(FlagCode, YearC, CatchUnit) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 0) %>%
  mutate(CatchUnit2 = ifelse(CatchUnit == "kg", "CPUE (kg/hook)", "CPUE (nr/hook)"))
LLsumry$FlagCode <- factor(LLsumry$FlagCode, levels = LLsumry$FlagCode %>% unique() %>% rev())

g <- ggplot(LLsumry, aes(YearC, FlagCode, fill = n)) + 
  geom_tile(colour = "black", width = 1) +
  facet_wrap(vars(CatchUnit2), ncol = 1) +
  theme_bw() +
  theme(panel.spacing = unit(0, "in")) + 
  scale_fill_viridis_b(trans = "log", breaks = LLsumry$n %>% log() %>% pretty() %>% exp() %>% round()) +
  labs(x = "Year", fill = "Number of records")
ggsave("Figures/LLCE/LLCE_summary.png", g, height = 8, width = 6)

# LL Nr - individual species CPUE > 0
LLsp <- LLCE %>% 
  filter(FlagCode %in% top10$FlagCode,
    #FlagCode %in% c("BRA", "JPN", "EU.ESP", "KOR", "TAI", "USA", "VEN"),
         TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")) %>%
  select(YearC, FlagCode, StrataID, BET, BSH, BUM, SMA, SWO, WHM, Eff1) %>%
  reshape2::melt(id.vars = c("YearC", "FlagCode", "StrataID", "Eff1"), 
                 variable.name = "Species", value.name = "Catch") %>%
  #filter(Catch > 0) %>%
  group_by(FlagCode, YearC, Species) %>%
  summarise(p = mean(Catch > 0))
LLsp$Species <- factor(LLsp$Species, levels = unique(LLsp$Species) %>% rev())

g <- ggplot(LLsp, aes(YearC, Species)) +
  geom_tile(aes(fill = p), colour = "black", width = 1) +
  facet_wrap(vars(FlagCode), ncol = 3) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.spacing = unit(0, "in")) + 
  scale_fill_viridis_c() + 
  #scale_fill_gradient2(high = scales::muted("red")) +
  labs(x = "Year", y = "Species", fill = "Proportion positive CPUE")
ggsave("Figures/LLCE/LLCEsp_summary.png", g, height = 6, width = 6)

# LL Nr spatial footprint
LLNRspace <- LLCE %>% 
  filter(#FlagCode %in% c("BRA", "JPN", "EU.ESP", "KOR", "TAI", "USA", "VEN"),
         TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")) %>%
  mutate(Decade = paste0(floor(YearC/10) * 10, "s")) %>%
  group_by(Decade, Latitude, Longitude) %>%
  summarise(Nhooks = sum(Eff1)) %>%
  group_by(Decade) %>%
  mutate(p = Nhooks/sum(Nhooks))

coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
#g <- ggplot(LLNRspace %>% filter(FlagCode == "JPN"), aes(Longitude, Latitude)) + 
#  geom_sf(data = coast, inherit.aes = FALSE) +
#  coord_sf(xlim = c(-100, 30), ylim = c(-45, 65)) + 
#  geom_point(shape = 21, aes(fill = n, colour = n), alpha = 0.85) +
#  facet_wrap(vars(Decade), ncol = 3) +
#  scale_color_viridis_c(trans = "log", name = "Number of records", breaks = c(1, 7, 54)) +
#  scale_fill_viridis_c(trans = "log", name = "Number of records", breaks = c(1, 7, 54)) +
#  theme(legend.position = "bottom",
#        panel.spacing = unit(0, "in"))
#ggsave("Figures/LLCE/LLNRspatial_JPN.png", g, height = 8, width = 6)

g <- ggplot(LLNRspace, aes(Longitude, Latitude)) + 
  geom_sf(data = coast, inherit.aes = FALSE) +
  coord_sf(xlim = c(-100, 30), ylim = c(-55, 65)) + 
  geom_point(shape = 21, aes(fill = Nhooks, colour = Nhooks), alpha = 0.85) +
  facet_wrap(vars(Decade), ncol = 3) +
  scale_color_viridis_c(trans = "log", name = "Number of hooks") +
  scale_fill_viridis_c(trans = "log", name = "Number of hooks") +
  theme(legend.position = "bottom",
        panel.spacing = unit(0, "in"))
ggsave("Figures/LLCE/LLNRspatial.png", g, height = 8, width = 6)

g <- ggplot(LLNRspace, aes(Longitude, Latitude)) + 
  geom_sf(data = coast, inherit.aes = FALSE) +
  coord_sf(xlim = c(-100, 30), ylim = c(-55, 65)) + 
  geom_point(shape = 21, aes(fill = p, colour = p), alpha = 0.85) +
  facet_wrap(vars(Decade), ncol = 3) +
  scale_color_viridis_c(trans = "log", name = "Distribution of hooks") +
  scale_fill_viridis_c(trans = "log", name = "Distribution of hooks") +
  theme(legend.position = "bottom",
        panel.spacing = unit(0, "in"))
ggsave("Figures/LLCE/LLNRspatial_p.png", g, height = 8, width = 6)

#g <- ggplot(LLNRspace %>% filter(FlagCode == "USA"), aes(Longitude, Latitude)) + 
#  geom_sf(data = coast, inherit.aes = FALSE) +
#  coord_sf(xlim = c(-100, 30), ylim = c(-45, 65)) + 
#  geom_point(shape = 21, aes(fill = n, colour = n), alpha = 0.85) +
#  facet_wrap(vars(Decade), ncol = 3) +
#  scale_color_viridis_c(trans = "log", name = "Number of records", breaks = c(1, 7, 54, 403)) +
#  scale_fill_viridis_c(trans = "log", name = "Number of records", breaks = c(1, 7, 54, 403)) +
#  theme(legend.position = "bottom",
#        panel.spacing = unit(0, "in"))
#ggsave("Figures/LLCE/LLNRspatial_USA.png", g, height = 8, width = 6)


LLNRspace <- LLCE %>% 
  filter(#FlagCode %in% c("JPN", "USA"),
         #TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         #SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")
         ) %>%
  group_by(YearC, Latitude, Longitude, QuadID) %>%
  summarise(Effort = sum(Eff1))


coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
g <- LLNRspace %>%
  filter(YearC >= 2010) %>%
  ggplot(aes(Longitude, Latitude)) + 
  geom_sf(data = coast, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_sf(xlim = range(LLNRspace$Longitude), ylim = range(LLNRspace$Latitude)) +
  geom_point(aes(colour = Effort), alpha = 0.85) + 
  facet_wrap(vars(YearC), ncol = 3) +
  scale_colour_viridis_c(trans = "log", breaks = log(LLNRspace$Effort) %>% pretty() %>% exp() %>% signif(2)) +
  theme_bw() + 
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45)) +
  labs(colour = "# Hooks")
ggsave("Figures/LLCE/spatial-effort-2010-2019.png", g, height = 8, width = 6)

g <- LLNRspace %>%
  filter(YearC %in% c(1956, seq(1960, 2010, 10), 2019)) %>%
  ggplot(aes(Longitude, Latitude)) + 
  geom_sf(data = coast, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_sf(xlim = range(LLNRspace$Longitude), ylim = range(LLNRspace$Latitude)) +
  geom_point(aes(colour = Effort)) + 
  facet_wrap(vars(YearC), ncol = 3) +
  scale_colour_viridis_c(trans = "log", breaks = log(LLNRspace$Effort) %>% pretty() %>% exp() %>% signif(2)) +
  theme_bw() + 
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45)) +
  labs(colour = "# Hooks")
ggsave("Figures/LLCE/spatial-effort-decadal.png", g, height = 8, width = 6)

g <- LLNRspace %>% 
  group_by(YearC, QuadID) %>%
  summarise(Effort = sum(Effort)) %>%
  ggplot(aes(YearC, Effort)) + 
  geom_line(aes(colour = as.factor(QuadID))) + 
  labs(x = "Year", y = "# Hooks", colour = "Quadrant")
ggsave("Figures/LLCE/spatial-effort-annual-quadrant.png", g, height = 3, width = 6)

g <- LLNRspace %>% 
  group_by(YearC, QuadID) %>%
  summarise(Effort = sum(Effort)) %>%
  ggplot(aes(YearC, Effort)) + 
  geom_col(colour = "black", aes(fill = as.factor(QuadID))) + 
  labs(x = "Year", y = "# Hooks", fill = "Quadrant")
ggsave("Figures/LLCE/spatial-effort-annual-quadrant-stack.png", g, height = 3, width = 6)

g <- LLNRspace %>% 
  group_by(YearC, QuadID) %>%
  summarise(Effort = sum(Effort)) %>%
  group_by(YearC) %>%
  mutate(p = Effort/sum(Effort)) %>%
  ggplot(aes(YearC, p)) + 
  geom_col(colour = "black", aes(fill = as.factor(QuadID))) + 
  labs(x = "Year", y = "Distribution of hooks", fill = "Quadrant")
ggsave("Figures/LLCE/spatial-effort-annual-quadrant-prop.png", g, height = 3, width = 6)




# Map longline effort annually

