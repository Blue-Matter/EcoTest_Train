
source("00_functions.R")


t1 <- get_t1("SAI")
t2sz <- read.csv(file.path(dir, "t2sz_SAI1950-19.csv"))

# Note there are "unkn" catches in SampAreaCode
W <- c("BIL91", "BIL92", "BIL93", "BIL94A", "BIL96")
E <- c("BIL94B", "BIL94C", "BIL95", "BIL97")

final_size_interval <- seq(5, 335, 5)

#### Catch
cat_E <- t1 %>% dplyr::filter(SampAreaCode %in% E)

# By gear - mix of gillnet and longline
cat_E %>% 
  group_by(YearC, GearGrp) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0) %>%
  ggplot(aes(YearC, catch, colour = GearGrp)) + 
  geom_point() + geom_line()

## By flag and gear
cat_E %>%
  group_by(YearC, GearGrp, FlagName) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0) %>%
  ggplot(aes(YearC, catch, colour = FlagName)) + facet_wrap(~GearGrp, scales = "free_y") + 
  geom_point() + geom_line()

## Longline by flag - predominantly Spain, Japan
cat_E %>% dplyr::filter(GearGrp == "LL") %>%
  group_by(YearC, FlagName) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0) %>%
  ggplot(aes(YearC, catch)) + facet_wrap(~FlagName) + 
  geom_point() + geom_line()

## Gillnet by flag - Ghana
cat_E %>% dplyr::filter(GearGrp == "GN") %>%
  group_by(YearC, FlagName) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0) %>%
  ggplot(aes(YearC, catch)) + facet_wrap(~FlagName) + 
  geom_point() + geom_line()


#### CAL
# Re-do length bins
cal_W <- t2sz %>% dplyr::filter(SampAreaCode %in% W, FreqTypeCode == "LJFL") %>% 
  mutate(ClassFrq2 = ifelse(FreqsGroup == "cm (ll)", ClassFrq + 0.5 * SzInterval, 
                            ifelse(FreqsGroup == "cm(ul)", ClassFrq - 0.5 * SzInterval, ClassFrq))) %>%
  mutate(ClassFrq3 = final_size_interval[findInterval(ClassFrq2, final_size_interval)]) %>%
  group_by(YearC, GearGrpCode, FlagName, ClassFrq3) %>% summarise(Nr = sum(Nr))

# By gear
cal_W_LL <- cal_W %>% dplyr::filter(GearGrpCode == "LL")
cal_W_RR <- cal_W %>% dplyr::filter(GearGrpCode == "RR")
cal_W_GN <- cal_W %>% dplyr::filter(GearGrpCode == "GN")
cal_W_SP <- cal_W %>% dplyr::filter(GearGrpCode == "SP")

ggplot(cal_W_LL, aes(ClassFrq3, Nr, colour = FlagName)) + facet_wrap(~ YearC, scales = "free_y") + 
  geom_point() + geom_line()

cal_W_LL <- cal_W %>% dplyr::filter(GearGrpCode == "LL") %>% 
  group_by(YearC, ClassFrq3) %>% summarise(Nr = sum(Nr))
ggplot(cal_W_LL, aes(ClassFrq3, Nr)) + facet_wrap(~ YearC, scales = "free_y") + 
  geom_point() + geom_line()

cal_W_OTH <- cal_W %>% dplyr::filter(GearGrpCode != "LL") %>% 
  group_by(YearC, ClassFrq3) %>% summarise(Nr = sum(Nr))
ggplot(cal_W_OTH, aes(ClassFrq3, Nr)) + facet_wrap(~ YearC, scales = "free_y") + 
  geom_point() + geom_line()


cal_matrix <- cal_W_LL %>% reshape2::acast(YearC ~ ClassFrq3, fill = 0)

plot_composition(as.numeric(rownames(cal_matrix)), cal_matrix, CAL_bins = as.numeric(colnames(cal_matrix)))



cat_E <- t1 %>% dplyr::filter(SampAreaCode %in% E) %>%
  group_by(YearC, GearGrp) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0)
ggplot(cat_E, aes(YearC, catch, colour = GearGrp)) + 
  geom_point() + geom_line()

cat <- t1 %>% group_by(YearC, FlagName, GearGrp) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0)

cat_gear <- t1 %>% group_by(YearC, FlagName, GearGrp) %>% summarise(catch = sum(Qty_t)) %>% dplyr::filter(catch > 0)
