
# Look at the fraction longline catch of BSH by country

dir <- "G:\\Shared drives\\BM shared\\1. Projects\\EcoTest\\Databases"

t1 <- readxl::read_excel(file.path(dir, "t1nc-20211002.xlsx"), skip = 3) %>% dplyr::filter(Species == "BSH")

# These are the fleet assignments in Dean Courtney's SS model
calc_F <- function(x = c("EU-Portugal", "Japan")) {
  ss_fleets <- c("EU-España" = "F1", "EU-Portugal" = "F1", "Japan" = "F2", "Chinese Taipei" = "F3", "USA" = "F4", "Venezuela" = "F5",
                 "Canada" = "F6", "China PR" = "F7", "Belize" = "F8")
  out <- ss_fleets[match(x, names(ss_fleets))]
  out[is.na(out)] <- "F9" # Other
  out
}

c <- t1 %>% dplyr::filter(Stock == "ATN") %>% mutate(FF = calc_F(FlagName)) %>% group_by(YearC, FF) %>%
  summarise(catch = sum(Qty_t))
cp <- t1 %>% dplyr::filter(Stock == "ATN") %>% mutate(FF = calc_F(FlagName)) %>% group_by(YearC, FF) %>%
  summarise(pLL = sum(Qty_t[GearGrp == "LL"]/sum(Qty_t)))

ggplot(c, aes(YearC, catch)) + geom_line() + facet_wrap(~FF) # Total catch by fleet, mostly F1
ggplot(cp, aes(YearC, pLL)) + geom_line() + facet_wrap(~FF) # Proportion longline
