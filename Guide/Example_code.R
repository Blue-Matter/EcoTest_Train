# Test Script

library(ETest)

dat = makerawdata_2(TD, sno=1, isBrel=F,
                   inc_Irel = F, inc_I = F, inc_CR = T, inc_CAL = T, inc_CAA = F,
                   stock_in = 1:3, fleet_in = 1:3, Bmin = 0.05)

out = train_NN(dat, c(4, 3), 20, T, NULL)



mody = load_mod(1)

# process data & depletion
# subset TD
# train NN & get standardization mu and sd
# predict depletion from data
# retros (strip processed data and calc depletion)

