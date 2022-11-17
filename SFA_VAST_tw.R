
# Spatial factor model using VAST (Tweedie distribution in likelihood)
source("SFA_VAST_fn.R")

LLCE <- local({
  dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases"
  read.csv(file.path(dir, "t2ce_LL.csv")) %>% 
    mutate(Latitude = ifelse(QuadID %in% 2:3, -Lat, Lat), 
           Longitude = ifelse(QuadID %in% 3:4, -Lon, Lon))
})

LLNRspace <- LLCE %>% 
  filter(FlagCode %in% c("JPN", "USA"),
         TimePeriodID <= 12, 
         Eff1Type == "NO.HOOKS", 
         SquareTypeCode == "5x5" | SquareTypeCode == "1x1",
         ifelse(DSetTypeID == "nw", CatchUnit == "nr", CatchUnit == "kg" | CatchUnit == "nr")) %>%
  mutate(Decade = paste0(floor(YearC/10) * 10, "s")) %>%
  dplyr::select(YearC, CatchUnit, SquareTypeCode, QuadID, Decade, FlagCode, 
                TimePeriodID, StrataID, Latitude, Longitude, BET, SWO, WHM, BUM, SMA, BSH, Eff1)

# 6 sp., 2 factors
DataFrame <- filter(LLNRspace, 
                    YearC >= 2010) %>%
  mutate(#Month = factor(TimePeriodID, levels = 1:12),
    Quarter = cut(TimePeriodID, 4, labels = 1:4),
    log_hook = log(Eff1)) %>%
  dplyr::select(YearC, FlagCode, StrataID, Quarter, Eff1, #log_hook, 
                SquareTypeCode, Latitude, Longitude, QuadID,
                BET, BUM, BSH, SMA, SWO, WHM) %>%
  reshape2::melt(id.vars = c(#"X", "Y", 
    "YearC", "FlagCode", "StrataID", "Quarter", "Eff1", #"log_hook", 
    "SquareTypeCode", "Latitude", "Longitude", "QuadID")) %>%
  mutate(Decade = floor(YearC/10) %>% `*`(10) %>% factor())

Return <- fit_VAST_model(DataFrame,
                         Aniso = TRUE,
                         Q1_formula = ~ FlagCode : variable,
                         Q2_formula = ~ FlagCode : variable + Quarter,
                         catchability_data = DataFrame %>% 
                           dplyr::select(FlagCode, Quarter, variable),
                         category_names = DataFrame$variable %>% table() %>% names(),
                         SD = TRUE,
                         run_model = TRUE,
                         Omega = c(0, 2),
                         Beta = c(0, 0),
                         ObsModel_ez = c(10, 2), 
                         CheckForErrors = FALSE)
saveRDS(Return, file = "SFA/SFA_VAST_6sp_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE,
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_tw", "/"))

plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_tw"))


# 6 sp, 3f
Return <- fit_VAST_model(DataFrame,
                         Aniso = TRUE,
                         Q1_formula = ~ FlagCode : variable,
                         Q2_formula = ~ FlagCode : variable + Quarter,
                         catchability_data = DataFrame %>% 
                           dplyr::select(FlagCode, Quarter, variable),
                         category_names = DataFrame$variable %>% table() %>% names(),
                         SD = TRUE,
                         run_model = TRUE,
                         Omega = c(0, 3),
                         Beta = c(0, 0),
                         ObsModel_ez = c(10, 2), 
                         CheckForErrors = FALSE)
saveRDS(Return, file = "SFA/SFA_VAST_6sp_3f_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE,
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_3f_tw", "/"))

plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_3f_tw"))

# 6 sp, 4f
Return <- fit_VAST_model(DataFrame,
                         Aniso = TRUE,
                         Q1_formula = ~ FlagCode : variable,
                         Q2_formula = ~ FlagCode : variable + Quarter,
                         catchability_data = DataFrame %>% 
                           dplyr::select(FlagCode, Quarter, variable),
                         category_names = DataFrame$variable %>% table() %>% names(),
                         SD = TRUE,
                         run_model = TRUE,
                         Omega = c(0, 4),
                         Beta = c(0, 0),
                         ObsModel_ez = c(10, 2), 
                         CheckForErrors = FALSE)
saveRDS(Return, file = "SFA/SFA_VAST_6sp_4f_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE,
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_4f_tw", "/"))

plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_4f_tw"))

# 6 sp, 5f
Return <- fit_VAST_model(DataFrame,
                         Aniso = TRUE,
                         Q1_formula = ~ FlagCode : variable,
                         Q2_formula = ~ FlagCode : variable + Quarter,
                         catchability_data = DataFrame %>% 
                           dplyr::select(FlagCode, Quarter, variable),
                         category_names = DataFrame$variable %>% table() %>% names(),
                         SD = TRUE,
                         run_model = TRUE,
                         Omega = c(0, 5),
                         Beta = c(0, 0),
                         ObsModel_ez = c(10, 2), 
                         CheckForErrors = FALSE)
saveRDS(Return, file = "SFA/SFA_VAST_6sp_5f_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE,
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_5f_tw", "/"))

plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_5f_tw"))

# 6 sp, 6f
Return <- fit_VAST_model(DataFrame,
                         Aniso = TRUE,
                         Q1_formula = ~ FlagCode : variable,
                         Q2_formula = ~ FlagCode : variable + Quarter,
                         catchability_data = DataFrame %>% 
                           dplyr::select(FlagCode, Quarter, variable),
                         category_names = DataFrame$variable %>% table() %>% names(),
                         SD = TRUE,
                         run_model = TRUE,
                         Omega = c(0, 6),
                         Beta = c(0, 0),
                         ObsModel_ez = c(10, 2), 
                         CheckForErrors = FALSE)
saveRDS(Return, file = "SFA/SFA_VAST_6sp_6f_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE, 
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_6f_tw", "/"))
#VAST_res <- list(dharmaRes = dharmaRes)
#dharmaRes = summary(Return, what = "residuals", #working_dir = working_dir, 
#                    type = 1)
#plot_quantile_residuals(dharmaRes = dharmaRes, fit = Return, 
#                        working_dir = working_dir, #year_labels = year_labels, 
#                        #years_to_plot = years_to_plot, 
#                        #n_cells_residuals = n_cells_residuals, 
#                        projargs = projargs, ...)
#DHARMa::plotResiduals(dharmaRes, form = DataFrame$variable)
#DHARMa::plotResiduals(dharmaRes, form = DataFrame$FlagCode)

plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_6f_tw"))



# Independent loadings (L_omega2_z is a diagonal matrix)
Return <- local({
  fit <- fit_VAST_model(DataFrame,
                        Aniso = TRUE,
                        Q1_formula = ~ FlagCode : variable,
                        Q2_formula = ~ FlagCode : variable + Quarter,
                        catchability_data = DataFrame %>% 
                          dplyr::select(FlagCode, Quarter, variable),
                        category_names = DataFrame$variable %>% table() %>% names(),
                        SD = TRUE,
                        run_model = FALSE,
                        Omega = c(0, 6),
                        Beta = c(0, 0),
                        ObsModel_ez = c(10, 2), 
                        CheckForErrors = FALSE)
  category_names <- DataFrame$variable %>% table() %>% names()
  
  
  L_omega2_z <- fit$tmb_list$Parameters$L_omega2_z
  
  map_L <- rep(NA, length(L_omega2_z))
  map_ind <- cumsum(1:length(category_names))
  map_L[map_ind] <- 1:length(category_names)
  
  L_omega2_z[-map_ind] <- 0
  
  fit$tmb_list$Parameters$L_omega2_z <- abs(L_omega2_z)
  fit$tmb_list$Map$L_omega2_z <- factor(map_L)
  
  data_list <- fit$tmb_list$Obj$env$data
  fit$tmb_list$Obj <- TMB::MakeADFun(data = data_list,
                                     parameters = fit$tmb_list$Parameters,
                                     map = fit$tmb_list$Map,
                                     random = fit$tmb_list$Random,
                                     DLL = "VAST_v13_1_0",
                                     silent = TRUE)
  
  fit$parameter_estimates <- TMBhelper::fit_tmb(obj = fit$tmb_list$Obj,
                                               #lower = tmb_list$Lower,
                                               #upper = tmb_list$Upper,
                                               control = list(eval.max = 10000, iter.max = 10000),
                                               loopnum = 0,
                                               quiet = TRUE,
                                               getsd = FALSE)
  fit$parameter_estimates$SD <- TMB::sdreport(fit$tmb_list$Obj)
  
  fit$Report <- fit$tmb_list$Obj$report(fit$tmb_list$Obj$env$last.par.best)
  
  fit
})
saveRDS(Return, file = "SFA/SFA_VAST_6sp_ind_tw.rds")

VAST_res <- plot_VAST_results(Return, RotationMethod = "Varimax", check_residuals = TRUE, 
                              working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_ind_tw", "/"))


plot_VAST_gg(Return, 
             country_list = DataFrame$FlagCode %>% unique() %>% new("list", .),
             country_index = DataFrame$FlagCode %>% factor() %>% as.integer() %>% `-`(1),
             strata_names = DataFrame$StrataID,
             dharma_resid = VAST_res$dharmaRes$scaledResiduals,
             working_dir = file.path(getwd(), "Figures", "SFA_VAST", "6sp_ind_tw"))


# AIC
model_dir <- paste0("SFA/SFA_VAST_6sp_", c("tw", "3f_tw", "4f_tw", "5f_tw", "6f_tw", "ind_tw"), ".rds")
AIC <- lapply(model_dir, function(x) {
  fit <- readRDS(x)
  
  data.frame(k = fit$Report$Omegainput2_gf %>% ncol(),
             nll = fit$parameter_estimates$objective,
             nll_marginal = fit$Report$jnll_comp %>% sum(),
             nfixed = length(fit$parameter_estimates$par),
             nrandom = fit$parameter_estimates$SD$par.random %>% length(),
             pdHess = fit$parameter_estimates$SD$pdHess) %>%
    mutate(AIC = 2 * nll + 2 * nfixed)
}) %>% bind_rows() %>%
  mutate(DAIC = AIC - min(AIC),
         AIC2 = 2 * nll + 2 * (nfixed + nrandom),
         DAIC2 = AIC2 - min(AIC2))
write.csv(AIC, file = 'SFA/AIC.csv')


# Loading matrix for 2-factor model
Return2 <- readRDS(file = "SFA/SFA_VAST_6sp_tw.rds")
g <- varimax(Return2$Report$L_omega2_cf)$loadings %>% structure(class = "matrix") %>% reshape2::melt() %>%
  reshape2::dcast(list("Category", "Var2")) %>%
  ggplot(aes(`1`, `2`, label = Category)) + 
  geom_point() +
  ggrepel::geom_text_repel() + 
  expand_limits(x = 0, y = 0) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Factor 1", y = "Factor 2")
ggsave("Figures/SFA_VAST/6sp_tw/Loadings_geometric.png", g, height = 4, width = 4)
