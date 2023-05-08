
library(VAST)
library(TMB)
library(tidyverse)
theme_set(theme_bw())

#dyn.load(system.file("executables", package = "VAST") %>% file.path("VAST_v13_1_0.dll"))

### TODO: Replace extrapolation grid with ecoregion shapefile
#ocean <- rnaturalearth::ne_download(scale = 110, type = "ocean", category = "physical")

fit_VAST_model <- function(DataFrame, 
                           Aniso = FALSE,
                           category_names = c("BET", "SWO"), 
                           SD = FALSE, 
                           X1_formula = ~ 0,
                           X2_formula = ~ 0,
                           covariate_data = NULL,
                           Q1_formula = ~ 0,
                           Q2_formula = ~ 0,
                           catchability_data = NULL,
                           Q1config_k = NULL,
                           Q2config_k = NULL,
                           n_factors = 2,
                           run_model = TRUE,
                           n_x = 150,
                           ObsModel_ez = c(10, 2),
                           CheckForErrors = TRUE,
                           Omega = c(0, 0),
                           Epsilon = c(0, 0),
                           Beta = c("IID", 0),
                           RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)) {
  
  settings <- FishStatsUtils::make_settings(n_x = n_x,
                                            Region = "user",
                                            purpose = "ordination",
                                            fine_scale = FALSE,
                                            n_categories = n_factors,
                                            max_cells = n_x)
  
  input_grid <- DataFrame %>%
    group_by(Latitude, Longitude) %>% 
    summarise(Area_km2 = 1) %>%
    rename(Lat = Latitude, Lon = Longitude)
  
  #xy_ind <- sf::st_intersects(input_grid %>% sf::st_as_sf(coords = c("Lon", "Lat"), crs = 4326),
  #                            ocean %>% as("sf")) %>%
  #  sapply(function(x) length(x) > 0)
  
  extrapolation_list <- make_extrapolation_info(Region = "user",
                                                max_cells = n_x,
                                                #input_grid = input_grid[xy_ind, ] %>% as.matrix(),
                                                input_grid = input_grid %>% as.matrix(),
                                                Save_results = FALSE,
                                                nstart = 1)
  
  spatial_list <- make_spatial_info(n_x = n_x,
                                    Lon_i = DataFrame$Longitude,
                                    Lat_i = DataFrame$Latitude,
                                    Extrapolation_List = extrapolation_list,
                                    Method = "Mesh",
                                    #Save_Results = TRUE,
                                    fine_scale = TRUE)
  
  settings$FieldConfig["Omega", ] <- Omega
  settings$FieldConfig["Epsilon", ] <- Epsilon
  settings$FieldConfig["Beta", ] <- Beta
  
  #### Year effect:
  # Default RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
  # 0 = each year as fixed effect (FE)
  # 1 = each year is IID as random effect (RE)
  # 2 = Random walk RE
  # 3 = single intercept FE
  # 4 = AR1 FE
  # 5 = AR1 FE by category
  settings$RhoConfig <- RhoConfig
  
  data_list <- make_data(b_i = as_units(DataFrame$value, "count"),                  # Response
                         a_i = as_units(DataFrame$Eff1, unitless),    # Offset
                         t_i = DataFrame$YearC,                                     # Time
                         c_iz = DataFrame$variable %>% as.integer() %>% `-`(1),     # Category (species)
                         e_i = model.matrix(~ 0 + variable : FlagCode, DataFrame) %>%  # SD per species x country combination
                           apply(1, function(x) which(x == 1)) %>% as.numeric() %>% `-`(1),
                         FieldConfig = settings$FieldConfig,
                         OverdispersionConfig = settings$OverdispersionConfig, 
                         RhoConfig = settings$RhoConfig,
                         VamConfig = settings$VamConfig, 
                         ObsModel_ez = ObsModel_ez,
                         covariate_data = covariate_data,
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         catchability_data = catchability_data,  
                         Q1_formula = Q1_formula,
                         Q2_formula = Q2_formula,
                         Q1config_k = Q1config_k,
                         Q2config_k = Q2config_k,
                         spatial_list = spatial_list,
                         Options = settings$Options,
                         Aniso = Aniso,
                         CheckForErrors = CheckForErrors)
  
  tmb_list <- make_model(TmbData = data_list,
                         Method = "Mesh",
                         Version = settings$Version,
                         RhoConfig = settings$RhoConfig,
                         loc_x = spatial_list$loc_x
  )
  tmb_list$Obj$env$beSilent()
  
  if (run_model) {
    
    message("Fitting model...")
    parameter_estimates1 <- TMBhelper::fit_tmb(obj = tmb_list$Obj,
                                               #lower = tmb_list$Lower,
                                               #upper = tmb_list$Upper,
                                               control = list(eval.max = 10000, iter.max = 10000),
                                               loopnum = 0,
                                               quiet = TRUE,
                                               getsd = FALSE)
    
    message("Model finished.")
    
    if(SD) {
      message("Calculating covariance matrix...")
      parameter_estimates1[["SD"]] <- TMB::sdreport(tmb_list$Obj)
      message("TMB::sdreport finished.")
    }
    ParHat = tmb_list$Obj$env$parList(parameter_estimates1$par)
  } else {
    parameter_estimates1 <- list()
    ParHat <- tmb_list$Obj$env$parList()
  }
  
  Report <- tmb_list$Obj$report(tmb_list$Obj$env$last.par.best)
  Report = amend_output(Report = Report, 
                        TmbData = data_list, 
                        Map = tmb_list$Map, 
                        Sdreport = NULL, 
                        year_labels = DataFrame$YearC %>% unique(), 
                        category_names = category_names, 
                        extrapolation_list = extrapolation_list)
  
  Return = list(data_frame = data.frame(Lat_i = DataFrame$Latitude, 
                                        Lon_i = DataFrame$Longitude,
                                        a_i = as_units(DataFrame$Eff1, unitless), 
                                        v_i = rep(0, nrow(DataFrame)),
                                        b_i = as_units(DataFrame$value, "count"),
                                        t_i = DataFrame$YearC,
                                        c_iz = DataFrame$variable %>% as.integer() %>% `-`(1)), 
                extrapolation_list = extrapolation_list, 
                spatial_list = spatial_list,
                data_list = data_list, 
                tmb_list = tmb_list, 
                parameter_estimates = parameter_estimates1, 
                Report = Report,
                ParHat = ParHat,
                year_labels = DataFrame$YearC %>% unique(),
                years_to_plot = 1:length(unique(DataFrame$YearC)),
                category_names = category_names,
                settings = settings) %>%
    structure(class = "fit_model")
  ##input_args = input_args, 
  #X1config_cp = X1config_cp, 
  #X2config_cp = X2config_cp, 
  #covariate_data = covariate_data, 
  #X1_formula = X1_formula, 
  #X2_formula = X2_formula, 
  #Q1config_k = Q1config_k, 
  #Q2config_k = Q1config_k, 
  #catchability_data = catchability_data, 
  #Q1_formula = Q1_formula, 
  #Q2_formula = Q2_formula)
  
  return(Return)
}

plot_VAST_results <- function (fit, settings = fit$settings, plot_set = 3, working_dir = paste0(getwd(), 
                                                                           "/"), year_labels = fit$year_labels, years_to_plot = fit$years_to_plot, 
          category_names = fit$category_names, strata_names = fit$strata_names, 
          use_biascorr = TRUE, map_list, check_residuals = TRUE, projargs = "+proj=longlat", 
          zrange, n_samples = 100, calculate_relative_to_average = FALSE, 
          type = 1, n_cells = NULL, n_cells_residuals = NULL, RotationMethod = "PCA", 
          quantiles = c(0.05, 0.5, 0.95), ...) 
{
  if (is.null(fit$Report)) 
    stop("`fit$Report` is missing, please check inputs")
  if (is.null(category_names)) 
    category_names = paste0("Category_", 1:fit$data_list$n_c)
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  message("\n### Creating plots in directory ", working_dir)
  message("\n### Making plots of data availability and knots")
  plot_data_args = list(...)
  plot_data_args = combine_lists(input = plot_data_args, args_to_use = formalArgs(plot_data), 
                                 default = list(Extrapolation_List = fit$extrapolation_list, 
                                                Spatial_List = fit$spatial_list, Lat_i = fit$data_frame[, 
                                                                                                        "Lat_i"], Lon_i = fit$data_frame[, "Lon_i"], 
                                                Year_i = fit$data_frame[, "t_i"], PlotDir = working_dir, 
                                                year_labels = year_labels, projargs = projargs))
  do.call(what = plot_data, args = plot_data_args)
  if (missing(map_list)) {
    message("\n### Obtaining default settings for plotting maps")
    map_list = make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, 
                             Extrapolation_List = fit$extrapolation_list)
  }
  message("\n### Making plot of anisotropy")
  plot_anisotropy(FileName = paste0(working_dir, "Aniso.png"), 
                  Obj = fit$tmb_list$Obj)
  plot_biomass_index_args = list(...)
  #if (!is.null(fit$parameter_estimates$SD)) {
  #  message("\n### Making plot of abundance index")
  #  plot_biomass_index_args = combine_lists(input = plot_biomass_index_args, 
  #                                          args_to_use = formalArgs(plot_biomass_index), default = list(DirName = working_dir, 
  #                                                                                                       fit = fit, year_labels = year_labels, years_to_plot = years_to_plot, 
  #                                                                                                       use_biascorr = use_biascorr, category_names = category_names, 
  #                                                                                                       strata_names = strata_names))
  #  Index = do.call(what = plot_biomass_index, args = plot_biomass_index_args)
  #}
  #else {
  #  Index = "Not run"
  #  message("\n### Skipping plot of abundance index; must re-run with standard errors to plot")
  #}
  #if (!is.null(fit$parameter_estimates$SD) & fit$data_list$n_c > 
  #    1) {
  #  message("\n### Making plot of composition data")
  #  Proportions = calculate_proportion(TmbData = fit$data_list, 
  #                                     Index = Index, year_labels = year_labels, years_to_plot = years_to_plot, 
  #                                     use_biascorr = use_biascorr, category_names = category_names, 
  #                                     DirName = working_dir)
  #}
  #else {
  #  Proportions = "Not run"
  #  message("\n### Skipping plot of composition data; must re-run with standard errors and multiple categories to plot")
  #}
  #if (!is.null(fit$parameter_estimates$SD)) {
  #  message("\n### Making plot of spatial indices")
  #  Range = plot_range_index(Report = fit$Report, TmbData = fit$data_list, 
  #                           Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), 
  #                           PlotDir = working_dir, year_labels = year_labels, 
  #                           years_to_plot = years_to_plot, use_biascorr = use_biascorr, 
  #                           category_names = category_names)
  #}
  #else {
  #  Range = "Not run"
  #  message("\n### Skipping plot of spatial indices; must re-run with standard errors to plot")
  #}
  #if ("jointPrecision" %in% names(fit$parameter_estimates$SD) & 
  #    n_samples > 0 & fit$data_list$Options_list$Options["Calculate_Range"] == 
  #    TRUE) {
  #  message("\n### Making plot of range edges")
  #  Edge = plot_range_edge(Obj = fit$tmb_list$Obj, Sdreport = fit$parameter_estimates$SD, 
  #                         working_dir = working_dir, year_labels = year_labels, 
  #                         years_to_plot = years_to_plot, category_names = category_names, 
  #                         n_samples = n_samples, quantiles = quantiles, calculate_relative_to_average = calculate_relative_to_average)
  #}
  #else {
  #  Edge = "Not run"
  #  message("\n### Skipping plot of range edge; only possible if `getJointPrecision=TRUE`, `Options['Calculate_Range']=TRUE`, and `n_samples`>0")
  #}
  message("\n### Making plots of spatial predictions")
  plot_maps_args = list(...)
  plot_maps_args = combine_lists(input = plot_maps_args, default = list(fit = fit, 
                                                                        plot_set = plot_set, category_names = category_names, 
                                                                        PlotDF = map_list[["PlotDF"]], MapSizeRatio = map_list[["MapSizeRatio"]], 
                                                                        working_dir = working_dir, year_labels = year_labels, 
                                                                        years_to_plot = years_to_plot, legend_x = map_list[["Legend"]]$x/100, 
                                                                        legend_y = map_list[["Legend"]]$y/100, projargs = projargs, 
                                                                        n_cells = n_cells))
  Dens_xt = do.call(what = plot_maps, args = plot_maps_args)
  message("\n### Making plots for factors (if present)")
  plot_factors_args = list(...)
  plot_factors_args = combine_lists(input = plot_factors_args, 
                                    default = list(fit = fit, mapdetails_list = map_list, 
                                                   projargs = projargs, n_cells = n_cells, RotationMethod = RotationMethod, 
                                                   plotdir = working_dir, category_names = category_names))
  Factors = do.call(what = plot_factors, args = plot_factors_args)
  if (check_residuals == TRUE) {
    message("\n### Making quantile residuals using conditional simulation and package DHARMa")
    dharmaRes = summary(fit, what = "residuals", working_dir = working_dir, 
                        type = type, ...)
    message("\n### Plotting quantile residuals ")
    plot_quantile_residuals(dharmaRes = dharmaRes, fit = fit, 
                            working_dir = working_dir, year_labels = year_labels, 
                            years_to_plot = years_to_plot, n_cells_residuals = n_cells_residuals, 
                            projargs = projargs, ...)
  }
  else {
    message("\n### Skipping quantile residuals using conditional simulation and package DHARMa")
    message("\n### Skipping plot of quantile residuals ")
    dharmaRes = NULL
  }
  Return = list(dharmaRes = dharmaRes, #Index = Index, Proportions = Proportions, 
                #Range = Range, 
                Dens_xt = Dens_xt, #Edge = Edge, 
                map_list = map_list, 
                plot_maps_args = plot_maps_args, 
                #plot_biomass_index_args = plot_biomass_index_args, 
                Factors = Factors)
  return(invisible(Return))
}


plot_VAST_gg <- function (fit, country_list,
                          country_index,
                          do_plot = c("data", "factor", "density"),
                          plot_set = 3, working_dir = paste0(getwd(), "/"),
                          year_labels = fit$year_labels, years_to_plot = fit$years_to_plot, 
                          category_names = fit$category_names, 
                          strata_names, 
                          dharma_resid = NA,
                          ...) {
  
                          #use_biascorr = TRUE, map_list, check_residuals = TRUE, projargs = "+proj=longlat", 
                          #zrange, n_samples = 100, calculate_relative_to_average = FALSE, 
                          #type = 1, n_cells = NULL, n_cells_residuals = NULL, RotationMethod = "PCA", 
                          #quantiles = c(0.05, 0.5, 0.95), ...) {
  
  #coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  ecoregion_shp <- local({
    dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/ICCAT_shp"
    sf::st_read(file.path(dir, "ICCAT_draft_ecoregions_May 2022.shp"))
  })
  
  sf::sf_use_s2(FALSE)
  # plot data
  if ("data" %in% do_plot) {
    dat <- fit$data_frame %>%
      mutate(b_i = as.numeric(b_i), a_i = as.numeric(a_i), 
             country = sapply(country_index + 1, function(x) country_list[[x]]),
             sp = category_names[c_iz + 1],
             D_i = fit$Report$D_i,
             strata = strata_names,
             Dobs_i = b_i/a_i, 
             log_Dobs_i = log(Dobs_i),
             dharma_resid = dharma_resid)
    
    dat_cc <- dat %>%
      filter(b_i > 0) %>%
      reshape2::dcast(strata + country ~ sp, value.var = "log_Dobs_i")
    g <- GGally::ggpairs(dat_cc,
                         columns = 3:ncol(dat_cc), 
                         xlab = "log CPUE", ylab = "log CPUE",
                         aes(colour = country, alpha = 0.2))
    ggsave(file.path(working_dir, "obs_pairs.png"), g, height = 6, width = 6)
    
    # Sample correlations
    #D_i_corr <- reshape2::acast(dat, list("strata", "sp"), value.var = "D_i")
    #samp_corr <- lapply(category_names, function(x) {
    #  corrdat %>% 
    #    filter(sp == x) %>%
    #    pull(Dobs_i)
    #  
    #})
    #
    #cov_j <- fit$Report[[covvec[i]]]
    #for(ii in 1:(nrow(cov_j)-1)) {
    #  for(j in (ii+1):nrow(cov_j)) cov_j[ii, j] <- NA 
    #}
    #
    #corr_j <- cov_j %>% 
    #  cov2cor() %>%
    #  structure(dimnames = list(Var1 = category_names, Var2 = category_names)) %>%
    #  signif(2) %>%
    #  reshape2::melt() %>%
    #  filter(!is.na(value)) %>% 
    #  mutate(Var1 = factor(Var1, levels = category_names %>% rev()),
    #         Var2 = factor(Var2, levels = category_names))
    #
    #g <- corr_j %>%
    #  ggplot(aes(Var2, Var1)) + 
    #  geom_tile(height = 1, width = 1, colour = "black", aes(fill = value), alpha = 0.5) + 
    #  geom_text(aes(label = value)) +
    #  coord_cartesian(expand = FALSE) +
    #  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) + 
    #  gfplot::theme_pbs() +
    #  theme(legend.position = "none") + 
    #  labs(x = "", y = "") + 
    #  guides(fill = "none")
    #ggsave(file.path(working_dir, paste0(lab[i], "_corr.png")), g, height = 3, width = 4)
    
    for(cc in unique(dat$c_iz)) { # Category (species)
      
      dat_cc <- dat %>% filter(c_iz == cc)
      
      g <- dat_cc %>% 
        ggplot(aes(b_i, D_i)) + 
        geom_abline(slope = 1, intercept = 0) + 
        geom_point(alpha = 0.2) + 
        facet_wrap(vars(country), scales = "free") +
        labs(x = "Observed", y = "Predicted") +
        ggtitle(category_names[cc + 1])
      ggsave(file.path(working_dir, paste0("obs_pred_", category_names[cc+1], ".png")), g, height = 6, width = 6)
      
      g <- dat_cc %>% 
        ggplot(aes(b_i, log(b_i/D_i))) + 
        geom_hline(yintercept = 0, linetype = 2) + 
        geom_point(alpha = 0.2) + 
        facet_wrap(vars(country), scales = "free") +
        labs(x = "Observed", y = "log Residual") +
        ggtitle(category_names[cc + 1])
      ggsave(file.path(working_dir, paste0("residual_", category_names[cc+1], ".png")), g, height = 6, width = 6)
      
      for(country in unique(country_index)) { # Country
        dat_out <- dat[country_index == country, ] %>%
          filter(c_iz == cc)
        
        # Data
        g <- dat_out %>%
          ggplot(aes(Lon_i, Lat_i)) + 
          #geom_sf(data = coast, inherit.aes = FALSE) +
          geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
          coord_sf(xlim = range(dat$Lon_i), ylim = range(dat$Lat_i)) +
          geom_jitter(aes(shape = b_i == 0, fill = log(b_i/a_i))) +
          facet_wrap(vars(t_i)) +
          #facet_grid(vars(t_i), vars(season)) + 
          labs(x = "Longitude", y = "Latitude", fill = "log CPUE") +
          theme_bw() +
          theme(panel.spacing = unit(0, "in")) + 
          scale_fill_viridis_c() + 
          scale_shape_manual(values = c(21, 4)) + 
          guides(shape = "none") + 
          ggtitle(category_names[cc + 1] %>% paste(country_list[[country + 1]]))
        ggsave(file.path(working_dir, paste0("data_map_", category_names[cc + 1], "_", country_list[[country + 1]], "_ecoregion.png")),
               g, height = 6, width = 8)
        
        # Residual
        g <- dat_out %>%
          ggplot(aes(Lon_i, Lat_i)) + 
          #geom_sf(data = coast, inherit.aes = FALSE) +
          geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
          coord_sf(xlim = range(dat$Lon_i), ylim = range(dat$Lat_i)) +
          geom_jitter(shape = 21, aes(fill = log(b_i/D_i))) +
          facet_wrap(vars(t_i)) +
          #facet_grid(vars(t_i), vars(season)) + 
          scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
          labs(x = "Longitude", y = "Latitude", fill = "Residual") +
          theme_bw() +
          theme(panel.spacing = unit(0, "in")) + 
          ggtitle(category_names[cc + 1] %>% paste(country_list[[country + 1]]))
        ggsave(file.path(working_dir, paste0("residual_map_", category_names[cc + 1], "_", country_list[[country + 1]], "_ecoregion.png")),
               g, height = 6, width = 8)
        
        # DHARMa
        #g <- dat_out %>%
        #  ggplot(aes(Lon_i, Lat_i)) + 
        #  #geom_sf(data = coast, inherit.aes = FALSE) +
        #  geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
        #  coord_sf(xlim = range(dat$Lon_i), ylim = range(dat$Lat_i)) +
        #  geom_jitter(shape = 21, aes(fill = dharma_resid)) +
        #  facet_wrap(vars(t_i)) +
        #  #facet_grid(vars(t_i), vars(season)) + 
        #  scale_fill_gradient2(high = scales::muted("red"), midpoint = 0.5, low = scales::muted("blue")) + 
        #  labs(x = "Longitude", y = "Latitude", fill = "DHARMa/nResidual") +
        #  theme_bw() +
        #  theme(panel.spacing = unit(0, "in")) + 
        #  ggtitle(category_names[cc + 1] %>% paste(country_list[[country + 1]]))
        #ggsave(file.path(working_dir, paste0("dharma_resid_map_", category_names[cc + 1], "_", country_list[[country + 1]], "_ecoregion.png")),
        #       g, height = 6, width = 8)
        
      }
    }
    
    message("Finished data and residuals.")
  }
  
  if("factor" %in% do_plot) {
    # Factor loadings - epsilon 2
    lab <- c("Epsilon1", "Epsilon2", "Omega1", "Omega2")
    Lvec <- c("L_epsilon1_cf", "L_epsilon2_cf", "L_omega1_cf", "L_omega2_cf")
    Om <- c("Epsiloninput1_sft", "Epsiloninput2_sft", "Omegainput1_sf", "Omegainput2_sf")
    covvec <- paste0("lowercov_uppercor_", c("epsilon1", "epsilon2", "omega1", "omega2"))
    
    for(i in 1:length(Lvec)) {
      Loadings <- fit$Report[[Lvec[i]]]
      
      if(length(Loadings)) {
        
        # Corr matrix
        cov_j <- fit$Report[[covvec[i]]]
        for(ii in 1:(nrow(cov_j)-1)) {
          for(j in (ii+1):nrow(cov_j)) cov_j[ii, j] <- NA 
        }
        
        corr_j <- cov_j %>% 
          cov2cor() %>%
          structure(dimnames = list(Var1 = category_names, Var2 = category_names)) %>%
          signif(2) %>%
          reshape2::melt() %>%
          filter(!is.na(value)) %>% 
          mutate(Var1 = factor(Var1, levels = category_names %>% rev()),
                 Var2 = factor(Var2, levels = category_names))
        
        g <- corr_j %>%
          ggplot(aes(Var2, Var1)) + 
          geom_tile(height = 1, width = 1, colour = "black", aes(fill = value), alpha = 0.5) + 
          geom_text(aes(label = value)) +
          coord_cartesian(expand = FALSE) +
          scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) + 
          gfplot::theme_pbs() +
          theme(legend.position = "none") + 
          labs(x = "", y = "") + 
          guides(fill = "none")
        ggsave(file.path(working_dir, paste0(lab[i], "_corr.png")), g, height = 3, width = 4)
        
        
        
        
        
        
        Psi_sjt <- fit$Report[[Om[i]]]
        logkappa = unlist(fit$ParHat[c("logkappa1", "logkappa2")])[c(1, 2, 1, 2)[i]]
        stopifnot(fit$data_list$Options_list$Options_vec[8] == 0)
        tau = 1/(exp(logkappa) * sqrt(4 * pi))
        
        Var_rot = rotate_factors(L_pj = Loadings, 
                                 Psi_sjt = Psi_sjt/tau, 
                                 RotationMethod = "Varimax", 
                                 testcutoff = 1e-4, 
                                 quiet = TRUE)
        
        # Loadings
        L_rot <- Var_rot$L_pj_rot %>% structure(dimnames = list(Category = category_names, Factor = 1:ncol(.))) %>%
          reshape2::melt() %>%
          mutate(Category = factor(Category, levels = rev(category_names)))
        
        g <- L_rot %>%
          ggplot(aes(factor(Factor), Category)) + 
          geom_tile(aes(fill = value), alpha = 0.7) +
          geom_text(aes(label = round(value, 2))) + 
          scale_fill_viridis_c() + 
          coord_cartesian(expand = FALSE) +
          guides(fill = "none") +
          theme(legend.position = "bottom") + 
          labs(x = "Factor", y = "Species", fill = expression("Loadings"~psi*minute))
        ggsave(file.path(working_dir, paste0(lab[i], "_Loadings.png")), g, height = 3, width = 4)
        
        # Factors
        Psi_rot <- Var_rot$Psi_rot %>% 
          structure(dimnames = list(Cell = 1:dim(.)[1], Factor = 1:dim(.)[2], 
                                    Year = if(i < 3) {
                                      fit$year_labels
                                    } else 1)) %>%
          reshape2::melt() %>%
          left_join(fit$spatial_list$latlon_s %>% as.data.frame() %>% mutate(Cell = 1:nrow(.)), by = "Cell")
        
        Psi_interp <- local({
          # Interpolation should be in UTM in the model
          coord_out <- expand.grid(xo = seq(min(Psi_rot$Lon), max(Psi_rot$Lon), by = 1),
                                   yo = seq(min(Psi_rot$Lat), max(Psi_rot$Lat), by = 1))
          full <- lapply(unique(Psi_rot$Factor), function(ff) {
            lapply(unique(Psi_rot$Year), function(yy) {
              mesh_x <- filter(Psi_rot, Factor == ff, Year == yy)
              interp::interp(x = mesh_x$Lon, 
                             y = mesh_x$Lat, 
                             z = mesh_x$value, 
                             xo = coord_out$xo,
                             yo = coord_out$yo,
                             output = "points") %>%
                bind_cols() %>% mutate(Factor = ff, Year = yy)
            }) %>% bind_rows()
          }) %>% bind_rows()
          
          xy_ind <- sf::st_intersects(full %>% sf::st_as_sf(coords = c("x", "y"), crs = 4326),
                                      ecoregion_shp %>% as("sf")) %>%
            sapply(function(x) length(x) > 0)
          full[xy_ind, ]
        })
        
        if(length(unique(Psi_rot$Year)) > 1) {
          for(ff in 1:max(Psi_rot$Factor)) {
            g <- Psi_rot %>%
              filter(Factor == ff) %>%
              ggplot(aes(Lon, Lat)) + 
              #geom_sf(data = coast, inherit.aes = FALSE) +
              geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
              coord_sf(xlim = range(Psi_rot$Lon), ylim = range(Psi_rot$Lat)) +
              geom_point(shape = 21, aes(fill = value)) +
              facet_wrap(vars(Year)) +
              scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
              labs(x = "Longitude", y = "Latitude", colour = paste("Factor", ff)) +
              theme_bw() +
              theme(panel.spacing = unit(0, "in")) 
            ggsave(file.path(working_dir, paste0(lab[i], "_Factor_", ff, "_ecoregion.png")), g, height = 6, width = 8)
            
            g <- Psi_interp %>%
              filter(Factor == ff) %>%
              ggplot(aes(x, y)) + 
              geom_tile(aes(fill = z, colour = z)) +
              #geom_sf(data = coast, inherit.aes = FALSE) +
              geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) + 
              coord_sf(xlim = range(Psi_rot$Lon), ylim = range(Psi_rot$Lat)) +
              facet_wrap(vars(Year)) +
              scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
              scale_colour_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
              labs(x = "Longitude", y = "Latitude", colour = paste("Factor", ff)) +
              theme_bw() +
              theme(panel.spacing = unit(0, "in")) + 
              guides(fill = "none")
            ggsave(file.path(working_dir, paste0(lab[i], "_Factor_", ff, "_interp_ecoregion.png")), g, height = 6, width = 8)
          }
        } else {
          g <- Psi_rot %>%
            ggplot(aes(Lon, Lat)) + 
            #geom_sf(data = coast, inherit.aes = FALSE) +
            geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) + 
            coord_sf(xlim = range(Psi_rot$Lon), ylim = range(Psi_rot$Lat)) +
            geom_point(shape = 21, aes(fill = value)) +
            facet_wrap(vars(paste("Factor", Factor))) +
            scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
            labs(x = "Longitude", y = "Latitude", colour = "Value") +
            theme_bw() +
            theme(panel.spacing = unit(0, "in")) 
          ggsave(file.path(working_dir, paste0(lab[i], "_Factors_ecoregion.png")), g, height = 6, width = 8)
          
          g <- Psi_interp %>%
            ggplot(aes(x, y)) + 
            geom_tile(aes(fill = z, colour = z)) +
            #geom_sf(data = coast, inherit.aes = FALSE) +
            geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
            coord_sf(xlim = range(Psi_rot$Lon), ylim = range(Psi_rot$Lat)) +
            facet_wrap(vars(paste("Factor", Factor))) +
            scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
            scale_colour_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
            labs(x = "Longitude", y = "Latitude", colour = "Value") +
            theme_bw() +
            theme(panel.spacing = unit(0, "in")) + 
            guides(fill = "none")
          ggsave(file.path(working_dir, paste0(lab[i], "_Factors_interp_ecoregion.png")), g, height = 6, width = 8)
        }
      }
    }
    
    message("Plotted factors.")
  }
  
  # Plot loadings x factors (by species, year)
  if("density" %in% do_plot) {
    
    D_gct <- local({
      x <- array(dim = dim(fit$Report$D_gct))
      x[] <- fit$Report$D_gct
      dimnames(x) <- list(Cell = 1:dim(x)[1], Category = dimnames(fit$Report$D_gct)[[2]], 
                          Year = dimnames(fit$Report$D_gct)[[3]])
      reshape2::melt(x) %>% 
        left_join(fit$spatial_list$latlon_g %>% as.data.frame() %>% mutate(Cell = 1:nrow(.)), by = "Cell") %>%
        left_join(fit$spatial_list$loc_g %>% as.data.frame() %>% mutate(Cell = 1:nrow(.)), by = "Cell")
    })
    
    
    g <- D_gct %>%
      dplyr::filter(Year == max(D_gct$Year)) %>%
      ggplot(aes(Lon, Lat)) + 
      #geom_sf(data = coast, inherit.aes = FALSE) +
      geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) +
      coord_sf(xlim = range(D_gct$Lon), ylim = range(D_gct$Lat)) +
      geom_point(shape = 21, aes(fill = log(value))) +
      facet_wrap(vars(Category)) +
      scale_fill_viridis_c() + 
      #scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
      labs(x = "Longitude", y = "Latitude", fill = "Log-Density") +
      theme_bw() +
      theme(panel.spacing = unit(0, "in")) 
    ggsave(file.path(working_dir, "Density_ecoregion.png"), g, height = 6, width = 8)
    
    D_gct_interp <- local({
      # Interpolation should be in UTM in the model
      coord_out <- expand.grid(xo = seq(min(D_gct$Lon), max(D_gct$Lon), by = 1),
                               yo = seq(min(D_gct$Lat), max(D_gct$Lat), by = 1))
      full <- lapply(unique(D_gct$Category), function(ss) {
        lapply(unique(D_gct$Year), function(yy) {
          mesh_x <- filter(D_gct, Category == ss, Year == yy)
          interp::interp(x = mesh_x$Lon, 
                         y = mesh_x$Lat, 
                         z = mesh_x$value, 
                         xo = coord_out$xo,
                         yo = coord_out$yo,
                         output = "points") %>%
            bind_cols() %>% mutate(Category = ss, Year = yy)
        }) %>% bind_rows()
      }) %>% bind_rows()
      
      xy_ind <- sf::st_intersects(full %>% sf::st_as_sf(coords = c("x", "y"), crs = 4326),
                                  ecoregion_shp %>% as("sf")) %>%
        sapply(function(x) length(x) > 0)
      full[xy_ind, ]
    })
    
    g <- D_gct_interp %>%
      dplyr::filter(Year == max(D_gct$Year), !is.na(z), z > 0) %>%
      ggplot(aes(x, y)) + 
      geom_tile(aes(fill = log(z), colour = log(z))) +
      #geom_sf(data = coast, inherit.aes = FALSE) +
      geom_sf(data = ecoregion_shp %>% as("sf"), colour = "black", fill = NA, inherit.aes = FALSE) + 
      coord_sf(xlim = range(D_gct$Lon), ylim = range(D_gct$Lat)) +
      facet_wrap(vars(Category)) +
      scale_fill_viridis_c() + 
      scale_colour_viridis_c() + 
      #scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
      #scale_colour_gradient2(high = scales::muted("red"), low = scales::muted("blue")) + 
      labs(x = "Longitude", y = "Latitude", colour = "Log Density", fill = "Log Density") +
      theme_bw() +
      theme(panel.spacing = unit(0, "in")) 
    ggsave(file.path(working_dir, "Density_interp_ecoregion.png"), g, height = 6, width = 8)
    
  }
  
  invisible()
}
