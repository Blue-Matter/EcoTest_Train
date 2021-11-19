
# What affects F on bycatch? 

# Fit autoregressive model

# From target
# F
# Abundance
# 

# From bycatch
# N

# Longline specific
library(MSEtool)


byc <- c("BSH", "SMA", "WHM", "BUM")
targ <- c("BET", "SWO")

multiHist_byc <- lapply(paste0('multiHist_', byc, '.rds'), readRDS)
multiHist_targ <- lapply(paste0('multiHist_', targ, '.rds'), readRDS)

MOM_byc <- lapply(paste0('MOM_', byc, '.rds'), readRDS)
MOM_targ <- lapply(paste0('MOM_', targ, '.rds'), readRDS)

##### SSB
fn_SB <- function(MOM, multiHist, sp, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      d <- multiHist[[s]][[p]][[1]]
      VB <- local({
        N <- d@AtAge$Number[1, , , ] %>% apply(1:2, sum)
        Fec <- d@SampPars$Stock$Fec_Age[1, , 1:ncol(N)]
        colSums(N * Fec)
      })
      
      #VB <- rowSums(d@TSdata$SBiomass[1, , ])
      data.frame(Year = d@OMPars$CurrentYr[1] - length(VB):1 + 1, value = VB,
                 Type = "Spawning biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:")) %>% dplyr::filter(Sex == "Female")
  
}

SB_byc <- fn_SB(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
SB_targ <- fn_SB(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)

rbind(SB_byc, SB_targ) %>% mutate(S = paste(Sp, Species)) %>%
  ggplot(aes(Year, value)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Female spawning biomass") + theme_bw()
ggsave("Figures/MOM/SB.png", width = 8, height = 4)

##### B
rbind(SB_byc, SB_targ) %>% mutate(S = paste(Sp, Species)) %>%
  ggplot(aes(Year, value)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Female spawning biomass") + theme_bw()
ggsave("Figures/MOM/SB.png", width = 8, height = 4)

##### Aggregate vulnerable biomass
fn_VBt <- function(MOM, multiHist, sp, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      d <- multiHist[[s]][[p]][[1]]
      VB <- rowSums(d@TSdata$VBiomass[1, , ])
      data.frame(Year = d@OMPars$CurrentYr[1] - length(VB):1 + 1, value = VB,
                 Type = "Aggregate vulnerable biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
  
}

VBt_byc <- fn_VBt(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
VBt_targ <- fn_VBt(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)

rbind(VBt_byc, VBt_targ) %>% mutate(S = paste(Sp, Species)) %>%
  ggplot(aes(Year, value, colour = Sex)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Aggregate vulnerable biomass") + theme_bw()
ggsave("Figures/MOM/VB_aggregate.png", width = 8, height = 4)


##### Fleet specific 
#fn_VBf <- function(MOM, multiHist, sp, isbycatch = TRUE) {
#  
#  lapply(1:length(sp), function(s) {
#    lapply(1:length(MOM[[s]]@Stocks), function(p) {
#      sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
#        d <- multiHist[[s]][[p]][[f]]
#        N <- d@AtAge$Number[1, , , ] %>% apply(1:2, sum)
#        V <- d@SampPars$Fleet$V[1, , 1:ncol(N)] 
#        if(!is.null(d@SampPars$Fleet$Wt_age_C)) {
#          Wt <- d@SampPars$Fleet$Wt_age_C[1, , 1:ncol(N)] 
#        } else {
#          Wt <- d@SampPars$Stock$Wt_age[1, , 1:ncol(N)] 
#        }
#        colSums(N * V * Wt) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
#      }) %>% structure(dimnames = list(Year = rownames(.), Fleet = names(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
#        dplyr::mutate(Type = "Vulnerable biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#    }) %>% bind_rows()
#  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
#  
#}
#
#VBf_byc <- fn_VBf(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
#VBf_targ <- fn_VBf(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)
#
#rbind(VBf_byc, VBf_targ) %>% mutate(S = paste(Sp, Species)) %>%
#  ggplot(aes(Year, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
#  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
#  labs(y = "Vulnerable biomass") + theme_bw()
#
#F_byc <- lapply(1:length(byc), function(s, MOM, multiHist, sp) {
#  lapply(1:length(MOM[[s]]@Stocks), function(p) {
#    sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
#      d <- multiHist[[s]][[p]][[f]]
#      apply(d@AtAge$F.Mortality[1, , , ], 2, max) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
#    }) %>% structure(dimnames = list(Year = rownames(.), Fleet = names(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
#      dplyr::mutate(Type = "Fishing mortality", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#  }) %>% bind_rows()
#}, MOM = MOM_byc, multiHist = multiHist_byc, sp = byc) %>% bind_rows() %>% mutate(Sp = "Bycatch:")
#
#F_targ <- lapply(1:length(targ), function(s, MOM, multiHist, sp) {
#  lapply(1:length(MOM[[s]]@Stocks), function(p) {
#    sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
#      d <- multiHist[[s]][[p]][[f]]
#      apply(d@AtAge$F.Mortality[1, , , ], 2, max) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
#    }) %>% structure(dimnames = list(Year = rownames(.), Fleet = names(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
#      dplyr::mutate(Type = "Fishing mortality", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#  }) %>% bind_rows()
#}, MOM = MOM_targ, multiHist = multiHist_targ, sp = targ) %>% bind_rows() %>% mutate(Sp = "Target:")
#
#rbind(F_byc, F_targ) %>% mutate(S = paste(Sp, Species)) %>%
#  ggplot(aes(Year, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
#  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
#  labs(y = "Fishing mortality") + theme_bw()
#
#Removals_byc <- lapply(1:length(byc), function(s, MOM, multiHist, sp) {
#  lapply(1:length(MOM[[s]]@Stocks), function(p) {
#    sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
#      d <- multiHist[[s]][[p]][[f]]
#      rowSums(d@TSdata$Removals[1, , ]) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
#    }) %>% structure(dimnames = list(Year = rownames(.), Fleet = names(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
#      dplyr::mutate(Type = "Fishery removals", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#  }) %>% bind_rows()
#}, MOM = MOM_byc, multiHist = multiHist_byc, sp = byc) %>% bind_rows() %>% mutate(Sp = "Bycatch:")
#
#Removals_targ <- lapply(1:length(targ), function(s, MOM, multiHist, sp) {
#  lapply(1:length(MOM[[s]]@Stocks), function(p) {
#    sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
#      d <- multiHist[[s]][[p]][[f]]
#      rowSums(d@TSdata$Removals[1, , ]) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
#    }) %>% structure(dimnames = list(Year = rownames(.), Fleet = names(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
#      dplyr::mutate(Type = "Fishery removals", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#  }) %>% bind_rows()
#}, MOM = MOM_targ, multiHist = multiHist_targ, sp = targ) %>% bind_rows() %>% mutate(Sp = "Target:")
#
#rbind(Removals_byc, Removals_targ) %>% mutate(S = paste(Sp, Species)) %>%
#  ggplot(aes(Year, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
#  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
#  labs(y = "Fishery removals") + theme_bw()

#### Separate to LL vs. other
LL_byc <- list(1:9, c(1:12)[-c(8, 10, 12)], 2, 2)
LL_targ <- list(10:18, 1:11)

fn_VLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE, biomass = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      F_LL <- local({
        Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
          d <- multiHist[[s]][[p]][[f]]
          apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
            structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
        }, simplify = "array")[, , LL[[s]], drop = FALSE]
        Fsum <- Fout %>% apply(1:2, sum)
        Fapical <- apply(Fsum, 2, max)
        V <- t(Fsum)/Fapical
        N <- multiHist[[s]][[p]][[1]]@AtAge$Number[1, , , ] %>% apply(1:2, sum)
        if(biomass) {
          Wt <- multiHist[[s]][[p]][[1]]@SampPars$Stock$Wt_age[1, , 1:ncol(N)]
        } else {
          Wt <- 1
        }
        data.frame(Year = dimnames(Fout)[[2]], Fleet = "Longline", value = colSums(N * t(V) * Wt))
      })
      
      if(length(LL[[s]]) < length(MOM[[s]]@Fleets[[p]])) {
        F_OTH <- local({
          Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
            d <- multiHist[[s]][[p]][[f]]
            apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
              structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
          }, simplify = "array")[, , -LL[[s]], drop = FALSE]
          Fsum <- Fout %>% apply(1:2, sum)
          Fapical <- apply(Fsum, 2, max)
          V <- t(Fsum)/Fapical
          N <- multiHist[[s]][[p]][[1]]@AtAge$Number[1, , , ] %>% apply(1:2, sum)
          if(biomass) {
            Wt <- multiHist[[s]][[p]][[1]]@SampPars$Stock$Wt_age[1, , 1:ncol(N)]
          } else {
            Wt <- 1
          }
          data.frame(Year = dimnames(Fout)[[2]], Fleet = "Other", value = colSums(N * t(V) * Wt))
        })
        
        F_LL <- rbind(F_LL, F_OTH)
      }
      F_LL %>% dplyr::mutate(Year = as.numeric(Year), Type = "Vulnerable biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
  
}

VBLL_byc <- fn_VLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc)
VBLL_targ <- fn_VLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE)
VBLL <- rbind(VBLL_byc, VBLL_targ) %>% mutate(S = paste(Sp, Species))

ggplot(VBLL, aes(Year, value, linetype = Fleet, colour = Sex)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Vulnerable biomass") + theme_bw()
ggsave("Figures/MOM/VB_LL.png", width = 8, height = 4)

VNLL_byc <- fn_VLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc, biomass = FALSE)
VNLL_targ <- fn_VLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE, biomass = FALSE)
VNLL <- rbind(VNLL_byc, VNLL_targ) %>% mutate(S = paste(Sp, Species))

ggplot(VNLL, aes(Year, value, linetype = Fleet, colour = Sex)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Vulnerable abundance") + theme_bw()
ggsave("Figures/MOM/VN_LL.png", width = 8, height = 4)


### 
fn_FLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      F_LL <- local({
        Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
          d <- multiHist[[s]][[p]][[f]]
          apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
            structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
        }, simplify = "array")[, , LL[[s]], drop = FALSE]
        Fsum <- Fout %>% apply(1:2, sum) %>% apply(2, max)
        data.frame(Year = dimnames(Fout)[[2]], Fleet = "Longline", value = Fsum)
      })
      
      if(length(LL[[s]]) < length(MOM[[s]]@Fleets[[p]])) {
        F_OTH <- local({
          Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
            d <- multiHist[[s]][[p]][[f]]
            apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
              structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
          }, simplify = "array")[, , -LL[[s]], drop = FALSE]
          Fsum <- Fout %>% apply(1:2, sum) %>% apply(2, max)
          data.frame(Year = dimnames(Fout)[[2]], Fleet = "Other", value = Fsum)
        })
        
        F_LL <- rbind(F_LL, F_OTH)
      }
      F_LL %>% dplyr::mutate(Year = as.numeric(Year), Type = "Fishing mortality", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
  
}

FLL_byc <- fn_FLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc)
FLL_targ <- fn_FLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, isbycatch = FALSE)
FLL <- rbind(FLL_byc, FLL_targ) %>% mutate(S = paste(Sp, Species)) %>% group_by(Year, Fleet, S, Type) %>% summarise(value = max(value))

ggplot(FLL, aes(Year, value, linetype = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Apical fishing mortality") + theme_bw()
ggsave("Figures/MOM/F_longline.png", width = 8, height = 4)


fn_CLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
        d <- multiHist[[s]][[p]][[f]]
        rowSums(d@TSdata$Removals[1, , ]) %>% structure(names = d@OMPars$CurrentYr[1] - length(.):1 + 1)
      }) %>% structure(dimnames = list(Year = rownames(.), Fleet = 1:length(MOM[[s]]@Fleets[[p]]))) %>% reshape2::melt() %>%
        dplyr::mutate(Type = "Fishery removals", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s], 
                      Fleet = ifelse(Fleet %in% LL[[s]], "Longline", "Other"))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    group_by(Year, Type, Species, Fleet) %>% summarise(value = sum(value)) %>% 
    mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
}

CLL_byc <- fn_CLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc)
CLL_targ <- fn_CLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE)
CLL <- rbind(CLL_byc, CLL_targ) %>% mutate(S = paste(Sp, Species))
ggplot(CLL, aes(Year, value, linetype = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
  labs(y = "Fishery removals") + theme_bw()
ggsave("Figures/MOM/Catch_longline.png", width = 8, height = 4)




## Correlations

# 1. F bycatch related to F target
FLL_df <- FLL %>% dplyr::filter(Fleet == "Longline") %>% reshape2::dcast(list("Year", "S"))
panel.cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  text(0.5, 0.5, txt, cex = 1)
}

png("Figures/MOM/F_longline_pairs.png", width = 6, height = 6, res = 400, units = "in")
FLL_df %>% dplyr::select(contains(":")) %>% 
  pairs(gap = 0, upper.panel = panel.cor, main = "Scatterplot and correlations in longline fishing mortality")
dev.off()

local({
  F1 <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(S1 = S, F1 = value)
  F2 <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(S2 = S, F2 = value)
  FF <- left_join(F1, F2, by = c("Year"))
  
  c <- FF %>% group_by(S1, S2) %>% summarise(c = cor(F1, F2, use = "complete.obs") %>% round(2))
  
  ggplot(FF, aes(F1, F2, colour = Year)) + facet_grid(S2 ~ S1, scales = "free", switch = "both") + 
    geom_path() + geom_point() + theme_bw() + 
    scale_colour_viridis_c() + labs(x = "F", y = "F") +
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    geom_text(data = c, inherit.aes = FALSE, vjust = 1.1, hjust = 0, x = 0, y = Inf, aes(label = c)) + 
    theme(strip.background = element_blank(), strip.placement	= "outside", strip.text = element_text(size = 10))
  ggsave("Figures/MOM/F_longline_pairs3.png", width = 10, height = 8)
})

## vector autoregression model
local({ # 1971 - 2013
  dat <- FLL_df %>% dplyr::filter(!is.na(rowSums(.)))
  out <- VAR(dat %>% dplyr::select(contains(":")), type = "both")
})


# 2. F bycatch related to N bycatch
local({
  FF <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(F = value)
  VB <- VBLL %>% dplyr::filter(Fleet == "Longline") %>% group_by(Year, Fleet, S) %>% summarise(VB = sum(value)) %>% 
    dplyr::filter(!is.na(VB))
  
  dat <- left_join(VB, FF, by = c("Year", "Fleet", "S"))
  
  ggplot(dat, aes(VB, F, colour = Year)) + geom_point() + geom_path() + facet_wrap(~S, scales = "free") + theme_bw() +
    scale_colour_viridis_c() + 
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    labs(x = "Longline vulnerable biomass", y = "Longline fishing mortality")
  ggsave("Figures/MOM/F_vs_VB_longline.png", width = 8, height = 4)
  
  ggplot(dat %>% dplyr::filter(grepl("Bycatch", S)), aes(VB, F, colour = Year)) + geom_point() + geom_path() + facet_wrap(~S, scales = "free") + theme_bw() +
    scale_colour_viridis_c() + 
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    labs(x = "Longline vulnerable biomass", y = "Longline fishing mortality")
  ggsave("Figures/MOM/F_vs_VB_longline2.png", width = 8, height = 4)
})


local({
  FF <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(F = value)
  VB <- VNLL %>% dplyr::filter(Fleet == "Longline") %>% group_by(Year, Fleet, S) %>% summarise(VB = sum(value)) %>% 
    dplyr::filter(!is.na(VB))
  
  dat <- left_join(VB, FF, by = c("Year", "Fleet", "S"))
  
  ggplot(dat, aes(VB, F, colour = Year)) + geom_point() + geom_path() + facet_wrap(~S, scales = "free") + theme_bw() +
    scale_colour_viridis_c() + 
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    labs(x = "Longline vulnerable abundance", y = "Longline fishing mortality")
  ggsave("Figures/MOM/F_vs_VN_longline.png", width = 8, height = 4)
  
  ggplot(dat %>% dplyr::filter(grepl("Bycatch", S)), aes(VB, F, colour = Year)) + geom_point() + geom_path() + facet_wrap(~S, scales = "free") + theme_bw() +
    scale_colour_viridis_c() + 
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    labs(x = "Longline vulnerable abundance", y = "Longline fishing mortality")
  ggsave("Figures/MOM/F_vs_VN_longline2.png", width = 8, height = 4)
})



# 2. F bycatch related to Removals
local({
  FF <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(F = value)
  RR <- CLL %>% dplyr::filter(Fleet == "Longline") %>% group_by(Year, Fleet, S) %>% summarise(Removals = sum(value)) %>% 
    dplyr::filter(!is.na(Removals))
  
  dat <- left_join(transmute(RR, S1 = paste("VBiomass", strsplit(S, " ") %>% sapply(function(x) x[2])), Removals = Removals), 
                   mutate(FF, S = paste("F ", strsplit(S, " ") %>% sapply(function(x) x[2]))), by = c("Year", "Fleet"))
  
  c <- dat %>% group_by(S, S1) %>% summarise(c = cor(Removals, F, use = "complete.obs") %>% round(2))
  
  ggplot(dat, aes(Removals, F, colour = Year)) + facet_grid(S ~ S1, scales = "free", switch = "both") + 
    geom_path() + geom_point() + theme_bw() + 
    scale_colour_viridis_c() + labs(x = NULL, y = NULL) +
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    geom_text(data = c, inherit.aes = FALSE, vjust = 1.1, hjust = 0, x = 0, y = Inf, aes(label = c)) + 
    theme(strip.background = element_blank(), strip.placement	= "outside", strip.text = element_text(size = 10))
  
  ggsave("Figures/MOM/F_vs_Catch_longline_all.png", width = 10, height = 8)
})

