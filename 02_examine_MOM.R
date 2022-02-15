
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
library(tidyverse)

byc <- c("BSH", "SMA", "WHM", "BUM")
targ <- c("BET", "SWO")

multiHist_byc <- lapply(paste0('MOM/multiHist_', byc, '.rds'), readRDS)
multiHist_targ <- lapply(paste0('MOM/multiHist_', targ, '.rds'), readRDS)

MOM_byc <- lapply(paste0('MOM/MOM_', byc, '.rds'), readRDS)
MOM_targ <- lapply(paste0('MOM/MOM_', targ, '.rds'), readRDS)

##### Spawning biomass
fn_SB <- function(MOM, multiHist, sp, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    lapply(1:length(MOM[[s]]@Stocks), function(p) {
      d <- multiHist[[s]][[p]][[1]]
      VB <- local({
        N <- d@AtAge$Number[1, , , ] %>% apply(1:2, sum)
        Fec <- d@SampPars$Stock$Fec_Age[1, , 1:ncol(N)]
        colSums(N * Fec)
      })
      
      data.frame(Year = d@OMPars$CurrentYr[1] - length(VB):1 + 1, value = VB,
                 Type = "Spawning biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:")) %>% dplyr::filter(Sex == "Female")
  
}

SB_byc <- fn_SB(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
SB_targ <- fn_SB(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)

rbind(SB_byc, SB_targ) %>% mutate(S = paste(Sp, Species)) %>%
  ggplot(aes(Year, value)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) +
  labs(y = "Female spawning biomass") + theme_bw()
ggsave("Figures/MOM/SB.png", width = 8, height = 4)

##### B
fn_B <- function(MOM, multiHist, sp, isbycatch = TRUE) {
  lapply(1:length(sp), function(s) {
    out <- lapply(1:length(MOM[[s]]@Stocks), function(p) {
      d <- multiHist[[s]][[p]][[1]]
      VB <- local({
        N <- d@AtAge$Number[1, , , ] %>% apply(1:2, sum)
        W <- d@SampPars$Stock$Wt_age[1, , 1:ncol(N)]
        colSums(N * W)
      })
      
      #VB <- rowSums(d@TSdata$SBiomass[1, , ])
      data.frame(Year = d@OMPars$CurrentYr[1] - length(VB):1 + 1, value = VB,
                 Type = "Total biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
    if(length(unique(out$Sex)) == 1) out$Sex <- "Unisex"
    return(out)
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
}

B_byc <- fn_B(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
B_targ <- fn_B(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)

rbind(B_byc, B_targ) %>% mutate(S = paste(Sp, Species)) %>%
  ggplot(aes(Year, value, linetype = Sex)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) + 
  scale_linetype_manual(values = c("Unisex" = 2, "Female" = 1, "Male" = 4)) +
  labs(y = "Total biomass") + theme_bw()
ggsave("Figures/MOM/B.png", width = 8, height = 4)

##### Aggregate vulnerable biomass
#fn_VBt <- function(MOM, multiHist, sp, isbycatch = TRUE) {
#  
#  lapply(1:length(sp), function(s) {
#    out <- lapply(1:length(MOM[[s]]@Stocks), function(p) {
#      d <- multiHist[[s]][[p]][[1]]
#      VB <- rowSums(d@TSdata$VBiomass[1, , ])
#      data.frame(Year = d@OMPars$CurrentYr[1] - length(VB):1 + 1, value = VB,
#                 Type = "Aggregate vulnerable biomass", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
#    }) %>% bind_rows()
#    if(length(unique(out$Sex)) == 1) out$Sex <- "Unisex"
#    return(out)
#  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
#  
#}
#
#VBt_byc <- fn_VBt(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc)
#VBt_targ <- fn_VBt(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, FALSE)
#
#rbind(VBt_byc, VBt_targ) %>% mutate(S = paste(Sp, Species)) %>%
#  ggplot(aes(Year, value, colour = Sex)) + facet_wrap(~S, scales = "free_y") + 
#  geom_line() + geom_hline(yintercept = 0, colour = NA) + 
#  labs(y = "Aggregate vulnerable biomass") + theme_bw()
#ggsave("Figures/MOM/VB_aggregate.png", width = 8, height = 4)


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

#### Calculate aggregate vulnerable biomass of longline vs. other fleet

## Index corresponding to longline fleets for each species
LL_byc <- list(1:9, c(1:12)[-c(8, 10, 12)], 2, 2) %>% structure(names = byc)
LL_targ <- list(10:18, 1:11) %>% structure(names = targ)

fn_VLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE, biomass = TRUE) {
  
  lapply(1:length(sp), function(s) {
    out <- lapply(1:length(MOM[[s]]@Stocks), function(p) {
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
    if(length(unique(out$Sex)) == 1) out$Sex <- "Unisex"
    return(out)
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
  
}

VBLL_byc <- fn_VLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc)
VBLL_targ <- fn_VLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE)
VBLL <- rbind(VBLL_byc, VBLL_targ) %>% mutate(S = paste(Sp, Species))

ggplot(VBLL, aes(Year, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) + 
  scale_linetype_manual(values = c("Unisex" = 2, "Female" = 1, "Male" = 4)) +
  labs(y = "Vulnerable biomass") + theme_bw()
ggsave("Figures/MOM/VB_LL.png", width = 8, height = 4)

VNLL_byc <- fn_VLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc, biomass = FALSE)
VNLL_targ <- fn_VLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE, biomass = FALSE)
VNLL <- rbind(VNLL_byc, VNLL_targ) %>% mutate(S = paste(Sp, Species))

ggplot(VNLL, aes(Year, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) + 
  scale_linetype_manual(values = c("Unisex" = 2, "Female" = 1, "Male" = 4)) +
  labs(y = "Vulnerable abundance") + theme_bw()
ggsave("Figures/MOM/VN_LL.png", width = 8, height = 4)

### Aggregate selectivity of longline vs. other fleet
fn_selLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE, biomass = TRUE) {
  
  lapply(1:length(sp), function(s) {
    out <- lapply(1:length(MOM[[s]]@Stocks), function(p) {
      sel_LL <- local({
        Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
          d <- multiHist[[s]][[p]][[f]]
          apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
            structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
        }, simplify = "array")[, , LL[[s]], drop = FALSE]
        Fsum <- Fout %>% apply(1:2, sum)
        Fapical <- apply(Fsum, 2, max)
        V <- t(Fsum)/Fapical
        data.frame(Fleet = "Longline", value = V[nrow(V), ], Age = 1:ncol(V) - 1)
      })
      
      if(length(LL[[s]]) < length(MOM[[s]]@Fleets[[p]])) {
        sel_OTH <- local({
          Fout <- sapply(1:length(MOM[[s]]@Fleets[[p]]), function(f) {
            d <- multiHist[[s]][[p]][[f]]
            apply(d@AtAge$F.Mortality[1, , , ], 1:2, max) %>% 
              structure(dimnames = list(Age = 1:nrow(.) - 1, Year = d@OMPars$CurrentYr[1] - ncol(.):1 + 1))
          }, simplify = "array")[, , -LL[[s]], drop = FALSE]
          Fsum <- Fout %>% apply(1:2, sum)
          Fapical <- apply(Fsum, 2, max)
          V <- t(Fsum)/Fapical
          
          data.frame(Fleet = "Other", value = V[nrow(V), ], Age = 1:ncol(V) - 1)
        })
        
        sel_LL <- rbind(sel_LL, sel_OTH)
      }
      sel_LL %>% dplyr::mutate(Type = "Selectivity", Sex = names(MOM[[s]]@Stocks)[p], Species = sp[s])
    }) %>% bind_rows()
    if(length(unique(out$Sex)) == 1) out$Sex <- "Unisex"
    return(out)
  }) %>% bind_rows() %>% mutate(Sp = ifelse(isbycatch, "Bycatch:", "Target:"))
  
}


selLL_byc <- fn_selLL(MOM = MOM_byc, multiHist = multiHist_byc, sp = byc, LL = LL_byc)
selLL_targ <- fn_selLL(MOM = MOM_targ, multiHist = multiHist_targ, sp = targ, LL = LL_targ, FALSE)
selLL <- rbind(selLL_byc, selLL_targ) %>% mutate(S = paste(Sp, Species))

ggplot(selLL, aes(Age, value, linetype = Sex, colour = Fleet)) + facet_wrap(~S, scales = "free_x") + 
  geom_line() + expand_limits(y = 0) + 
  scale_linetype_manual(values = c("Unisex" = 2, "Female" = 1, "Male" = 4)) +
  labs(y = "Selectivity") + theme_bw()
ggsave("Figures/MOM/sel_LL.png", width = 8, height = 4)



### Apical F of (aggregate) longline vs. other fleet
fn_FLL <- function(MOM, multiHist, sp, LL, isbycatch = TRUE) {
  
  lapply(1:length(sp), function(s) {
    out <- lapply(1:length(MOM[[s]]@Stocks), function(p) {
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
FLL <- rbind(FLL_byc, FLL_targ) %>% mutate(S = paste(Sp, Species)) %>% group_by(Year, Fleet, S, Species, Type) %>% summarise(value = max(value))

ggplot(FLL, aes(Year, value, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) + 
  labs(y = "Apical fishing mortality") + theme_bw()
ggsave("Figures/MOM/F_longline.png", width = 8, height = 4)

### Fishery removals of longline vs. other fleet
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
ggplot(CLL, aes(Year, value, colour = Fleet)) + facet_wrap(~S, scales = "free_y") + 
  geom_line() + expand_limits(y = 0) + 
  labs(y = "Fishery removals") + theme_bw()
ggsave("Figures/MOM/Catch_longline.png", width = 8, height = 4)




###### Predicting F bycatch - looking through correlations and doing some regressions

# 1. F bycatch as a function of F target
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
  
  cc <- FF %>% group_by(S1, S2) %>% summarise(c = cor(F1, F2, use = "complete.obs") %>% round(2))
  
  ggplot(FF, aes(F1, F2, colour = Year)) + facet_grid(S2 ~ S1, scales = "free", switch = "both") + 
    geom_path() + geom_point() + theme_bw() + 
    scale_colour_viridis_c() + labs(x = "F", y = "F") +
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    geom_text(data = cc, inherit.aes = FALSE, vjust = 1.1, hjust = 0, x = 0, y = Inf, aes(label = c)) + 
    theme(strip.background = element_blank(), strip.placement	= "outside", strip.text = element_text(size = 10))
  ggsave("Figures/MOM/F_longline_pairs3.png", width = 10, height = 8)
  
  dplyr::filter(FF, grepl("Target", S1), grepl("Bycatch", S2)) %>%
    ggplot(aes(F1, F2, colour = Year)) + facet_grid(S2 ~ S1, scales = "free", switch = "both") + 
    geom_path() + geom_point() + theme_bw() + 
    scale_colour_viridis_c() + labs(x = "F", y = "F") +
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    geom_text(data = dplyr::filter(cc, grepl("Target", S1), grepl("Bycatch", S2)), 
              inherit.aes = FALSE, vjust = 1.1, hjust = 0, x = 0, y = Inf, aes(label = c)) + 
    theme(strip.background = element_blank(), strip.placement	= "outside", strip.text = element_text(size = 10))
  ggsave("Figures/MOM/F_longline_pairs4.png", width = 5, height = 6)
  
})

# Multivariate regression to predict F bycatch from F_target
#1971-2013
FLL_df <- FLL %>% dplyr::filter(Fleet == "Longline") %>% reshape2::dcast(list("Year", "Species"))
FLL_df <- FLL_df[22:64, ]
mod <- lm(cbind(BSH, BUM, SMA, WHM) ~ 1, FLL_df)
mod2 <- lm(cbind(BSH, BUM, SMA, WHM) ~ BET + SWO + 0, FLL_df)
mod3 <- lm(cbind(BSH, BUM, SMA, WHM) ~ BET + SWO + 1, FLL_df)

# Generate predictive surface
F_predict <- expand.grid(SWO = seq(0, 0.5, 0.01), BET = seq(0, 0.5, 0.01))
F_out <- predict(mod2, newdata = F_predict)

cbind(F_predict, F_out) %>% reshape2::melt(id.vars = c("SWO", "BET")) %>%
  ggplot(aes(SWO, BET)) + 
  #geom_contour_filled(aes(z = value), binwidth = 0.05) + 
  geom_contour(aes(z = value), binwidth = 0.1) + 
  geom_label_contour(aes(z = value), binwidth = 0.1) +
  facet_wrap(~paste(variable, "Longline F")) + theme_bw() +
  labs(x = "SWO Longline F", y = "BET Longline F")
ggsave("Figures/MOM/F_lm.png", width = 5, height = 5)
               
saveRDS(mod2, "Frel/F_lm.rds")

## vector autoregression model
#local({ # 1971 - 2013
#  dat <- FLL_df %>% dplyr::filter(!is.na(rowSums(.)))
#  out <- VAR(dat %>% dplyr::select(contains(":")), type = "both")
#})


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



# 2. F bycatch is predicted by Removals
local({
  FF <- FLL %>% dplyr::filter(Fleet == "Longline") %>% transmute(F = value)
  RR <- CLL %>% dplyr::filter(Fleet == "Longline") %>% group_by(Year, Fleet, S) %>% summarise(Removals = sum(value)) %>% 
    dplyr::filter(!is.na(Removals))
  
  dat <- left_join(transmute(RR, S1 = paste("Removals", strsplit(S, " ") %>% sapply(function(x) x[2])), Removals = Removals), 
                   mutate(FF, S = paste("F ", strsplit(S, " ") %>% sapply(function(x) x[2]))), by = c("Year", "Fleet"))
  
  cc <- dat %>% group_by(S, S1) %>% summarise(c = cor(Removals, F, use = "complete.obs") %>% round(2))
  
  ggplot(dat, aes(Removals, F, colour = Year)) + facet_grid(S ~ S1, scales = "free", switch = "both") + 
    geom_path() + geom_point() + theme_bw() + 
    scale_colour_viridis_c() + labs(x = NULL, y = NULL) +
    geom_hline(yintercept = 0, colour = NA) + geom_vline(xintercept = 0, colour = NA) + 
    geom_text(data = cc, inherit.aes = FALSE, vjust = 1.1, hjust = 0, x = 0, y = Inf, aes(label = c)) + 
    theme(strip.background = element_blank(), strip.placement	= "outside", strip.text = element_text(size = 10))
  ggsave("Figures/MOM/F_vs_Catch_longline_all.png", width = 10, height = 8)
})

head(FLL)

FVB <- local({
  a <- FLL %>% mutate(F = value) %>% select(c("Year", "Fleet", "Species", "F", "S"))
  b <- VBLL %>% group_by(Year, Fleet, Species) %>% summarise(VB = sum(value, na.rm = TRUE)) #CLL %>% mutate(Removals = value) %>% select(c("Year", "Fleet", "Species", "Removals", "S"))
  c <- dplyr::left_join(a, b, by = c("Year", "Fleet", "Species")) %>% dplyr::filter(Year >= 1971 & Year <= 2013, Fleet == "Longline")
  reshape2::melt(c[, c(1, 3, 4, 6)], id.vars = c("Year", "Species")) %>%
    reshape2::dcast(Year ~ Species + variable)
})

#mod <-  lm(cbind(BSH_F, BUM_F, SMA_F, WHM_F) ~ 1, FLL_df)
mod2 <- lm(cbind(BSH_F, BUM_F, SMA_F, WHM_F) ~ BET_F + SWO_F + BET_VB + SWO_VB + 0, FVB)
mod3 <- lm(cbind(BSH_F, BUM_F, SMA_F, WHM_F) ~ BET_F + SWO_F + BET_VB + SWO_VB + 1, FVB)

