library(noisemaker)
library(fungible)
library(tidyverse)
library(pbmcapply)

mc_cores <- 3

# Create a matrix of all of the simulation conditions
condition_matrix <- tidyr::expand_grid(
  factors = c(1, 3, 5),
  items_per_factor = c(5, 15),
  factor_cor = c(0, .3, .6),
  loading_type = c("SS", "MM", "WW", "SW", "SM", "MW"),
  target_rmsea = c(0.025, 0.065, 0.090)
) %>% 
  filter(!(factors == 1 & factor_cor != 0)) %>% 
  mutate(condition_num = 1:n())

condition_matrix <- mutate(
  condition_matrix,
  target_cfi = case_when(target_rmsea == 0.025 ~ 0.99,
                         target_rmsea == 0.065 ~ 0.95,
                         target_rmsea == 0.090 ~ 0.90)
) %>%
  select(condition_num, factors:loading_type, target_rmsea, target_cfi)

# TODO: Add documentation to this section.
reps <- 3
muffled_noisemaker <- possibly(noisemaker, otherwise = NA)

generate_loading_range <- function(loading_type, items) {
  switch(loading_type,
         "SS" = c(.8, .8),
         "MM" = c(.6, .6),
         "WW" = c(.4, .4),
         "SM" = c(.8, .6),
         "SW" = c(.8, .4),
         "MW" = c(.6, .4))
}

# Simulation loop
results_list <- pbmcapply::pbmclapply(
  X = seq_along(condition_matrix$condition_num),
  FUN = function(condition) {
    
    factors <- condition_matrix$factors[condition]
    items_per_factor <- condition_matrix$items_per_factor[condition]
    factor_cor <- condition_matrix$factor_cor[condition]
    loading_type <- condition_matrix$loading_type[condition]
    target_rmsea <- condition_matrix$target_rmsea[condition]
    target_cfi <- condition_matrix$target_cfi[condition]
    
    FacLoadRange <- generate_loading_range(loading_type)
    
    # Generate factor model
    mod <- simFA(Model = list(NFac = factors,
                              NItemPerFac = items_per_factor,
                              Model = "oblique"),
                 Loadings = list(FacLoadDist = "sequential",
                                 FacLoadRange = FacLoadRange),
                 Phi = list(PhiType = "fixed",
                            MaxAbsPhi = factor_cor))
    
    wb_mod <- get_wb_mod(mod, n = 100, values = 15)
    sigma_cb <- muffled_noisemaker(mod, 
                                   method = "CB", 
                                   target_rmsea = target_rmsea)
    
    sigma_list <- purrr::map(
      .x = seq_len(reps), 
      .f = function(i, mod, target_rmsea, target_cfi) {
        sigma_tkl_rmsea <- muffled_noisemaker(mod, method = "TKL", 
                                              target_rmsea = target_rmsea)
        sigma_tkl_rmsea_cfi <- muffled_noisemaker(mod, method = "TKL",
                                                  target_rmsea = target_rmsea,
                                                  target_cfi = target_cfi)
        sigma_wb  <- muffled_noisemaker(mod, method = "WB", 
                                        target_rmsea = target_rmsea, 
                                        wb_mod = wb_mod)
        
        list(sigma_tkl_rmsea = sigma_tkl_rmsea,
             sigma_tkl_rmsea_cfi = sigma_tkl_rmsea_cfi,
             sigma_cb = NA,
             sigma_wb = sigma_wb)
      }, 
      mod = mod, 
      target_rmsea = target_rmsea, 
      target_cfi = target_cfi
    )
    sigma_list[[1]]$sigma_cb <- sigma_cb
    sigma_list
  }, mc.cores = mc_cores
)
