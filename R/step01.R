if (!require(noisemaker)) {
  devtools::install_github("JustinKracht/noisemaker",
                           dependencies = TRUE)
}

library(noisemaker)
library(fungible)
library(tidyverse)
library(pbmcapply)
library(tictoc)

# Set variable values -----------------------------------------------------

# Number of simulation reps
reps <- 5

# Create a matrix of all of the simulation conditions
condition_matrix <- tidyr::expand_grid(
  factors = c(1, 3, 5, 10),
  items_per_factor = c(5, 15),
  factor_cor = c(0, .3, .6),
  loading = c("weak", "moderate", "strong", "random"),
  target_rmsea = c(0.025, 0.065, 0.090)
) %>% filter(
  !(factors == 1 & factor_cor != 0)
) %>% mutate(
  condition_num = 1:n()
) %>% dplyr::relocate(condition_num)

# Add target CFI values
condition_matrix <- condition_matrix %>% mutate(
  target_cfi = case_when(target_rmsea == 0.025 ~ 0.99,
                         target_rmsea == 0.065 ~ 0.95,
                         target_rmsea == 0.090 ~ 0.90)
)

# Make matrix of optimization method choices for TKL
method_matrix <- tidyr::expand_grid(
  rmsea_weight = c(0,1),
  cfi_weight = c(0,1)
) %>% filter(rmsea_weight != 0 | cfi_weight != 0)

# Make a safe noisemaker function
quiet_noisemaker <- purrr::possibly(noisemaker, otherwise = NA, quiet = TRUE)

# Simulation loop ---------------------------------------------------------

set.seed(666)
seed_list <- sample(1e7, size = nrow(condition_matrix), replace = FALSE)

tic() # Start timing

results_list <- pbmcapply::pbmclapply(
  X = seq_along(condition_matrix$condition_num),
  FUN = function(condition) {
    
    set.seed(seed_list[condition])
    
    factors <- condition_matrix$factors[condition]
    items_per_factor <- condition_matrix$items_per_factor[condition]
    factor_cor <- condition_matrix$factor_cor[condition]
    loading <- condition_matrix$loading[condition]
    target_rmsea <- condition_matrix$target_rmsea[condition]
    target_cfi <- condition_matrix$target_cfi[condition]
    
    # Set factor loading range given the loading condition
    if (loading == "random") {
      FacLoadDist <- "sequential"
    } else {
      FacLoadDist <- "fixed"
    }
    
    FacLoadRange <- switch(loading,
                           "weak" = 0.4,
                           "moderate" = 0.6,
                           "strong" = 0.8,
                           "random" = c(0.4, 0.8))
    
    # Generate factor model
    mod <- simFA(Model = list(NFac = factors,
                              NItemPerFac = items_per_factor,
                              Model = "oblique"),
                 Loadings = list(FacLoadDist = FacLoadDist,
                                 FacLoadRange = FacLoadRange),
                 Phi = list(PhiType = "fixed",
                            MaxAbsPhi = factor_cor))
    
    wb_mod <- get_wb_mod(mod, n = 100, values = 15)
    
    sigma_list <- purrr::map(
      .x = seq_len(reps), 
      .f = function(i, mod, target_rmsea) {
        
        # TKL method variations
        sigma_tkl_rmsea_constraints <- quiet_noisemaker(
          mod,
          method = "TKL",
          target_rmsea = target_rmsea,
          target_cfi = NULL,
          tkl_ctrl = list(WmaxLoading = .3,
                          NWmaxLoading = 2,
                          max_tries = 100,
                          maxit = 5000)
        )
        
        sigma_tkl_cfi_constraints <- quiet_noisemaker(
          mod,
          method = "TKL",
          target_rmsea = target_rmsea,
          target_cfi = NULL,
          tkl_ctrl = list(WmaxLoading = .3,
                          NWmaxLoading = 2,
                          max_tries = 100,
                          maxit = 5000)
        )
        
        sigma_tkl_rmsea_cfi_constraints <- quiet_noisemaker(
          mod,
          method = "TKL",
          target_rmsea = target_rmsea,
          target_cfi = target_cfi,
          tkl_ctrl = list(WmaxLoading = .3,
                          NWmaxLoading = 2,
                          max_tries = 100,
                          maxit = 5000)
        )
        
        # Other methods
        sigma_cb  <- quiet_noisemaker(mod, method = "CB", 
                                target_rmsea = target_rmsea)
        sigma_wb  <- quiet_noisemaker(mod, method = "WB", 
                                target_rmsea = target_rmsea, wb_mod = wb_mod)
        
        list(sigma_tkl_rmsea_constraints = sigma_tkl_rmsea_constraints,
             sigma_tkl_cfi_constraints = sigma_tkl_cfi_constraints,
             sigma_tkl_rmsea_cfi_constraints = sigma_tkl_rmsea_cfi_constraints,
             sigma_cb  = sigma_cb,
             sigma_wb  = sigma_wb)
      }, 
      mod = mod,
      target_rmsea = target_rmsea
    )
    
    sigma_list
  }, mc.cores = 1
)

toc() # End timing; how much time elapsed?
