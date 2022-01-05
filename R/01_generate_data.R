if (!require(noisemaker)) {
  devtools::install_github("JustinKracht/noisemaker",
                           dependencies = TRUE)
}

library(noisemaker)
library(here)
library(fungible)
library(tidyverse)
library(pbmcapply)

# Set variable values -----------------------------------------------------

DATA_DIR <- here("data")

# Number of simulation reps
reps <- 500

# Create a matrix of all of the simulation conditions
condition_matrix <- tidyr::expand_grid(
  factors = c(1, 3, 5, 10),
  items_per_factor = c(5, 15),
  factor_cor = c(0, .3, .6),
  loading = c("weak", "moderate", "strong"),
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

saveRDS(condition_matrix, here("data", "condition_matrix.RDS"))

# Make matrix of optimization method choices for TKL
method_matrix <- tidyr::expand_grid(
  rmsea_weight = c(0,1),
  cfi_weight = c(0,1)
) %>% filter(rmsea_weight != 0 | cfi_weight != 0)

# tryCatch function to wrap noisemaker so that warnings and errors are captured
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# Simulation loop ---------------------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(666)
seed_list <- sample(1e7, size = nrow(condition_matrix), replace = FALSE)

results_list <- pbmclapply(
  X = seq_along(condition_matrix$condition_num),
  # X = 154:nrow(condition_matrix),
  FUN = function(condition) {
    
    set.seed(seed_list[condition])
    
    factors <- condition_matrix$factors[condition]
    items_per_factor <- condition_matrix$items_per_factor[condition]
    factor_cor <- condition_matrix$factor_cor[condition]
    loading <- condition_matrix$loading[condition]
    target_rmsea <- condition_matrix$target_rmsea[condition]
    target_cfi <- condition_matrix$target_cfi[condition]
    
    FacLoadRange <- switch(loading,
                           "weak" = 0.4,
                           "moderate" = 0.6,
                           "strong" = 0.8)
    
    # Generate factor model
    mod <- simFA(Model = list(NFac = factors,
                              NItemPerFac = items_per_factor,
                              Model = "oblique"),
                 Loadings = list(FacLoadDist = "fixed",
                                 FacLoadRange = FacLoadRange),
                 Phi = list(PhiType = "fixed",
                            MaxAbsPhi = factor_cor))
    
    wb_mod <- get_wb_mod(mod, n = 100, values = 15)
    
    sigma_list <- purrr::map(
      .x = seq_len(reps), 
      .f = function(i, mod, target_rmsea) {
        
        # TKL method variations
        sigma_tkl_rmsea <- myTryCatch(
          noisemaker(
            mod,
            method = "TKL",
            target_rmsea = target_rmsea,
            target_cfi = NULL,
            tkl_ctrl = list(WmaxLoading = .3,
                            NWmaxLoading = 2,
                            max_tries = 100,
                            maxit = 1000)
          )
        )
        
        sigma_tkl_cfi <- myTryCatch(
          noisemaker(
            mod,
            method = "TKL",
            target_rmsea = NULL,
            target_cfi = target_cfi,
            tkl_ctrl = list(WmaxLoading = .3,
                            NWmaxLoading = 2,
                            max_tries = 100,
                            maxit = 1000)
          )
        )
        
        sigma_tkl_rmsea_cfi <- myTryCatch(
          noisemaker(
            mod,
            method = "TKL",
            target_rmsea = target_rmsea,
            target_cfi = target_cfi,
            tkl_ctrl = list(WmaxLoading = .3,
                            NWmaxLoading = 2,
                            max_tries = 100,
                            maxit = 1000)
          )
        )
        
        # Other methods
        # Skip CB for the largest conditions; takes way too long.
        if (condition %in% 154:nrow(condition_matrix)) {
          sigma_cb <- list("value" = NA,
                           "warning" = NULL,
                           "error" = NULL)
        } else {
          sigma_cb  <- myTryCatch(noisemaker(mod, method = "CB",
                                             target_rmsea = target_rmsea))
        }
        sigma_wb  <- myTryCatch(noisemaker(mod, method = "WB", 
                                           target_rmsea = target_rmsea, 
                                           wb_mod = wb_mod))
        
        list(sigma_tkl_rmsea = sigma_tkl_rmsea,
             sigma_tkl_cfi = sigma_tkl_cfi,
             sigma_tkl_rmsea_cfi = sigma_tkl_rmsea_cfi,
             sigma_cb  = sigma_cb,
             sigma_wb  = sigma_wb)
      }, 
      mod = mod,
      target_rmsea = target_rmsea
    )
    
    # Save condition results
    saveRDS(
      sigma_list, 
      file = here(
        "data",
        paste0("results_", 
               formatC(condition, width = 3, flag = 0),
               ".RDS"
        )
      )
    )
    
  }, 
  mc.preschedule = FALSE, 
  mc.cores = 18
)
