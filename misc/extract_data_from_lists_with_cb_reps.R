# Extract simulation data from lists and calculate alternative fit statistics

library(dplyr)
library(tidyr)
library(purrr)
library(pbmcapply)
library(stringr)
library(fungible)
library(here)
library(purrrgress)

condition_matrix <- readRDS(here("data", "condition_matrix.RDS"))

# Create results matrix ---------------------------------------------------

reps <- 500

# Create a matrix to hold results
results <- condition_matrix[rep(1:nrow(condition_matrix), each = reps),]

# For each rep in each condition, record the rep number
results$rep_num <- rep(1:reps, times = nrow(condition_matrix))

# Create model error method variable
results <- expand_grid(
  results,
  error_method = c("CB")
)

# Function definitions ----------------------------------------------------

# Compute the fit statistics that were not computed in the main simulation
# loop

compute_fit_stats <- function(error_method_results, Rpop, k) {
  
  # Check if value is NULL; if so, return empty list
  if (length(error_method_results$value) < 2) {
    if (is.null(error_method_results$value)) {
      out = data.frame(RMSEA_thetahat = NA,
                       CFI_thetahat = NA,
                       CRMR_theta = NA,
                       CRMR_thetahat = NA,
                       TLI_theta = NA,
                       TLI_thetahat = NA)
      return(out)
    } else if (is.na(error_method_results$value)) {
      out = data.frame(RMSEA_thetahat = NA,
                       CFI_thetahat = NA,
                       CRMR_theta = NA,
                       CRMR_thetahat = NA,
                       TLI_theta = NA,
                       TLI_thetahat = NA)
      return(out)
    }
  }
  
  RpopME <- pluck(error_method_results, "value", "Sigma")
  
  p <- ncol(RpopME)
  tp <- p * (p + 1)/2
  
  fout <- tryCatch({
    factanal(covmat = RpopME, factors = k, 
             rotation = "none", n.obs = 1000, 
             control = list(nstart = 100))
  }, error = function(e) NA)
  
  Tr <- function(X) sum(diag(X))
  
  # Compute degrees of freedom
  DF <- (p * (p - 1)/2) - (p * k) + (k * (k - 1)/2)
  DF_B <- (p * (p - 1)/2) - (p * p) + (p * (p - 1)/2)
  
  # If fout isn't NA, compute factor loadings and correlation matrix
  if (length(fout) > 1) {
    Fhat <- fout$loadings
    Sigma_k <- Fhat %*% t(Fhat)
    diag(Sigma_k) <- 1 
    
    num_thetahat <- log(det(Sigma_k)) - log(det(RpopME)) + 
      Tr(RpopME %*% solve(Sigma_k)) - p
    
    F_T_thetahat <- log(det(Sigma_k)) - log(det(RpopME)) + 
      Tr(RpopME %*% solve(Sigma_k)) - p
    F_B_thetahat <- -log(det(RpopME))
    
    # Calculate RMSEA thetahat
    RMSEA_thetahat <- sqrt(num_thetahat/DF)
    
    # Calculate CFI thetahat
    CFI_thetahat <- 1 - F_T_thetahat/F_B_thetahat
    
    # Calculate CRMR thetahat
    CRMR_thetahat <- sqrt(
      sum(((RpopME - Sigma_k)[upper.tri(Rpop, diag = FALSE)]^2)) / (tp - p)
    )
    
    TLI_thetahat <- 1 - ( (F_T_thetahat / DF) / (F_B_thetahat / DF_B) )
  } else {
    RMSEA_thetahat <- NA
    CFI_thetahat <- NA
    CRMR_thetahat <- NA
    TLI_thetahat <- NA
  }
  
  F_T_theta <- log(det(Rpop)) - log(det(RpopME)) + 
    Tr(RpopME %*% solve(Rpop)) - p
  F_B_theta <- -log(det(RpopME))
  
  # Calculate CRMR values
  CRMR_theta <- sqrt(
    sum(((RpopME - Rpop)[upper.tri(Rpop, diag = FALSE)]^2)) / (tp - p)
  )
  
  # Calculate TLI values (from Xia & Yang, 2019, p. 412)
  TLI_theta <- 1 - ( (F_T_theta / DF) / (F_B_theta / DF_B) )
  
  data.frame(RMSEA_thetahat = RMSEA_thetahat,
             CFI_thetahat = CFI_thetahat,
             CRMR_theta = CRMR_theta,
             CRMR_thetahat = CRMR_thetahat,
             TLI_theta = TLI_theta,
             TLI_thetahat = TLI_thetahat)
}

# Compute d values
compute_d_values <- function(error_method_data, target_rmsea, target_cfi) {
  
  if (length(error_method_data$value) < 2) {
    if (is.null(error_method_data$value)) {
      out = data.frame(d1 = NA,
                       d2 = NA,
                       d3 = NA)
      return(out)
    } else if (is.na(error_method_data$value)) {
      out = data.frame(d1 = NA,
                       d2 = NA,
                       d3 = NA)
      return(out)
    }
  }
  
  
  rmsea <- pluck(error_method_data, "value", "rmsea")
  cfi <- pluck(error_method_data, "value", "cfi")
  
  d1 <- abs(rmsea - target_rmsea)
  d2 <- abs(cfi - target_cfi)
  d3 <- d1 + d2 
  
  data.frame(
    d1 = d1,
    d2 = d2,
    d3 = d3
  )
}

# Check if there are any major factors in W
check_w_major_factors <- function(W) {
  if (!is.matrix(W)) {
    NA
  } else {
    sum(W[,1] >= .3) > 2
  }
}

# Create a vector of result file paths ------------------------------------

results_files <- list.files(
  path = here("data", "cb_results"),
  pattern = "cb_results_[0-9]+\\.RDS",
  full.names = TRUE
)

# Read data from each condition and calculate statistics ------------------

# When partially complete;
complete_results_files <- list.files(
  path = here("data", "cb_results"),
  pattern = "cb_results_matrix_[0-9]+\\.RDS"
)

complete_conditions <- str_extract(complete_results_files, "[0-9]+") %>%
  as.numeric()
incomplete_conditions <- setdiff(1:153, complete_conditions)

if (length(incomplete_conditions > 0)) {
  results_files <- results_files[incomplete_conditions]
}

# Calculate alternative fit statistics

pbmcapply::pbmclapply(
  X = seq_along(results_files), 
  FUN = function(i, results_files, condition_matrix) {
    
    # cat("/nWorking on condition:", i)
    
    # Read in loading matrix
    condition_results <- readRDS(results_files[i])
    
    # Extract the condition number from the results file name
    condition_num <- as.numeric(
      str_extract(results_files[i], pattern = "(?<=results_)[0-9]+")
    )
    
    # Extract condition information (loading strength, num. factors, num. items)
    j <- which(condition_matrix$condition_num == condition_num)
    loading_condition <- condition_matrix$loading[j]
    k <- condition_matrix$factors[j]
    p <- condition_matrix$items_per_factor[j] * k
    
    # Extract target RMSEA and CFI values
    target_rmsea <- condition_matrix$target_rmsea[j]
    target_cfi <- condition_matrix$target_cfi[j]
    
    # Create factor loading matrix
    loadings <- switch(loading_condition,
                       "weak" = .4,
                       "moderate" = .6,
                       "strong" = .8)
    
    # Extract Rpop
    mod <- fungible::simFA(Model = list(NFac = condition_matrix$factors[j],
                                        NItemPerFac = condition_matrix$items_per_factor[j],
                                        Model = "oblique"),
                           Loadings = list(FacLoadDist = "fixed",
                                           FacLoadRange = loadings),
                           Phi = list(PhiType = "fixed",
                                      MaxAbsPhi = condition_matrix$factor_cor[j]))
    Rpop <- mod$Rpop
    
    condition_results_matrix <- purrr::map_dfr(
      .x = seq_along(condition_results),
      .f = function(z, Rpop, p, k, target_rmsea, target_cfi,
                    condition_num) {
        
        rep_results <- condition_results[[z]]
        other_fit_stats <- compute_fit_stats(rep_results, Rpop = Rpop, k = k)
        d_values <- compute_d_values(rep_results, 
                                               target_rmsea = target_rmsea, 
                                               target_cfi = target_cfi)
        
        # PICK UP HERE; MODIFY DATA FRAME TO ONLY INLCUDE CB RESULTS.
        out <- data.frame(
          condition_num = condition_num,
          rep_num = z,
          cfi = pluck(rep_results, "value", "cfi", .default = NA),
          cfi_thetahat = other_fit_stats$CFI_thetahat,
          rmsea = pluck(rep_results, "value", "rmsea", .default = NA),
          rmsea_thetahat = other_fit_stats$RMSEA_thetahat,
          crmr = other_fit_stats$CRMR_theta,
          crmr_thetahat = other_fit_stats$CRMR_thetahat,
          tli = other_fit_stats$TLI_theta,
          tli_thetahat = other_fit_stats$TLI_thetahat,
          m = NA,
          v = NA,
          eps = NA,
          fn_value = pluck(rep_results, "value", "fn_value", .default = NA),
          error_method = c("cb"),
          w_has_major_factors = NA,
          d1 = d_values$d1,
          d2 = d_values$d2,
          d3 = d_values$d3,
          warning = ifelse(is.null(rep_results$warning), NA, as.character(rep_results$warning)),
          error = ifelse(is.null(rep_results$error), NA, as.character(rep_results$error))
        )
        
        rownames(out) <- NULL
        out
      },
      Rpop = Rpop, p = p, k = k, target_rmsea = target_rmsea, 
      target_cfi = target_cfi, condition_num
    )
    
    saveRDS(condition_results_matrix, 
            file = paste0("data/cb_results/cb_results_matrix_", 
                          formatC(condition_num, width = 3, flag = 0), 
                          ".RDS"))
    
  }, results_files = results_files, condition_matrix = condition_matrix, 
  mc.cores = 8
)
