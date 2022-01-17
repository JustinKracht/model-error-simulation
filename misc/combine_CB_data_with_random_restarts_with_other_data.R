# The first simulation run used the MBESS::Sigma.2.SigmaStar function, which
# is deterministic (all inputs give the same output). The cb_results files
# have CB results from a modified function that is not deterministic so that
# the 500 reps are actually doing something useful.

library(here)

for (condition in 1:153) {
  cat("\nWorking on condition:", condition)
  
  # Read in results files
  RESULTS_LIST_PATH <- here(
    "data", 
    paste0("results_", formatC(condition, width = 3, flag = "0"), ".RDS")
  )
  
  CB_RESULTS_LIST_PATH <- here(
    "data",
    "cb_results",
    paste0("cb_results_", formatC(condition, width = 3, flag = "0"), ".RDS")
  )
  
  results_list <- readRDS(RESULTS_LIST_PATH)
  cb_results_list <- readRDS(CB_RESULTS_LIST_PATH)
  
  # Replace the CB elements of the results list
  for (rep in 1:500) {
    results_list[[rep]]$sigma_cb <- cb_results_list[[rep]]
  }
  
  saveRDS(results_list, RESULTS_LIST_PATH)
}
