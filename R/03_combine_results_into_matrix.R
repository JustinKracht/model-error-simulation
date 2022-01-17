# Combine results files into a results matrix
library(purrrgress)
library(tidyverse)
library(here)

results_files <- list.files(path = here("data"), 
                            pattern = "results_matrix_[0-9]+\\.RDS",
                            full.names = TRUE)
results_files_cb <- list.files(path = here("data", "cb_results"),
                               pattern = "cb_results_matrix_[0-9]+\\.RDS",
                               full.names = TRUE)

results_matrix <- pro_map(
  .x = results_files,
  .f = function(x) {
    result <- readRDS(x)
    condition <- as.numeric(str_extract(x, pattern = "(?<=results_matrix_)[0-9]+"))
    if (condition %in% 1:153) {
      result <- result %>% filter(error_method != "cb")
      result_cb <- readRDS(here("data", "cb_results",
                                paste0("cb_results_matrix_",
                                       formatC(condition, width = 3, flag = "0"),
                                       ".RDS")))
      result <- result %>% 
        bind_rows(result_cb) %>%
        arrange(condition, rep_num)
    } 
    result
  }
)

results_matrix <- do.call(bind_rows, results_matrix)
saveRDS(results_matrix, here("data", "results_matrix.RDS"))
