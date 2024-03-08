getwd()
source('R/libs.R')
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
source("R/helper_f.R")
load("data/lung_boot_lst.rdata")
load("data/sim_boot_lst.rdata")
#load("output/sim_boot_lst.rdata")
boot_lst <- lung_boot_lst

cov_interact <- FALSE
method <- c('lasso',"ridge")

## check if HDSI_model_f works it does work!!
# HDSI_model_f(boot_df_train = lung_boot_lst[[1]],
#              cov_interact = FALSE,
#              method= 'lasso'
#               )


# inter_perf_f(boot_lst=lung_boot_lst[c(1,2)], ## list input
#                      method = c("ridge"),
#                      k = 2,
#                      cov_interact = FALSE,
#                      qtl = 0.02,
#                      Rf = 0)


# Define hyperparameter grid
hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.05, by = 0.01),
                                   Rf = seq(0,2,by=0.1))

hyperparameter_grid <- expand.grid(qtl = 0.05,
                                   Rf = 0.3)

library(parallel)
# parallelize hyperparameter tuning
parallel_tuning <- function(row) {
  qtl_value <- row['qtl']
  Rf_value <- row['Rf']

  # Call function with current hyperparameter values
  performance_metric_value <- inter_perf_f(boot_lst=boot_lst,
                                           cov_interact = cov_interact,
                                           k = 2,
                                           method = method,
                                           qtl = qtl_value, Rf = Rf_value)

  # Return results for each row
  c(qtl = qtl_value, Rf = Rf_value, performance_metric = performance_metric_value)
}

# Set up a parallel cluster
cl <- makeCluster(detectCores())

# Export necessary functions to the worker nodes
clusterExport(cl, c("inter_perf_f","boot_lst","cov_interact",'method',"HDSI_model_f"))

#
clusterEvalQ(cl, {
  library(caret)
  library(survival)
  library(dplyr)
  library(stringr)
})

# Use pblapply for parallel computing with progress bar
results_list <- pbapply::pblapply(split(hyperparameter_grid, 1:nrow(hyperparameter_grid)), parallel_tuning, cl = cl)


# Stop the parallel cluster
stopCluster(cl)

### find the best perfromance hyperparameters and



m <- method
real_dat_internal_res <- lapply(1:length(m),function(x){
  lapply(1:nrow(hyperparameter_grid),function(y){
    results_list[[y]][[x+2]] %>% bind_rows() ## 2 is the number of hyperparamters in the grid df
  }) %>% bind_rows()
}) %>%  bind_rows()

save(real_dat_internal_res,file = 'data/real_dat_internal_res.rdata')

