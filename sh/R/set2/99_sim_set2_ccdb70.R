getwd()
source('R/libs.R')
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
source("R/helper_f.R")
#load("data/lung_boot_lst.rdata")
#load("data/sim_boot_lst.rdata")
#load("output/sim_boot_lst.rdata")

getwd()
source('R/libs.R')
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
source("R/helper_f.R")
#load("data/lung_boot_lst.rdata")
#load("data/sim_boot_lst.rdata")
#load("output/sim_boot_lst.rdata")

# load("data/sim_set1.rdata")

# Set the path to the folder containing .rdata files
folder_path <- "data/sim_set"

# Get a list of all .rdata files in the folder
rdata_files <- list.files(folder_path, pattern = "\\.rdata$", full.names = TRUE)

# Load all .rdata files into the global environment
lapply(rdata_files, function(file) {
  load(file, envir = .GlobalEnv)
})



sim_data <- TRUE

cov_interact <- TRUE
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
# hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.05, by = 0.01),
#                                    Rf = seq(0,2,by=0.1))

# hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.15, by = 0.01),
#                                    Rf = seq(-2.5,2.5,by=0.1))

hyperparameter_grid_tmp <- expand.grid(qtl = seq(0.01, 0.5, by = 0.1),
                                       Rf = seq(-2.5,4,by=0.5))

num_holder <- seq(1,nrow(hyperparameter_grid_tmp),4)+3

i <- 18
hyperparameter_grid <-hyperparameter_grid_tmp[c( (num_holder[i]-3):nrow(hyperparameter_grid_tmp)),]


library(parallel)

sim_set2_hyper_out <- NULL
for(i in c(1:length(sim_set2))){
  print(paste("i=",i))
  boot_lst <- sim_set2[[i]]
  parallel_tuning <- function(row) {
    qtl_value <- row['qtl']
    Rf_value <- row['Rf']

    # Call function with current hyperparameter values
    performance_metric_value <- inter_perf_f(boot_lst=boot_lst,
                                             cov_interact = cov_interact,
                                             k = 5,
                                             method = method,
                                             qtl = qtl_value, Rf = Rf_value)

    # Return results for each row
    c(qtl = qtl_value, Rf = Rf_value, performance_metric = performance_metric_value)
  }
  # parallelize hyperparameter tuning


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
  internal_res_tmp <- lapply(1:length(m),function(x){
    lapply(1:nrow(hyperparameter_grid),function(y){
      results_list[[y]][[x+2]] %>% bind_rows() ## 2 is the number of hyperparamters in the grid df
    }) %>% bind_rows()
  }) %>%  bind_rows()

  sim_set2_hyper_out <- rbind(sim_set2_hyper_out,internal_res_tmp)
}

## in case overwrite
randomid <- 70
save(sim_set2_hyper_out,file = paste0('data/sim_set2_hyper_out_',Sys.Date(),randomid,".rdata"))


