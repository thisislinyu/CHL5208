# Install and load required packages
install.packages(c("survival", "foreach", "doParallel"))
library(survival)
library(foreach)
library(doParallel)
library(dplyr)

# Set the number of cores
num_cores <- detectCores()

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Assuming your data is in a data frame called 'your_data' with a column 'time' for survival times and 'status' for censoring
# and columns V1, V2, ..., V20000 for your independent variables

# Create a survival object

work_dat1 <- work_dat1 %>%
  mutate(OS.time_month = OS.time/30) %>%
  filter(!is.na(age_at_diagnosis)& !is.na(ajcc_pathologic_stage))

surv_object <- Surv(time = cox_uni_dat$OS.time_month, event = cox_uni_dat$OS)

# Create a data frame with independent variables
# ncol(cox_uni_dat)
gene_vars <- work_dat1[, 26:ncol(work_dat1)] %>%
  apply(.,2,as.numeric)
colnames(gene_vars) <- rlang::set_names(paste0('X_',
                                               make.names(colnames(gene_vars))))# Assuming your independent variables start from the third column




