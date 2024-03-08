

# adjp_dat$BH %>% summary()
# adjp_dat$BH %>% hist()

## apply our algorithm---
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")




lung_boot_lst <- boot_f1(
  df = reg_dat,
  nboot = 100,
  boot_ft = 20,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)

reg_dat <-cox_dat %>%
  select("time", "status", "age_at_diagnosis", "cancer_stage",
         (uni_cox_top200$var_name))

lung_boot_lst <- boot_f1(
  df = reg_dat,
  nboot = 1,
  boot_ft = 20,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)

save(lung_boot_lst, file = "output/lung_boot_lst.rdata")

# tmp <- lapply(1:length(lung_boot_lst), function(x){
#   trainIndex <- createDataPartition(cox_dat$status,
#                                     p = 0.8, list = FALSE)
#
#   trainingSet <- cox_dat[trainIndex,]
#   testSet <- cox_dat[-trainIndex,]
#
#   return(list(trainingSet,testSet))
# })

#
tmp1 <- tmp[[1]]
#
# tmp2 <- boot_lst[[1]]

test <- HDSI_model_f(boot_df_train = lung_boot_lst[[1]],
                         # boot_df_test = NA,
                         method = "lasso")

## till now we have the results for one bootstrap set
## we need to get n bootstrap results, using lapply

## and finally, we select features by averaging the bootstrap results

dat$outcome %>% table()
trainingSet$outcome %>% table()
testSet$outcome %>% table()



