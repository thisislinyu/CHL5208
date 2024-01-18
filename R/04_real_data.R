## UVA: cox apply to real-data set-------
load('data/cox_dat.rdata')
library(caret)

## define univariate cox function
uni_cox_f <- function(data=NA, ## outcome variable names must be time and status
                      gen_var_index = 5, ## from which the gene data is started
                      cov_vars = 'age_at_diagnosis+cancer_stage',
                      ...
){
  surv_object <- Surv(time = data$time,
                      event = data$status
  )
  vars = colnames(data[-c(1:gen_var_index-1)])
  results <- lapply(vars, function(var) {
    #print(vars)
    var = make.names(var)
    # formula <- as.formula(paste("surv_object ~", var))
    # cox_model <- coxph(formula, data = data)
    # result_summary <- summary(cox_model)
    # result_summary$coefficients %>% data.frame() %>% bind_rows()

    tryCatch({
      formula <- as.formula(paste("surv_object ~", paste(cov_vars,var,sep = "+")))
      cox_model <- coxph(formula, data = data)
      result_summary <- summary(cox_model)
      result_summary$coefficients %>% data.frame() %>% bind_rows()
    }, error = function(err) {
      # 在发生错误时返回当前结果
      message(paste("Error for variable:", var, " - ", conditionMessage(err)))
      return(data.frame(Variable = var, Error = conditionMessage(err)))
    })
  })

  return(results)
}

t1 = Sys.time()
uni_cox_res <- uni_cox_f(data = cox_dat,
                         gen_var_index = 5,
                         cov_vars = 'age_at_diagnosis+cancer_stage')
t2 = Sys.time()
t2-t1
uni_cox_res1 <- uni_cox_res %>% bind_rows()

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("multtest")
library(rio)
#install.packages("multtest")
library(multtest)


adjusted_p_values <- mt.rawp2adjp(rawp = uni_cox_res1$Pr...z.., proc= "BH")

adjp_dat <- adjusted_p_values$adjp %>% data.frame()

uni_cox_res2 <- cbind(uni_cox_res1,adjp_dat)
uni_cox_res2$var_name <- rownames(uni_cox_res2)
save(uni_cox_res2,file = 'output/uni_cox_res2.rdata')

# adjp_dat$BH %>% summary()
# adjp_dat$BH %>% hist()

## apply our algorithm---
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
## boot_f1 for real data
boot_f1 <- function(df = NA,
                   nboot = 100,
                   boot_ft = 5,
                   gene_index = 5,
                   cov_vars = c('age_at_diagnosis','cancer_stage'),
                   seed = 20231106,
                   ...) {
  ## get the feature names in df total features
  tot_ft <- df[,-c(1:gene_index-1)] %>%
    # starts_with("X") # ^ start sign;
    # $ end of a string, d+ one or more digital
    ## This may be problematic if df is not generated from sim_surv_f
    colnames()

  set.seed(seed)

  boot_lst <- lapply(1:nboot, function(x) {
    ### sample boot_ft from tot_ft
    boot_feature <- sample(tot_ft, boot_ft, replace = FALSE) %>%
      naturalsort::naturalsort()

    ### sample nrow(df) people from df with replacement;

    boot_sample <- df %>%
      #### sample status = 1
      filter(status == 1) %>%
      sample_n(size = nrow(.), replace = TRUE) %>%
      bind_rows(df %>% filter(status == 0) %>%
                  #### sample status = 0
                  sample_n(size = nrow(.), replace = TRUE)) %>%
      select(time, status,all_of(cov_vars) ,all_of(boot_feature))
  })

  return(boot_lst)
}
trainIndex <- createDataPartition(cox_dat$status,
                                  p = 0.8, list = FALSE)

trainingSet <- cox_dat[trainIndex,]
testSet <- cox_dat[-trainIndex,]


lung_boot_lst <- boot_f1(
  df = trainingSet,
  nboot = 13,
  boot_ft = 5,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)



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

test <- HDSI_model_f(boot_df_train = tmp1,
                         # boot_df_test = NA,
                         method = "lasso")

## till now we have the results for one bootstrap set
## we need to get n bootstrap results, using lapply

## and finally, we select features by averaging the bootstrap results

dat$outcome %>% table()
trainingSet$outcome %>% table()
testSet$outcome %>% table()



