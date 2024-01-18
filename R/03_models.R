library(survival)
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
dat <- sim_survdat_f(nsample = 1000, varnum = 25, dist = "g", lambda = 0.01, rho = 1, beta = c(1.8,0.5,0.4,-0.4,0.45,0.6,-0.6), crate = 0.001, cor = TRUE, seed = 20231106)
fit <- survival::coxph(Surv(time, status) ~ X1 + X2 + X3 + X4+ X5+X1*X2+X3*X4, data = dat)
load('data/cox_dat.rdata')
library(caret)

pacman::p_load(glmnet)

pacman::p_load(caret)

HDSI_model_f <- function(boot_df_train = NA,
                         # boot_df_test = NA,
                         method = "lasso") {
  library(pacman)
  pacman::p_load(glmnet)

  ### x input in cv.glmnet
  y <- "survival::Surv(time,status)"
  ## f=stats::as.formula(paste(y," ~ .*.*."))
  f <- stats::as.formula(paste(y, " ~ .*.")) ## pair-wise interaction
  X <- stats::model.matrix(f, boot_df_train)[, -1]
  #X_test <- stats::model.matrix(f, boot_df_test)[, -1]
  ### y input in cv.glmnet
  time <- boot_df_train$time
  status <- boot_df_train$status
  Y <- Surv(boot_df_train$time, boot_df_train$status)
  #Y_test <- Surv(boot_df_test$time, boot_df_test$status)

  if (method == "lasso") {
    cv_fit <- glmnet::cv.glmnet(
      x = X, y = Y, alpha = 1,
      standardize = F, family = "cox", nfolds = 5
    )
    lambda.1se <- cv_fit$lambda.1se
    lambda.min <- cv_fit$lambda.min

    model_fit <- glmnet::glmnet(x = X, y = Y, lambda = lambda.1se, alpha = 0, standardize = F, family = "cox")
    model_coef <- as.matrix(coef(model_fit, s = lambda.1se)) %>%
      data.frame() %>%
      mutate(varname = rownames(.))

    pred <- predict(model_fit, s = lambda.1se, newx = X, type = "link")

    model_cindex <- apply(pred, 2, Cindex, y = Y)

    model_out <- model_coef %>% mutate(model_cindex = model_cindex)

    # y1 = cbind(time = time, status = status)
    #
    # apply(risk_scores, 2, glmnet::Cindex, y=y1)
  }
  if (method == "ridge") {
    cv_fit <- glmnet::cv.glmnet(
      x = X, y = Y, alpha = 0,
      standardize = F, family = "cox", nfolds = 5
    )

    lambda.1se <- cv_fit$lambda.1se
    lambda.min <- cv_fit$lambda.min

    model_fit <- glmnet::glmnet(x = X, y = Y, lambda = lambda.1se, alpha = 1, standardize = F, family = "cox")

    model_coef <- as.matrix(coef(model_fit, s = lambda.1se)) %>%
      data.frame() %>%
      mutate(varname = rownames(.))

    pred <- predict(model_fit, s = lambda.1se, newx = X, type = "link")

    model_cindex <- apply(pred, 2, Cindex, y = Y)

    model_out <- model_coef %>% mutate(model_cindex = model_cindex)
  }

  return(model_out)
}

uni_cox_res2 <- uni_cox_res2 %>%
  mutate(index = grepl("X_",var_name))

save(uni_cox_res2,file = "output/uni_cox_res2.rdata")

uni_gene_res <- uni_cox_res2 %>% filter(index==TRUE) %>%
  filter(BH<=0.05)

cox_sig_dat <- cox_dat[,colnames(cox_dat) %in% c("time","status","age_at_diagnosis","cancer_stage",uni_gene_res$var_name)]

t1 = Sys.time()
boot_lst <- boot_f1(
  df = cox_sig_dat,
  nboot = 100,
  boot_ft = 500,
  seed = 20231106
)


# trainIndex <- createDataPartition(boot_df$status, p = .8,
#                                   list = FALSE,
#                                   k =5,
#                                   times = 1)

#folds <- createFolds(boot_df$status, k = 5, list = TRUE, returnTrain = FALSE)




m <- c("lasso")
inter_perf_f <- function(boot_lst=NA,
                         method = "lasso",
                         qtl = NA,
                         Rf = NA){
  m <- method
  model_out_all <- lapply(1:length(m), function(x) {
    lapply(1:length(boot_lst), function(y) {
      folds <- createFolds(boot_lst[[y]]$status, k = 5, list = TRUE, returnTrain = FALSE)

      ## the validation performance
      lapply(1:length(folds),function(z){

        res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][-folds[[z]],],
                             # boot_df_test = boot_lst[[y]][folds[[z]],],
                             method = m[[x]])
        res$boot = y
        res$fold = -z
        res
      })
      # %>% bind_rows() %>%
      #   group_by(varname) %>%
      #   summarise(X1 = mean(X1),
      #             model_cindex = mean(model_cindex))

    })
  })

  fs_out5 <- lapply(1:length(m), function(x) {
    ## bind B boot sample results into a dataframe
    perf_tmp <- model_out_all[[x]] %>% bind_rows()
    lapply(1:5,function(y){
      perf_onefold <-  perf_tmp[perf_tmp$fold==-y,]

      fe_star <-  perf_onefold %>%  group_by(model_cindex) %>% summarise(n=n())
      ### compute min_cindex
      min_cindex <- perf_onefold %>%
        group_by(varname) %>%
        summarise(min_cindex = min(model_cindex)) %>%
        mutate(
          miu_min_cindex = mean(min_cindex),
          # sigma_min_cindex = sqrt(mean((min_cindex - miu_min_cindex)^2) / (fe_star$n %>% unique() - 1)),
          sd_min_cindex = sd(min_cindex),
          Rf = Rf, ## hyperparameter
          include_yn = min_cindex > (miu_min_cindex + Rf * sd_min_cindex)
        )

      ### compute coef estimate and CI
      beta_quantile <- perf_onefold %>%
        group_by(varname) %>%
        mutate(quantile = qtl) %>% ## quantile is a hyperparameter
        summarise(
          qtl = qtl,
          beta_hat = mean(X1),
          qtl_lower = quantile(X1, probs = qtl/2),
          qtl_upper = quantile(X1, probs = 1 - (qtl /2) ),
          include0_yn = !(((qtl_lower <= 0) & (0 <= qtl_upper)) | ((qtl_upper <= 0) & (0 <= qtl_lower)))
        )

      ### join beta and cindex in a dataframe; filter selected features
      perf <- full_join(beta_quantile, min_cindex, by = "varname") %>%
        filter(include0_yn == TRUE & include_yn == TRUE)

      ### adding back main effect if only interaction terms are selected
      included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

      ## detect main effect terms from interaction
      included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

      ## final selected features
      included_fe <- c(perf$varname, included_inter1) %>% unique()
      #included_fe
      list(included_fe,perf)

    })

  })

  ## calculate averaged performance over 5 datasets; internal validation
  inter_cindex_tmp <- lapply(1:length(m),function(x){
    lapply(1:5,function(y){
      fs_out5[[x]][[y]][[2]] %>% bind_rows()
    })
  })

  inter_val_hyper <- lapply(1:length(m),function(x){
    avg_cindex <- inter_cindex_tmp[[x]] %>% bind_rows()%>%
      summarise(avg_cindex = mean(min_cindex))
    hyper_table <- data.frame(method = m[x],
                              avg_cindex = avg_cindex$avg_cindex,
                              hyper_quantile = qtl,
                              Rf = Rf
    )
  })


  return(inter_val_hyper)

}

t3 = Sys.time()
test <- inter_perf_f(boot_lst=boot_lst,
                         method = c("lasso"),
                         qtl = 0.02,
                         Rf = 0.6)
t4 = Sys.time()
#
#
# test <- inter_perf_f(boot_lst=boot_lst,
#                      method = c("lasso","ridge"),
#                      qtl = 0.05,
#                      Rf = 1)
## hyperparameter tuning--------------

# Define hyperparameter grid
hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.02, by = 0.01),
                                   Rf = seq(0.1, 1, by = 0.5))

library(parallel)
# parallelize hyperparameter tuning
parallel_tuning <- function(row) {
  qtl_value <- row['qtl']
  Rf_value <- row['Rf']

  # Call function with current hyperparameter values
  performance_metric_value <- inter_perf_f(boot_lst=boot_lst, method = c("lasso", "ridge"), qtl = qtl_value, Rf = Rf_value)

  # Return results for each row
  c(qtl = qtl_value, Rf = Rf_value, performance_metric = performance_metric_value)
}

# Set up a parallel cluster
cl <- makeCluster(detectCores())

# Export necessary functions to the worker nodes
clusterExport(cl, c("inter_perf_f","boot_lst","HDSI_model_f"))

#
clusterEvalQ(cl, {
  library(caret)
  library(survival)
  library(dplyr)
  library(stringr)
})

# Use pblapply for parallel computing with progress bar
results_list <- pbapply::pblapply(split(hyperparameter_grid, 1:nrow(hyperparameter_grid)), parallel_tuning, cl = cl)

results_list[[2]]
# Stop the parallel cluster
stopCluster(cl)

### find the best perfromance hyperparameters and




internal_res <- lapply(1:length(m),function(x){
    lapply(1:nrow(hyperparameter_grid),function(y){
      results_list[[y]][[x+2]] %>% bind_rows() ## 2 is the number of hyperparamters in the grid df
    }) %>% bind_rows()
  }) %>%  bind_rows()

opt_dat <-left_join(internal_res,
  internal_res %>% group_by(method) %>%
  summarise(opt = max(avg_cindex)),by ='method') %>%
  filter(avg_cindex == opt)


### until now, we have the best hyper parameters, the we select the data
opt_dat_tmp <- opt_dat %>% filter(method=='lasso')
opt_qtl <- opt_dat_tmp$avg_cindex

t2 = Sys.time()

trainIndex <- createDataPartition(cox_dat$status,
                                  p = 0.8, list = FALSE)

trainingSet <- cox_dat[trainIndex,]
testSet <- cox_dat[-trainIndex,]


lung_boot_lst <- boot_f1(
  df = trainingSet,
  nboot = 100,
  boot_ft = 500,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)

boot_lst <- lung_boot_lst
m = "lasso"
opt_qtl =
model_out_all <- lapply(1:length(m), function(x) {
  lapply(1:length(boot_lst), function(y) {
    HDSI_model_f(boot_df = boot_lst[[y]], method = m[[x]])
  })
})

Sys.time()
selected_dat <- lapply(1:length(m), function(x) {
  ## bind B boot sample results into a dataframe
  perf_tmp <- model_out_all[[x]] %>% bind_rows()
  fe_star <- perf_tmp %>% group_by(model_cindex) %>% summarise(n=n())
  ### compute min_cindex
  min_cindex <- perf_tmp %>%
    group_by(varname) %>%
    summarise(min_cindex = min(model_cindex)) %>%
    mutate(
      miu_min_cindex = mean(min_cindex),
      # sigma_min_cindex = sqrt(mean((min_cindex - miu_min_cindex)^2) / (fe_star$n %>% unique() - 1)),
      sd_min_cindex = sd(min_cindex),
      Rf = 2, ## hyperparameter
      include_yn = min_cindex > (miu_min_cindex + Rf * sd_min_cindex)
    )

  ### compute coef estimate and CI
  beta_quantile <- perf_tmp %>%
    group_by(varname) %>%
    mutate(quantile = opt_qtl) %>% ## quantile is a hyperparameter
    summarise(
      qtl = opt_qtl,
      beta_hat = mean(X1),
      qtl_lower = quantile(X1, probs = qtl/2),
      qtl_upper = quantile(X1, probs = 1 - (qtl /2) ),
      include0_yn = !(((qtl_lower <= 0) & (0 <= qtl_upper)) | ((qtl_upper <= 0) & (0 <= qtl_lower)))
    )

  ### join beta and cindex in a dataframe; filter selected features
  perf <- full_join(beta_quantile, min_cindex, by = "varname") %>%
    filter(include0_yn == TRUE & include_yn == TRUE)

  ### adding back main effect if only interaction terms are selected
  included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

  ## detect main effect terms from interaction
  included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

  ## final selected features
  included_fe <- c(perf$varname, included_inter1) %>% unique()

  res <- list(perf,included_fe)

  res
})

## [[1]]代表lasso方法；
selected_vars_dat <- data.frame(var_name =selected_dat[[1]][[2]]) %>%
mutate(new_var = str_replace_all(var_name,":","\\*"),
       index = grepl("age",selected_vars_dat$new_var)) %>%
  filter(index==FALSE)

select_gene <- selected_vars_dat$new_var





selected_genes <- paste(select_gene, collapse = "+")
gene_string <- capture.output(cat(selected_genes)) %>% paste(collapse = "\n")

surv_object <- Surv(time = cox_dat$time,
                    event = cox_dat$status)
formula <- as.formula(paste("surv_object ~",
                            paste(gene_string,
                                  'age_at_diagnosis+cancer_stage',sep = "+")))
cox_model <- coxph(formula, data = cox_dat)
result_summary <- summary(cox_model)
 <- result_summary$coefficients %>% data.frame() %>% bind_rows()

c_index <- concordance(cox_model)
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




#####supporting codes---------------
### internal performance function


fs_out5 <- lapply(1:length(m), function(x) {
  ## bind B boot sample results into a dataframe
  perf_tmp <- model_out_all[[x]] %>% bind_rows()
  lapply(1:5,function(y){
   perf_onefold <-  perf_tmp[perf_tmp$fold==-y,]

    fe_star <-  perf_onefold %>%  group_by(model_cindex) %>% summarise(n=n())
    ### compute min_cindex
    min_cindex <- perf_onefold %>%
      group_by(varname) %>%
      summarise(min_cindex = min(model_cindex)) %>%
      mutate(
        miu_min_cindex = mean(min_cindex),
        # sigma_min_cindex = sqrt(mean((min_cindex - miu_min_cindex)^2) / (fe_star$n %>% unique() - 1)),
        sd_min_cindex = sd(min_cindex),
        Rf = 1, ## hyperparameter
        include_yn = min_cindex > (miu_min_cindex + Rf * sd_min_cindex)
      )

    ### compute coef estimate and CI
    beta_quantile <- perf_onefold %>%
      group_by(varname) %>%
      mutate(quantile = qtl) %>% ## quantile is a hyperparameter
      summarise(
        qtl = qtl,
        beta_hat = mean(X1),
        qtl_lower = quantile(X1, probs = qtl/2),
        qtl_upper = quantile(X1, probs = 1 - (qtl /2) ),
        include0_yn = !(((qtl_lower <= 0) & (0 <= qtl_upper)) | ((qtl_upper <= 0) & (0 <= qtl_lower)))
      )

    ### join beta and cindex in a dataframe; filter selected features
    perf <- full_join(beta_quantile, min_cindex, by = "varname") %>%
      filter(include0_yn == TRUE & include_yn == TRUE)

    ### adding back main effect if only interaction terms are selected
    included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

    ## detect main effect terms from interaction
    included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

    ## final selected features
    included_fe <- c(perf$varname, included_inter1) %>% unique()
    #included_fe
    list(included_fe,perf)

  })

})

## calculate averaged performance over 5 datasets; internal validation
inter_cindex_tmp <- lapply(1:length(m),function(x){
  lapply(1:5,function(y){
    fs_out5[[x]][[y]][[2]] %>% bind_rows()
  })
})

inter_cindex <- lapply(1:length(m),function(x){
  avg_cindex <- inter_cindex_tmp[[x]] %>% bind_rows()%>%
    summarise(avg_cindex = mean(min_cindex))
  hyper_table <- data.frame(method = m[x],avg_cindex = avg_cindex$avg_cindex,
                            hyper_quantile = qtl
  )
})






## fit model on whole training set



  ### compute coef estimate and CI
  beta_quantile <- perf_onefold %>%
    group_by(varname) %>%
    mutate(quantile = 0.05) %>% ## quantile is a hyperparameter
    summarise(
      qtl = 0.05,
      beta_hat = mean(X1),
      qtl_lower = quantile(X1, probs = qtl/2),
      qtl_upper = quantile(X1, probs = 1 - (qtl /2) ),
      include0_yn = !(((qtl_lower <= 0) & (0 <= qtl_upper)) | ((qtl_upper <= 0) & (0 <= qtl_lower)))
    )

  ### join beta and cindex in a dataframe; filter selected features
  perf <- full_join(beta_quantile, min_cindex, by = "varname") %>%
    filter(include0_yn == TRUE & include_yn == TRUE)

  ### adding back main effect if only interaction terms are selected
  included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

  ## detect main effect terms from interaction
  included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

  ## final selected features
  included_fe <- c(perf$varname, included_inter1) %>% unique()

  res <- list(included_fe,perf)


})

Sys.time()

tmp21 <- tmp[[2]][[1]] %>% data.frame()
# var_freq <- stringr::str_split(min_cindex$varname,":") %>% unlist() %>% table() %>% data.frame()
