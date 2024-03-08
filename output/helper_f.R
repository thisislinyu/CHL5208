## boot_f1 for real data

# > colnames(reg_dat) %>% head()
# [1] "time"
# [2] "status"
# [3] "age_at_diagnosis"
# [4] "cancer_stage"
# [5] "X_CIDEC"
# [6] "X_LOC441869"

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

HDSI_model_f <- function(boot_df_train = NA,
                         # boot_df_test = NA,
                         cov_interact = FALSE,
                         method = "lasso") {

  library(glmnet)

  ### x input in cv.glmnet
  y <- "survival::Surv(time,status)"
  ## f=stats::as.formula(paste(y," ~ .*.*."))
  f <- stats::as.formula(paste(y, " ~ .*.")) ## pair-wise interaction
  f1 <- stats::as.formula(paste(y, " ~ ."))

  if(cov_interact == TRUE){
    X <- stats::model.matrix(f, boot_df_train)[, -1]
  }else{
    X_tmp <-  boot_df_train %>%
      select(time, status,age_at_diagnosis,cancer_stage) %>%
      bind_cols(
        stats::model.matrix(f, boot_df_train %>%
                              select(-age_at_diagnosis,-cancer_stage))[, -1] )
    X <- stats::model.matrix(f1, X_tmp)[,-1]
  }

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
    print('so far so good 1')
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
  print("so far so good ")
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


inter_perf_f <- function(boot_lst=NA,
                         cov_interact = FALSE,
                         method = "lasso",
                         k = 5, ## cv fold number
                         qtl = NA,
                         Rf = NA){
  library(caret)
  m <- method
  model_out_all <- lapply(1:length(m), function(x) {
    lapply(1:length(boot_lst), function(y) {
      folds <- createFolds(boot_lst[[y]]$status, k = k, list = TRUE, returnTrain = FALSE)

      ## the validation performance
      lapply(1:length(folds),function(z){

        res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][-folds[[z]],],
                             cov_interact = cov_interact,
                             # boot_df_test = boot_lst[[y]][folds[[z]],],
                             method = m[[x]])
        res$boot = y
        res$fold = -z
        res
      })
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
