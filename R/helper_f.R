## boot_f1 for real data

# > colnames(reg_dat) %>% head()
# [1] "time"
# [2] "status"
# [3] "age_at_diagnosis"
# [4] "cancer_stage"
# [5] "X_CIDEC"
# [6] "X_LOC441869"

boot_f <- function(df = NA,
                   nboot = 100,
                   boot_ft = 5,
                   seed = 20231106) {
  ## get the feature names in df total features
  tot_ft <- df %>% select(matches("^X\\d+$")) %>%
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
      select(time, status, all_of(boot_feature))
  })

  return(boot_lst)
}



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
      standardize = T, family = "cox", nfolds = 3
    )
  #  print('so far so good lasso')
    lambda.1se <- cv_fit$lambda.1se
    lambda.min <- cv_fit$lambda.min

    model_fit <- glmnet::glmnet(x = X, y = Y, lambda = lambda.1se, alpha = 1, standardize = T, family = "cox")
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
  #  print("so far so good (ridge) ")
    cv_fit <- glmnet::cv.glmnet(
      x = X, y = Y, alpha = 0,
      standardize = T, family = "cox", nfolds = 3
    )

    lambda.1se <- cv_fit$lambda.1se
    lambda.min <- cv_fit$lambda.min

    model_fit <- glmnet::glmnet(x = X, y = Y, lambda = lambda.1se, alpha = 0, standardize = T, family = "cox")

    model_coef <- as.matrix(coef(model_fit, s = lambda.1se)) %>%
      data.frame() %>%
      mutate(varname = rownames(.))

    pred <- predict(model_fit, s = lambda.1se, newx = X, type = "link")

    model_cindex <- apply(pred, 2, Cindex, y = Y)

    model_out <- model_coef %>% mutate(model_cindex = model_cindex)
  }

  return(model_out)
}

#boot_lst <- sim_set1[[11]]
inter_perf_f <- function(boot_lst=NA,
                         cov_interact = NA,
                         method = "lasso",
                         k = 2, ## cv fold number
                         qtl = NA,
                         Rf = NA){
  library(caret)
  library(dplyr)
  library(survival)
  library(stringr)
  m <- method
  model_out_all <- lapply(1:length(m), function(x) {
    lapply(1:length(boot_lst), function(y) {
      print(y)
      folds <- createFolds(boot_lst[[y]]$status, k = k, list = TRUE, returnTrain = FALSE)

      ## I want to know how many true effects are
      ## selected in each boot dataset
      selected_true <- c("X1","X2","X3","X4","X5") %in%
        (boot_lst[[y]] %>% colnames())

      X1X2 <- ifelse((selected_true[1]==TRUE & selected_true[2]==TRUE),
                     TRUE, FALSE )
      X3X4 <- ifelse((selected_true[3]==TRUE & selected_true[4]==TRUE),
                     TRUE, FALSE )
      selected_X1 <- c("X1") %in%
        (boot_lst[[y]] %>% colnames())
      selected_X2 <- c("X2") %in%
        (boot_lst[[y]] %>% colnames())
      selected_X3 <- c("X3") %in%
        (boot_lst[[y]] %>% colnames())
      selected_X4 <- c("X4") %in%
        (boot_lst[[y]] %>% colnames())
      selected_X5 <- c("X5") %in%
        (boot_lst[[y]] %>% colnames())

      boot_true_margin <- selected_true %>% sum()
      boot_true_inter <- sum(X1X2,X3X4)


      ## the validation performance
      if(k!=1){
      res <- lapply(1:length(folds),function(z){
       # print(paste0('This is the ',y,'th boot sample ', z,'th fold'))
        res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][-folds[[z]],],
                             cov_interact = cov_interact,
                             # boot_df_test = boot_lst[[y]][folds[[z]],],
                             method = m[[x]])
        res$boot = y
        res$fold = -z
        res$boot_true_margin = boot_true_margin
        res$boot_true_inter = boot_true_inter
        res$boot_trueX1= selected_X1
        res$boot_trueX2= selected_X2
        res$boot_trueX3= selected_X3
        res$boot_trueX4= selected_X4
        res$boot_trueX5= selected_X5

      })
      }
      if(k==1){
        res <-  lapply(1:length(folds),function(z){
          print(paste0('This is the ',y,'th boot sample ', z,'th fold'))
          res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][folds[[z]],],
                               cov_interact = cov_interact,
                               # boot_df_test = boot_lst[[y]][folds[[z]],],
                               method = m[[x]])
          res$boot = y
          res$fold = -z
          res$boot_true_margin = boot_true_margin
          res$boot_true_inter = boot_true_inter
          res

        })
      }

      res

    })
  })

  print('model_out_all computation done!! good job!!!')
  fs_out5 <- lapply(1:length(m), function(x) {
    ## bind B boot sample results into a dataframe
    perf_tmp <- model_out_all[[x]] %>% bind_rows()
    lapply(1:k,function(y){
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
      perf_tmp1  <- full_join(beta_quantile, min_cindex, by = "varname")

      perf <-  perf_tmp1 %>%
        filter(include0_yn == TRUE & include_yn == TRUE)

      ### adding back main effect if only interaction terms are selected
      included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

      ## detect main effect terms from interaction
      included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

      ## final selected features
      included_fe <- c(perf$varname, included_inter1) %>% unique()
      #included_fe
      perf2  <- perf_tmp1 %>%
        filter(varname %in% included_fe)
      list(included_fe,perf2)

    })

  })

  ###
  inter_vars_tmp <- lapply(1:length(m),function(x){
    lapply(1:k,function(y){
      fs_out5[[x]][[y]][[1]]
    })
  })

  ## calculate averaged performance over 5 datasets; internal validation
  inter_cindex_tmp <- lapply(1:length(m),function(x){
    lapply(1:k,function(y){
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


  return(list(model_out_all,fs_out5,inter_val_hyper))

}

# inter_perf_f(boot_lst=lung_boot_lst[c(1,2)],
#                          cov_interact = FALSE,
#                          method = "ridge",
#                          k = 2, ## cv fold number
#                          qtl = 0.01,
#                          Rf = 0)

#
# tmp <- HDSI_model_f(boot_df_train = sim_boot_lst[[1]],
#                          # boot_df_test = NA,
#                          cov_interact = TRUE,
#                          method = "ridge")
#
# inter_perf_f(boot_lst=lung_boot_lst,
#              cov_interact = TRUE, ## for simulation data, always true
#              method = "ridge",
#              k = 2, ## cv fold number
#              qtl = 0.01,
#              Rf = 0.1)


inter_perf_f1 <- function(boot_lst=NA,
                         cov_interact = NA,
                         method = "ridge",
                         k = 2, ## cv fold number
                         qtl = NA,
                         Rf = NA){
  library(caret)
  library(dplyr)
  library(survival)
  library(stringr)
  m <- method
  model_out_all <- lapply(1:length(m), function(x) {
    lapply(1:length(boot_lst), function(y) {
      print(y)
      folds <- createFolds(boot_lst[[y]]$status, k = k, list = TRUE, returnTrain = FALSE)



      ## the validation performance
      if(k!=1){
        res <- lapply(1:length(folds),function(z){
          # print(paste0('This is the ',y,'th boot sample ', z,'th fold'))
          res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][-folds[[z]],],
                               cov_interact = cov_interact,
                               # boot_df_test = boot_lst[[y]][folds[[z]],],
                               method = m[[x]])
          res$boot = y
          res$fold = -z

          res
        })
      }
      if(k==1){
        res <-  lapply(1:length(folds),function(z){
          print(paste0('This is the ',y,'th boot sample ', z,'th fold'))
          res <-  HDSI_model_f(boot_df_train = boot_lst[[y]][folds[[z]],],
                               cov_interact = cov_interact,
                               # boot_df_test = boot_lst[[y]][folds[[z]],],
                               method = m[[x]])
          res$boot = y
          res$fold = -z

          res
        })
      }

      res

    })
  })

  print('model_out_all output!! good job!!!')
  fs_out5 <- lapply(1:length(m), function(x) {
    ## bind B boot sample results into a dataframe
    perf_tmp <- model_out_all[[x]] %>% bind_rows()
    lapply(1:k,function(y){
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
      perf_tmp1  <- full_join(beta_quantile, min_cindex, by = "varname")

      perf <-  perf_tmp1 %>%
        filter(include0_yn == TRUE & include_yn == TRUE)

      ### adding back main effect if only interaction terms are selected
      included_inter <- perf$varname[str_detect(perf$varname, ":")] ## detect interaction terms

      ## detect main effect terms from interaction
      included_inter1 <- stringr::str_split(included_inter, ":") %>% unlist()

      ## final selected features
      included_fe <- c(perf$varname, included_inter1) %>% unique()
      #included_fe
      perf2  <- perf_tmp1 %>%
        filter(varname %in% included_fe)
      list(included_fe,perf2)

    })

  })

  ###
  inter_vars_tmp <- lapply(1:length(m),function(x){
    lapply(1:k,function(y){
      fs_out5[[x]][[y]][[1]]
    })
  })

  ## calculate averaged performance over 5 datasets; internal validation
  inter_cindex_tmp <- lapply(1:length(m),function(x){
    lapply(1:k,function(y){
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


  #return(list(model_out_all,fs_out5,inter_val_hyper))
  inter_val_hyper
}

set.seed(20240307)
test <- inter_perf_f1(boot_lst=boot_lst,
                          cov_interact =FALSE, ## for real data
                          method = "ridge",
                          k = 1, ## cv fold number
                          qtl = 0.03,
                          Rf = 2)

test

