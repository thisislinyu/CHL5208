library(survival)
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
dat <- sim_survdat_f(nsample = 1000, varnum = 25, dist = "g", lambda = 0.01, rho = 1, beta = c(2.3, 0.3, 0.4), crate = 0.001, cor = TRUE, seed = 20231106)
fit <- survival::coxph(Surv(time, status) ~ X1 + X3 + X1 * X3, data = dat)
boot_lst <- boot_f(
  df = dat,
  nboot = 100,
  boot_ft = 5,
  seed = 20231106
)
boot_df = boot_lst[[1]]

pacman::p_load(glmnet)

HDSI_model_f <- function(boot_df = NA,
                         method = "lasso") {
  library(pacman)
  pacman::p_load(glmnet)


  y="survival::Surv(time,status)"
  ## f=stats::as.formula(paste(y," ~ .*.*."))
  f=stats::as.formula(paste(y," ~ .*.")) ## pair-wise interaction
  X =stats::model.matrix(f,boot_df)[,-1]

  time <- boot_df$time
  status <- boot_df$status
  Y=Surv(boot_df$time,boot_df$status)

  if(method=='lasso'){
    cv_fit <- glmnet::cv.glmnet(x = X ,y = Y, alpha=0,
                                 standardize=F, family="cox",nfolds = 3)
    lambda.1se=cv_fit$lambda.1se
    lambda.min=cv_fit$lambda.min

    model_fit = glmnet::glmnet(x = X ,y = Y, lambda = lambda.1se, alpha=0, standardize=F, family="cox")
    Coef=as.matrix(coef(model_fit,s=lambda.min))

    risk_scores <- predict(model_fit,s = lambda.min,newx = X, type='link')
    cindex <- concordance.index(Surv(time, status), risk_scores)

    # y1 = cbind(time = time, status = status)
    #
    # apply(risk_scores, 2, glmnet::Cindex, y=y1)

  }
  if(method == 'ridge'){
    cv_fit <- glmnet::cv.glmnet(x = X ,y = Y, alpha=1,
                                 standardize=F, family="cox",nfolds = 3)

    lambda.1se=cv_fit$lambda.1se
    lambda.min=cv_fit$lambda.min

    model_fit = glmnet::glmnet(x = X ,y = Y, lambda = lambda.1se, alpha=1, standardize=F, family="cox")
  }

  return(model_fit)

}

cv_fit <- glmnet::cv.glmnet(x = X ,y = Y, alpha=1,
                            standardize=F, family="cox",nfolds = 3)

lambda.1se=cv_fit$lambda.1se
lambda.min=cv_fit$lambda.min

model_fit = glmnet::glmnet(x = X ,y = boot_df$time, lambda = lambda.1se, alpha=1, standardize=F, family="gaussian")
risk_scores <- predict(model_fit,newx = X, type='link')

auc.ord<-unlist(auc.perf@y.values)

m = c('lasso','ridge')

tmp_res12 <- lapply(1:length(m),function(x){

  lapply(1:length(boot_lst),function(y){

    HDSI_model_f(boot_df = boot_lst[[y]],method = m[[x]])

  })


})


tmp13 <-  lapply(1:length(boot_lst),function(y){

  HDSI_model_f(boot_df = boot_lst[[y]],method = 'lasso')

})
















