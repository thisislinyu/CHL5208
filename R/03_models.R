library(survival)
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
dat <- sim_survdat_f(nsample = 1000, varnum = 25, dist = "g", lambda = 0.01, rho = 1, beta = c(1.8,0.5,0.4,-0.4,0.45,0.6,-0.6), crate = 0.001, cor = TRUE, seed = 20231106)
fit <- survival::coxph(Surv(time, status) ~ X1 + X2 + X3 + X4+ X5+X1*X2+X3*X4, data = dat)

boot_lst <- boot_f(
  df = dat,
  nboot = 100,
  boot_ft = 5,
  seed = 20231106
)
#boot_df <- boot_lst[[1]]

pacman::p_load(glmnet)

HDSI_model_f <- function(boot_df = NA,
                         method = "lasso") {
  library(pacman)
  pacman::p_load(glmnet)

  ### x input in cv.glmnet
  y <- "survival::Surv(time,status)"
  ## f=stats::as.formula(paste(y," ~ .*.*."))
  f <- stats::as.formula(paste(y, " ~ .*.")) ## pair-wise interaction
  X <- stats::model.matrix(f, boot_df)[, -1]
  ### y input in cv.glmnet
  time <- boot_df$time
  status <- boot_df$status
  Y <- Surv(boot_df$time, boot_df$status)

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



m <- c("lasso", "ridge")

model_out_all <- lapply(1:length(m), function(x) {
  lapply(1:length(boot_lst), function(y) {
    HDSI_model_f(boot_df = boot_lst[[y]], method = m[[x]])
  })
})

Sys.time()
tmp <- lapply(1:length(m), function(x) {
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

  res <- list(perf,included_fe)

  res
})

Sys.time()

tmp21 <- tmp[[2]][[1]] %>% data.frame()
# var_freq <- stringr::str_split(min_cindex$varname,":") %>% unlist() %>% table() %>% data.frame()
