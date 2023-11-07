#' Generating bootstrap sample datasets
#'
#' @param df working dataset,simulated or real data collected
#' @param nboot number of bootstrap datasets
#' @param boot_ft number of features in each bootstrap dataset
#' @param seed random seed
#'
#' @return a list, size = nboot, each element is a dataframe
#' @export
#'
#' @examples
#'library(survival)
#'dat <- sim_survdat_f(nsample =1000, varnum = 25,dist='g',lambda=0.01, rho=1, beta=c(2.3,0.3,0.4), crate=0.001,cor=TRUE,seed=20231106)
#'fit <- survival::coxph(Surv(time, status) ~ X1+X3+X1*X3, data=dat)
#' boot_lst <- boot_f(df = ,nboot = 100,boot_ft = 5,seed = 20231106)

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
