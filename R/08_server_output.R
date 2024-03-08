load('server_output/sim_internal_res.rdata') ## -->internal_res


#


boot_lst <- sim_boot_lst_external
m = c("lasso",'ridge')
m = "ridge"
m = 'lasso'
cov_interact = TRUE

model_out_all <- lapply(1:length(m), function(x) {
  lapply(1:length(boot_lst), function(y) {
    HDSI_model_f(boot_df = boot_lst[[y]],
                 cov_interact = cov_interact,
                 method = m[[x]])
  })
})

sim_hyper_dat <- internal_res

opt_qtl = 0.04
Rf = 2

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
      Rf = Rf, ## hyperparameter
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
         index = grepl("age",new_var)) %>%
  filter(index==FALSE)

select_gene <- selected_vars_dat$new_var

length(select_gene)
