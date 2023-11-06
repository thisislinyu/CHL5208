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
