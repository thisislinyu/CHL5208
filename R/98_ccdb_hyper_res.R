load("server_output/sim_set1_hyper_out.rdata")


hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.15, by = 0.01),
                                   Rf = seq(-2.5,2.5,by=0.1))


sim_set1_hyper_dat <- sim_set1_hyper_out %>%
  mutate(nboot_ft = rep(c(10:20),
                        8
                       # nrow(hyperparameter_grid)*2
                        ) %>%
           sort())

sim_set2_hyper_dat <- sim_set2_hyper_out %>%
  mutate(nboot_ft = rep(c(10:20),
                        4*2
                        # nrow(hyperparameter_grid)*2
  ) %>%
    sort())

tmp <- sim_set2_hyper_dat %>%
  select(-qtl,-Rf) %>%group_by(nboot_ft) %>%
  mutate(avg_nboot_perf = mean(avg_cindex)) %>%  ## performance of different nboot,
 select(-avg_cindex) %>%
  unique()

tmp %>%
  ggplot(aes(x = as.factor(nboot_ft), y = avg_nboot_perf, color = method)) +
  geom_line(aes(group = method)) +
  geom_point()+
  theme_bw()+
  xlab('Number of bootstrap features')+
  ylab('Averaged performance (Cindex)')


opt_sim_set2_hyper <- sim_set2_hyper_dat %>% filter(!is.na(avg_cindex)) %>%
  group_by(method) %>%
  mutate(max_cindex = max(avg_cindex)) %>%
  filter(max_cindex==avg_cindex)




opt_sim_hyper <- sim_set1_hyper_dat %>% filter(!is.na(avg_cindex)) %>%
  group_by(method) %>%
  mutate(max_cindex = max(avg_cindex)) %>%
  filter(max_cindex==avg_cindex)

# method avg_cindex   qtl    Rf nboot_ft max_cindex
# <chr>       <dbl> <dbl> <dbl>    <int>      <dbl>
#   1 lasso       0.797  0.09   2.5        8      0.797
# 2 ridge       0.808  0.02   2.4        8      0.808


## generate 10 test samples
res_hyper10 <- list()
beta_coef <-   c(1, 1, 0.75, -0.75, 0.75, 1, -1)

## the best bootstrap features
opt_boot_ft <- 20

for(i in c(1:5)){

sim_ext_dat <- sim_survdat_f(nsample = 500,
                     varnum = 25,
                     dist = "g",
                     lambda = 0.01,
                     rho = 1,
                     beta = beta_coef,
                     crate = 0.001, cor = TRUE, seed = i)



sim_ext_lst <- boot_f(df = sim_ext_dat,
                   nboot = 25,
                   boot_ft = opt_boot_ft,
                   seed = i)

m <- c('lasso','ridge')
sim_ext_out <- lapply(1:length(m), function(x) {
  lapply(1:length(sim_ext_lst), function(y) {
    HDSI_model_f(boot_df = sim_ext_lst[[y]], method = m[[x]],cov_interact = TRUE)
  })
})

model_out_all <- sim_ext_out

 opt_qtl <- c(0.09,0.02)
 opt_Rf <- c(2,2)

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
      Rf_tuned = opt_Rf[x], ## hyperparameter
      include_yn = min_cindex > (miu_min_cindex + Rf_tuned * sd_min_cindex)
    )

  ### compute coef estimate and CI
  beta_quantile <- perf_tmp %>%
    group_by(varname) %>%
    mutate(quantile = opt_qtl[x]) %>% ## quantile is a hyperparameter
    summarise(
      qtl = opt_qtl[x],
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

res_hyper10[[i]] <- selected_dat
}

#hyper10_res <- NULL

#lasso results
for(i in c(1:length(res_hyper10))){

 print(res_hyper10[[i]][[1]][[2]] )  ## 1--ith; 2 -- ridge ；2--features
}

#lridge results
for(i in c(1:length(res_hyper10))){

  print(res_hyper10[[i]][[2]][[2]] )  ## 1--ith; 2 -- ridge ；2--features
}

# [1] "X1:X2"  "X1:X9"  "X4:X14" "X9:X13" "X9:X14" "X1"     "X2"
# [8] "X9"     "X4"     "X14"    "X13"
# [1] "X1"     "X1:X11" "X1:X15" "X1:X16" "X1:X2"  "X1:X22" "X1:X24"
# [8] "X1:X5"  "X1:X6"  "X11"    "X15"    "X16"    "X2"     "X22"
# [15] "X24"    "X5"     "X6"
# [1] "X1:X16" "X1:X17" "X1:X2"  "X1:X23" "X1"     "X16"    "X17"
# [8] "X2"     "X23"
# [1] "X1:X2"  "X1:X25" "X2:X5"  "X3:X5"  "X5:X11" "X5:X17" "X1"
# [8] "X2"     "X25"    "X5"     "X3"     "X11"    "X17"
# [1] "X1:X12" "X1:X2"  "X1:X20" "X1:X22" "X1:X24" "X1:X5"  "X1:X8"
# [8] "X2:X5"  "X1"     "X12"    "X2"     "X20"    "X22"    "X24"
# [15] "X5"     "X8"
# [1] "X1:X11" "X1:X2"  "X1:X20" "X1:X21" "X1:X23" "X1:X5"  "X2:X5"
# [8] "X3:X15" "X3:X5"  "X1"     "X11"    "X2"     "X20"    "X21"
# [15] "X23"    "X5"     "X3"     "X15"
# [1] "X1:X2" "X2:X5" "X1"    "X2"    "X5"
# [1] "X1:X2"  "X1:X5"  "X2:X5"  "X5:X10" "X5:X11" "X5:X17" "X5:X23"
# [8] "X5:X24" "X5:X6"  "X5:X9"  "X1"     "X2"     "X5"     "X10"
# [15] "X11"    "X17"    "X23"    "X24"    "X6"     "X9"
# [1] "X1:X2"  "X2:X21" "X2:X5"  "X1"     "X2"     "X21"    "X5"
# [1] "X1:X18" "X1:X2"  "X1:X23" "X1:X9"  "X2:X5"  "X3:X5"  "X1"
# [8] "X18"    "X2"     "X23"    "X9"     "X5"     "X3"

## [[1]]代表lasso方法；
selected_vars_dat <- data.frame(var_name =selected_dat[[2]][[2]]) %>%
  mutate(new_var = str_replace_all(var_name,":","\\*"),
         index = grepl("age",new_var)) %>%
  filter(index==FALSE)



# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'nboot_ft', 'performance_metric', 'method', etc.

# Install and load ggplot2 if not already installed
# install.packages("ggplot2")
library(ggplot2)

# Example: Plotting performance metrics against 'qtl' for different 'Rf' values, faceted by 'nboot_ft' and colored by 'method'
sim_set1_hyper_dat %>%
  ggplot(aes(x = qtl, y = avg_cindex, color = as.factor(Rf), linetype = method)) +
  geom_line() +
  labs(x = "qtl", y = "Performance Metric", title = "Hyperparameter Tuning Results") +
  facet_wrap(~nboot_ft, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(legend.position = "top")


# Install and load scatterplot3d if not already installed
# install.packages("scatterplot3d")
library(scatterplot3d)

# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'nboot_ft', 'avg_cindex', 'method', etc.

# Convert 'method' to a factor and get numeric codes
method_codes <- as.numeric(factor(results_list$method))

# Create a 3D scatter plot
scatterplot3d(results_list$qtl, results_list$Rf, results_list$avg_cindex,
              color = method_codes, pch = 16, main = "3D Scatter Plot",
              xlab = "qtl", ylab = "Rf", zlab = "Average C-Index")

# Add planes for each method
for (method in levels(results_list$method)) {
  subset_data <- subset(results_list, method == method)
  # Fit a plane using lm
  model <- lm(avg_cindex ~ qtl + Rf, data = subset_data)
  # Extract coefficients
  coef <- coef(model)
  # Add the plane to the scatter plot
  planes(coef[1], coef[2], coef[3], col = "gray", alpha = 0.5)
}

# Add legend
legend("topright", legend = levels(results_list$method), col = 1:2, pch = 16, title = "Method")


# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'nboot_ft', 'performance_metric', 'method', etc.

# Create a 3D scatter plot
method_codes <- as.numeric(factor(results_list$method))
results_list <- sim_set1_hyper_dat
scatterplot3d(results_list$qtl, results_list$Rf, results_list$avg_cindex,
              color = method_codes, pch = 16, main = "3D Scatter Plot",
              xlab = "qtl", ylab = "Rf", zlab = "Performance Metric")

# Add legend
legend("topright", legend = levels(results_list$method), col = 1:2, pch = 16, title = "Method")

############----

# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'avg_cindex', 'method', etc.

# Install and load the 'rgl' package
install.packages("rgl")
library(rgl)

# Create a 3D scatter plot
scatter3d <- plot3d(results_list$qtl, results_list$Rf, results_list$avg_cindex,
                    col = method_codes, type = 's', size = 3,
                    xlab = "qtl", ylab = "Rf", zlab = "Average C-Index",
                    main = "3D Scatter Plot")

# Add planes for each method
for (method in levels(results_list$method)) {
  subset_data <- subset(results_list, method == method)
  # Fit a plane using lm
  model <- lm(avg_cindex ~ qtl + Rf, data = subset_data)
  # Extract coefficients
  coef <- coef(model)
  # Add the plane to the scatter plot
  planes3d(coef[1], coef[2], coef[3], col = "gray", alpha = 0.5, add = TRUE)
}

# Add legend
legend("topright", legend = levels(results_list$method), col = 1:2, pch = 16, title = "Method")

# Interactively rotate the plot
rglwidget(scatter3d)

######

# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'avg_cindex', 'method', etc.

# Install and load the 'plotly' package
# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'avg_cindex', 'method', etc.

# Install and load the required packages
# Assuming 'results_list' is your dataframe containing tuning results
# It should have columns like 'qtl', 'Rf', 'avg_cindex', 'method', etc.

# Install and load the required packages
install.packages("plotly")

# Convert avg_cindex to numeric
results_list$avg_cindex <- as.numeric(results_list$avg_cindex)

# Create 3D scatter plot
plot_3d <- plot_ly(results_list, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                   marker = list(size = 1)) %>%
  layout(scene = list(aspectmode = "cube", aspectratio = list(x = 1, y = 1, z = 1)),
         margin = list(l = 0, r = 0, b = 0, t = 0))

# Display the plot
plot_3d




install.packages("plotly")
library(plotly)

# Assuming you have a results_list dataframe with columns qtl, Rf, and avg_cindex

# Filter the results_list
results_list <- sim %>%
  filter(method == 'lasso' & nboot_ft == 5)

# Create a 3D scatter plot
plot_3d <- plot_ly(results_list, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                   marker = list(size = 1))

# Add a fitted plane
fitted_plane <- lm(avg_cindex ~ qtl + Rf, data = results_list)
plane <- predict(fitted_plane, newdata = results_list)
plot_3d <- add_surface(plot_3d, z = ~plane, colorscale = "Viridis")

# Customize the layout
plot_3d <- layout(plot_3d, scene = list(aspectmode = "cube", aspectratio = list(x = 1, y = 1, z = 1)),
                  margin = list(l = 0, r = 0, b = 0, t = 0))

# Display the plot
plot_3d



# Optionally, you can save the ggplot object and use ggplotly
# gg <- ggplot(results_list, aes(x = qtl, y = avg_cindex, color = method)) +
#   geom_point() +
#   facet_wrap(~method, scales = "free_y") +
#   labs(title = "Tuning Results", x = "qtl", y = "Average C-Index") +
#   theme_minimal()
# plotly::ggplotly(gg)


# Ensure that res_hyper10[[i]][[2]] is not empty
if (length(res_hyper10[[i]][[2]]) == 0) {
  stop("Error: The list of variables is empty.")
}

# Construct formula using paste
formula_str <- paste(res_hyper10[[i]][[2]], collapse = "+")

# Print out the constructed formula
cat("Constructed formula:", formula_str, "\n")

# Create formula using as.formula
formula_obj <- as.formula(paste("Surv(time, status) ~", formula_str))

# Alternatively, you can use reformulate
# formula_obj <- reformulate(response = "Surv(time, status)", termlabels = res_hyper10[[i]][[2]])

# Print out the final formula


# Fit the Cox model
fit <- survival::coxph(formula_obj, data = sim_ext_dat)



