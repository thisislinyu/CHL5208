load("server_output/sim_set1_hyper_out.rdata")


hyperparameter_grid <- expand.grid(qtl = seq(0.01, 0.15, by = 0.01),
                                   Rf = seq(-2.5,2.5,by=0.1))


sim_set1_hyper_dat <- sim_set1_hyper_out %>%
  mutate(nboot_ft = rep(c(5:10),nrow(hyperparameter_grid)*2) %>% sort())

set1 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 5)

# Create a 3D scatter plot
plot_3d_set1 <- plot_ly(set1, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                   marker = list(size = 5))

set2 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 6)

plot_3d_set2 <- plot_ly(set2, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))

set3 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 7)

plot_3d_set3 <- plot_ly(set3, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))


set4 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 8)

plot_3d_set4 <- plot_ly(set4, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))

set5 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 9)

plot_3d_set5 <- plot_ly(set5, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))

set6 <- sim_set1_hyper_dat %>%
  filter(method == 'lasso' & nboot_ft == 10)

plot_3d_set6 <- plot_ly(set6, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))


subplot <- subplot(plot_3d_set1, plot_3d_set2, plot_3d_set3,nrows = 1)
# Create a 3D scatter plot
plot_3d_set1 <- plot_ly(set1, x = ~qtl, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                        marker = list(size = 5))

# Add a fitted plane


# Customize the layout
plot_3d <- layout(plot_3d, scene = list(aspectmode = "cube", aspectratio = list(x = 1, y = 1, z = 1)),
                  margin = list(l = 0, r = 0, b = 0, t = 0))

# Display the plot
plot_3d
