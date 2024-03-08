# Set the directory where your .rdata files are located
directory <- "server_output/sim_set1_ccdb_output"

# List all .rdata files in the directory
rdata_files <- list.files(directory, pattern = "\\.rdata$", full.names = TRUE)

# Read each .rdata file
data_list <- lapply(rdata_files, function(file) {
  load(file)
  # Return the loaded object(s)
  return(get(load(file)))
})

# Optionally, you can combine all loaded objects into a single list or data frame
combined_data <- do.call(rbind, data_list)  # Example for combining into a data frame


## 2*4*11*16+ 2*2*11
combined_data1 <- combined_data %>%
  mutate(nboot_ft = rep(c(10:20),
                        nrow(combined_data)/11
                        # nrow(hyperparameter_grid)*2
  ) %>%
    sort())

combined_data1 %>%
  ggplot(aes(x = qtl, y = avg_cindex, color = as.factor(Rf), linetype = method)) +
  geom_line() +
  labs(x = "qtl", y = "Performance Metric", title = "Hyperparameter Tuning Results") +
  facet_wrap(~nboot_ft, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(legend.position = "top")


plot_3d <- plot_ly(combined_data1 %>% filter(method=='lasso'), x = ~nboot_ft, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
                   marker = list(size = 1))

plot_ly(combined_data1 %>% filter(method=='ridge'), x = ~nboot_ft, y = ~Rf, z = ~avg_cindex, type = "scatter3d",
        marker = list(size = 1))


tmp <- combined_data1 %>%
  group_by(nboot_ft) %>%
  mutate(avg_nboot_perf = mean(na.omit(avg_cindex))) %>%  ## performance of different nboot,
  select(-avg_cindex) %>%
  unique()

tmp %>%
  ggplot(aes(x = as.factor(nboot_ft), y = avg_nboot_perf, color = method)) +
  geom_line(aes(group = method)) +
  geom_point()+
  theme_bw()+
  xlab('Number of bootstrap features')+
  ylab('Averaged performance (Cindex)')


