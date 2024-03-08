# Install and load required packages
install.packages(c("survival", "foreach", "doParallel"))
library(survival)
library(foreach)
library(doParallel)
library(dplyr)

load('data/work_dat1.rdata')

# Set the number of cores
num_cores <- detectCores()

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Assuming your data is in a data frame called 'your_data' with a column 'time' for survival times and 'status' for censoring
# and columns V1, V2, ..., V20000 for your independent variables

# Create a survival object

work_dat1 <- work_dat1 %>%
  mutate(OS.time_month = OS.time/30) %>%
  filter(!is.na(age_at_diagnosis)& !is.na(ajcc_pathologic_stage))

surv_object <- Surv(time = cox_uni_dat$OS.time_month, event = cox_uni_dat$OS)

# Create a data frame with independent variables
# ncol(cox_uni_dat)
gene_vars <- work_dat1[, 26:ncol(work_dat1)] %>%
  apply(.,2,as.numeric)
colnames(gene_vars) <- rlang::set_names(paste0('X_',
                                               make.names(colnames(gene_vars))))# Assuming your independent variables start from the third column


outcome_vars <- work_dat1 %>% select(OS.time_month,OS) %>%
  rlang::set_names("time","status")

cox_dat <- cbind(outcome_vars,gene_vars)

t1 = Sys.time()
uni_cox_res <- uni_cox_f(data = cox_dat,vars = colnames(cox_dat)[-c(1:2)])
t2 = Sys.time()
t2-t1
uni_cox_f <- function(data=NA, ## outcome variable names must be time and status
                      vars = NA,
                      ...
){
 surv_object <- Surv(time = data$time, event = data$status)
 results <- lapply(vars, function(var) {
   #print(vars)
   var = make.names(var)
   # formula <- as.formula(paste("surv_object ~", var))
   # cox_model <- coxph(formula, data = data)
   # result_summary <- summary(cox_model)
   # result_summary$coefficients %>% data.frame() %>% bind_rows()

   tryCatch({
     formula <- as.formula(paste("surv_object ~", var))
     cox_model <- coxph(formula, data = data)
     result_summary <- summary(cox_model)
     result_summary$coefficients %>% data.frame() %>% bind_rows()
   }, error = function(err) {
     # 在发生错误时返回当前结果
     message(paste("Error for variable:", var, " - ", conditionMessage(err)))
     return(data.frame(Variable = var, Error = conditionMessage(err)))
   })
 })

 return(results)
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("multtest")
uni_cox_res1 <- uni_cox_res %>% bind_rows()
library(rio)
install.packages("multtest")
library(multtest)
uni_cox_res2 <- uni_cox_res1 %>%
  mutate(
    var_name = rownames(.),
    index = grepl("age", var_name)
  ) %>%
  filter(index==FALSE) %>%
  select(-index) %>%
  mutate(index = row_number())



uni_cox_res2$Pr...z.. %>% na.omit() %>% hist()

(uni_cox_res2$Pr...z..<=0.05) %>% table()

adjusted_p_values <- mt.rawp2adjp(rawp = uni_cox_res2$Pr...z..,
                                  proc= c("Bonferroni","BH"))

adjp_dat <- adjusted_p_values$adjp %>% data.frame() %>%
  mutate(index = adjusted_p_values$index)

uni_cox_res2 <- left_join(uni_cox_res2,adjp_dat,by='index')
save(uni_cox_res2,file="output/uni_cox_res2.rdata")

data_long <- tidyr::gather(adjp_dat, key = "Variable", value = "Value") %>%
  arrange((Variable))

# Create faceted histogram using facet_grid
uni_p_dat <- ggplot(data_long, aes(x = Value)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "white", alpha = 0.7) +
  facet_grid(Variable ~ ., scales = "free") +
  labs(title = "",
       x = "P-value",
       y = "Frequency")+
  theme_bw()

ggsave(uni_p_dat,file = 'output/uni_p_dat.png',dpi=300)

ggplot(adjp_dat, aes(x = rawp)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~BH, scales = "free") +
  labs(title = "",
       x = "Variable 1",
       y = "Frequency")



adjp_dat$BH %>% hist()
adjp_dat$rawp %>% hist()

table(adjp_dat<0.05)

reg_vars <- uni_cox_res2 %>% arrange(BH) %>%
  mutate(p_order =row_number()) %>%
  filter(p_order<=50 | (p_order>100 & p_order<200 ) | p_order> nrow(.)-50)


load('data/cox_dat.rdata')

reg_dat <-cox_dat %>%
  select("time", "status", "age_at_diagnosis", "cancer_stage",
         (reg_vars$var_name))


