## 对每一个gene进行uni variate data analysis adjusting for age and cancer stage
## UVA: cox apply to real-data set-------
load('data/cox_dat.rdata')
library(caret)

## define univariate cox function
uni_cox_f <- function(data=NA, ## outcome variable names must be time and status
                      gen_var_index = 5, ## from which the gene data is started
                      cov_vars = 'age_at_diagnosis+cancer_stage',
                      ...
){
  surv_object <- Surv(time = data$time,
                      event = data$status
  )
  vars = colnames(data[-c(1:gen_var_index-1)])
  results <- lapply(vars, function(var) {
    #print(vars)
    var = make.names(var)

    tryCatch({
      formula <- as.formula(paste("surv_object ~", paste(cov_vars,var,sep = "+")))
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

t1 = Sys.time()
uni_cox_res <- uni_cox_f(data = cox_dat,
                         gen_var_index = 5,
                         cov_vars = 'age_at_diagnosis+cancer_stage')
t2 = Sys.time()
t2-t1
uni_cox_res1 <- uni_cox_res %>% bind_rows()

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("multtest")
library(rio)
#install.packages("multtest")
library(multtest)


adjusted_p_values <- mt.rawp2adjp(rawp = uni_cox_res1$Pr...z.., proc= "BH")

adjp_dat <- adjusted_p_values$adjp %>% data.frame()

#############
uni_cox_res2 <- uni_cox_res1 %>%
  mutate(
    var_name = rownames(.),
    index = grepl("age", var_name)
  ) %>%
  filter(index==FALSE) %>%
  select(-index) %>%
  mutate(index = row_number())

adjusted_p_values <- mt.rawp2adjp(rawp = uni_cox_res2$Pr...z..,
                                  proc= c("Bonferroni","BH"))

adjp_dat <- adjusted_p_values$adjp %>% data.frame() %>%
  mutate(index = adjusted_p_values$index)

uni_cox_res2 <- left_join(uni_cox_res2,adjp_dat,by='index')

uni_cox_top200 <- uni_cox_res2%>%
  arrange(BH) %>%
  filter(BH <= BH[200])

save(uni_cox_res2,file="output/uni_cox_res2.rdata")
save(uni_cox_top200,file="output/uni_cox_top200.rdata")

