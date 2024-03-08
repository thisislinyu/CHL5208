load('output/uni_cox_top200.rdata')
load('data/cox_dat.rdata')
source('R/helper_f.R')

reg_dat <-cox_dat %>%
  select("time", "status", "age_at_diagnosis", "cancer_stage",
         (uni_cox_top200$var_name))

reg_dat50 <- reg_dat[,c(1:104)]
reg_dat <- reg_dat50
set.seed(20240123)
trainIndex <- createDataPartition(reg_dat$status,
                                  p = 0.8, list = FALSE)

trainingSet <- reg_dat[trainIndex,]
testSet <- reg_dat[-trainIndex,]

#save(trainingSet,file = 'output/trainingSet.rdata')
#save(testSet,file = 'output/testSet.rdata')

lung_boot_lst100 <- boot_f1(
  df = trainingSet,
  nboot = 100,
  boot_ft = 15,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)

lung_boot_lst1000 <- boot_f1(
  df = trainingSet,
  nboot = 1000,
  boot_ft = 15,
  gene_index = 5,
  cov_vars = c('age_at_diagnosis','cancer_stage'),
  seed = 20231106
)


set.seed(20240307)
test100 <- inter_perf_f1(boot_lst=lung_boot_lst100,
                      cov_interact =FALSE, ## for real data
                      method = "ridge",
                      k = 1, ## cv fold number
                      qtl = 0.05,
                      Rf = 2)


