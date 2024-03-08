
#### simulation bootstrap dataset preparation----------

### prepare different nfeature of bootstrap sets
beta_coef <-c(1,1,0.75,-0.75,0.75,1,-1)
source("R/01_simulate_survival_data.R")
source("R/02_bootstrap_datasets.R")
source("R/helper_f.R")

dat <- sim_survdat_f(nsample = 1000,
                     varnum = 25,
                     dist = "g",
                     lambda = 0.01,
                     rho = 1,
                     beta = beta_coef,
                     crate = 0.001, cor = TRUE, seed = 20231106)
library(survival)
fit <- survival::coxph(Surv(time, status) ~ X1 + X2 + X3 + X4+ X5+X1*X2+X3*X4, data = dat)

varnum_vals <- c(25,50,100)
varnum_vals <- c(25)
sim_set <- lapply(varnum_vals,function(varnum){
  sim_survdat_f(nsample = 1000,
                varnum = varnum,
                dist = "g",
                lambda = 0.01,
                rho = 1,
                beta = beta_coef,
                crate = 0.001, cor = TRUE, seed = 20230131)
})

#[1] "event rate is 0.437"
#[1] "event rate is 0.425"
#[1] "event rate is 0.433"


boot_ft1 <- c(10:20)
boot_ft2 <- c(10:20)
boot_ft3 <- c(15:20)
boot_ft1 <- 15
sim_set1 <- lapply(boot_ft1,function(boot_ft){
  boot_f (df = sim_set[[1]],
          nboot = 100,
          boot_ft = boot_ft,
          seed = seed)

})

# sim_set2 <- lapply(boot_ft2,function(boot_ft){
#   boot_f (df = sim_set[[2]],
#           nboot = 100,
#           boot_ft = boot_ft,
#           seed = 20231106)
#
# })
#
# sim_set3 <- lapply(boot_ft3,function(boot_ft){
#   boot_f (df = sim_set[[3]],
#           nboot = 100,
#           boot_ft = boot_ft,
#           seed = 20231106)
#
# })
#
#
# save(sim_set1,file = paste0('data/sim_set/sim_set1_',Sys.Date(),".rdata"))
# save(sim_set2,file = paste0('data/sim_set/sim_set2_',Sys.Date(),".rdata"))
# save(sim_set3,file = paste0('data/sim_set/sim_set3_',Sys.Date(),".rdata"))
