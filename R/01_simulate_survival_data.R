install.packages("pacman")
library(pacman)
pacman::p_load(coxed,tidyverse)

simdata <- sim.survdata(N=1000, T=250, xvars=25, censor=.2,num.data.frames = 1)
summary(simdata$data)


# baseline hazard: Weibull

# N = sample size
# lambda = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C

sim_survdat_f <- function(nsample = 1000,
                          varnum = 25, ## total number of covariates
                          dist = 'w',# c("exponential",'weibull','gompertz')
                          alpha = 1, ## shape parameter  v for weibull
                          lambda = 0.01, ## scale parameter
                          beta = c(0.2,0.3,0.4), ## pre-specified effect size
                          seed = 20231105,
                          cor = TRUE
                          ){
  Sigma=matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)

  for(i in 1:varnum){Sigma[i,i]=10}
  # Correlation Settings
  if(cor){
    Sigma[1,2]=3;Sigma[1,3]=3;Sigma[1,4]=6;Sigma[1,5]=6
    Sigma[2,1]=3;Sigma[3,1]=3;Sigma[4,1]=6;Sigma[5,1]=6
    Sigma[2,3]=3;Sigma[2,4]=2;Sigma[2,5]=1
    Sigma[3,2]=3;Sigma[4,2]=2;Sigma[5,2]=1
    Sigma[3,4]=2;Sigma[3,5]=1
    Sigma[4,3]=2;Sigma[5,3]=1
    Sigma[4,5]=1
    Sigma[5,4]=1
  }

  set.seed(seed)
  covar_df <-  data.frame(MASS::mvrnorm(n = nsample, rep(0, varnum), Sigma/10))

  # beta_df <- matrix(rep(beta,each = nsample),ncol = length(beta)) %>% data.frame()
  #
  # modelvar_df <- covar_df %>% select(X1,X3) %>% mutate(X1_X3 = X1*X3)
  # linear_pred <- apply((beta_df*modelvar_df), 1,sum) %>% data.frame()## element-wise product of two data frames
  #
  X1 <- covar_df$X1
  X3 <- covar_df$X3
  Z = beta[1]*X1+beta[2]*X3+beta[3]*X1*X3

  U = runif(nsample,min=0,max=1)
  if(dist == 'w'){
    Time = (- (log(U))/(lambda*exp(Z)))^(1/alpha) ## weibull
    #(- log(v) / (lambda * exp(Z)))^(1 / rho)
  }

  if(dist == 'e'){
    Time =(- (log(U))/( lambda*exp(linear_pred)))
  }

  if(dist == 'g'){
    Time = (1/alpha)*(log(1- ( (alpha*log(U))/(lambda*exp(linear_pred)) ) ))
  }

  censor_time = rexp(n=nsample, rate=0.001)

  event_time <- pmin( as.vector(Time) %>% unlist(), censor_time) %>% data.frame()

  status <- as.numeric(Time <= censor_time)

  survdat <- data.frame(id = 1:nsample,
             event_time = event_time$.,
             status = status,
             covar_values)

  return(survdat)
}




tmpdat <- sim_survdat_f(nsample = 1000,
                          varnum = 25, ## total number of covariates
                          dist = 'w',# c("exponential",'weibull','gompertz')
                          alpha = 1, ## shape parameter  v for weibull
                          lambda = 0.01, ## scale parameter
                          beta = c(0.6,0.3,0.4), ## pre-specified effect size
                          seed = 20231105,
                          cor = TRUE
)


fit <- coxph(Surv(event_time, status) ~ X1+X3+X1*X3, data=survdat)

set.seed(1234)
betaHat <- rep(NA, 1e3)
for(k in 1:1e3)
{
  tmpdat <- sim_survdat_f(nsample = 1000,
                          varnum = 25, ## total number of covariates
                          dist = 'w',# c("exponential",'weibull','gompertz')
                          alpha = 1, ## shape parameter  v for weibull
                          lambda = 1, ## scale parameter
                          beta = c(0.2,0.3,0.4), ## pre-specified effect size
                          seed = 20231105,
                          cor = TRUE)
  fit <- coxph(Surv(event_time, status) ~ X1+X3+X1*X3, data=tmpdat)
  betaHat[k] <- fit$coef[[1]]
}
mean(betaHat)

simulWeib <- function(varnum =25, N, lambda, rho, beta, rateC)
{
  # covariate --> N Bernoulli trials

  Sigma=matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)

  for(i in 1:varnum){Sigma[i,i]=10}
  # Correlation Settings
  if(cor){
    Sigma[1,2]=3;Sigma[1,3]=3;Sigma[1,4]=6;Sigma[1,5]=6
    Sigma[2,1]=3;Sigma[3,1]=3;Sigma[4,1]=6;Sigma[5,1]=6
    Sigma[2,3]=3;Sigma[2,4]=2;Sigma[2,5]=1
    Sigma[3,2]=3;Sigma[4,2]=2;Sigma[5,2]=1
    Sigma[3,4]=2;Sigma[3,5]=1
    Sigma[4,3]=2;Sigma[5,3]=1
    Sigma[4,5]=1
    Sigma[5,4]=1
  }

  set.seed(seed)
  covar_df <-  data.frame(MASS::mvrnorm(n = nsample, rep(0, varnum), Sigma/10))

  # beta_df <- matrix(rep(beta,each = nsample),ncol = length(beta)) %>% data.frame()
  #
  # modelvar_df <- covar_df %>% select(X1,X3) %>% mutate(X1_X3 = X1*X3)
  # linear_pred <- apply((beta_df*modelvar_df), 1,sum) %>% data.frame()## element-wise product of two data frames
  #
  X1 <- covar_df$X1
  X3 <- covar_df$X3
  Z = beta[1]*X1+beta[2]*X3+beta[3]*X1*X3


  # x3 <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
  # x1 <- rnorm(n = N,mean = 0,sd = 1)
  # x2 <- rnorm(n = N,mean = 0,sd = 1)
  # x = x1*x2


  # Z = x1*0.5+x2*0.3
  # Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(Z)))^(1 / rho)

  # censoring times
  C <- rexp(n=N, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(id=1:N,
             time=time,
             status=status,
             covar_df)
}




set.seed(1234)
betaHat <- rep(NA, 1e3)
for(k in 1:1e3)
{
  dat <- simulWeib(varnum = 25,N=100, lambda=0.01, rho=1, beta=c(2.3,0.3,0.4), rateC=0.001)
  fit <- coxph(Surv(time, status) ~ X1+X3+X1*X3, data=dat)
  betaHat[k] <- fit$coef[[2]]
}
mean(betaHat)

