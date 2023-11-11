#' Simulate a survival data set
#'
#' @param nsample total number of observations
#' @param varnum  total number of covariates
#' @param dist distribution of baseline hazard function
#' @param lambda scale parameter of the distribution
#' @param rho shape parameter of the distribution
#' @param beta pre-specified effect size of covariates
#' @param crate rate parameter of the exponential distribution
#' @param cor logical; whether consider the correlation between covariates
#' @param seed
#'
#' @return a dataframe; time, status, covariates
#' @export
#' @examples
#'library(survival)
#'dat <- sim_survdat_f(nsample =1000, varnum = 25,dist='g',lambda=0.01, rho=1, beta=c(2.3,0.3,0.4), crate=0.001,cor=TRUE,seed=20231106)
#'fit <- survival::coxph(Surv(time, status) ~ X1+X3+X1*X3, data=dat)

sim_survdat_f <- function(nsample = 100,
                          varnum =25,
                          dist = 'w',
                          lambda = 0.01,
                          rho = 1,
                          beta=NA,
                          crate = 0.001,
                          cor=TRUE,
                          seed = 20231106,
                          ...)
{
  library(pacman)
  pacman::p_load(tidyverse)

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

  beta_df <- matrix(rep(beta,each = nsample),ncol = length(beta)) %>% data.frame()

  modelvar_df <- covar_df %>% select(X1,X2,X3,X4,X5) %>% mutate(X1_X2 = X1*X2,X3_X4=X3*X4)
  Z <- apply((beta_df*modelvar_df), 1,sum) %>% unlist()## element-wise product of two data frames

  # X1 <- covar_df$X1
  # X3 <- covar_df$X3
  # Z = beta[1]*X1+beta[2]*X3+beta[3]*X1*X3

  U <- runif(n=nsample)

  if(dist == 'e'){
    time <- (- log(U) / (lambda * exp(Z)))
  }

  if(dist == 'w'){
    time <- (- log(U) / (lambda * exp(Z)))^(1 / rho)
  }

  if(dist =='g'){
    nd_term = ( (rho * log(U)) / (lambda * exp(Z)))
    time <-(1/rho) * log(1- nd_term)
  }
  # censoring times
  C <- rexp(n=nsample, rate=crate)

  # follow-up times and event indicators
  status <- as.numeric(time <= C)
  time <- pmin(time, C) ## over-write time variable, compare between time and censor time


  # data set
  data.frame(time=time,
             status=status,
             covar_df)
}


#
# set.seed(1234)
# betaHat <- rep(NA, 1e3)
# for(k in 1:1e3)
# {
#   dat <- sim_survdat_f(varnum = 25,nsample =1000, lambda=0.01, rho=1, beta=c(2.3,0.3,0.4), crate=0.001,cor=TRUE,seed=20231106)
#   fit <- coxph(Surv(time, status) ~ X1+X3+X1*X3, data=dat)
#   betaHat[k] <- fit$coef[[2]]
# }
# mean(betaHat)
#dat <- sim_survdat_f(nsample =1000, varnum = 25,dist='w',lambda=0.01, rho=1, beta=c(0.1,0.2,0.3), crate=0.001,cor=TRUE,seed=20231106)
