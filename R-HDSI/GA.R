#Set working directory
setwd("/gpfs/fs1/home/w/wxu/qianzzq/HDSI2022")
#Load packages and functions
library(memoise)
library(glmnet)
library(ROCR)
library(Rmisc)
library(GA)
library(mosaic)
library(parallel)
library(doParallel)
library(caret)
library(msgps)
library(pROC)
library(plyr)
#Read parameter values from the bash script
args=commandArgs(T)
scenario=args[1]
folder=args[2]
method=args[3]
k=as.numeric(args[4])
sapply(list.files("HDSI-code",pattern="*.R$"), function(x) source(paste0("HDSI-code/",x)))
#read training dataset and test dataset
if (folder=="realdata") {
  sample<-readRDS(file= paste0("realdata",scenario,"/sample1.RData"))
  traindf<-sample$train
} else {
  traindf<-readRDS(file= paste0("scenario",scenario,"/sample.RData"))
}

traindf<-sample$train
testdf<-sample$test

#Generate formula & design matrix
int_term=2
x<-c(rep(".*",int_term-1),".")
x<-Reduce('paste0',x)
f=as.formula(paste("y ~",x))
Matrix=model.matrix(f,traindf)[,-1] #model.matrix: create a design matrix
outvar="y"

##################################
##      Genetic Algorithm       ##
##################################
GA_auc<-ga(type="real-valued",fitness = fitness,train=traindf,effectsize="large",rule="auc",method=method,lower=c(2,0.9,1),upper=c(k,1,2.5),popSize=20,run=10,maxiter=200,parallel=T,pcrossover=0.8,pmutation=0.4,nfold=5)
summary(GA_auc)
capture.output(summary=summary(GA_auc),file=paste0(folder,scenario,"/results/GA_",method,"_auc.txt"),append=T)