# Set working directory
setwd("/scratch/w/wxu/qianzzq/simulation")
#load packages and functions
library(memoise)
library(glmnet)
library(ROCR)
library(Rmisc)
library(GA)
library(mosaic)
library(parallel)

#file folder name
scenario=8
sapply(list.files("HDSI-code",pattern="*.R$"), function(x) source(paste0("HDSI-code/",x)))
#read training dataset and test dataset
sample<-readRDS(file= paste0("output/scenario",scenario,"/sample.RData"))
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
## Ordinary logistic regression ##
##################################

glm <- glm(formula=f,family = "binomial", data = traindf,maxit = 200)
capture.output(summary(glm),file=paste0("output/scenario",scenario,"/results/standardlogistic.txt"),append=T)
interaction<-sum(grepl(":",names(glm$coefficients)[-1]))
marginal<-length(names(glm$coefficients))-interaction-1
if(length(names(glm$coefficients))==1) {
  res<-data.frame(marginal=0,interaction=0,auc=0)
  }else{
  auc.fin<-CalAUC(glm,testdf,testdf[,outvar])
  res<-data.frame(marginal=marginal,interaction=interaction,auc=auc.fin)
}
capture.output(res,file=paste0("output/scenario",scenario,"/results/standardlogistic.txt"),append=T)


##################################
##      Genetic Algorithm       ##
##################################
GA_glm_auc<-ga(type="real-valued",fitness = fitness,method="glm",lower=c(1,0,-2),upper=c(10,1,2),popSize=50,run=100,maxiter=200,parallel=T,pcrossover=0.8,pmutation=0.4)
#GA_lasso_2<-ga(type="real-valued",fitness = fitness,scenario=2,method="lasso",min=c(1,0,-2),max=c(10,1,2),popSize=50,run=3)
capture.output(summary(GA_glm_auc),file=paste0("output/scenario",scenario,"/results/GA_glm_auc.txt"),append=T)
GA_glm_aic<-ga(type="real-valued",fitness = fitness2,method="glm",lower=c(1,0,-2),upper=c(10,1,2),popSize=50,run=100,maxiter=200,parallel=T,pcrossover=0.8,pmutation=0.4)
capture.output(summary(GA_glm_auc),file=paste0("output/scenario",scenario,"/results/GA_glm_aic.txt"),append=T)
##################################
##            HDSI              ##
##################################

HDSI_binary(df=sample,int=2,method="glm",scenario=2,q=5,a=0.40,b=1.70)
HDSI_binary(df=sample,int=2,method="forward",scenario=2,q=5,a=0.40,b=1.77)
HDSI_binary(df=sample,int=2,method="forward",scenario=2,q=6,a=0.66,b=0.35)
HDSI_binary(df=sample,int=2,method="lasso",scenario=2,q=5,a=0.51,b=1.40)

HDSI_binary(df=sample,int=2,method="glm",scenario=3,q=4,a=0.68,b=1.36)
HDSI_binary(df=sample,int=2,method="lasso",scenario=3,q=4,a=0.60,b=1.23)
HDSI_binary(df=sample,int=2,method="forward",scenario=3,q=4,a=0.68,b=1.36)




########################
##    Delong Test     ##
########################
library(pROC)
V1<-testdf[,"y"]
V2<-as.vector(lasso.pred.1se)
roc_1<-roc(V1,V2)
glm.2 <- glm(y~X2+X6+X20+X9+X11+X17+X21+X3+X2:X6+X2:X20+X9:X11+X17:X21+X20:X3+X2:X17,family = "binomial", data = traindf,maxit = 200)
pred=predict(glm.2,newdata = testdf,type='response')
roc_2<-roc(V1,pred)
roc.test(roc_1,roc_2,method="delong")
