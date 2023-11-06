#Calculate AUC 
CalAUC =function(model,testdf,real){
  pred=predict(model,newdata = testdf,type='response')
  rocr.pred=prediction(pred,real)
  rocr.perf=performance(rocr.pred,'auc')
  as.numeric(rocr.perf@y.values)
}

#Feature selection
feature_selection=function(coef,auc,bic,a,b,selection=c("auc","bic")){
  ind<-apply(coef,2,function(x){
    x<-x[!is.na(x)]
    mean_x<-mean(x)
    sd_x<-sd(x)
    q_lower<-mean_x-qnorm(1-(1-a)/2)*sd_x
    q_upper<-mean_x+qnorm(1-(1-a)/2)*sd_x
    #q_lower<-quantile(x,(1-a)/2)
    #q_upper<-quantile(x,1-(1-a)/2)
    a<-ifelse((q_lower<0 & q_upper>0)| is.na(sd_x),0,1)})
  coef_name<-names(ind)[ind!=0]
  #auc threshold,b=
  #max.auc<-apply(auc,2,function(x){
  #  x<-x[!is.na(x)]
  #  if(length(x)==0) return(NA)
  #  else
  #  return(max(x,na.rm=T))})
  #max.auc<-max.auc[-1]
  #max.auc<-max.auc[!is.na(max.auc)]
  #mean_auc<-mean(max.auc)
  #sd_auc<-sd(max.auc)
  #auc_upper<-mean_auc+1.64*sd_auc
  #y<-names(max.auc)[max.auc>=quantile(max.auc,0.9)]
  if (selection=="auc") {
    min.auc<-apply(auc,2,function(x){
      x<-x[!is.na(x)]
      if(length(x)==0) return(NA)
      else
        return(min(x,na.rm=T))})
    min.auc<-min.auc[-1]
    min.auc<-min.auc[!is.na(min.auc)]
    mean.min_auc<-mean(min.auc)
    sd.min_auc<-sd(min.auc)
    auc_lower<-mean.min_auc+b*sd.min_auc
    auc_name<-names(min.auc)[min.auc>auc_lower]
    comb<-intersect(coef_name,auc_name)
    return(comb)
  }else{
    max.bic<-apply(bic,2,function(x){
      x<-x[!is.na(x)]
      if(length(x)==0) return(NA)
      else
        return(max(x,na.rm=T))})
    max.bic<-max.bic[-1]
    max.bic<-max.bic[!is.na(max.bic)]
    mean.max.bic<-mean(max.bic)
    sd.max.bic<-sd(max.bic)
    bic_upper<-mean.max.bic+b*sd.max.bic
    bic_name<-names(max.bic)[max.bic<bic_upper]
    comb<-intersect(coef_name,bic_name)
    return(comb)
  }
}

#Generate formula
gen_formula=function(int,outvar){
  x<-c(rep(".*",int-1),".")
  x<-Reduce('paste0',x)
  f<-as.formula(paste(outvar," ~ ",x))
  return(f)
}

#Create result list
make_result = function(type=c("coef","auc"),method=method,n=length(boots),s=scenario){
  if(!require(plyr)) install.packages("plyr")
  library(plyr)
  res<-lapply(1:n,function(x) readRDS(file= paste0("output/scenario",s,"/",method,x,".RData"))[[type]])
  res<-rbind.fill(res,all=T)
  return(res)
}

#Make result list - when no RData is saved
fun1 <- function(lst, n){
  sapply(lst, `[`, n)
}

#BIC
bic_cal<-function(model){
  tLL <-  - deviance(model)
  k <- model$df
  n <- model$nobs
  #AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  #AICc
  
  BIC<-log(n)*k - tLL
  return(BIC)
}

#sampling data for simulation scenarios
gendata<-function(ta){
  index_0<-which(ta$y==0)
  index_1<-which(ta$y==1)
  index_0<-sample(index_0,150,replace=F)
  index_1<-sample(index_1,150,replace=F)
  index<-c(index_0,index_1)
  return(ta[index,])
}

#$summary
#marginal interaction       auc
#        4           2 0.7381333

#features
#"X1"    "X6"    "X2"    "X4"    "X1:X6" "X2:X4"