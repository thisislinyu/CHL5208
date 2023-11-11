#fitness function
fitness = function(x,nfold=3,train=traindf,interaction_numb=2,outvar="y",method="glm",effectsize="large",seed=1,rule=c("auc","bic")){
  q<-ceiling(x[1])
  datasize<-nrow(train)
  index<-CVgroup(K=nfold,datasize=datasize,seed=seed)
  #cl<-makeCluster(39,type="FORK")
  res<-sapply(1:nfold,K_HDSI,data=train,index=index,interaction_numb=interaction_numb,method=method,effectsize=effectsize,
              seed=seed,q=q,a=x[2],b=x[3],selection=rule,outvar=outvar)
  #stopCluster(cl)
  mean_auc<-mean(res)
  return(mean_auc)
}

#Split K groups
CVgroup <- function(K,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:K,ceiling(datasize/K))[1:datasize]
  temp <- sample(n,datasize)
  x <- 1:K
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}

#K-fold Cross-validation
K_HDSI = function(x,data=train,index=index,interaction_numb=interaction_numb,method=method,effectsize=effectsize,seed=seed,q=q,a=a,b=b,selection,outvar=outvar){
  traindf<-data[-index[[x]],]
  testdf<-data[index[[x]],]
  df<-list(train=traindf,test=testdf)
  return(HDSI_binary(df=df,int=interaction_numb,outvar=outvar,method=method,effectsize=effectsize,seed=seed,q=q,a=a,b=b,selection=selection)$summary$auc)
}
