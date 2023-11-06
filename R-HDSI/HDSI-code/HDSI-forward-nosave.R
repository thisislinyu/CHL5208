HDSI_forward=function(traindf = df, outvar = outvar, f = f,boot){
  x=traindf[,names(traindf)!=outvar] 
  y=traindf[,outvar] 
  
  glm <- glm(f,family = binomial(link = logit), data = traindf,maxit = 200)
  #Forward Selection
  initial<-glm(y~1,family=binomial(link=logit),data=traindf,maxit=200)
  logit.forward <- step(initial,direction = "forward",scope=list(upper=glm,lower=initial))
  #summary(logit.step)
  coef.forward<-as.data.frame(t(logit.forward$coefficients))
  #auc.forward
  pred<-predict(logit.forward,newdata = traindf,type='response')
  pred.new <- prediction(pred,y)
  auc.perf<-performance(pred.new,measure="auc")
  auc.ord<-auc.perf@y.values
  aucdata.for<-as.data.frame(t(logit.forward$coefficients))
  aucdata.for[]<-unlist(auc.ord)
  bic<-BIC(logit.forward)
  bic.for<-as.data.frame(t(logit.forward$coefficients))
  bic.for[]<-bic
  res<-list(coef=coef.forward,auc=aucdata.for,bic=bic.for)
  #saveRDS(res, file= paste0("output/scenario",scenario,"/forward",boot,".RData"))
}