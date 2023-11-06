HDSI_regression=function(traindf = df, outvar = outvar, f = f,boot){
  x=traindf[,names(traindf)!=outvar] 
  y=traindf[,outvar] 
  
  glm <- glm(f,family = binomial(link = logit), data = traindf,maxit = 200)
  #summary(glm)
  coef.glm<-as.data.frame(t(glm$coefficients))
  #auc
  #pred<-predict(glm,newdata = traindf,type='response')
  #pred.new <- prediction(pred,y)
  #auc.perf<-performance(pred.new,measure="auc")
  #capture.output(predict(glm,traindf,type="response"),file=paste0("output/realdata6/results/test.txt"),append=T)
  auc.ord<-auc(y,predict(glm,traindf,type="response"))
  #auc.ord<-auc.perf@y.values
  aucdata.glm<-as.data.frame(t(glm$coefficients))
  aucdata.glm[]<-auc.ord
  bic<-BIC(glm)
  bic.glm<-as.data.frame(t(glm$coefficients))
  bic.glm[]<-bic
  res<-list(coef=coef.glm,auc=aucdata.glm,bic=bic.glm)
  #saveRDS(res, file= paste0("output/scenario",scenario,"/glm",boot,".RData"))
}