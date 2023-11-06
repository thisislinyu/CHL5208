#Regression for HDSI-BO
HDSI_regression=function(traindf = df, outvar = outvar, f = f,boot){
  x=traindf[,names(traindf)!=outvar] 
  y=traindf[,outvar] 
  glm <- glm(f,family = binomial(link = logit), data = traindf,maxit = 200)
  coef.glm<-as.data.frame(t(glm$coefficients))
  auc.ord<-auc(y,predict(glm,traindf,type="response"))
  aucdata.glm<-as.data.frame(t(glm$coefficients))
  aucdata.glm[]<-auc.ord
  bic<-BIC(glm)
  bic.glm<-as.data.frame(t(glm$coefficients))
  bic.glm[]<-bic
  res<-list(coef=coef.glm,auc=aucdata.glm,bic=bic.glm)
}
#LASSO for HDSI-BO
HDSI_lasso=function(traindf=df,outvar=outvar,f=f,boot){
  Matrix=stats::model.matrix(f,traindf)[,-1]
  set.seed(1)
  cvfit=cv.glmnet(Matrix,traindf[,outvar],family="binomial",alpha=1,type.measure="class") 
  lambda.1se=cvfit$lambda.1se
  lambda.min=cvfit$lambda.min
  model=glmnet(Matrix, traindf[,outvar], lambda = lambda.1se, alpha=1, standardize=F, family="binomial")
  coef=as.data.frame(t(as.matrix(coef(model,s=lambda.1se))))
  if(sum(coef!=0)<=1){ #test if the model with lambda.1se is too simple (only no parameters)
    model = glmnet::glmnet(Matrix, traindf[,outvar], lambda = lambda.min, alpha=1, standardize=F, family="binomial")
    coef=as.data.frame(t(as.matrix(coef(model,s=lambda.min))))
  }
  lasso.pred=predict(model,newx=Matrix,type="response")   
  pred<- prediction(lasso.pred,traindf[,outvar])
  auc.perf<-performance(pred,measure="auc")
  auc.ord<-unlist(auc.perf@y.values)
  aucdata.lasso<-coef
  aucdata.lasso[aucdata.lasso!=0]<-unlist(auc.ord)
  aucdata.lasso[aucdata.lasso==0]<-NA
  coef[coef==0]<-NA
  bic<-bic_cal(model)
  bic.lasso<-coef
  bic.lasso[bic.lasso!=0]<-bic
  bic.lasso[bic.lasso==0]<-NA
  res<-list(coef=coef,auc=aucdata.lasso,bic=bic.lasso)
  #res<-list(coef=coef,auc=aucdata.lasso)
  #saveRDS(res, file= paste0("output/scenario",scenario,"/lasso",boot,".RData"))
}
#Elastic net for HDSI-BO
HDSI_elastic=function(traindf=df,outvar=outvar,f=f,boot){
set.seed(1)
  Matrix=stats::model.matrix(f,traindf)[,-1]
  elastic<-train(f,data=traindf,method="glmnet",trControl=trainControl("cv",number=5),tuneLength=10)
  alpha<-elastic$bestTune[1]
  lambda<-elastic$bestTune[2]
  model=glmnet(Matrix, traindf[,outvar], lambda = lambda, alpha=alpha, standardize=F, family="binomial")
  coef=as.data.frame(t(as.matrix(coef(model,s=lambda))))
  coef[coef==0]<-NA
  elastic.pred=predict(model,newx=Matrix,type="response")   
  pred<- prediction(elastic.pred,traindf[,outvar])
  auc.perf<-performance(pred,measure="auc")
  auc.ord<-unlist(auc.perf@y.values)
  aucdata.elastic<-coef
  aucdata.elastic[!is.na(aucdata.elastic)]<-unlist(auc.ord)
  bic<-bic_cal(model)
  bic.elastic<-coef
  bic.elastic[!is.na(coef)]<-bic
  res<-list(coef=coef,auc=aucdata.elastic,bic=bic.elastic)
}
#Ridge regression for HDSI-BO
HDSI_ridge=function(traindf=df,outvar=outvar,f=f,boot){
set.seed(1)
  Matrix=stats::model.matrix(f,traindf)[,-1]
  cvfit=cv.glmnet(Matrix,traindf[,outvar],family="binomial",alpha=0,type.measure="class") 
  lambda.1se=cvfit$lambda.1se
  lambda.min=cvfit$lambda.min
  model=glmnet(Matrix, traindf[,outvar], lambda = lambda.1se, alpha=0, standardize=F, family="binomial")
  coef=as.data.frame(t(as.matrix(coef(model,s=lambda.1se))))
  if(sum(coef!=0)<=1){ #test if the model with lambda.1se is too simple (only no parameters)
    model = glmnet::glmnet(Matrix, traindf[,outvar], lambda = lambda.min, alpha=0, standardize=F, family="binomial")
    coef=as.data.frame(t(as.matrix(coef(model,s=lambda.min))))
  }
  lasso.pred=predict(model,newx=Matrix,type="response")   
  pred<- prediction(lasso.pred,traindf[,outvar])
  auc.perf<-performance(pred,measure="auc")
  auc.ord<-unlist(auc.perf@y.values)
  aucdata.lasso<-coef
  aucdata.lasso[aucdata.lasso!=0]<-unlist(auc.ord)
  aucdata.lasso[aucdata.lasso==0]<-NA
  coef[coef==0]<-NA
  bic<-bic_cal(model)
  bic.lasso<-coef
  bic.lasso[!is.na(bic.lasso)]<-bic
  res<-list(coef=coef,auc=aucdata.lasso,bic=bic.lasso)
  #res<-list(coef=coef,auc=aucdata.lasso)
  #saveRDS(res, file= paste0("output/scenario",scenario,"/lasso",boot,".RData"))
}
#Adaptive lasso for HDSI-BO
HDSI_Alasso=function(traindf=df,outvar=outvar,f=f,boot){
set.seed(1)
  Matrix=stats::model.matrix(f,traindf)[,-1]
  Alasso<-cv.glmnet(Matrix,traindf[,outvar],family="binomial",alpha=0,type.measure="class")
  lambda.min<-Alasso$lambda.min
  best_ridge_coef<-as.numeric(coef(Alasso,s=lambda.min))[-1]
  alasso_cv<-cv.glmnet(Matrix,traindf[,outvar],family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef),type.measure="class")
  alasso<-glmnet(Matrix,traindf[,outvar],family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef),lambda=alasso_cv$lambda.1se)
  Coef.alasso=as.data.frame(t(as.matrix(coef(alasso_cv,s=alasso_cv$lambda.1se))))
  if(sum(Coef.alasso!=0)<=1){ #test if the model with lambda.1se is too simple (only no parameters)
    alasso<-glmnet(Matrix,traindf[,outvar],family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef),lambda=alasso_cv$lambda.min)
    Coef.alasso=as.data.frame(t(as.matrix(coef(alasso_cv,s=alasso_cv$lambda.min))))
  }
  Coef.alasso[Coef.alasso==0]<-NA
  alasso.pred=predict(alasso,newx=Matrix,type="response")  
  auc<-auc(traindf[,outvar],as.vector(alasso.pred))
  #rownames(Coef.alasso)[Coef.alasso!=0]
  aucdata.alasso<-Coef.alasso
  aucdata.alasso[!is.na(aucdata.alasso)]<-auc
  bic<-bic_cal(alasso)
  bic.alasso<-Coef.alasso
  bic.alasso[!is.na(bic.alasso)]<-bic
  res<-list(coef=Coef.alasso,auc=aucdata.alasso,bic=bic.alasso)
}

HDSI_svm=function(traindf=df,outvar=outvar,f=f,boot){
  fit.svm<-svm(formula=f,data=traindf,type='C',kernel="linear")
  #extract auc
  auc<-auc(roc(traindf[,outvar],as.numeric(fit.svm$fitted)))
  #extract the weights and constant from the SVM model:
  w <- t(fit.svm$coefs) %*% fit.svm$SV
  #b <- -1 * fit_svm$rho; #(sometimes called w0)
  w<-as.data.frame(w)
  aucdata.svm<-w
  aucdata.svm[]<-auc
  res<-list(weight=w,auc=aucdata.svm)
}