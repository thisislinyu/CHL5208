HDSI_binary<-function(df,int=interaction_numb,outvar="y",method=c("glm","lasso","ridge","elastic","Alasso"),effectsize="large",seed=1,q=NA,a,b,selection=c("auc","bic")){
  #read training dataset and testdataset
  traindf<-df$train
  testdf<-df$test
  traindf[,outvar]=as.factor(traindf[,outvar])
  testdf[,outvar]=as.factor(testdf[,outvar])
  #Generate formula & design matrix
  f<-gen_formula(int=int,outvar=outvar)
  Matrix=model.matrix(f,traindf)[,-1] #model.matrix: create a design matrix
  #Bootstraps
  boots<-mbootsample(k=q,interaction_numb=int, effectsize=effectsize, inputdf=traindf, outvar=outvar,seed_multiplier=seed)
  
  #Fit models with HDSI methods
  methodlist=list(glm = HDSI_regression, lasso = HDSI_lasso, ridge=HDSI_ridge,elastic=HDSI_elastic,Alasso=HDSI_Alasso,svm=HDSI_svm)
  #op <-pbapply::pboptions(nout=9000) #the maximum number of times the progress bar is updated
  result=lapply(1:length(boots), function(x){
    rows=boots[[x]][[2]]     #samples
    columns=boots[[x]][[1]]  #features
    df=traindf[rows, columns]
    y<-traindf[rows,outvar]
    df<-cbind(df,y)
    # Run the model
    return(methodlist[[method]](traindf = df, outvar = outvar, f = f,boot=x))
  })
  #pbapply::pboptions(op)
  
  #Summarize the final results 
  res<-rbind.fill(fun1(result,1),all=T)
  auc<-rbind.fill(fun1(result,2),all=T)
  if (!is.null(dim(res[3]))) bic<-rbind.fill(fun1(result,3),all=T) else bic=list()
  #Feature selection
  feature<-feature_selection(res,auc,bic,a=a,b=b,selection=selection)
  if (method=="svm") feature=gsub(pattern="[.]",replacement=":",feature)
  split_var = unlist(strsplit(feature, ":"))
  feature.new<-unique(c(split_var,feature))
  interaction<-sum(grepl(":",feature.new))
  marginal<-length(feature.new)-interaction
  if(length(feature.new)==0) return(list(summary=data.frame(marginal=0,interaction=0,auc=0)))
  
  #Fit new models
  f.new=as.formula(paste(outvar,"~", paste(feature.new, collapse = "+")))
  glm.new <- glm(f.new,family = binomial(link = logit), data = traindf,maxit = 200)
  rocdata<-roc(testdf[,outvar],as.vector(predict(glm.new,testdf,type="response")))
  saveRDS(rocdata,file=paste0("output/scenario",scenario,"/roc_",method,".RData"))
  auc.fin<-auc(testdf[,outvar],as.vector(predict(glm.new,testdf,type="response")))
  fin<-list(summary=data.frame(marginal=marginal,interaction=interaction,auc=auc.fin),features=feature.new)
  return(fin)
}