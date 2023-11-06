HDSI_binary<-function(df,int=interaction_numb,outvar="y",method=c("glm","forward","lasso","ridge","elastic","Alasso"),effectsize="large",seed=1,q=NA,a,b,selection=c("auc","bic")){
  #read training daset and testdataset
  traindf<-df$train
  testdf<-df$test
  #Generate formula & design matrix
  f<-gen_formula(int=int,outvar=outvar)
  Matrix=model.matrix(f,traindf)[,-1] #model.matrix: create a design matrix
  #Bootstraps
  boots<-mbootsample(k=q,interaction_numb=int, effectsize=effectsize, inputdf=traindf, outvar=outvar,seed_multiplier=seed)
  
  #Fit models with HDSI methods
  methodlist=list(glm = HDSI_regression, lasso = HDSI_lasso, forward =  HDSI_forward,ridge=HDSI_ridge,elastic=HDSI_elastic,Alasso=HDSI_Alasso)
  #op <-pbapply::pboptions(nout=9000) #the maximum number of times the progress bar is updated
  result=lapply(1:length(boots), function(x){
    rows=boots[[x]][[2]]     #samples
    columns=boots[[x]][[1]]  #features
    df=traindf[rows, columns]
    y<-traindf[rows,outvar]
    df<-cbind(df,y)
    #cat("df", " ")
    # Run the model
    return(methodlist[[method]](traindf = df, outvar = outvar, f = f,boot=x))
    #result<-readRDS(file= paste0("output/scenario",scenario,"/modelR",x,".RData"))
    #auc<-readRDS(file= paste0("output/scenario",scenario,"/aucR",x,".RData"))
    #cat("res", " ")
    # if(any(grepl("X1_:X2_", res[[2]]$Variable)) | any(grepl("X2_:X1_", res[[2]]$Variable))){
    #   print(res[[2]])
    # }
    
    #result=res[[2]]
  })
  #pbapply::pboptions(op)
  
  #make coefficients and AUC lists
  #res<-make_result(type="coef",method=method,n=length(boots),s=scenario)
  #auc<-make_result(type="auc",method=method,n=length(boots),s=scenario)
  res<-rbind.fill(fun1(result,1),all=T)
  #saveRDS(res,file= paste0("output/scenario9/res.RData"))
  auc<-rbind.fill(fun1(result,2),all=T)
  #saveRDS(auc,file= paste0("output/scenario9/auc.RData"))
  bic<-rbind.fill(fun1(result,3),all=T)
  #saveRDS(bic,file= paste0("output/scenario9/bic.RData"))
  #feature selection
  feature<-feature_selection(res,auc,bic,a=a,b=b,selection=selection)
  split_var = unlist(strsplit(feature, ":"))
  feature.new<-unique(c(split_var,feature))
  interaction<-sum(grepl(":",feature.new))
  marginal<-length(feature.new)-interaction
  if(length(feature.new)==0) return(list(summary=data.frame(marginal=0,interaction=0,auc=0)))
  
  #fit new models
  f.new=as.formula(paste(outvar,"~", paste(feature.new, collapse = "+")))
  glm.new <- glm(f.new,family = binomial(link = logit), data = traindf,maxit = 200)
  #capture.output(summary(glm.new),file=paste0("output/scenario",scenario,"/results/HDSI",method,".txt"),append=T)
  auc.fin<-auc(testdf[,outvar],as.vector(predict(glm.new,testdf,type="response")))
  fin<-list(summary=data.frame(marginal=marginal,interaction=interaction,auc=auc.fin),features=feature.new)
  #capture.output(fin,file=paste0("output/scenario",scenario,"/results/HDSI",method,".txt"),append=T)
  return(fin)
}