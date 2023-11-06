#' Generate the artificial dataset for testing
#'
#' @param varnum creates teh number of variable needed in the model determines the HDSI models that need to be run
#' @param setting defines to create crorelation or no correlation between the input variables
#' @param var defines the type of relationship between the outcome and input variables
#' @param seed is a random seed generator
#' @return The performance of the different models alongwith feature selection by the models
#' @export
# Create simulated linear dataset
dataset=function(varnum, setting="No_Correlation", var=c("Mar", "No_Mar", "No_Var"), seed=2){
  Sigma=matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)
  for(i in 1:varnum){Sigma[i,i]=10}
  # Correlation Settings
  if(setting=="Correlation"){
    Sigma[1,2]=3;Sigma[1,3]=3;Sigma[1,4]=6;Sigma[1,5]=6
    Sigma[2,1]=3;Sigma[3,1]=3;Sigma[4,1]=6;Sigma[5,1]=6
    Sigma[2,3]=3;Sigma[2,4]=2;Sigma[2,5]=1
    Sigma[3,2]=3;Sigma[4,2]=2;Sigma[5,2]=1
    Sigma[3,4]=2;Sigma[3,5]=1
    Sigma[4,3]=2;Sigma[5,3]=1
    Sigma[4,5]=1
    Sigma[5,4]=1
  }
  #print(Sigma)
  set.seed(seed)
  ta=data.frame(MASS::mvrnorm(n = 1000, rep(0, varnum), Sigma/10))
  c=1
  if(var=="Mar"){beta_a=1; beta_b=1}
  else if(var=="No_Mar"){beta_a=0; beta_b=1}
  else{beta_a=0; beta_b=0}
  b1=0.2*beta_a
  b2=0.3*beta_a
  b3=0.4*beta_b
  b4=0.3*beta_b
  variablelist=list()
  for(i in 1:varnum){
    variablelist[[i]]=gsub(" ", "",paste("X",i))
    ta[,i]=mosaic::zscore(ta[,i])
  }
  variablelist=unlist(variablelist)
  ta$y= c + (b1*ta$X1) + (b2*ta$X2) + (b3*ta$X3) + (b4*ta$X1*ta$X2) + rnorm(n=1000, mean=0, sd=0.25)
  index=sample(1:nrow(ta), nrow(ta)/2, replace = F)
  traindf=ta[index,]
  validationdf=ta[-index,]
  return(list(train=traindf,test=validationdf))
}


dataset
