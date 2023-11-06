#Estimate the number of boostrap samples B
B_est=function(p,k=NA,rows=NA, interaction_numb=2, effectsize= c("large", "medium", "small")){
  #print(c(p,k))
  denominator=choose(p,interaction_numb)
  numerator=choose(k,interaction_numb) # where, 2 is the order of interaction eg of 2 order interaction is X1_X2
  #print(c(numerator,denominator))
  init_prob=numerator/denominator
  
  # Minimum bootstraps which will give minimum 13 occurrences (consider large effect) of an interaction variable with 99% confidence
  if(effectsize=="large"){prefer_occurence = 13}
  else if (effectsize == "medium"){prefer_occurence = 32}
  else if(effectsize == "small"){prefer_occurence = 200}
  else {prefer_occurence = effectsize}
  # print(c(prefer_occurence,init_prob))
  min_boot= ceiling(prefer_occurence/init_prob) # where 32 is minimum number of occurence desired for an interaction variable during bootstrapping.
  #print(min_boot)
  max_boot= 8*min_boot
  #print(max_boot)
  f_opt=function(x,y=init_prob,max_x=max_boot,prefer=prefer_occurence){
    Val=qbinom(0.05, floor(x), y) #P(X>Val)>=95%
    value=(prefer-Val)^2+(x/max_x)
    return(value)
  }
  bootvalue = optimize(f_opt,interval = c(min_boot,max_boot))
  # print(ceiling(bootvalue$minimum))
  return(ceiling(bootvalue$minimum))
}

#Bootstrapping
bootsample=function(k=NA,interaction_numb=2, effectsize="large", inputdf, outvar,seed_multiplier=1){
  p<-dim(inputdf)[2]-length(outvar)
  rows<-dim(inputdf)[1]
  # Find k if not provided by user
  if(is.na(k)==T){
    optimal_k=function(x, r=rows){
      int=choose(x,interaction_numb)
      value=abs((int+x)-r)/r
      return(value)
    }
    k=stats::optimize(optimal_k,interval = c(interaction_numb,p))$minimum
    k=floor(k)
  }
  
  #print(p)
  # Find the number of sample that need to be created
  boots=B_est(p=p,k=k, rows=rows, interaction_numb=interaction_numb, effectsize=effectsize)
  #cat(c("Bootstraps:",boots))
  
  ## Create Variable list for each bootstrap
  res=bootsample_binary(boots=boots, p=p,rows=rows, k=k, seed_multiplier = seed_multiplier,inputdf=inputdf,outvar=outvar)
  
  return(res)
}

bootsample_binary=function(boots=boots, p, rows, k, seed_multiplier=1, inputdf,outvar){
  if(!require(pbapply)) install.packages("pbapply")
  library(pbapply)
  ## Create Variable list for each bootstrap
  samples=1:rows
  op <- pbapply::pboptions(type = "timer") #adds progress bar to vectorized R functions
  res=pblapply(1:boots, function(x) {
    set.seed(x*seed_multiplier)
    ### Create variable combinations
    features <- sample(p, k, replace = FALSE) # replaced samples with features/3
    features<- sort(features)
    ## Create Sample list for each bootstrap
    CAT_1 <- which(inputdf[,outvar]==1)
    CAT_0 <- which(inputdf[,outvar]==0)
    S_1<-sample(CAT_1,length(CAT_1),replace=T)
    S_0<-sample(CAT_0,length(CAT_0),replace=T)
    samples <- c(S_1,S_0)
    list(features, samples)
  })
  pbapply::pboptions(op)
  return(res)
}
mbootsample = memoise::memoise(bootsample)