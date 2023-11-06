#' Get the performance of HDSI model on any given dataset
#'
#' The functions takes a dataframe list as an input and evalaute the performance of different HDSI modeling techniques
#' like LASSO, ALASSO, RIDGE, ARIDGE, Forward, and Regression
#'
#' @param model determines the HDSI models that need to be run
#' @param inputdf is the dataframe containing the training and the test datasets
#' @param covariate takes the covariate which needs to be controlled during the feature selection process.
#' This is functional only for forward selection approach.
#' @param outvar is the outcome variable in the dataset.
#' @param para are the hyperparameters used to tweak the model settings.
#' It includes interactions, interaction level, number of latent factors for PLS, intercept term and datatype.
#' @return The performance of the different models alongwith feature selection by the models
#' @export
HDSI_model=function(model="all", inputdf=df , covariate=c(1, "X8_"), outvar="y", seed=1, bootstrap=T, effectsize="large",
                    k=7, min_max= c("min", "ci", "quartile"), cint = 0.95, sd_level=1, model_tech="reg",
                    para=HDSI_para_control(interactions=T, int_term=2, intercept=T, out_type="continuous", perf_metric=c("beta", "mp", "mp_beta"))){

  out_type = para['out_type']
  perf_metric = para['perf_metric']
  #str(inputdf)
  # Create Fakenames
  {
    if(out_type=="survival"){
      fake_name=fakename_gen(datafile = inputdf, train_test=T, outcome_var=c("time", outvar), num_outcome=0)
      fake_df=fake_name$modified_file
    }
    else{
      fake_name=fakename_gen(datafile = inputdf, train_test=T, outcome_var=c(outvar), num_outcome=0)
      #str(fake_name)
      fake_df=fake_name$modified_file
      #str(fake_df)
    }
    #str(fake_df)
    #cat("Fakename", " ")
    traindf=fake_df[[1]]
    testdf=fake_df[[2]]
  }
  #str(fake_df)
  # Select the models that need to run
  {
    if(model == "all"){ model = c("reg", "lasso", "alasso", "ridge", "aridge", "Forward")}
    else{model}

    # Prepare the formula for the models
    f=HDSI::HDSI_formula_gen(other_para=para)
    #print(f)
    methodlist=list(reg = HDSI_Regression, lasso = HDSI_Lasso, alasso = HDSI_Alasso, ridge = HDSI_Ridge, aridge = HDSI_Aridge,
                    Forward =  HDSI_Forward)
    # methodlist=list(reg = m_HDSI_Regression, lasso = m_HDSI_Lasso, alasso = m_HDSI_Alasso, ridge = m_HDSI_Ridge,
    #                 aridge = m_HDSI_Aridge, Forward =  m_HDSI_Forward)
  }
  #str(traindf[1:5,1])
  #cat("Modelname", " ")
  # Prepare the bootstraps
  {
    feature_names=setdiff(names(traindf), union(covariate, c("time", outvar)))
    #print(feature_names)
    boots= HDSI::mbootsample(p = length(feature_names), k=k, rows=nrow(traindf), interaction_numb=2, effectsize=effectsize,
                      feature_name = feature_names, inputdf=traindf, type=out_type, seed_multiplier=seed,
                      bootstrap=bootstrap)
    #print(length(boots))
  }
  # cat("Bootstraps", " ")
  # Run the model

  model_fit = lapply(model, function(x){
    #print(x)

    op <-pbapply::pboptions(nout=9000)
    res=lapply(1:length(boots), function(y){
      rows=boots[[y]][[2]]
      columns=boots[[y]][[1]]
      if(out_type == "survival"){ Columns=union(columns, union(covariate[-1], c("time",outvar)))}
      else{Columns =  union(columns, union(covariate[-1], outvar))}

      #print(Columns)
      #print(names(traindf))

      df=traindf[rows, Columns]
      #cat("df", " ")
      # Run the model
      res=methodlist[[x]](df = df, outvar = outvar, f = f, other_para = para, covariate = covariate)
      #cat("res", " ")
      # if(any(grepl("X1_:X2_", res[[2]]$Variable)) | any(grepl("X2_:X1_", res[[2]]$Variable))){
      #   print(res[[2]])
      # }

      result=res[[2]]
    })
    pbapply::pboptions(op)
    #cat("boots", " ")
    full_res=do.call(rbind, res)
    full_res<<-full_res
    #str(full_res)

    # # Perform the variable selection
    res_summary=HDSI::HDSI_feature_selection(approach = perf_metric, df=full_res, min_max= min_max, cint = cint, sd_level=sd_level)
    #cat("summary", " ")
    #print(res_summary[[1]])
    #str(res_summary)

    # Evaluate the model performance
    raw_list=res_summary$raw_list
    #print(res_summary$sel_feature_list)
    #if(!any(grep("X1_X2", res_summary$sel_feature_list))){print(res_summary$sel_feature_list)}

    perf = HDSI::HDSI_performance(raw_list = raw_list, df = fake_df, outvar = outvar, model_tech = model_tech, covariate = covariate, output=out_type)

    result=list(performance = perf, feature = res_summary$sel_feature_list, feature_number = res_summary$feature_numb,
                technique = x,
                MV_list= res_summary$MV_list, MV_numb = res_summary$MV_numb,
                IV_list = res_summary$IV_list, IV_numb = res_summary$IV_numb)
    performance = data.frame(result$performance, MV = result$MV_numb, IV = result$IV_numb,
                             technique = result$technique, stringsAsFactors = F)
    list(performance, result)
  })
  #cat("Model", " ")
  performance = lapply(model_fit, function(x) x[[1]])
  performance = do.call(rbind, performance)
  #cat("Performance", " ")
  outcome=list(performance = performance, fulldata = model_fit, realname=fake_name$realcolnames, fakename=names(fake_name$modified_file[[1]]))
  #print("Result")
  return(outcome)
}

