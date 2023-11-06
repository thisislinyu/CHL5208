# HDSI
High Dimensional Selection with Interactions (HDSI) Algorithm perform featue selection. It selects both the marginal and interaction terms. This repository contains various codes to help in feature selection.

In order to use HDSI, it needs to be installed in R. 
```R
# install.packages("githubinstall") # In case githubinstall is installed
# library(githubinstall)
# githubinstall("HDSI")
library(HDSI)
```
Once the HDSI is installed. The first step is to generate a dataset. The `dataset()` function generate the artificial dataset. It needs four inputs.<br />
`varnum` needs to know the number of input features in the dataset. In current example, a dataset with 25 input features is created.<br /> `setting` needs to know if correlation among the features is needed or not. The function has a pre-defined correlation structure which will be used to generate the sample data for the inpur features. Default is "Correlation".<br />`var` defines the type of target features to be added in the dataset. The function has three options: `Mar` Target features contain three marginal variables and one interaction term, `No_Mar` Target features contain one marginal variables and one interaction term, `No_Var` No target feature are added into the model.<br />
`seed` is an integer value to create replicable randomness.<br />In case, any dataset is already available, this step could be skipped.
```R
# Generate the dataset
df = dataset(varnum = 25, setting = "Correlation", var = "Mar", seed=2)
# Target Features are "X1", "X2", "X3", "X1_X2". "X1_X2" is the interaction term.
```
`df` is a list containing training data and test data. Now, `df` is used for feature selection using HDSI.
```R
sm_para = HDSI_para_control(interactions=T, # Should model consider interaction terms. Currently, it is fixed to TRUE. 
                                            # Plan is to make it flexible in future updates.
                            int_term=2, # Level of interactions a model should consider. Currently, it is fixed to 2.
                                        # Plan is to make it flexible in future updates.
                            intercept=T, # Should model consider intercept term during statistical model preparation. 
                                         # Currently, it is fixed to TRUE. Plan is to make it flexible in future updates.
                            out_type="continuous", # Should model consider "continuous" outcome or "survival" outcome. 
                                                   # Currently, it is fixed to "continuous".
                                                   # Plan is to make it flexible in future updates.
                            perf_metric="mp_beta") # The preformance metric of a statistical model. 
                                                   # Currently, it is fixed at "mp_beta" to consider both the model and 
                                                   # feature performance.
                                                   # Plan is to incoporate only model performance ("mp") option and 
                                                   # only feature performance ("beta") option. 

res = HDSI_model(model="lasso", # To select the base statistical model. Options are c("lasso", "alasso", "reg") 
                 inputdf=df,  # Input the data dataset containing input features and outcome variable. 
                 covariate=c(1), # Fixes the covariate in the model (Control). 
                                 # The selected covariate interaction terms are not considered and is present in all 
                                 # statistical models.
                 outvar="y", # Name of the outcome variable in the dataset. Currently, outcome variable name must be "y". 
                             # Plan is to make it flexible in future updates.
                 seed=1, # an integer value to create replicable randomness.
                 bootstrap=T, # Fixed. Will be made flexible in future iteration
                 effectsize= "large", # The minimum number of times a feature should be sampled. "large" means 13,
                                      # "medium" means 32 and "small" means 200.
                 k=5, # Hyperparameter q. Ideally, should be obtained from optimisation.
                 min_max= "min", # Fixed. It determine the metric on which model performance cut-off value is applied.
                 cint = 0.95, # Hyperparameter Qi. Ideally, should be obtained from optimisation.
                 sd_level= 1.0, # Hyperparameter Rf. Ideally, should be obtained from optimisation.
                 model_tech = "reg", # The model used to test the predictive performance of the selected features. 
                                     # "reg" is for linear regression and "ridge" penalised L2- regularisation. 
                 para= sm_para) # Some additional parameters.
```
`res` will provide the lot of results. Some of the main results are as follows:
1) use `res$performance` to check for model performance in training and testing data
2) use `res$fulldata[[2]]$feature` to check the selected features
3) use `res$realname` to get the real names of the features
4) use `res$fakename` to get the names used by the model. They are the names which are displayed in `res$fulldata[[2]]$feature`
