p_est(p = 25,k=6,rows=1000, interaction_numb=2, effectsize= c("large"))


model="lasso"# To select the base statistical model. Options are c("lasso", "alasso", "reg")
inputdf=df # Input the data dataset containing input features and outcome variable.
covariate=c(1) # Fixes the covariate in the model (Control).
# The selected covariate interaction terms are not considered and is present in all
# statistical models.
outvar="y"# Name of the outcome variable in the dataset. Currently, outcome variable name must be "y".
# Plan is to make it flexible in future updates.
seed=1# an integer value to create replicable randomness.
bootstrap=T # Fixed. Will be made flexible in future iteration
effectsize= "large" # The minimum number of times a feature should be sampled. "large" means 13,
# "medium" means 32 and "small" means 200.
k=5 # Hyperparameter q. Ideally, should be obtained from optimisation.
min_max= "min" # Fixed. It determine the metric on which model performance cut-off value is applied.
cint = 0.95 # Hyperparameter Qi. Ideally, should be obtained from optimisation.
sd_level= 1.0 # Hyperparameter Rf. Ideally, should be obtained from optimisation.
model_tech = "reg"

out_type = para['out_type']
perf_metric = para['perf_metric']
