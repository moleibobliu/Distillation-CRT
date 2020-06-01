# This scripts run some examples using the dCRT function for different needs.

############### Require necessary packages ##########

source('functions.R')
library(glmnet)
library(MASS)
library(knockoff)
library(matrixStats)
library(randomForest)


######################## Use for Gaussian X ###################


###########################################################
##################### List of Choices #####################
###########################################################


# X: observed covariates; Y: observed response.
# FDR: desirable FDR level.
# d.interaction: 'T' for using dICRT; F for using 'd0CRT'
# k: dimension of the control variables in distilled data, only useful when "d.interaction = T"

# model = c('Gaussian_lasso', 'Binomial_lasso', 'RF'), for linear model, logistics model (binary Y) and random forest (dICRT (RF)).
#
# RF.num.trees = c(100, 30): Number of trees used for distilling (y, Z) and constructing importance measure.
#
# Gen_X: NULL (default) for Gaussian X with mean 0 and Sigma_X specifying its covariance structure.
#        Can be a function Gen_X(X, indx, num) where X represents the observed covariates, indx is 
#        the indx for the variable needs to be resampled and num is the number of samples one want.
#        It should output a matrix with each column represent a conditional samples with size n.
#        Please see our example
#
# Sigma_X: covariance matrix of X. When Gen_X=NULL, i.e. X is gaussian, it is used for generating
#          X and x-distillation. When Gen_X is not NULL, it is used for solely x-distillation. 
#          Default as NULL but needs to be specified by the user when k1 > 0 or Gen_X=NULL.
#
# candidate_set: indices of the covariates to be tested. Default is all.
#
# mean_X_mat: matrix with j-th column being the conditional mean of covariate j.

# CDR: 'Consist w model' (default) for using the same choice as model; 'No' for not doing CDR;
#      can be a function CDR(Y, X) outputing a subset of 1:p.

# M = 50000: Number of resampling if needed
# MC_free: whether to use resampling-free d0CRT or dICRT when X is not gaussian, i.e. Gen_X is not NULL.
# central = T/F: whether to centralize data marginally before running dCRT.
# lambda.seq: list of candidate tuning parameters for LASSO, default = NULL.
# eps_dI: parameter to control how precise the p-values of resampling-free dICRT be.
# X_quantile_mat: matrix with j-th column being the conditional quantile of covariate j,
#                 useful when X is not gaussian and one use MC_free = T.
# X_cond_var_mat: matrix with j-th column being the conditional variance of covariate j,
#                 useful when X is not gaussian and one use MC_free = T.
# return.pvalues: If or not to return the estimated p-values.

######################## Generate data (linear Y~X and Gaussian X) ###################

n <- 500
p <- 300
s <- 30

Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = 'linear', Y_dist = 'Gaussian',
                          Sigma = 'AR', r = 0.5, magn = 0.3,
                          X_dist = 'gauss')
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma
support <- which(Gen_data$beta_true != 0)


############### d0CRT (resampling-free) ###############

d0CRT_results <- dCRT(X, Y, Sigma_X = Sigma_true,
                      FDR = 0.1, model = 'Gaussian_lasso') 

############### d0CRT (natural form) ###############

d0CRT_results <- dCRT(X, Y, Sigma_X = Sigma_true,
                      FDR = 0.1, MC_free = F, model = 'Gaussian_lasso') 


############### dICRT (resampling-free) ###############

dICRT_results <- dCRT(X, Y, Sigma_X = Sigma_true, 
                      d.interaction = T, k = as.integer(2 * log(p)),
                      FDR = 0.1, model = 'Gaussian_lasso') 


############### dICRT (natural form) ###############

dICRT_results <- dCRT(X, Y, Sigma_X = Sigma_true, 
                      d.interaction = T, k = as.integer(2 * log(p)),
                      FDR = 0.1, MC_free = F, model = 'Gaussian_lasso') 

############### dICRT (RF) ###############

dkCRT_results <- dCRT(X, Y, Sigma_X = Sigma_true, candidate_set = c(1),
                      d.interaction = T, k = as.integer(2 * log(p)),
                      CDR = 'No', model = 'RF', M = 1000) 


##################### dCRT for non gaussian X ###################

## Take X_ij's are independently gamma (normalized) as an example.
## Generate independent X with shape = 3 and rate = 0.5

Gen_data <- Generate_data(n, p, s = 30,
                          model = 'linear', Y_dist = 'Gaussian', magn = 0.3,
                          X_dist = 'indept_gamma') 
X <- Gen_data$X
Y <- Gen_data$Y


############### d0CRT (resampling-free) ###############

## In this case, one needs to specify X_quantile_mat and X_cond_var_mat. For example:

gamma_quantile <- function(X){
  return(pgamma(6 + X * sqrt(12), shape = 3, rate = 0.5))
}
X_quantile_mat <- gamma_quantile(X)
X_con_var_mat <- matrix(1, n, p)

d0CRT_results <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                      Gen_X = Gen_X_example, CDR = 'Consist w model',
                      X_quantile_mat = X_quantile_mat, X_cond_var_mat = X_con_var_mat,
                      model = 'Gaussian_lasso', k = 0, MC_free = T)

############### d0CRT (natural form) ###############

## One needs to specify the Gen_X function to generate X, for example:

Gen_X_example <- function(X, indx, num = 1){
  n <- length(X[,1])
  X_sample <- matrix(rgamma(n * num, shape = 3, rate = 0.5), n, num)
  return((X_sample - 6) / sqrt(12))
}

d0CRT_results <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                      Gen_X = Gen_X_example, CDR = 'Consist w model',
                      model = 'Gaussian_lasso', k = 0, MC_free = F)


