
### We summariz the setup for generating data and running the methods of the simulation studies
### in our paper. We recommand using multiple computing clusters for parallel and replicate each setting
### for many times (says 300), and summarizing the paralleled results. 

## Variable to save in all the cases for summary on power, FDR and computation time: 
## results, time_lst


library(CVglasso)
library(glmnet)
library(MASS)
library(knockoff)
library(matrixStats)
library(CompQuadForm)
library(ncvreg)
library(randomForest)

source('functions.R')



FDR_rate <- 0.1

##########################################################################
                ###### Moderate size simulations ######
##########################################################################


# liner model:

Y_dist = 'Gaussian'
magn <- 0.075 # the magnitude of the signals, in our implementation: magn in 0.075 * c(1:5)
model = 'gaussian'
model.dCRT = 'Gaussian_lasso'

# logistic model:

Y_dist = 'binom'
magn <- 0.5 # magn in 0.5 * c(1:5)
model = 'binomial'
model.dCRT = 'Binomial_lasso'

n <- 300
p <- 300
s <- 30
r <- 0.5

## Generate data: 

Gen_data <- Generate_data(n, p, s, intercept = 0, model = 'linear', Y_dist = Y_dist,
                          Sigma = 'AR', r = r, magn = magn, X_dist = 'gauss', support = 'equal')

# support = 'equal' for equally spaced support
# support = 'first' for adjacent support


X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma
support <- which(Gen_data$beta_true != 0)

results <- vector('list', 8)
time_lst <- c()

## Run the approaches:

# d0CRT:

ptm <- proc.time()[1]
results[[1]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model.dCRT) 
time_lst <- c(time_lst, proc.time()[1] - ptm) # Record the consumption time

# knockoffs:

ptm <- proc.time()[1]
knockoffs = function(X) create.gaussian(X, rep(0, p), Sigma_true, method = "sdp",
                                        diag_s = create.solve_sdp(Sigma_true))
foo = stat.glmnet_coefdiff
k_stat = function(X, X_k, y) foo(X, X_k, y, family = model)
Knockoff_Result <- knockoff.filter(X, Y, statistic = k_stat, knockoffs = knockoffs,
                                   fdr = FDR_rate, offset = 1)
selection_set_KO <- Knockoff_Result$selected
a <- list()
a[['select_set']] <- as.vector(selection_set_KO)
results[[2]] <- a
time_lst <- c(time_lst, proc.time()[1] - ptm)

# HRT: 

ptm <- proc.time()[1]
results[[3]] <- HRT(Y, X, Sigma = Sigma_true, FDR = FDR_rate, N = 50000, 
                    model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)

# DML: 

ptm <- proc.time()[1]
results[[4]] <- DML(X, Y, Sigma_true, FDR = FDR_rate, K = 8, model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)

# dICRT: 

ptm <- proc.time()[1]
results[[5]] <- dCRT(X, Y, Sigma_X = Sigma_true, 
                     d.interaction = T, k = as.integer(2 * log(p)),
                     FDR = 0.1, MC_free = T, model = model.dCRT) 
time_lst <- c(time_lst, proc.time()[1] - ptm)


## oCRT (LASSO):
ptm <- proc.time()[1]
results[[6]] <- CRT_sMC(Y, X, Sigma_X = Sigma_true, m = 15000, 
                        type = 'LASSO', FDR = FDR_rate, model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)

## oCRT (Elastic net):

ptm <- proc.time()[1]
results[[7]] <- CRT_sMC(Y, X, Sigma_X = Sigma_true, m = 15000, 
                        type = 'Elasnet', FDR = FDR_rate, model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)

## oCRT (AdaLASSO):

ptm <- proc.time()[1]
results[[8]] <- CRT_sMC(Y, X, Sigma_X = Sigma_true, m = 15000, 
                        type = 'AdaLASSO', FDR = FDR_rate, model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)




##########################################################################
                ###### Large size simulations ######
##########################################################################

n <- 800
p <- 800
s <- 50
r <- 0.5
FDR_rate <- 0.1

# AR correlation design:

sig_type = 'AR' 

# Constant correlation design:

sig_type = 'Cons'

# generation of Y:

# logistic model:

Y_dist = 'binom'
model = 'Binomial_lasso'
model_hrt = 'binomial'
model_type = 'linear' 
magn <- 0.5

# Poisson model:

Y_dist = 'Poisson'
model = 'Gaussian_lasso'
model_hrt = 'gaussian'
model_type = 'exp'
magn <- 0.175


# Polynomial model:

Y_dist = 'Gaussian'
model = 'Gaussian_lasso'
model_hrt = 'gaussian'
model_type = 'poly' 
magn <- 0.105

# linear model:

Y_dist = 'Gaussian'
model = 'Gaussian_lasso'
model_hrt = 'gaussian'
model_type = 'linear' 
magn <- 0.175 


Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = model_type, Y_dist = Y_dist,
                          Sigma = sig_type, r = r, magn = magn,
                          X_dist = 'gauss', support = 'equal')

# support = 'equal' for equally spaced support
# support = 'first' for adjacent support

X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma
support <- which(Gen_data$beta_true != 0)

## Run the approaches: d0CRT, dICRT, DML, knockoffs and HRT.
## Their implemetation are the same as in "Moderate size simulations".



##########################################################################
               ###### Settings with interaction ######
##########################################################################

n <- 800
p <- 800
s <- 50
r <- 0.5
FDR_rate <- 0.1

### parametric model
# linear model 

model <- 'Gaussian_lasso'
model.hrt <- 'gaussian'
Y_dist <- 'Gaussian'
magn <- 0.15

# logistic model

model <- 'Binomial_lasso'
model.hrt <- 'binomial'
Y_dist <- 'binom'
magn_lst <- 0.3

# Generate data:
Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = 'non_linear_single', Y_dist = Y_dist,
                          Sigma = 'AR', r = r, magn = magn,
                          prop = 0.1, support = 'random')
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma
results <- vector('list', 4)

## Run the approaches (only test for the first variable):

# d0CRT:

results[[1]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', MC_free = T) 

# dICRT:

results[[2]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', d.interaction = T,
                     k = as.integer(2 * log(p)), MC_free = T) 


### Random forest

n <- 800
p <- 800
r <- 0.5
magn <- 0.18
  
# Generate data:
Gen_data <- Gen_forest(n, p, s = 5, r = r, intercept = 0, 
                       Sigma = 'AR', X_dist = 'gauss', prop = 0, 
                       magn = magn, Y_dist = 'Gaussian')
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma


# d0CRT (lasso):

results[[1]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', MC_free = T) 

# dICRT (lasso):

results[[2]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', d.interaction = T,
                     k = as.integer(2 * log(p)), MC_free = T) 

# dICRT (random forest):

results[[3]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = 'RF',
                     candidate_set = 1, CDR = 'No', d.interaction = T,
                     RF.num.trees = c(100, 30), k = as.integer(2 * log(p)),
                     M = 1000, MC_free = F) 



##########################################################################
                 ###### Robustness evaluation ######
##########################################################################

# In-sample estimate of X's model

n <- 800
p <- 800
s <- 50
r <- 0.5

# linear model:

Y_dist = 'Gaussian'
model <- 'gaussian'
model.dCRT <- 'Gaussian_lasso'
magn <- 0.175 

# logistic model:

Y_dist = 'binom'
model <- 'binomial'
model.dCRT <- 'Binomial_lasso'
magn <- 0.5

support <- 'first' # support <- 'equal'
Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = 'linear', Y_dist = Y_dist,
                          Sigma = 'AR', r = r, magn = magn, support = support)
X <- Gen_data$X
Y <- Gen_data$Y

# Consider three approaches to estimate the covariance structure:

## 1. Ledoit–Wolf (optimal shrinkage) estimator: 

Sigma_est <- linshrink_cov(X, k = 1, normalize = T)

## 2. Graphic lasso estimator: 

Sigma_est <- glasso_cov(X)

## 3. Nodewise lasso for conditional mean estimation 
# (just set the input Sigma_est as NULL and dCRT() will use Nodewise lasso by default):

Sigma_est <- NULL

results <- vector('list', 5)
time_lst <- c()

## Run the algorithm with fitted Sigma_est:

## d0CRT

ptm <- proc.time()[1]
results[[1]] <-  dCRT(X, Y, Sigma_X = Sigma_est, FDR = 0.1, model = model.dCRT) 
time_lst <- c(time_lst, proc.time()[1] - ptm)

## knockoff

if (is.null(Sigma_est)){
  
}else{
  ptm <- proc.time()[1]
  knockoffs = function(X) create.gaussian(X, rep(0, p), Sigma_est, method = "sdp",
                                          diag_s = create.solve_sdp(Sigma_est))
  foo = stat.glmnet_coefdiff
  k_stat = function(X, X_k, y) foo(X, X_k, y, family = model)
  Knockoff_Result <- knockoff.filter(X, Y, statistic = k_stat, knockoffs = knockoffs,
                                     fdr = FDR_rate, offset = 1)
  selection_set_KO <- Knockoff_Result$selected
  a <- list()
  a[['select_set']] <- as.vector(selection_set_KO)
  results[[2]] <- a
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
}

## HRT

ptm <- proc.time()[1]
results[[3]] <- HRT(Y, X, Sigma = Sigma_est, FDR = FDR_rate, N = 50000, 
                    model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)


### DML

ptm <- proc.time()[1]
results[[4]] <- DML(X, Y, Sigma_est, FDR = FDR_rate, K = 8, model = model)
time_lst <- c(time_lst, proc.time()[1] - ptm)


## dCRT

ptm <- proc.time()[1]
results[[5]] <- dCRT(X, Y, Sigma_X = Sigma_est, 
                     d.interaction = T, k = as.integer(2 * log(p)),
                     FDR = 0.1, MC_free_dI = T, model = model.dCRT) 
time_lst <- c(time_lst, proc.time()[1] - ptm)



# Gaussian approximation to non-gaussian X (known first and second moments)

# linear model:

Y_dist = 'Gaussian'
model <- 'gaussian'
model.dCRT <- 'Gaussian_lasso'
#magn <- 0.1
magn <- 0

# logistic model:

Y_dist = 'binom'
model <- 'binomial'
model.dCRT <- 'Binomial_lasso'
#magn <- 0.2
magn <- 0

n <- 800
p <- 800
s <- 50
r <- 0.5
level <- 0.05
para_lst <- c(0.5, 2, 8, 32, 64)
para_x = para_lst[3]
power_lst <- vector('list', 2)
FDR_rate <- 0.1
k <- 2 * log(p)

Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = 'linear', Y_dist = Y_dist,
                          Sigma = 'AR', r = 0.5, magn = magn,
                          X_dist = 'poisson_single', para_x = para_x)
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma


## Implementation of the algorithms are the same as in 

# d0CRT:

results[[1]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', MC_free = T) 

# dICRT:

results[[2]] <- dCRT(X, Y, Sigma_X = Sigma_true, FDR = 0.1, model = model,
                     candidate_set = 1, CDR = 'No', d.interaction = T,
                     k = as.integer(2 * log(p)), MC_free = T) 
# HRT:

results[[3]] <- HRT(Y, X, Sigma = Sigma_true, FDR = 0.1, N = 1000, 
                    model = model, test_group = c(1))

# DML:
results[[4]] <- DML(X, Y, Sigma_true, FDR = FDR_rate, K = 8, model = model,
                    test_set = 1)


##########################################################################
          ###### Sensitiviy of dICRT to the choice of k ######
##########################################################################

# The setup of data generation and approaches could be found from the "Large size simulations" part.
# One just need to vary the setup of k for implementing dICRT with different k's, for example:

results <- dCRT(X, Y, Sigma_X = Sigma_true, 
                d.interaction = T, k = 2, FDR = 0.1, MC_free = T, 
                model = model.dCRT) 


##########################################################################
 ### dCRT implementation with no computation dimensionality reduction ###
       ############ (assess the impact of screening) ############ 
##########################################################################

# The setup of data generation and approaches could be found from the "Large size simulations" part.
# One just need to set: CDR='no' when implementing dCRT approaches without screening, for example:

results <- dCRT(X, Y, Sigma_X = Sigma_true, CDR = 'No',
                d.interaction = T, k = as.integer(2 * log(p)),
                FDR = 0.1, MC_free = T, model = model.dCRT) 

##########################################################################
         ###### Compare oCRT with dCRT and variation of oCRT ######
##########################################################################

n <- 800
p <- 800
s <- 50
r <- 0 # r <- 0.5
space <- 1 # space <- 16
model <- 'gaussian'
model.dCRT <- 'Gaussian_lasso'
magn <- sqrt(0.5 / (1 - r)) * 3 * 0.04

# Generate data:

Gen_data <- Generate_data(n, p, s, intercept = 0, 
                          model = 'linear', Y_dist = 'Gaussian',
                          Sigma = 'AR', r = r, magn = magn,
                          X_dist = 'gauss', sign_design = 'half', support_dist = space)
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- Gen_data$Sigma
beta_true <- Gen_data$beta_true

results <- vector('list', 2)

# d0CRT

results[[1]] <- dCRT(X, Y, Sigma_X = Sigma_true, candidate_set = 1, CDR = 'No', 
                     model = model.dCRT, k = 0, M = 2000, return.pvalues = T)
results[[1]]$p_values[1]

# oCRT
results[[2]] <- CRT_modified(Y, X, Sigma_X = Sigma_true, m = 2000, 
                             candidate_set = 1)

# p-value of original CRT with lasso:

results[[2]]$p_value[1]

# p-value of original CRT with lasso and no soft-thresholding:

results[[2]]$p_value[2]

# p-value of original CRT with lasso, no soft-thresholding and centering (fixing the non-centering issue):

results[[2]]$p_value[5]


##########################################################################
              ##### Assess algorithmic variability #####
##########################################################################


# The setup of data generation and approaches could be found from the "Large size simulations" part.
# One needs to generate 100 data sets independently for each setting (with different magnitudes),
# and then run the algorithms on each for (says) 50 times independently.
# Finally, he/she can implement the following code for calculating the Jaccard index:


jac_mat <- c()

for (magn in 0.025 + 0.05 * c(1:5)){
  jac_magn_vec <- 0
  for (data_ID in 1:100) {
    num <- 0 
    reject_set_lst <- vector('list', 4)
    for (t in 1:4) {
      reject_set_lst[[t]] <- vector('list', 50)
    }
    
    for (set in 1:50){
      
      ## The file name of the results for each repetition on each data set:
      
      filename <-  paste('stab_magn', as.character(magn), '_dataID', as.character(data_ID), 
                         '_seed', as.character(set), '.rda', sep = '')
      res <- try(load(filename))
      if(inherits(res, "try-error")){
        #error handling code, maybe just skip this iteration using
        next
      }
      
      for (t in 1:4) {
        reject_set_lst[[t]][[set]] <- results[[t]]$select_set
      }
    }
    
    num <- 0
    jac_index <- 0
    for (set1 in 1:49){
      for (set2 in (set1 + 1):50){
        if ((! is.null(reject_set_lst[[1]][[set1]])) & (! is.null(reject_set_lst[[1]][[set2]]))){
          num <- num + 1
          JJJ <- c()
          for (t in 1:4){
            JJJ <- c(JJJ, (length(intersect(reject_set_lst[[t]][[set1]], reject_set_lst[[t]][[set2]])) + 1e-20) /
                       (length(union(reject_set_lst[[t]][[set1]], reject_set_lst[[t]][[set2]])) + 1e-20))
          }
          jac_index <- jac_index + JJJ
        }
      }
    }
    jac_index <- jac_index / num
    jac_magn_vec <- jac_magn_vec + jac_index
  }
  jac_magn_vec <- jac_magn_vec / 100
  jac_mat <- rbind(jac_mat, jac_magn_vec)
}




##########################################################################
  ##### dCRT: resample-free version v.s. non-resample-free version #####
##########################################################################

n <- 800
p <- 800
s <- 50
r <- 0.5
magn <- 0.175
results <- vector('list', 7)

### Gamma X

# Generate data:

X_dist <- 'indept_gamma'
Gen_data <- Generate_data(n, p, s, model = 'linear', 
                          Y_dist = 'Gaussian', magn = magn,
                          X_dist = X_dist)
X <- Gen_data$X
Y <- Gen_data$Y

# Convert gamma X to gaussian according to their quantiles for resample-free dCRT:

gamma_quantile <- function(X){
  return(pgamma(6 + X * sqrt(12), shape = 3, rate = 0.5))
}

X_quantile_mat <- gamma_quantile(X)
X_con_var_mat <- matrix(1, n, p)

# Define some function for data generation of gamma X

Gen_X <- function(X, indx, num = 1){
  n <- length(X[,1])
  X_sample <- matrix(rgamma(n * num, shape = 3, rate = 0.5), n, num)
  return((X_sample - 6) / sqrt(12))
}

### Binary X

# Generate data:

X_dist <- 'indept_binary'
Gen_data <- Generate_data(n, p, s, model = 'linear', 
                          Y_dist = 'Gaussian', magn = magn,
                          X_dist = X_dist)
X <- Gen_data$X
Y <- Gen_data$Y

# Convert gamma X to gaussian according to their quantiles for resample-free dCRT:

bernoulli_quantile <- function(X){
  X <- (X + 1) / 2
  n <- length(X[,1])
  p <- length(X[1,])
  X_bar <- 0.5
  X_var <- X_bar * (1 - X_bar)
  set_0 <- which(X == 0)
  set_1 <- which(X == 1)
  X_unif <- rep(0, n * p)
  
  # perturb to make the quantile uniform:
  X_unif[set_1] <- runif(length(set_1), X_bar, 1)
  X_unif[set_0] <- runif(length(set_0), 0, X_bar)
  
  return(matrix(X_unif, n, p))
}

X_quantile_mat <- bernoulli_quantile(X)
X_con_var_mat <- matrix(1, n, p)

# Define some function for data generation of gamma X

Gen_X <- function(X, indx, num = 1){
  n <- length(X[,1])
  X_sample <- matrix(rbinom(n * num, 1, 0.5), n, num)
  return((X_sample - 0.5) * 2)
}


### Run the approaches:

# Resample-free dCRT:

results[[1]] <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                     Gen_X = Gen_X, X_quantile_mat = X_quantile_mat, 
                     X_cond_var_mat = X_con_var_mat, model = 'Gaussian_lasso',
                     k = 0, MC_free = T)

results[[2]] <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                     Gen_X = Gen_X, X_quantile_mat = X_quantile_mat, 
                     X_cond_var_mat = X_con_var_mat, model = 'Gaussian_lasso',
                     d.interaction = T, MC_free = T)

# non-resample-free dCRT (natural version):

results[[3]] <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                     Gen_X = Gen_X, X_quantile_mat = X_quantile_mat, 
                     X_cond_var_mat = X_con_var_mat, model = 'Gaussian_lasso',
                     k = 0, MC_free = F)

results[[4]] <- dCRT(X, Y,FDR = 0.1, mean_X_mat = matrix(6, n, p), 
                     Gen_X = Gen_X, X_quantile_mat = X_quantile_mat, 
                     X_cond_var_mat = X_con_var_mat, model = 'Gaussian_lasso',
                     d.interaction = T, MC_free = F)

# HRT:

results[[5]] <- HRT(Y, X, Sigma = Sigma_true, FDR = FDR_rate, N = 50000, 
                     model = 'gaussian', x_type = X_dist)

# knockoffs: 
  
  
if (X_dist == 'indept_binary'){
  knockoffs = function(X){
    n <- nrow(X)
    p <- ncol(X)
    X_knock <- rbinom(n * p, 1, 0.5)
    X_knock <- matrix(X_knock, n, p)
    X_knock <- (X_knock - 0.5) * 2
    return(X_knock)
  }
}

if (X_dist == 'indept_gamma'){
  knockoffs = function(X){
    n <- nrow(X)
    p <- ncol(X)
    X_knock <- rgamma(n * p, shape = 3, rate = 0.5)
    X_knock <- matrix(X_knock, n, p)
    X_knock <- (X_knock - 6) / sqrt(12)
    return(X_knock)
  }
}

foo = stat.glmnet_coefdiff
k_stat = function(X, X_k, y) foo(X, X_k, y, family = 'gaussian')
Knockoff_Result <- knockoff.filter(X, Y, statistic = k_stat, knockoffs = knockoffs,
                                   fdr = FDR_rate, offset = 1)
selection_set_KO <- Knockoff_Result$selected
a <- list()
a[['select_set']] <- as.vector(selection_set_KO)
results[[6]] <- a
  
# DML:

results[[7]] <- DML(X, Y, Sigma_true, FDR = FDR_rate, K = 8, model = 'gaussian')


##########################################################################
###### Compare with DML/GCM with small size and laplace tailed data ######
##########################################################################

n <- 30
p <- 150
s <- 10
magn <- 2 # magn in 2 * c(1:5)
para_x <- 3

## Generate data

Gen_data <- Generate_data(n, p, s, intercept = 0, Y_dist = 'Laplace',
                          Sigma = 'AR', r = 0.4, magn = magn,
                          X_dist = 'Laplace', para_x = para_x)
X <- Gen_data$X
Y <- Gen_data$Y
Sigma_true <- diag(rep(2 / 9, p))
results <- vector('list', 5)

## Run the approaches:

Gen_X <- function(X, indx, num = 1){
  n <- nrow(X)
  sign_X <- 2 * rbinom(num * n, 1, 0.5) - 1
  X_sample <- sign_X * rexp(num * n, 3) 
  X_sample <- matrix(X_sample, n, num)
  return(X_sample)
}

results[[1]] <- dCRT(X, Y, Sigma_X = NULL, FDR = 0.1,
                     mean_X_mat = matrix(0, n, p), Gen_X = Gen_X, CDR = 'No',
                     model = 'Gaussian_lasso', d.interaction = F, 
                     k = 0, M = 30000, MC_free = F)

results[[2]] <- dCRT(X, Y, Sigma_X = NULL, FDR = 0.1,
                     mean_X_mat = matrix(0, n, p), Gen_X = Gen_X, CDR = 'No',
                     model = 'Gaussian_lasso', d.interaction = T,  M = 30000, MC_free = F)

results[[3]] <- HRT(Y, X, FDR = 0.1, N = 30000, 
                    model = 'gaussian', x_type = 'laplace')

results[[4]] <- DML(X, Y, Sigma_true, FDR = 0.1, K = 8, model = 'gaussian')

results[[5]] <- GCM(X, Y, Sigma_true, FDR = 0.1, model = 'gaussian')




