source("generate-data.R");
source("candidate-models.R");
source("ARM-linear-regression.R");
source("BIC-linear-regression.R");
source("ARM-logistic-regression.R");
source("BIC-logistic-regression.R");

# @params n: number of rows of the design matrix
# @params p: number of covariates of the design matrix
# @params beta_star: vector of betas that generates the data
# @params sd: the standard deviation udes to generate the design matrix
# @params rho: parameter used to generate the covariance matrix for sampling 
#         from a multivariate normal distribution
# @params phi: a positive number to control the improvement of the prior weight (defaults to 1)
# @params nsim: number of simulations when weight_type ARM is used
# @params L: number of regularization parameters to coose the candidate models
# @params family: "gaussian" for linear regression or "binomial" for logistic regression
# @params weight_type: either "ARM" for adaprive regression by mixing or "BIC" for
#         using BIC-p weighting
soil <- function(n, p, beta_star, sd, rho, nsim = 100, phi = 1, L = 100, family = "gaussian", weight_type = "ARM") {
  
  if ( p != length(beta_star) ) {
    stop("The number of independent variables must match the length of the vector 'betas'!");
  }
  if ( !(family %in% c("gaussian", "binomial")) ) {
    stop("The parameter 'family' has to be 'gaussian' or 'binomial'!");
  }
  if (!(weight_type %in% c("ARM", "BIC")) ) {
    stop("The parameter 'weight_type' has to be 'ARM' or 'BIC'!");
  }
  data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta_star, family = family);
  X <- data$X;
  y <- data$y;
  candidate_models <- candidateModels(X = X, y = y, L = L, family = family);
  if (family == "gaussian") {
    if (weight_type == "ARM") {
      soil <- linearRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
    } else {
      soil <- linearRegressionBIC(X = X, y = y, candidate_models = candidate_models);
    }
  } else {
    if (weight_type == "ARM") {
      soil <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
    } else {
      soil <- logisticRegressionBIC(X = X, y = y, candidate_models = candidate_models);
    }
  }
  soil$candidate_models <- candidate_models;
  return(soil);
  
}

soil(n = 100, p = 5, beta_star = c(4,0,0,0,0), sd = 0.01, 
     rho = 0.9, nsim = 100, phi = 1, L = 100, family = "gaussian", weight_type = "ARM");
soil(n = 100, p = 5, beta_star = c(4,0,0,0,0), sd = 0.01, 
     rho = 0.9, nsim = 100, phi = 1, L = 100, family = "gaussian", weight_type = "BIC");
soil(n = 100, p = 5, beta_star = c(4,0,0,0,0), sd = 0.01, 
     rho = 0.9, nsim = 100, phi = 1, L = 100, family = "binomial", weight_type = "ARM");
soil(n = 100, p = 5, beta_star = c(4,0,0,0,0), sd = 0.01, 
     rho = 0.9, nsim = 100, phi = 1, L = 100, family = "binomial", weight_type = "BIC")
