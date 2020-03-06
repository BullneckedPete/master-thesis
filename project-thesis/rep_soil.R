source("candidate-models.R");
source("ARM-linear-regression.R");
source("BIC-linear-regression.R");
source("ARM-logistic-regression.R");
source("BIC-logistic-regression.R");

# @params psi: a positive number to control the improvement of the prior weight (defaults to 1)
# @params nsim: number of simulations when weight_type ARM is used
# @params L: number of regularization parameters to coose the candidate models
# @params family: "gaussian" for linear regression or "binomial" for logistic regression
# @params weight_type: either "ARM" for adaptive regression by mixing or "BIC" for
#         using BIC-p weighting
rep_SOIL <- function(X, y, nsim = 100, psi = 1, L = 100, family = "gaussian", 
                     weight_type = "ARM", customize = FALSE, candidate_models = matrix(0)) {
  
  if ( !(family %in% c("gaussian", "binomial")) ) {
    stop("The parameter 'family' has to be 'gaussian' or 'binomial'!");
  }
  if (!(weight_type %in% c("ARM", "BIC")) ) {
    stop("The parameter 'weight_type' has to be 'ARM' or 'BIC'!");
  }
  
  if (!customize) {
    candidate_models <- candidateModels(X = X, y = y, L = L, family = family);
  } else {
    candidate_models <- candidate_models;
  }
  
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



