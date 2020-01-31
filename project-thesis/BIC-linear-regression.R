source("helpers.R");
linearRegressionBIC <- function(X, y, candidate_models, phi = 1) {
  s_k <- apply(candidate_models, 2, sum);
  k <- ncol(candidate_models);
  I_k <- rep(0, ncol(candidate_models));
  n <- nrow(X);
  p <- ncol(X);
  w_k_numerator <- C_k <- numeric(k);
  for (i in 1:k) {
    # get the indices of the model to select the right positions of the design matrix
    indices <- as.vector(which(candidate_models[,i] != 0));
    # extract the right indices form the candidate model in the current iteration
    # out of the design matrix
    Xs_k <- X[ ,indices];
    # prepare the data for the linear regression
    reg_data <- as.data.frame(cbind(y, Xs_k));
    colnames(reg_data)[1] <- "y";
    # fit linear regression of y on Xs_k
    fit_lm_k <- lm( y ~ . , data = reg_data );
    # for more infos, check: https://en.wikipedia.org/wiki/Residual_sum_of_squares
    residualSumOfSquares <- sum(resid(fit_lm_k)^2);
    # definition of I_k_BIC, see 'Gaussian special case' in https://en.wikipedia.org/wiki/Bayesian_information_criterion
    I_k[i] <- BIC(fit_lm_k);
    C_k[i] <- calculateCk(s_k[i], p);
    w_k_numerator[i] <- exp(-I_k[i] / 2 - phi * C_k[i]);
  }
  w_k_numerator[is.nan(w_k_numerator)] <- 0;
  weight_vector <- w_k_numerator / sum(w_k_numerator);
  soil_importance <- weight_vector %*% t(candidate_models);
  return (
    list(
      weight_vector = round(weight_vector, 6), 
      soil_importance =  soil_importance
    )
  );
}
