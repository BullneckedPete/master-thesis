linearRegressionBIC <- function(X, y, candidate_models) {
  if (nrow(candidate_models) > 1) {
    s_k <- apply(candidate_models, 2, sum);
  } else {
    s_k <- candidate_models;
  }
  k <- ncol(candidate_models);
  I_k <- rep(0, ncol(candidate_models));
  n <- nrow(X);
  p <- ncol(X);
  w_k_numerator <- C_k <- numeric(k);
  phi <- 1;
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
    # calculate the residual sum of squares, defined as t(y)%*%y + t(y)%*% X %*% solve(t(X)%*% X) %*% t(X)%*% y
    # for more infos, check: https://en.wikipedia.org/wiki/Residual_sum_of_squares
    residualSumOfSquares <- sum(resid(fit_lm_k)^2);
    # definition of I_k_BIC, see 'Gaussian special case' in https://en.wikipedia.org/wiki/Bayesian_information_criterion
    I_k[i] <- n*log(residualSumOfSquares / n) + p*log(n);
    C_k[i] <- s_k[i] * log( exp(1) * p / s_k[i] ) + 2 * log(s_k[i] + 2);
    w_k_numerator[i] <- exp(-I_k[i] / 2 - phi * C_k[i]);
  }
  w_k_denominator <- sum( exp(-I_k/2 - phi*C_k) );
  weight_vector <- w_k_numerator / w_k_denominator;
  soil_importance <- weight_vector %*% t(candidate_models);
  return (
    list(
      weight_vector = round(weight_vector, 6), 
      soil_importance =  soil_importance
    )
  );
}
