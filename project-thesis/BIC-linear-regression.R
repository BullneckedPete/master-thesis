linRegBIC <- function(X, y, candidate_models, likelihood) {
  
  s_k <- apply(candidate_models, 2, sum);
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
    Xs_k <- X[,indices];
    
    l_k <- 0; ## maximized likelihood of model k
    I_k[i] <- -2*log(lk) + sk[i] * log(n); ## definition of I_k_BIC
    C_k[i] <- s_k[i] * log( exp(1) * p / s_k[i] ) + 2 * log(s_k[i] + 2);
    w_k_numerator[i] <- exp(-I_k[i] / 2 - phi * C_k);
    
  }
  w_k_denominator <- sum( exp(I_k/2 - phi*C_k) );
  weight_vector <- w_k_numerator / w_k_denominator;
  soil_importance <- weight_vector %*% t(candidate_models);
  return (
    list(
      weight_vector = round(weight_vector, 2), 
      soil_importance =  soil_importance
    )
  );
}