source("helpers.R");
logisticRegressionBIC <- function(X, y, candidate_models, psi = 1) {
  I_k <- rep(0, ncol(candidate_models));
  w_k_numerator <- C_k <- numeric(ncol(candidate_models));
  s_k <- apply(candidate_models, 2, sum);
  k <- ncol(candidate_models);
  p <- ncol(X);
  for (i in 1:k) {
    indices <- as.vector(which(candidate_models[,i] != 0));
    if (length(indices) == 0) {
      fit_logReg_k <- glm( y ~ 1, family = "binomial" );
    } else {
      Xs_k <- X[ ,indices];
      reg_data <- as.data.frame(cbind(y, Xs_k));
      colnames(reg_data)[1] <- "y";
      fit_logReg_k <- glm( y ~ . , data = reg_data, family = "binomial" );
    }
    I_k[i] <- BIC(fit_logReg_k);
    C_k[i] <- calculateCk(s_k[i], p);
    if(is.infinite(I_k[i])) {
      w_k_numerator[i] <- rep(0, times = k);
    } else {
      w_k_numerator[i] <- exp(-I_k[i]/2 - psi*C_k[i]);
    }
   
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
