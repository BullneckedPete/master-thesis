source("helpers.R");
linearRegressionARM <- function(X, y, nsim, candidate_models, psi = 1) {
  weight_vectors <- matrix(0, ncol = ncol(candidate_models), nrow = nsim);
  n <- nrow(X);
  p <- ncol(X);
  for (j in 1:nsim) {
    sampleSize <- ceiling((nrow(X)/2));
    randomSample <- selectionRejectionSample(U = X, n = sampleSize);
    D1 <- X[randomSample, ]; 
    y1 <- y[randomSample];
    D2 <- X[-randomSample, ]; 
    y2 <- y[-randomSample];
    s_k <- apply(candidate_models, 2, sum);
    k <- ncol(candidate_models);
    w_k_numerator <- numeric(k);
    for (i in 1:k) {
      indices <- as.vector(which(candidate_models[,i] != 0));
      if (length(indices) == 0) {
        fit_lm_k <- lm( y1 ~ 1 );
      } else {
        Xs_k <- D1[ ,indices];
        reg_data <- as.data.frame(cbind(y1, Xs_k));
        fit_lm_k <- lm( y1 ~ . , data = reg_data );
      }
      if (any(is.na(fit_lm_k$coefficients[-1]))) {
        w_k_numerator[i] <- rep(0, times = length(fit_lm_k$coefficients[-1]));
      } else {
        sigma_hat <- summary(fit_lm_k)$sigma;
        coeff <-  fit_lm_k$coefficients;
        bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
        prediction_data <- as.matrix(cbind(1,D2[,indices]));
        y2_prediction <- as.vector(prediction_data %*% bs_hat_k);
        C_k <- calculateCk(s_k = s_k[i], p = p);
        w_k_numerator[i] <- exp(1)^(-psi*C_k)*(sigma_hat^(-n/2))*prod(exp(-(sigma_hat^-2)*((y2 - y2_prediction)^2)/2));
      }
    }
    w_k_numerator[is.nan(w_k_numerator)] <- 0;
    weight_vectors[j, ] <- w_k_numerator / sum(w_k_numerator);
  }
  weight_vector <- colMeans(checkIfWeightsValid(weight_vectors));
  soil_importance <- weight_vector %*% t(candidate_models);
  return (
    list(
      weight_vector = round(weight_vector, 6), 
      soil_importance =  soil_importance
    )
  );
}
