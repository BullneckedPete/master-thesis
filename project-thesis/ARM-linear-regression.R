source("helpers.R");
linearRegressionARM <- function(X, y, nsim, candidate_models, phi = 1) {
  weight_vectors <- matrix(0, ncol = ncol(candidate_models), nrow = nsim);
  n <- nrow(X);
  p <- ncol(X);
  for (j in 1:nsim) {
    # create a training and a test set of equal size
    # round up, if nrow(X) is odd
    sampleSize <- ceiling((nrow(X)/2));
    # randomly change rows
    randomSample <- selectionRejectionSample(U = X, n = sampleSize);
    # training sets
    D1 <- X[randomSample, ]; 
    y1 <- y[randomSample];
    # test sets
    D2 <- X[-randomSample, ]; 
    y2 <- y[-randomSample];
    # s_k: number of non constant predictors for model k, that is |A^k|
    s_k <- apply(candidate_models, 2, sum);
    k <- ncol(candidate_models);
    w_k_numerator <- numeric(k);
    for (i in 1:k) {
      # get the indices of the model to select the right positions of the design matrix
      indices <- as.vector(which(candidate_models[,i] != 0));
      # extract the right indices form the candidate model in the current iteration
      # out of the design matrix
      Xs_k <- D1[,indices];
      # prepare the data for the linear regression
      reg_data <- as.data.frame(cbind(y1, Xs_k));
      # fit standard linear regression of y on Xs_k using the training set D1
      fit_lm_k <- lm( y1 ~ . , data = reg_data );
      if (any(is.na(fit_lm_k$coefficients[-1]))) {
        # not enough information to estimate the parameters
        w_k_numerator[i] <- rep(0, times = length(fit_lm_k$coefficients[-1]));
      } else {
        # get estimated standard deviation and coefficients
        sigma_hat <- summary(fit_lm_k)$sigma;
        coeff <-  fit_lm_k$coefficients[-1];
        bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
        # prepare the prediction data and compute the prediction t(Xs_k) %*% bs_hat on the test set D2
        prediction_data <- as.matrix(D2[,indices]);
        y2_prediction <- as.vector(prediction_data %*% bs_hat_k);
        # compute the constant c_k used to compute the weight w_k later
        C_k <- calculateCk(s_k = s_k[i], p = p);
        # compute the nominator of the weight vector w_k for each candidate model
        w_k_numerator[i] <- exp(1)^(-phi*C_k) * (sigma_hat^(-n/2)) * 
                              prod(exp(-(sigma_hat^-2) * ((y2 - y2_prediction)^2)/2));
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
