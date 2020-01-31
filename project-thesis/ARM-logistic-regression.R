source("helpers.R");
logisticRegressionARM <- function(X, y, nsim, candidate_models, phi = 1) {
  weight_vectors <- matrix(NA, ncol = ncol(candidate_models), nrow = nsim);
  n <- nrow(X);
  p <- ncol(X);
  for (j in 1:nsim) {
    # create a training and a test set of equal size
    # round up, if nrow(X) is odd
    sampleSize <- ceiling((nrow(X)/2));
    # randomly change rows
    randomSample <- selectionRejectionSample(U = X, n = sampleSize);
    # training set
    D1 <- X[randomSample, ];
    y1 <- y[randomSample];
    # test set
    D2 <- X[-randomSample, ]; 
    y2 <- y[-randomSample];
    # s_k: number of non constant predictors for model k, that is |A^k|
    s_k <- apply(candidate_models, 2, sum);
    k <- ncol(candidate_models);
    w_k_nominator <- numeric(k);
    for (i in 1:k) {
      # get the indices of the model to select the right positions of the design matrix
      indices <- as.vector(which(candidate_models[,i] != 0));
      # extract the right indices form the candidate model in the current iteration
      # out of the design matrix
      Xs_k <- D1[ ,indices];
      # prepare the data for the logistic regression
      reg_data <- as.data.frame(cbind(y1, Xs_k));
      # fit logistic regression of y on Xs_k using the training set D1
      fit_logreg_k <- glm( y1 ~ . , data = reg_data, family = "binomial" );
      if (any(is.na(fit_logreg_k))) {
        # not enough information to estimate the parameters
        w_k_nominator[i] <- rep(0, times = length(fit_logreg_k$coef));
      } else {
        # get estimated coefficients
        coeff <-  fit_logreg_k$coef;
        bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
        # responding function of the predicted conditional probability:
        # compute the predicted probability on the test set {i | i E D2}
        pred_p_hat_k <-  1 / (1 + exp(-((cbind(1, D2[ ,indices]) %*% bs_hat_k))));
        # compute the constant c_k used to compute the weight w_k later
        C_k <- calculateCk(s_k = s_k[i], p = p);
        # compute the nominator of the weight vector w_k for each candidate model
        w_k_nominator[i] <- exp(1)^(-phi*C_k) * 
          prod((pred_p_hat_k)^y2 * (1-pred_p_hat_k)^(1-y2));
      }
    }
    w_k_nominator[is.nan(w_k_nominator)] <- 0;
    weight_vectors[j, ] <- w_k_nominator / sum(w_k_nominator);
  }
  weight_vector <- colMeans(checkIfWeightsValid(weight_vectors));
  soil_importance <- weight_vector %*% t(candidate_models);
  return(
    list(
      weight_vector = round(weight_vector, 6), 
      soil_importance = soil_importance
    )
  );
}
