logisticRegressionARM <- function(X, y, nsim, candidate_models) {
  sim_soil_importance <- matrix(NA, ncol = nrow(candidate_models), nrow = nsim);
  weight_vectors <- matrix(NA, ncol = ncol(candidate_models), nrow = nsim);
  n <- nrow(X);
  p <- ncol(X);
  for (j in 1:nsim) {
    # create a training and a test set of equal size
    # round up, if nrow(X) is odd
    cut <- ceiling((nrow(X)/2));
    # randomly change rows
    randomSample <- sample(1:nrow(X), nrow(X), replace = FALSE);
    X <- X[randomSample, ];
    y <- y[randomSample];
    # training set
    D1 <- X[1:cut, ];
    y1 <- y[1:cut];
    # test set
    D2 <- X[(cut+1):nrow(X), ]; 
    y2 <- y[(cut+1):length(y)];
    # s_k: number of non constant predictors for model k, that is |A^k|
    if (nrow(candidate_models) > 1) {
      s_k <- apply(candidate_models, 2, sum);
    } else {
      s_k <- candidate_models;
    }
    # phi: a positive number to control the improvement of the prior weight (defaults to 1)
    phi <- 1;
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
      # get estimated coefficients
      coeff <-  fit_logreg_k$coef;
      bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
      # responding function of the predicted conditional probability:
      # compute the predicted probability on the test set {i | i E D2}
      pred_p_hat_k <-  1 / (1 + exp(-((cbind(1, D2[ ,indices]) %*% bs_hat_k))));
      # compute the constant c_k used to compute the weight w_k later
      if (s_k[i] != 0) {
        c_k <-  2 * log(s_k[i] + 2) + s_k[i] * log( exp(1) * p / s_k[i] );
      } else {
        c_k <-  2 * log(s_k[i] + 2);
      }
      # compute the nominator of the weight vector w_k for each candidate model
      w_k_nominator[i] <- exp(1)^(-phi*c_k) * prod( (pred_p_hat_k)^y2 * (1-pred_p_hat_k)^(1-y2));
      #w_k_nominator[i] <- abs(w_k_nominator[i]);
      if (!is.finite(w_k_nominator[i]))
         w_k_nominator[i] <- 1;
    }
    #w_k_nominator[!is.finite(w_k_nominator)] <- 0;
    weight_vectors[j, ] <- w_k_nominator / sum(w_k_nominator);
    # sim_soil_importance[j, ] <- soil_importance;
  }
  weight_vector <- colMeans(weight_vectors);
  soil_importance <- weight_vector %*% t(candidate_models);
  return(
    list(
      weight_vector = round(weight_vector, 6), 
      soil_importance = soil_importance
      #soil_importance =  colSums(sim_soil_importance)/nsim
    )
  );
}
