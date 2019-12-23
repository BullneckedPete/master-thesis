linRegARM <- function(X, y, nsim, candidate_models) {
  sim_soil_importance <- matrix(0, ncol = nrow(candidate_models), nrow = nsim);
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
    s_k <- apply(candidate_models, 2, sum);
    # phi: a positive number to control the improvement of the prior weight (defaults to 1)
    phi <- 1;
    k <- ncol(candidate_models);
    w_k_nominator <- numeric(k);
    for (i in 1:k) {
      
      # get the indices of the model to select the right positions of the design matrix
      indices <- as.vector(which(candidate_models[,i] != 0));
      
      # extract the right indices form the candidate model in the current iteration
      # out of the design matrix
      Xs_k <- D1[,indices];
      
      # prepare the data for the linear regression
      reg_data <- as.data.frame(cbind(y1, Xs_k));
      
      # fit standard linear regression of y on Xs_k using the training set D1
      fit_lm_k <- lm( y1 ~ . , reg_data );
      
      # get estimated standard deviation and coefficients
      sigma_hat <- summary(fit_lm_k)$sigma;
      coeff <-  fit_lm_k$coefficients[-1];
      bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
      
      # prepare the prediction data and compute the prediction t(Xs_k) %*% bs_hat on the test set D2
      prediction_data <- as.matrix(D2[,indices]);
      y2_prediction <- prediction_data %*% bs_hat_k;
      
      # compute the constant c_k used to compute the weight w_k later
      c_k <- s_k[i] * log( exp(1) * p / s_k[i] ) + 2 * log(s_k[i] + 2);
      
      # compute the nominator of the weight vector w_k for each candidate model
      w_k_nominator[i] <- exp(1) ^ (-phi*c_k) * (sigma_hat^(-n/2)) * prod(exp ( -sigma_hat^-2 * (y2 - as.vector(y2_prediction))^2 / 2 ));
    }
    w_k_nominator[!is.finite(w_k_nominator)] <- 0;
    weight_vector <- w_k_nominator / sum(w_k_nominator);
    soil_importance <- weight_vector %*% t(candidate_models);
    sim_soil_importance[j, ] <- soil_importance;
  }
  return(
    list(
      weight_vector = round(weight_vector, 2), 
      soil_importance =  colSums(sim_soil_importance)/nsim
  ));
}

