source("generate-data.R");
source("candidate-models.R");
source("ARM-linear-regression.R");
source("BIC-linear-regression.R");
source("ARM-logistic-regression.R");
source("BIC-logistic-regression.R");
require("SOIL");

betas_linear <- list(
  c(4, rep(0, times = 4)),
  c(4, 2, 3, rep(0, times = 5)),
  c(rep(1, times = 10)),
  c(4, 4, 4, -6*sqrt(2), 4/3, rep(0, times = 5))
);

n <- 100;
rho <- 0.9;
sd <- 0.1;
L <- 100;
nsim <- 100;

par(mfrow=c(2,2));
b <- "Beta's: ";

# repeating the simulation 100 times to compute the average variable importance measures
nsimu <- 100;
# Comparison of ARM linear regression
for(beta_star in betas_linear) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  soilArmLinearRegression <- matrix(0, nrow = nsimu, ncol = p);
  soilArmLinearRegression_author <- matrix(0, nrow = nsimu, ncol = p);
  for (i in 1:nsimu) {
    data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "gaussian");
    X <- data$X;
    y <- data$y;
    candidate_models <- candidateModels(X = X, y = y, L = L, family = "gaussian");
    v_ARM <- SOIL(X, y, family = "gaussian", weight_type = "ARM", prior = TRUE);
    soilArmLinearRegression_author[i, ] <- as.vector(v_ARM$importance);
    linear_regression_arm <- linearRegressionARM(X = X, y = y, nsim = nsim, candidate_models = t(v_ARM$candidate_models_cleaned[-1,]));
    soilArmLinearRegression[i, ] <- as.vector(linear_regression_arm$soil_importance);
    
  }
  print(colMeans(soilArmLinearRegression));
  print(colMeans(soilArmLinearRegression_author));
  betas_string <- paste(round(beta_star,2), collapse = ", ");
  plot(colMeans(soilArmLinearRegression), type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "Importance", xlab = "Variable Index", 
       lwd = 2, sub = paste(b, betas_string));
  points(colMeans(soilArmLinearRegression_author), type = "o", col = "red");
  legend("left", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) ); 
}


# Comparison of BIC linear regression
for(beta_star in betas_linear) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  soilBicLinearRegression <-  matrix(0, nrow = nsimu, ncol = p);
  soilBicLinearRegression_author <-  matrix(0, nrow = nsimu, ncol = p);
  for(i in 1:nsimu) {
    data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "gaussian");
    X <- data$X;
    y <- data$y;
    candidate_models <- candidateModels(X = X, y = y, L = L, family = "gaussian");
    v_BIC <- SOIL(X, y, family = "gaussian", weight_type = "BIC", prior = TRUE);
    soilBicLinearRegression_author[i, ] <- as.vector(v_BIC$importance);
    linear_regression_bic <- linearRegressionBIC(X = X, y = y, candidate_models = t(v_BIC$candidate_models_cleaned[-1,]));
    soilBicLinearRegression[i, ] <- as.vector(linear_regression_bic$soil_importance);
    
  }
  print(colMeans(soilBicLinearRegression));
  print(colMeans(soilBicLinearRegression_author));
  betas_string <- paste(round(beta_star,2), collapse = ", ");
  plot(colMeans(soilBicLinearRegression), type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "Importance", xlab = "Variable Index", 
      lwd = 2, sub = paste(b, betas_string));
  points(colMeans(soilBicLinearRegression_author), type = "o", col = "red");
  legend("left", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}




betas_logistic <- list(
  #c(4, rep(0, times = 4)),
  #c(4, 2, 3, rep(0, times = 5)),
  c(1, 1/2, 1/3, 1/4, 1/5, 1/6, 0),
  c(4, 4, 4, 6*sqrt(2), 4/3, rep(0, times = 5))
);

# Comparison of ARM logistic regression
for(beta_star in betas_logistic) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  soilArmLogisticRegression <- matrix(0, nrow = nsimu, ncol = p);
  soilArmLogisticRegression_author <- matrix(0, nrow = nsimu, ncol = p);
  for (i in 1:nsimu) {
    data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "binomial");
    X <- data$X;
    y <- data$y;
    candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
    v_ARM_logistic <- SOIL(X, y, family = "binomial", weight_type = "ARM", prior = TRUE);
    soilArmLogisticRegression_author[i, ] <- as.vector(v_ARM_logistic$importance);
    #logistic_regression_arm <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
    logistic_regression_arm <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = t( v_ARM_logistic$candidate_models_cleaned));
    
    soilArmLogisticRegression[i, ] <- as.vector(logistic_regression_arm$soil_importance);
  }
  print(colMeans(soilArmLogisticRegression));
  print(colMeans(soilArmLogisticRegression_author));
  betas_string <- paste(round(beta_star,2), collapse = ", ");
  plot(colMeans(soilArmLogisticRegression), type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "Importance", xlab = "Variable Index", 
      lwd = 2, sub = paste(b, betas_string));
  points(colMeans(soilArmLogisticRegression_author), type = "o", col = "red");
  legend("topright", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}

# Comparison of BIC logistic regression
for(beta_star in betas_logistic) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  soilBicLogisticRegression <- matrix(0, nrow = nsimu, ncol = p);
  soilBicLogisticRegression_author <- matrix(0, nrow = nsimu, ncol = p);
  for (i in 1:nsimu) {
    data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "binomial");
    X <- data$X;
    y <- data$y;
    candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
    v_BIC_logistic <- SOIL(X, y, family = "binomial", weight_type = "BIC", prior = TRUE);
    soilBicLogisticRegression_author[i, ] <- as.vector(v_BIC_logistic$importance);
    logistic_regression_bic <- logisticRegressionBIC(X = X, y = y, candidate_models = t(v_BIC_logistic$candidate_models_cleaned));
    soilBicLogisticRegression[i, ] <- as.vector(logistic_regression_bic$soil_importance);
  }
  print(colMeans(soilBicLogisticRegression));
  print(colMeans(soilBicLogisticRegression_author));
  betas_string <- paste(round(beta_star,2), collapse = ", ");
  plot(colMeans(soilBicLogisticRegression), type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "Importance", xlab = "Variable Index", 
      lwd = 2, sub = paste(b, betas_string));
  points(colMeans(soilBicLogisticRegression_author), type = "o", col = "red");
  legend("topright", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}
