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
  c(4, 4, 4, 6*sqrt(2), 4/3, rep(0, times = 5))
);

n <- 100;
rho <- 0.9;
sd <- 0.1;
L <- 100;
nsim <- 100;

par(mfrow=c(2,2));
b <- "Beta's: ";

# Comparison of ARM linear regression 
for(beta_star in betas_linear) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "gaussian");
  X <- data$X;
  y <- data$y;
  candidate_models <- candidateModels(X = X, y = y, L = L, family = "gaussian");
  linear_regression_arm <- linearRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
  soilArmLinearRegression <- as.vector(linear_regression_arm$soil_importance);
  print(soilArmLinearRegression);
  v_ARM <- SOIL(X, y, family = "gaussian", weight_type = "ARM", prior = TRUE);
  soilArmLinearRegression_author <- as.vector(v_ARM$importance);
  betas_string <- paste(beta_star, collapse = ", ");
  plot(soilArmLinearRegression, type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "SOIL importance", xlab = "variable number", 
       main="Comparison ARM linear regression",lwd = 2, sub = paste(b, betas_string));
  points(soilArmLinearRegression_author, type = "o", col = "red");
  legend("left", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}

# Comparison of BIC linear regression
for(beta_star in betas_linear) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "gaussian");
  X <- data$X;
  y <- data$y;
  candidate_models <- candidateModels(X = X, y = y, L = L, family = "gaussian");
  linear_regression_bic <- linearRegressionBIC(X = X, y = y, candidate_models = candidate_models);
  soilBicLinearRegression <- as.vector(linear_regression_bic$soil_importance);
  print(soilBicLinearRegression)
  v_BIC <- SOIL(X, y, family = "gaussian", weight_type = "BIC", prior = TRUE);
  soilBicLinearRegression_author <- as.vector(v_BIC$importance);
  plot(soilBicLinearRegression, type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "SOIL importance", xlab = "variable number", 
       main="Comparison BIC linear regression",lwd = 2, sub = paste(b, betas_string));
  points(soilBicLinearRegression_author, type = "o", col = "red");
  legend("left", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
  
}

betas_logistic <- list(
  c(4, rep(0, times = 4)),
  c(4, 2, 3, rep(0, times = 5)),
  c(1, 1/2, 1/3, 1/4, 1/5, 1/6, 0),
  c(4, 4, 4, 6*sqrt(2), 4/3, rep(0, times = 5))
);

# Comparison of ARM logistic regression
for(beta_star in betas_logistic) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "binomial");
  X <- data$X;
  y <- data$y;
  candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
  logistic_regression_arm <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
  soilArmLogisticRegression <- as.vector(logistic_regression_arm$soil_importance);
  print(soilArmLogisticRegression);
  v_ARM_logistic <- SOIL(X, y, family = "binomial", weight_type = "ARM", prior = TRUE);
  soilArmLogisticRegression_author <- as.vector(v_ARM_logistic$importance);
  betas_string <- paste(beta_star, collapse = ", ");
  plot(soilArmLogisticRegression, type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "SOIL importance", xlab = "variable number", 
       main="Comparison ARM logistic regression",lwd = 2, sub = paste(b, betas_string));
  points(soilArmLogisticRegression_author, type = "o", col = "red");
  legend("topright", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}

# Comparison of BIC logistic regression
for(beta_star in betas_logistic) {
  beta <- matrix(beta_star, ncol = 1);
  p <- nrow(beta);
  data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta, family = "binomial");
  X <- data$X;
  y <- data$y;
  candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
  logistic_regression_bic <- logisticRegressionBIC(X = X, y = y, candidate_models = candidate_models);
  soilBicLogisticRegression <- as.vector(logistic_regression_bic$soil_importance);
  print(soilBicLogisticRegression);
  v_BIC_logistic <- SOIL(X, y, family = "binomial", weight_type = "BIC", prior = TRUE);
  soilBicLogisticRegression_author <- as.vector(v_ARM_logistic$importance);
  betas_string <- paste(beta_star, collapse = ", ");
  plot(soilBicLogisticRegression, type = "o", ylim = c(0,1), xlim = c(1,p), ylab = "SOIL importance", xlab = "variable number", 
       main="Comparison BIC logistic regression",lwd = 2, sub = paste(b, betas_string));
  points(soilBicLogisticRegression_author, type = "o", col = "red");
  legend("topright", legend=c("Jonas", "Ye et al."), col=c("black", "red"), lwd = c(2,1) );
}
