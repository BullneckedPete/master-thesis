source("generate-data.R");
source("candidate-models.R");
source("ARM-linear-regression.R");
source("BIC-linear-regression.R");
require("SOIL");

plotSimulationLinearRegression <- function(betas, n, rho, sd, nreps, limitAxis, noise = FALSE, quad = FALSE) {
  par(mfrow=c(2,2));
  nsim <- 100;
  results <- list();
  if (length(rho) != length(sd)) {
    stop("rho and sd need to be vectors of same length!");
  }
  for(j in 1:length(rho)) {
    beta <- matrix(betas, ncol = 1);
    if (!noise) {
      p <- nrow(beta);
    } else {
      p <- nrow(beta) + 1;
    }
    soilArmLinearRegression <- matrix(0, nrow = nreps, ncol = p);
    soilArmLinearRegression_author <- matrix(0, nrow = nreps, ncol = p);
    soilBicLinearRegression <-  matrix(0, nrow = nreps, ncol = p);
    soilBicLinearRegression_author <-  matrix(0, nrow = nreps, ncol = p);
    for (i in 1:nreps) {
      # generate the data
      if (!noise && !quad) {
        data <- generateData(n = n, p = p, rho = rho[j], sd = sd[j], beta_star = betas, family = "gaussian");
      } else if (quad) {
        data <- generateData(n = n, p = p-6, rho = rho[j], sd = sd[j], beta_star = betas, family = "gaussian", quad = TRUE);
      } else if (noise) {
        data <- generateData(n = n, p = p-1, rho = rho[j], sd = sd[j], beta_star = betas, family = "gaussian");
      }
      
      X <- data$X;
      y <- data$y;
      
      # Example 2: specific case when adding arbitrary noise
      if (noise) {
        X <- cbind(X, 0.5*X[,1]+ 2*X[,4]+ rnorm(n, mean = 0, sd = 0.01));
      }
      
      # obtain the candidate models
      candidate_models <- candidateModels(X = X, y = y, family = "gaussian");
      
      # ARM
      v_ARM <- SOIL(X, y, family = "gaussian", weight_type = "ARM", prior = TRUE);
      soilArmLinearRegression_author[i, ] <- as.vector(v_ARM$importance);
      linear_regression_arm <- linearRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
      soilArmLinearRegression[i, ] <- as.vector(linear_regression_arm$soil_importance);
      
      # BIC 
      v_BIC <- SOIL(X, y, family = "gaussian", weight_type = "BIC", prior = TRUE);
      soilBicLinearRegression_author[i, ] <- as.vector(v_BIC$importance);
      linear_regression_bic <- linearRegressionBIC(X = X, y = y, candidate_models = candidate_models);
      soilBicLinearRegression[i, ] <- as.vector(linear_regression_bic$soil_importance);
    }
  
    rep_mean_arm <- colMeans(soilArmLinearRegression);
    org_mean_arm <- colMeans(soilArmLinearRegression_author);
    rep_mean_bic <- colMeans(soilBicLinearRegression);
    org_mean_bic <- colMeans(soilBicLinearRegression_author);
    
    # create the line charts
    title <- bquote( rho == .(rho[j]) ~ ", " ~ sigma == .(sd[j]) );
    plot(rep_mean_arm, type = "b", ylim = c(0,1.4), xlim = c(1,limitAxis), ylab = "Importance", xlab = "Variable Index", 
         pch = 15, main = title);
    points(org_mean_arm, type = "b", col = "red", pch = 16);
    points(rep_mean_bic, type = "b", col = "blue", pch = 17);
    points(org_mean_bic, type = "b", col = "forestgreen", pch = 18);
    legend("topright", legend=c("Replication ARM", "SOIL ARM", "Replication BIC", "SOIL BIC"), 
           col=c("black", "red", "blue", "forestgreen"), pch = 15:18, ncol=2); 
    results[[paste("Plot ",j)]] <- list(
        "Replication ARM" = rep_mean_arm,
        "SOIL ARM" = org_mean_arm,
        "Replication BIC" = rep_mean_bic,
        "SOIL BIC" = org_mean_bic
      )
  }
  return(results);
}


