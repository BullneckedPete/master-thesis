source("generate-data.R");
source("candidate-models.R");
source("ARM-logistic-regression.R");
source("BIC-logistic-regression.R");
require("SOIL");

plotSimulationLogisticRegression <- function(betas, n, rho, sd, nreps, limitAxis) {
  par(mfrow=c(1,3));
  nsim <- 100;
  results <- list();
  for(j in 1:length(rho)) {
    beta <- matrix(betas, ncol = 1);
    p <- nrow(beta);
    soilArmLogisticRegression <- matrix(0, nrow = nreps, ncol = p);
    soilArmLogisticRegression_author <- matrix(0, nrow = nreps, ncol = p);
    soilBicLogisticRegression <- matrix(0, nrow = nreps, ncol = p);
    soilBicLogisticRegression_author <- matrix(0, nrow = nreps, ncol = p);
    
    for (i in 1:nreps) {
      # generate the data
      data <- generateData(n = n, p = p, rho = rho[j], sd = sd[j], beta_star = beta, family = "binomial");
      X <- data$X;
      y <- data$y;
      
      #obtain the candidate models
      candidate_models <- candidateModels(X = X, y = y, family = "binomial");
      
      # ARM
      v_ARM_logistic <- SOIL(X, y, family = "binomial", weight_type = "ARM", method = "customize", candidate_models = t(candidate_models));
      soilArmLogisticRegression_author[i, ] <- as.vector(v_ARM_logistic$importance);
      logistic_regression_arm <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
      soilArmLogisticRegression[i, ] <- as.vector(logistic_regression_arm$soil_importance);
      
      #BIC
      v_BIC_logistic <- SOIL(X, y, family = "binomial", weight_type = "BIC", method = "customize", candidate_models = t(candidate_models));
      soilBicLogisticRegression_author[i, ] <- as.vector(v_BIC_logistic$importance);
      logistic_regression_bic <- logisticRegressionBIC(X = X, y = y, candidate_models = candidate_models);
      soilBicLogisticRegression[i, ] <- as.vector(logistic_regression_bic$soil_importance);
    }
    
    rep_mean_arm <- colMeans(soilArmLogisticRegression);
    org_mean_arm <- colMeans(soilArmLogisticRegression_author);
    rep_mean_bic <- colMeans(soilBicLogisticRegression);
    org_mean_bic <- colMeans(soilBicLogisticRegression_author);
    
    # create the line charts
    title <- bquote( rho == .(rho[j]) );
    plot(rep_mean_arm, type = "b", ylim = c(0,1.2), xlim = c(1,limitAxis), 
         ylab = "Importance", xlab = "Variable Index",pch = 15, main = title);
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




