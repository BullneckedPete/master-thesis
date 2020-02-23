source("rep_soil.R");
require(SOIL);
load("Bardet.rda");

methods <- c("rep_arm", "rep_bic", "soil_arm", "soil_bic");
reps <- 100;


final_rep_arm <- 

for (method in methods) {
  
  if (method == "rep_arm") {
    soil_importance <- as.vector(rep_SOIL(X = X, y = y, weight_type = "ARM")$soil_importance);
  } else if (method == "rep_bic") {
    soil_importance <- as.vector(rep_SOIL(X = X, y = y, weight_type = "BIC")$soil_importance);
  } else if (method == "soil_arm") {
    soil_importance <- as.vector(SOIL(X, y, weight_type = "ARM", prior = TRUE)$importance);
  } else if (method == "soil_bic") {
    soil_importance <- as.vector(SOIL(X, y, weight_type = "BIC", prior = TRUE)$importance);
  }
  
  max_importance1 <- which(soil_importance == max(soil_importance));
  max_importance2 <- which(soil_importance == max(soil_importance[-max(soil_importance)]));
  
  rep_arm <- matrix(NA, ncol = nreps);
  rep_bic <- matrix(NA, ncol = nreps);
  soil_arm <- matrix(NA, ncol = nreps);
  soil_bic <- matrix(NA, ncol = nreps);
  
  for (i in 1:reps) {
    new_X <- X[ ,c(max_importance1, max_importance2)];
    reg_data <- as.data.frame(cbind(y, new_X));
    fit <- lm( y ~ . , data = reg_data );
    sigma_hat <- summary(fit)$sigma;
    coeff <-  fit$coefficients[-1];
    beta_hat <- matrix(coeff, ncol = 1 , nrow = length(coeff));
    y_new <- new_X %*% beta_hat + sigma_hat*rnorm(mean = 0, sd = 1);
    
    rep_arm <- as.vector(rep_SOIL(X = X, y = y, weight_type = "ARM")$soil_importance);
    rep_bic <- as.vector(rep_SOIL(X = X, y = y, weight_type = "BIC")$soil_importance);
    soil_arm <- as.vector(SOIL(X, y, weight_type = "ARM", prior = TRUE)$importance);
    soil_bic <- as.vector(SOIL(X, y, weight_type = "BIC", prior = TRUE)$importance);
  }
  
  
  
}
