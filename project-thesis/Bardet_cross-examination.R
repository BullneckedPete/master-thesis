source("candidate-models.R");
require(randomForest);
require(SOIL);
set.seed(456);
load("data/Bardet.rda");
Bardet <- as.data.frame(cbind(y,x));
attach(Bardet);
methods <- c("soil_arm", "soil_bic", "rf1", "rf2");
candidate_models <- candidateModels(x, y, family = "gaussian"); # A = { A_lasso, A_scad, A_mcp }
reps <- 100;
soil_arm <- matrix(NA, ncol = ncol(x), nrow = reps);
soil_bic <- matrix(NA, ncol = ncol(x), nrow = reps);
rf1 <- matrix(NA, ncol = ncol(x), nrow = reps);
rf2 <- matrix(NA, ncol = ncol(x), nrow = reps);
methods_top_variables_Bardet <- list("soil_arm"=c(), "soil_bic"=c(), "rf1"=c(), "rf2"=c());
importance_results_Bardet <- list();
for (method in methods) {
  print("BASE METHOD: ");
  print(method);
  if (method == "soil_arm") {
    imp_arm <- SOIL(x, y, weight_type = "ARM", method = "customize", candidate_models = t(candidate_models))$importance;
    top_importance <- imp_arm[ ,order(imp_arm[1,], decreasing = TRUE)][1:4];
    methods_top_variables_Bardet$soil_arm <- names(top_importance);
  } else if (method == "soil_bic") {
    imp_bic <- SOIL(x, y, weight_type = "BIC", method = "customize", candidate_models = t(candidate_models))$importance;
    top_importance <- imp_bic[ ,order(imp_bic[1,], decreasing = TRUE)][1:4];
    methods_top_variables_Bardet$soil_bic <- names(top_importance);
  } else if (method == "rf1") {
    rf_Bardet <- randomForest(y ~ ., data=Bardet, ntree=1000, keep.forest=FALSE, importance=TRUE);
    imp_rf1 <- importance(rf_Bardet, type = 1);
    top_importance <- imp_rf1[order(imp_rf1[,1], decreasing = TRUE), ][1:10];
    methods_top_variables_Bardet$rf1 <- names(top_importance);
  } else if (method == "rf2") {
    rf_Bardet <- randomForest(y ~ ., data=Bardet, ntree=1000, keep.forest=FALSE, importance=TRUE);
    imp_rf2 <- importance(rf_Bardet, type = 2);
    top_importance <- imp_rf2[order(imp_rf2[,1], decreasing = TRUE), ][1:10];
    methods_top_variables_Bardet$rf2 <- names(top_importance);
  }
  
  top_variables <- names(top_importance);
  for(i in 1:reps) {
    new_x <- x[ , top_variables];
    reg_data <- as.data.frame(cbind(y, new_x));
    fit <- lm( y ~ . , data = reg_data );
    sigma_hat <- summary(fit)$sigma;
    beta_hat <- matrix(fit$coefficients, ncol = 1 , nrow = length(fit$coefficients));
    y_new <- cbind(1,new_x) %*% beta_hat + sigma_hat*rnorm(n = 1, mean = 0, sd = 1);
    new_Bardet <- as.data.frame(cbind(y_new, x));
    colnames(new_Bardet)[1] <- "y_new"; 
    soil_arm[i, ] <- as.vector(SOIL(x, y_new, weight_type = "ARM", method = "customize", candidate_models = t(candidate_models))$importance);
    soil_bic[i, ] <- as.vector(SOIL(x, y_new, weight_type = "BIC", method = "customize", candidate_models = t(candidate_models))$importance);
    rf1[i, ] <- as.vector(importance(randomForest(y_new ~ ., data=new_Bardet, ntree=1000, keep.forest=FALSE, importance=TRUE), type = 1));
    rf2[i, ] <- as.vector(importance(randomForest(y_new ~ ., data=new_Bardet, ntree=1000, keep.forest=FALSE, importance=TRUE), type = 2));
  }
  mean_soil_arm <- colMeans(soil_arm);
  mean_soil_bic <- colMeans(soil_bic);
  mean_rf1 <- colMeans(rf1);
  mean_rf2 <- colMeans(rf2);
  names(mean_soil_arm) <- names(mean_soil_bic) <- names(mean_rf1) <- names(mean_rf2) <- colnames(x);
  
  importance_results_Bardet[[paste("Base measure: ", method)]] <- list (
    "soil_arm" = sort(mean_soil_arm, decreasing = TRUE),
    "soil_bic" = sort(mean_soil_bic, decreasing = TRUE),
    "rf1" = sort(mean_rf1, decreasing = TRUE),
    "rf2" = sort(mean_rf2, decreasing = TRUE)
  );
  
}
save(importance_results_Bardet, file = "cross-examination_results/importance_results_Bardet.RData")
save(methods_top_variables_Bardet, file = "cross-examination_results/methods_top_variables_Bardet.RData")