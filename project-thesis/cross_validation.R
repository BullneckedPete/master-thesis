source("rep_soil.R");
require(SOIL);
#load("Bardet.rda");
load("BGSboys.rda")
BGSboys
y <- as.vector(BGSboys$HT18)
x <- as.matrix(BGSboys[, -c(8)])


# ARM
a <- round(as.matrix(rep_SOIL(X = x, y = y, weight_type = "ARM")$soil_importance, ncol=ncol(x)), 3)
c <- round(SOIL(x, y, weight_type = "ARM", prior = TRUE)$importance, 3)

# BIC
b <- round(rep_SOIL(X = x, y = y, weight_type = "BIC")$soil_importance, 3)
d <- round(SOIL(x, y, weight_type = "BIC", prior = TRUE)$importance, 3)

a <- matrix(a, nrow = 1)
c <- matrix(c, nrow = 1)
b <- matrix(b, nrow = 1)
d <- matrix(d, nrow = 1)
colnames(a) <- colnames(b) <- colnames(c) <- colnames(d) <- colnames(x)
sort(a[ ,which(as.vector(a) > 0)], decreasing = TRUE)
sort(c[ ,which(as.vector(c) > 0)], decreasing = TRUE)
sort(b[ ,which(as.vector(b) > 0)], decreasing = TRUE)
sort(d[ ,which(as.vector(d) > 0)], decreasing = TRUE)


## cross-examination
methods <- c("rep_arm", "soil_arm", "rep_bic", "soil_bic");
reps <- incr_reps <- 10;
rep_arm <- matrix(NA, ncol = ncol(x), nrow = reps);
rep_bic <- matrix(NA, ncol = ncol(x), nrow = reps);
soil_arm <- matrix(NA, ncol = ncol(x), nrow = reps);
soil_bic <- matrix(NA, ncol = ncol(x), nrow = reps);
results <- list();
for (method in methods) {
    print("BASE METHOD: ");
    print(method);
    if (method == "rep_arm") {
      soil_importance <- as.vector(rep_SOIL(X = x, y = y, weight_type = "ARM")$soil_importance);
    } else if (method == "rep_bic") {
      soil_importance <- as.vector(rep_SOIL(X = x, y = y, weight_type = "BIC")$soil_importance);
    } else if (method == "soil_arm") {
      soil_importance <- as.vector(SOIL(x, y, weight_type = "ARM", prior = TRUE)$importance);
    } else if (method == "soil_bic") {
      soil_importance <- as.vector(SOIL(x, y, weight_type = "BIC", prior = TRUE)$importance);
    }
  
    max_importance1 <- which(soil_importance == max(soil_importance));
    if(length(max_importance1) > 1) {
      max_importance2 <- max_importance1[2];
      max_importance1 <- max_importance1[1];
    } else {
      max_importance2 <- which(soil_importance == max(soil_importance[-max_importance1]))[1];
    }
    variable1 <- colnames(x)[max_importance1];
    variable2 <- colnames(x)[max_importance2];
    print(variable1); print(variable2);
    for(i in 1:reps) {
      new_x <- x[ ,c(max_importance1, max_importance2)];
      reg_data <- as.data.frame(cbind(y, new_x));
      fit <- lm( y ~ . , data = reg_data );
      sigma_hat <- summary(fit)$sigma;
      coeff <-  fit$coefficients;
      beta_hat <- matrix(coeff, ncol = 1 , nrow = length(coeff));
      y_new <- cbind(1,new_x) %*% beta_hat + sigma_hat*rnorm(n = 1, mean = 0, sd = 1);
      
      rep_arm[i, ] <- as.vector(rep_SOIL(X = x, y = y_new, weight_type = "ARM")$soil_importance);
      rep_bic[i, ] <- as.vector(rep_SOIL(X = x, y = y_new, weight_type = "BIC")$soil_importance);
      soil_arm[i, ] <- as.vector(SOIL(x, y_new, weight_type = "ARM", prior = TRUE)$importance);
      soil_bic[i, ] <- as.vector(SOIL(x, y_new, weight_type = "BIC", prior = TRUE)$importance);
  }
  results[[paste("Base measure: ", method)]] <- list (
    "rep arm" = colMeans(rep_arm),
    "soil arm" = colMeans(soil_arm),
    "rep bic" = colMeans(rep_bic),
    "soil bic" = colMeans(soil_bic)
  )
  
}
results



