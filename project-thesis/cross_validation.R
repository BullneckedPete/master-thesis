source("candidate-models.R");
source("rep_soil.R");
require(SOIL);

set.seed(123);

#### Linear Regression: BGS boys ####
load("data/BGSboys.rda");
y <- as.vector(BGSboys$HT18);
x <- as.matrix(BGSboys[, -c(8)]);
candidate_models <- candidateModels(x, y, family = "gaussian");

# ARM
BGS_rep_arm <- round(rep_SOIL(X = x, y = y, weight_type = "ARM")$soil_importance, 3);
BGS_soil_arm <- round(SOIL(x, y, weight_type = "ARM", method = "customize", candidate_models = t(candidate_models))$importance, 3);

# BIC
BGS_rep_bic <- round(rep_SOIL(X = x, y = y, weight_type = "BIC")$soil_importance, 3);
BGS_soil_bic <- round(SOIL(x, y, weight_type = "BIC", method = "customize", candidate_models = t(candidate_models))$importance, 3);

BGS_rep_arm <- matrix(BGS_rep_arm, nrow = 1);
BGS_soil_arm<- matrix(BGS_soil_arm, nrow = 1);
BGS_rep_bic <- matrix(BGS_rep_bic, nrow = 1);
BGS_soil_bic <- matrix(BGS_soil_bic, nrow = 1);
colnames(BGS_rep_arm) <- colnames(BGS_soil_arm) <- colnames(BGS_rep_bic) <- colnames(BGS_soil_bic) <- colnames(x);
sort(BGS_rep_arm[ ,which(as.vector(BGS_rep_arm) > 0)], decreasing = TRUE);
sort(BGS_soil_arm[ ,which(as.vector(BGS_soil_arm) > 0)], decreasing = TRUE);
sort(BGS_rep_bic[ ,which(as.vector(BGS_rep_bic) > 0)], decreasing = TRUE);
sort(BGS_soil_bic[ ,which(as.vector(BGS_soil_bic) > 0)], decreasing = TRUE);



#### Linear Regression: Bardet ####
load("data/Bardet.rda");
candidate_models <- candidateModels(x, y, family = "gaussian");
# ARM
BARDET_rep_arm <- round(rep_SOIL(X = x, y = y, weight_type = "ARM")$soil_importance, 3);
BARDET_soil_arm <- round(SOIL(x, y, weight_type = "ARM", method = "customize", candidate_models = t(candidate_models))$importance, 3);

# BIC
BARDET_rep_bic <- round(rep_SOIL(X = x, y = y, weight_type = "BIC")$soil_importance, 3);
BARDET_soil_bic <- round(SOIL(x, y, weight_type = "BIC", method = "customize", candidate_models = t(candidate_models))$importance, 3);

BARDET_rep_arm <- matrix(BARDET_rep_arm, nrow = 1);
BARDET_soil_arm <- matrix(BARDET_soil_arm, nrow = 1);
BARDET_rep_bic <- matrix(BARDET_rep_bic, nrow = 1);
BARDET_soil_bic <- matrix(BARDET_soil_bic, nrow = 1);
colnames(BARDET_rep_arm) <- colnames(BARDET_soil_arm) <- colnames(BARDET_rep_bic) <- colnames(BARDET_soil_bic) <- colnames(x);
sort(BARDET_rep_arm[ ,which(as.vector(BARDET_rep_arm) > 0)], decreasing = TRUE);
sort(BARDET_soil_arm[ ,which(as.vector(BARDET_soil_arm) > 0)], decreasing = TRUE);
sort(BARDET_rep_bic[ ,which(as.vector(BARDET_rep_bic) > 0)], decreasing = TRUE);
sort(BARDET_soil_bic[ ,which(as.vector(BARDET_soil_bic) > 0)], decreasing = TRUE);



#### Logistic regression: Lung cancer ####
load("data/Lung.rda");
x <- t(Lung_Boston_x[,-c(1,2)]);
gene_names <- Lung_Boston_x$NAME;
colnames(x) <- gene_names;
y <- numeric(length(Lung_Boston_y));
for (i in 1:length(y)) {
  ifelse(Lung_Boston_y[i] == "A", y[i] <- 0, y[i] <- 1);
}
candidate_models <- candidateModels(x, y, family = "binomial");

#ARM
LUNG_rep_arm <- round(rep_SOIL(X = x, y = y, weight_type = "ARM", family = "binomial")$soil_importance, 3);
LUNG_soil_arm <- round(SOIL(x, y, weight_type = "ARM", method = "customize", family = "binomial", candidate_models = t(candidate_models))$importance, 3);

# BIC
LUNG_rep_bic <- round(rep_SOIL(X = x, y = y, weight_type = "BIC", family = "binomial")$soil_importance, 3);
LUNG_soil_bic <- round(SOIL(x, y, weight_type = "BIC", method = "customize", family = "binomial", candidate_models = t(candidate_models))$importance, 3);

LUNG_rep_arm <- matrix(LUNG_rep_arm, nrow = 1);
LUNG_soil_arm <- matrix(LUNG_soil_arm, nrow = 1);
LUNG_rep_bic <- matrix(LUNG_rep_bic, nrow = 1);
LUNG_soil_bic <- matrix(LUNG_soil_bic, nrow = 1);
colnames(LUNG_rep_arm) <- colnames(LUNG_soil_arm) <- colnames(LUNG_rep_bic) <- colnames(LUNG_soil_bic) <- colnames(x);
sort(LUNG_rep_arm[ ,which(as.vector(LUNG_rep_arm) > 0)], decreasing = TRUE);
sort(LUNG_soil_arm[ ,which(as.vector(LUNG_soil_arm) > 0)], decreasing = TRUE);
sort(LUNG_rep_bic[ ,which(as.vector(LUNG_rep_bic) > 0)], decreasing = TRUE);
sort(LUNG_soil_bic[ ,which(as.vector(LUNG_soil_bic) > 0)], decreasing = TRUE);



#### Cross examination procedure for linear regression ####
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
    soil_importance <- as.vector(SOIL(x, y, weight_type = "ARM", method = "union")$importance);
  } else if (method == "soil_bic") {
    soil_importance <- as.vector(SOIL(x, y, weight_type = "BIC", method = "union")$importance);
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
    soil_arm[i, ] <- as.vector(SOIL(x, y_new, weight_type = "ARM", method = "union")$importance);
    soil_bic[i, ] <- as.vector(SOIL(x, y_new, weight_type = "BIC", method = "union")$importance);
  }
  results[[paste("Base measure: ", method)]] <- list (
    "rep arm" = colMeans(rep_arm),
    "soil arm" = colMeans(soil_arm),
    "rep bic" = colMeans(rep_bic),
    "soil bic" = colMeans(soil_bic)
  )
  
}
results

