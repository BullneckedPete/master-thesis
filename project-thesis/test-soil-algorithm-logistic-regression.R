source("generate-data.R");
source("candidate-models.R");
source("ARM-logistic-regression.R");
source("BIC-logistic-regression.R");
require("SOIL");

#### Generating data ####
# number of observations
n <- 100;
# number of dimensions
p <- 10;
rho <- 0.9;
sd <- 0.1;
beta_star <- matrix(c(4, 4, 4, -6*sqrt(2), 4/3, rep(0, times = 5)), ncol = 1);
data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta_star, family = "binomial");
X <- data$X;
y <- data$y;

#### Generating the candidate models ####
L <- 100 # number of different regularization parameters
candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
print(candidate_models);


v_ARM <- SOIL(X, y, family = "binomial", weight_type = "ARM", prior = TRUE);
t(v_ARM$candidate_models_cleaned)[,-1]

#### Weighting using ARM with non uniform priors ####
nsim <- 100;
logistic_regression_arm <- logisticRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
logistic_regression_arm



#### Weighting using information criteria with nonuniform priors ####
logistic_regression_bic <- logisticRegressionBIC(X = X, y = y, candidate_models = candidate_models);
logistic_regression_bic

v_BIC <- SOIL(X, y, family = "binomial", weight_type = "BIC", prior = TRUE);
v_BIC
