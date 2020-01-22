source("generate-data.R");
source("candidate-models.R");
source("ARM-linear-regression.R");
source("BIC-linear-regression.R");
require("SOIL");

#### Generating data ####
# number of observations
n <- 1000;
# number of dimensions
p <- 5;
rho <- 0.9;
sd <- 0.1;
beta_star <- matrix(c(4,-4,0,0,0), ncol = 1);
data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta_star, family = "gaussian");
X <- data$X;
y <- data$y;

#### Generating the candidate models ####
L <- 100 # number of different regularization parameters
candidate_models <- candidateModels(X = X, y = y, L = L, family = "gaussian");
print(candidate_models);

#### Weighting using ARM with non uniform priors ####
nsim <- 100;
linear_regression_arm <- linearRegressionARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
linear_regression_arm

# Compare with the implementation from the autor, see example from official package manual
v_ARM <- SOIL(X, y, family = "gaussian", weight_type = "ARM", prior = TRUE);
v_ARM


#### Weighting using information criteria with nonuniform priors ####
linear_regression_bic <- linearRegressionBIC(X = X, y = y, candidate_models = candidate_models);
linear_regression_bic

v_BIC <- SOIL(X, y, family = "gaussian", weight_type = "BIC", prior = TRUE);
v_BIC
