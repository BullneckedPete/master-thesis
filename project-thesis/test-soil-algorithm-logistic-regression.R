source("generate-data.R");
source("candidate-models.R");
source("ARM-logistic-regression.R");

#### Generating data ####
# number of observations
n <- 100;
# number of dimensions
p <- 5;
rho <- 0.9;
sd <- 0.1;
beta_star <- matrix(c(4,0,0,0,0), ncol = 1);
data <- generateData(n = n, p = p, rho = rho, sd = sd, beta_star = beta_star, type = "logistic");
X <- data$X;
y <- data$y;
y
#### Generating the candidate models ####
L = 100 # number of different regularization parameters
candidate_models <- candidateModels(X = X, y = y, L = L, family = "binomial");
print(candidate_models);

#### Weighting using ARM with non uniform priors ####
nsim = 1000;
logRegARM <- logRegARM(X = X, y = y, nsim = nsim, candidate_models = candidate_models);
logRegARM

library(SOIL);
v_ARM <- SOIL(X, y, family = "binomial", weight_type = "ARM", prior = TRUE);
v_ARM