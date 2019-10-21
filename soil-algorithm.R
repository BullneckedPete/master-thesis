#### Generating data ####

# number of observations
n = 100
# number of dimensions
p = 5

# generate the (p x p ) covariance matrix
rho <- 0.9
genCovMAt <- function(rho, p) {
  Sigma <- matrix(rep(0, times = p*p), nrow = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i,j] = rho ^ abs(i-j)
    }
  }
  return(Sigma)
}
Sigma <- genCovMAt(rho = rho, p = p)
Sigma

# generate sample data from multivariate normal distribution using cholesky decomposition
rMvNorm <- function(n, p, Sigma) {
  mu <- rep(0, times = p)
  Z <- matrix(rnorm(p*n, mean = 0, sd = 1), nrow = p, ncol = n)
  L <- t(chol(Sigma))
  X <- mu + L %*% Z
  return(X)
}
X <- rMvNorm(n = n, p = p, Sigma = Sigma) 
beta_star <- matrix(c(2,1,0,0,0), ncol = 1) 
beta_star

genY <- function(X, beta_star, sd, n) {
  e <- rnorm(n, mean = 0, sd = sd)
  return( t(X) %*% beta_star + e )
}
sd <- 0.1
y <- genY(X = X, beta_star = beta_star, sd = sd, n = n)
y


#### candidate models ####
X = t(X)
X
dim(X) 
dim(y)

require(glmnet)
require(ncvreg)
L = 100 # number of different regularization parameters
# SCAD
A_scad <- ncvreg(X = X, y = y, family = "gaussian", alpha = 1, penalty = "SCAD", nlambda = L )
# get rid of the intercept
A_scad_beta <- as.matrix(A_scad$beta[-1,])
A_scad_beta 
# MCP
A_mcp <- ncvreg(X = X, y = y, family = "gaussian", alpha = 1, penalty = "MCP", nlambda = L )
# get rid of the intercept
A_mcp_beta <- as.matrix(A_scad$beta[-1,])
A_mcp_beta
# LASSO
A_lasso <- glmnet(x = X, y = y, family = "gaussian", alpha = 1, nlambda = L)
# the betas are in format  ("CsparseMatrix")., so we need to transform them
A_lasso_beta <- as.matrix(A_lasso$beta)
A_lasso_beta


# get final models: with A_lambda_l = supp(b_hat^lambda_l) = {j: beta_hat_j^lambda_l != 0}
A_combined <- cbind(A_scad_beta, A_mcp_beta, A_lasso_beta)
A_combined
colnames(A_combined) <- NULL
# transform A_combined to a binary vector since later we have to compute the SOIL importance as:
# S_j = S(j;w,A) = sum(w_k*I(j E A^k)) for k = 1, ... , K and j = 1, ... , p where K is the number
# of candidate models and p is the number of variables
A_combined <- apply(A_combined, 1:2, function(x) ifelse(x == 0, 0, 1))
# remove duplicated columns
candidate_models <- A_combined[ , !duplicated(t(A_combined))]
candidate_models <- candidate_models[,-1]
print(candidate_models)


#### weighting using ARM with non uniform priors ####
nsim <- 1000; #number of simulations
sim_soil_importance <- matrix(0, ncol = nrow(candidate_models), nrow = nsim);
for (j in 1:nsim) {
  # create a training and a test set of equal size
  # round up, if nrow(X) is odd
  cut <- ceiling((nrow(X)/2));
  # randomly change rows
  randomSample <- sample(1:nrow(X), nrow(X), replace = FALSE);
  X <- X[randomSample, ];
  y <- y[randomSample];
  # training set
  D1 <- X[1:cut, ];
  y1 <- y[1:cut];
  # test set
  D2 <- X[(cut+1):nrow(X), ]; 
  y2 <- y[(cut+1):length(y)];
  
  # s_k: number of non constant predictors for model k, that is |A^k|
  s_k <- colSums(candidate_models);
  # phi: a positive number to control the improvement of the prior weight (defaults to 1)
  phi <- 1;
  k <- ncol(candidate_models);
  w_k_nominator <- numeric(k);
  for (i in 1:k) {
    
    # get the indices of the model to select the right positions of the design matrix
    indices <- as.vector(which(candidate_models[,i] != 0));
    
    # extract the right indices form the candidate model in the current iteration
    # out of the design matrix
    Xs_k <- D1[,indices];
    
    # prepare the data for the linear regression
    reg_data <- as.data.frame(cbind(y1, Xs_k));
    
    # fit standard linear regression of y on Xs_k using the training set D1
    fit_lm_k <- lm( y1 ~ . , reg_data );
    
    # get estimated standard deviation and coefficients
    sigma_hat <- summary(fit_lm_k)$sigma;
    coeff <-  fit_lm_k$coefficients[-1];
    bs_hat_k <- matrix(coeff, ncol = 1 , nrow = length(coeff));
    
    # prepare the prediction data and compute the prediction t(Xs_k) %*% bs_hat on the test set D2
    prediction_data <- as.matrix(D2[,indices]);
    y2_prediction <- prediction_data %*% bs_hat_k;
    
    # compute the constant c_k used to compute the weight w_k later
    c_k <- s_k[i] * log( exp(1) * p / s_k[i] ) + 2 * log(s_k[i] + 2);
    
    # compute the nominator of the weight vector w_k for each candidate model
    w_k_nominator[i] <- exp(1) ^ (-phi*c_k) * (sigma_hat^(-n/2)) * prod(exp ( -sigma_hat^-2 * (y2 - as.vector(y2_prediction))^2 / 2 ));
  }
  w_k_nominator[!is.finite(w_k_nominator)] <- 0;
  weight_vector <- w_k_nominator / sum(w_k_nominator);
  soil_importance <- weight_vector %*% t(candidate_models);
  sim_soil_importance[j, ] <- soil_importance;
}

final_weight_vector <- round(weight_vector, 2);
print(final_weight_vector);
final_soil_importance <- colSums(sim_soil_importance)/nsim;
print(final_soil_importance);
print(candidate_models);


# compare with the implementation from the autor, see example from official package manual
library(SOIL)
v_ARM <- SOIL(X, y, family = "gaussian",
              weight_type = "ARM", prior = TRUE)
v_ARM

















