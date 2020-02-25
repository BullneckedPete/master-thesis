generateData <- function(n, p, rho, sd, beta_star, family, quad = FALSE) {
  
  generateCovarianceMatrix <- function(rho, p) {
    Sigma <- matrix(rep(0, times = p*p), nrow = p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i,j] = rho ^ abs(i-j);
      }
    }
    return(Sigma);
  }
  
  randomMvNorm <- function(n, p, Sigma) {
    mu <- rep(0, times = p);
    Z <- matrix(rnorm(p*n, mean = 0, sd = 1), nrow = p, ncol = n);
    L <- t(chol(Sigma));
    X <- mu + L %*% Z;
    return(X);
  }
  
  generateLinearY <- function(X, beta_star, sd, n) {
    e <- rnorm(n, mean = 0, sd = sd);
    return( t(X) %*% beta_star + e );
  }
  
  generateLogisticY <- function(X, beta_star) {
    y <- 1 / (1 + exp(-(t(X) %*% beta_star)));
    return ( sapply(y, function(x) ifelse(x < 0.5, 0, 1)) );
  }
  
  Sigma <- generateCovarianceMatrix(rho = rho, p = p);
  
  if (!quad){
    X <- randomMvNorm(n = n, p = p, Sigma = Sigma);
  } else {
    X <- randomMvNorm(n = n, p = 6, Sigma = Sigma);
    X <- rbind(X, X[1,]^2, X[2,]^2, X[3,]^2, X[4,]^2, X[5,]^2, X[6,]^2);
  }
  
  if (family == "gaussian") {
    y <- generateLinearY(X = X, beta_star = beta_star, sd = sd, n = n);
  } else if (family == "binomial") {
    y <- generateLogisticY(X =X, beta_star = beta_star);
  }
  
  return ( list(Sigma = Sigma, X = t(X), y = y) );
}
