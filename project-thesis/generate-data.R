
# generate the (p x p) covariance matrix
generateCovarianceMatrix <- function(rho, p) {
  Sigma <- matrix(rep(0, times = p*p), nrow = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i,j] = rho ^ abs(i-j);
    }
  }
  return(Sigma);
}

# generate sample data from multivariate normal distribution using cholesky decomposition
rMvNorm <- function(n, p, Sigma) {
  mu <- rep(0, times = p);
  Z <- matrix(rnorm(p*n, mean = 0, sd = 1), nrow = p, ncol = n);
  L <- t(chol(Sigma));
  X <- mu + L %*% Z;
  return(X);
}

genY <- function(X, beta_star, sd, n) {
  e <- rnorm(n, mean = 0, sd = sd);
  return( t(X) %*% beta_star + e );
}

generateData <- function(n, p, rho, sd, beta_star) {
  Sigma <- generateCovarianceMatrix(rho = rho, p = p);
  X <- rMvNorm(n = n, p = p, Sigma = Sigma);
  y <- genY(X = X, beta_star = beta_star, sd = sd, n = n);
  return ( list(Sigma = Sigma, X = t(X), y = y) );
}
