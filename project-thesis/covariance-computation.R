covariance.matrix <- function(designMatrix) {
  n <- nrow(M) 
  
  # create a matrix of column means to center the data
  M_mean <- matrix(data=1, nrow=n) %*% colMeans(M)
  
  # center the data
  D <- as.matrix(M - M_mean)
  
  #compute the covariance function
  C <- (n-1)^-1 * t(D) %*% D
  
  return(C)
}

M <- matrix(c(1,2,3,4), ncol = 2)

#  same function as the built in cov() function
covariance.matrix(M) == cov(M)
