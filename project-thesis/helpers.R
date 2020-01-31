selectionRejectionSample <- function(U, n) {
  N <- nrow(U);
  j <- 0;
  S <- numeric();
  for(i in 0:(N-1)) {
    u_i <- runif(n = 1, min = 0, max = 1);
    if (u_i < (n-j)/(N-i)) {
      S <- c(S,i);
      j = j+1;
    }
  }
  return(S+1);
}

calculateCk <- function(s_k, p) {
  if (s_k != 0) {
    C_k <-  2 * log(s_k+2) + s_k * log(exp(1) * p / s_k);
  } else {
    C_k <-  2 * log(s_k+2);
  }
  return(C_k);
}

checkIfWeightsValid <- function(weight_vectors) {
  return(
    as.matrix(weight_vectors[!rowSums(!is.finite(weight_vectors)), ])
  );
}
