candidateModels <- function(X, y, L = 100, family) {
  require(glmnet);
  require(ncvreg);
  A_scad <- ncvreg(X = X, y = y, family = family, 
                     warn = FALSE, penalty = "SCAD",  max.iter = 1e+05);
  A_scad_beta <- as.matrix(A_scad$beta[-1,]);
  A_mcp <- ncvreg(X = X, y = y, family = family, 
                    warn = FALSE, penalty = "MCP",  max.iter = 1e+05);
  A_mcp_beta <- as.matrix(A_mcp$beta[-1,]);
  A_lasso <- glmnet(x = X, y = y, family = family, 
                      alpha = 1, nlambda = L, maxit = 1e+05);
  A_lasso_beta <- as.matrix(A_lasso$beta);
  A_combined <- cbind(A_scad_beta, A_mcp_beta, A_lasso_beta);
  colnames(A_combined) <- NULL;
  # transform A_combined to a binary vector since later we have to compute
  # the SOIL importance as:
  # S_j = S(j;w,A) = sum(w_k*I(j E A^k)) for k = 1, ... , K and j = 1, ... , p 
  # where K is the number of candidate models and p is the number of variables
  eps <- 0;
  A_combined <- apply(A_combined, 1:2, function(x) ifelse(abs(x) > 0.001, 1, 0));
  candidate_models <- A_combined[ , !duplicated(t(A_combined))];
  candidate_models <- candidate_models[,-1];
  if (is.vector(candidate_models)){
    candidate_models <- as.matrix(candidate_models);
  }
  return( candidate_models );
}
