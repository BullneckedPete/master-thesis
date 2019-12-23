candidateModels <- function(X, y, L) {
  require(glmnet);
  require(ncvreg);
  
  # SCAD
  A_scad <- ncvreg(X = X, y = y, family = "gaussian", alpha = 1, penalty = "SCAD", nlambda = L);
  # get rid of the intercept
  A_scad_beta <- as.matrix(A_scad$beta[-1,]);
  
  # MCP
  A_mcp <- ncvreg(X = X, y = y, family = "gaussian", alpha = 1, penalty = "MCP", nlambda = L);
  # get rid of the intercept
  A_mcp_beta <- as.matrix(A_scad$beta[-1,]);
  
  # LASSO
  A_lasso <- glmnet(x = X, y = y, family = "gaussian", alpha = 1, nlambda = L);
  # the betas are in format  ("CsparseMatrix")., so we need to transform them
  A_lasso_beta <- as.matrix(A_lasso$beta);
  
  # get final models: with A_lambda_l = supp(b_hat^lambda_l) = {j: beta_hat_j^lambda_l != 0}
  A_combined <- cbind(A_scad_beta, A_mcp_beta, A_lasso_beta);
  colnames(A_combined) <- NULL;
  
  # transform A_combined to a binary vector since later we have to compute the SOIL importance as:
  # S_j = S(j;w,A) = sum(w_k*I(j E A^k)) for k = 1, ... , K and j = 1, ... , p where K is the number
  # of candidate models and p is the number of variables
  A_combined <- apply(A_combined, 1:2, function(x) ifelse(x == 0, 0, 1));
  # remove duplicated columns
  candidate_models <- A_combined[ , !duplicated(t(A_combined))];
  candidate_models <- candidate_models[,-1];
  return( candidate_models );
}
