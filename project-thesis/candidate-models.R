candidateModels <- function(X, y, L, family) {
  require(glmnet);
  require(ncvreg);
  # SCAD
  A_scad <- ncvreg(X = X, y = y, family = family, warn = FALSE, penalty = "SCAD",  max.iter = 1e+04);
  # get rid of the intercept
  A_scad_beta <- as.matrix(A_scad$beta[-1,]);
  # MCP
  A_mcp <- ncvreg(X = X, y = y, family = family, warn = FALSE, penalty = "MCP",  max.iter = 1e+04);
  # get rid of the intercept
  A_mcp_beta <- as.matrix(A_mcp$beta[-1,]);
  # LASSO
  A_lasso <- glmnet(x = X, y = y, family = family, alpha = 1, nlambda = L, maxit = 1e+06);
  # the betas are in format  ("CsparseMatrix")., so we need to transform them
  A_lasso_beta <- as.matrix(A_lasso$beta);
  thelasso.cv<-cv.glmnet(X, y,family = family, alpha=1) ## first stage ridge
  bhat<-as.matrix(coef(thelasso.cv, s="lambda.1se"))[-1,1] ## coef() is a sparseMatrix
  if( all(bhat==0) ){
    ## if bhat is all zero then assign very close to zero weight to all.
    ## Amounts to penalizing all of the second stage to zero.
    bhat<-rep(0.0000000001,length(bhat))
  }
  adpen<-(1/pmax(abs(bhat))) ## the adaptive lasso weight
  ## Second stage lasso (the adaptive lasso)
  A_adaptive_lasso_beta <- glmnet(X, y,family = family, alpha=1, exclude=which(bhat==0), penalty.factor=adpen)
  A_adaptive_lasso_beta <- as.matrix(A_adaptive_lasso_beta$beta);
  # get final models: with A_lambda_l = supp(b_hat^lambda_l) = {j: beta_hat_j^lambda_l != 0}
  A_combined <- cbind(A_scad_beta, A_mcp_beta, A_lasso_beta, A_adaptive_lasso_beta);
  colnames(A_combined) <- NULL;
  # define an epsilon for the params
  # transform A_combined to a binary vector since later we have to compute the SOIL importance as:
  # S_j = S(j;w,A) = sum(w_k*I(j E A^k)) for k = 1, ... , K and j = 1, ... , p where K is the number
  # of candidate models and p is the number of variables
  eps <- 10e-2;
  A_combined <- apply(A_combined, 1:2, function(x) ifelse(abs(x) > eps, 1, 0));
  # remove duplicated columns
  candidate_models <- A_combined[ , !duplicated(t(A_combined))];
  candidate_models <- candidate_models[,-1];
  if (is.vector(candidate_models)){
    candidate_models <- as.matrix(candidate_models);
  }
  return( candidate_models );
}
