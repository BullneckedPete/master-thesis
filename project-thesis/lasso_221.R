source("generate-data.R");
require(glmnet);
require(ncvreg);
require(plotmo);
set.seed(123)
n <- 1000;
p <- 7;
rho <- 0.9;
sd <- 0.1;
beta_star <- matrix(c(4,-4*sqrt(2),1, 1/2, 1/10, 0, 0), ncol = 1);
data <- generateData(n = n, p = p, rho = rho, sd = sd, 
            beta_star = beta_star, family = "gaussian");
X <- data$X;
y <- data$y;

# least absolute shrinkage and selection operator (LASSO)
lasso <- glmnet(x = X, y = y, family = "gaussian", alpha = 1);
plot_glmnet(lasso,  xvar = "norm");
lasso_coeff <- as.matrix(lasso$beta)
lasso_coeff_cleaned <- apply(lasso_coeff, 1:2, 
                          function(x) ifelse(abs(x) > 0, 1, 0));
lasso_models <- lasso_coeff_cleaned[ , !duplicated(t(lasso_coeff_cleaned))];
lasso_models


# smoothly clipped absolute deviation penalty (SCAD)
scad <- ncvreg(X = X, y = y, family = "gaussian", penalty = "SCAD");
scad_coeff <- as.matrix(scad$$beta[-1,])

##c <- citation(package = "glmnet");
##print(c, bibtex = TRUE)
