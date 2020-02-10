source("simulation_linear-regression_charts.R");
source("simulation_logistic-regression_charts.R");

#### Example 1: Linear regression ####
n <- 100;
betas <-  c(4, 4, 4, -6*sqrt(2), 4/3, rep(0, times = 100));
rho <- c(0, 0.9, 0, 0.9);
sd <- c(0.1, 0.1, 5, 5);
nreps <- 10;
plotSimulationLinearRegression(betas = betas, n = n, 
                      rho = rho, sd = sd, nreps = nreps, limitAxis = 20);


#### Example 2: Linear Regression ####
n <- 150;
betas <-  c(4, 4, 4, -6*sqrt(2), 4/3, rep(0, times = 9));
rho <- c(0, 0.9, 0, 0.9);
sd <- c(0.1, 0.1, 5, 5);
nreps <- 10;
plotSimulationLinearRegression(betas = betas, n = n, 
                               rho = rho, sd = sd, nreps = nreps, limitAxis = 15, noise = TRUE);

#### Example 3: Linear Regression ####
n <- 150;
betas <-  c(rep(1, times = 8));
rho <- c(0, 0.9, 0, 0.9);
sd <- c(0.1, 0.1, 5, 5);
nreps <- 10;

plotSimulationLinearRegression(betas = betas, n = n, 
                               rho = rho, sd = sd, nreps = nreps, limitAxis = 8);


#### Example 4: Logistic regression ####
n <- 80;
betas <- c(1, 1/2, 1/3, 1/4, 1/5, 1/6, 0);
rho <- c(0, 0.5, 0.9);
sd <- c(NA, NA, NA);
nreps <- 10;
plotSimulationLogisticRegression(betas = betas, n = n, rho = rho, sd = sd, nreps= nreps, limitAxis = 7) 


#### Example 5: Logistic regression ####
n <- 100;
betas <- c(4, 4, 4, -6*sqrt(2), 4/3, rep(0, times = 195));
rho <- c(0, 0.5, 0.9);
sd <- c(NA, NA, NA);
nreps <- 10;
plotSimulationLogisticRegression(betas = betas, n = n, rho = rho, sd = sd, nreps= nreps, limitAxis = 20) 
