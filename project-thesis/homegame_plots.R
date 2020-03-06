load("cross-examination_results/importance_results_Bardet.RData");
load("cross-examination_results/methods_top_variables_Bardet.RData");
load("cross-examination_results/importance_results_lung_cancer.RData");
load("cross-examination_results/methods_top_variables_lung_cancer.RData");

plot_homgame_importance <- function(baseimp, topvariables, title) {
  base_arm <- baseimp$soil_arm[c(topvariables)];
  max_base_arm <- tail(baseimp$soil_arm, n = 1);
  base_arm <- c(base_arm, max_base_arm);
  
  base_bic <- baseimp$soil_bic[c(topvariables)];
  max_base_bic <- tail(baseimp$soil_bic, n = 1);
  base_bic <- c(base_bic, max_base_bic);
  
  base_rf1 <- baseimp$rf1[c(topvariables)];
  max_base_rf1 <- tail(baseimp$rf1, n = 1);
  base_rf1 <- c(base_rf1, max_base_rf1);
  base_rf1 <- base_rf1 / sqrt(sum(base_rf1^2));
  
  base_rf2 <- baseimp$rf2[c(topvariables)];
  max_base_rf2 <- tail(baseimp$rf2, n = 1);
  base_rf2 <- c(base_rf2, max_base_rf2);
  base_rf2 <- base_rf2 / sqrt(sum(base_rf2^2));
  
  names(base_arm)[length(base_arm)] <- names(base_bic)[length(base_bic)] <- names(base_rf1)[length(base_rf1)] <- names(base_rf2)[length(base_rf2)] <- "max";
  max <- -1e10;
  min <- 1e10;
  for(i in 1:length(base_arm)) {
    if (base_arm[i] > max) max <- base_arm[i];
    if (base_arm[i] < min) min <- base_arm[i];
    if (base_bic[i] > max) max <- base_bic[i];
    if (base_bic[i] < min) min <- base_bic[i];
    if (base_rf1[i] > max) max <- base_rf1[i];
    if (base_rf1[i] < min) min <- base_rf1[i];
    if (base_rf2[i] > max) max <- base_rf2[i];
    if (base_rf2[i] < min) min <- base_rf2[i];
  }
  plot(base_arm, type = "b", ylim = c(min,max+0.5), xaxt="n", ylab = "Importance", xlab = "Variable Name",pch = 15, main = title);
  points(base_bic, type = "b", col = "red", pch = 16);
  points(base_rf1, type = "b", col = "blue", pch = 17);
  points(base_rf2, type = "b", col = "forestgreen", pch = 18);
  axis(1, at=1:length(base_arm), labels=names(base_arm)); 
  legend("topright", legend=c("SOIL-ARM", "SOIL-BIC", "RF1", "RF2"), col=c("black", "red", "blue", "forestgreen"), pch = 15:18, ncol=2); 
}

#### Bardet data ####
## Base method: SOIL ARM
topvariables_arm <- methods_top_variables_Bardet$soil_arm
importance_base_arm <- importance_results_Bardet$`Base measure:  soil_arm`
plot_homgame_importance(baseimp = importance_base_arm, topvariables = topvariables_arm, title="Homegame ARM")

# Base method: SOIL BIC
topvariables_bic <- methods_top_variables_Bardet$soil_bic
importance_base_bic <- importance_results_Bardet$`Base measure:  soil_bic`
plot_homgame_importance(baseimp = importance_base_bic, topvariables = topvariables_bic, title="Homegame BIC")

## Base method: RF1
topvariables_rf1 <- methods_top_variables_Bardet$rf1
importance_base_rf1 <- importance_results_Bardet$`Base measure:  rf1`
plot_homgame_importance(baseimp = importance_base_rf1, topvariables = topvariables_rf1, title="Homegame RF1")

## Base method: RF2
topvariables_rf2 <- methods_top_variables_Bardet$rf2
importance_base_rf2 <- importance_results_Bardet$`Base measure:  rf2`
plot_homgame_importance(baseimp = importance_base_rf2, topvariables = topvariables_rf2, title="Homegame RF2")

#### Lung cancer data ####
## Base method: SOIL ARM
topvariables_arm <- methods_top_variables_lung_cancer$soil_arm
importance_base_arm <- importance_results_lung_cancer$`Base measure:  soil_arm`
plot_homgame_importance(baseimp = importance_base_arm, topvariables = topvariables_arm, title="Homegame ARM")

# Base method: SOIL BIC
topvariables_bic <- methods_top_variables_lung_cancer$soil_bic
importance_base_bic <- importance_results_lung_cancer$`Base measure:  soil_bic`
plot_homgame_importance(baseimp = importance_base_bic, topvariables = topvariables_bic, title="Homegame BIC")

## Base method: RF1
topvariables_rf1 <- methods_top_variables_lung_cancer$rf1
importance_base_rf1 <- importance_results_lung_cancer$`Base measure:  rf1`
plot_homgame_importance(baseimp = importance_base_rf1, topvariables = topvariables_rf1, title="Homegame RF1")

## Base method: RF2
topvariables_rf2 <- methods_top_variables_lung_cancer$rf2
importance_base_rf2 <- importance_results_lung_cancer$`Base measure:  rf2`
plot_homgame_importance(baseimp = importance_base_rf2, topvariables = topvariables_rf2, title="Homegame RF2")

