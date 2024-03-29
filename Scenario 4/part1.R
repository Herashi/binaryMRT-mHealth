setwd("~/binary_mrt/Scenario 4")

source("estimators.R")
source("estimators_robust_adhocery.R") #didn't use the function in this file

source("data_generating.R")
data_generating_process <- dgm_binary_categorical_covariate

library(tidyverse)

# library(foreach)
# library(doMC)
# library(doRNG)

compute_result_beta <- function(beta_true, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level,
                                na.rm = FALSE) {
  
  beta_true_array <- matrix(NA, nrow = nrow(beta),ncol = length(beta_true))
  for (ind1 in 1:dim(beta_true_array)[1]) {
    beta_true_array[ind1,] <- beta_true
  }
  
  p <- length(moderator_vars) + 1
  q <- length(control_vars) + 1
  
  bias <- mean(beta - beta_true_array,na.rm = na.rm)
  sd <- sd(beta, na.rm = na.rm)
  rmse <- sqrt(mean((beta - beta_true_array)^2))
  
  critical_factor <- qnorm(1 - significance_level/2)
  ci_left <- beta - critical_factor * beta_se
  ci_right <- beta + critical_factor * beta_se
  beta_se_mean = mean(beta_se)
  coverage_prob <- mean((ci_left < beta_true_array) & (ci_right > beta_true_array))
  
  critical_factor_adj <- qt(1 - significance_level/2, df = sample_size - 1 - q)
  ci_left_adj <- beta - critical_factor_adj * beta_se_adjusted
  ci_right_adj <- beta + critical_factor_adj * beta_se_adjusted
  beta_se_adjusted_mean = mean(beta_se_adjusted)
  coverage_prob_adj <- mean((ci_left_adj < beta_true_array) & (ci_right_adj > beta_true_array))
  
  return(list(bias = bias,se=beta_se_mean, se_adjusted=beta_se_adjusted_mean,
              sd = sd, rmse = rmse, coverage_prob = coverage_prob, coverage_prob_adjusted = coverage_prob_adj))
}


total_T <- as.numeric(Sys.getenv("tmax"))

nsim <- 1

control_vars <- "S"
moderator_vars <- c()

i_ss = 1

group_ls = group_all[[i_ss]]
sample_size <- group_ls[["n"]] 

dta <- data_generating_process(sample_size, total_T,group_ls)

fit_wcls <- weighted_centered_least_square(
  dta = dta,
  group_ls=group_ls,
  id_varname = "pairid",
  decision_time_varname = "day",
  treatment_varname = "IE",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.2,
  estimator_initial_value = NULL
)
result  <- list(fit_wcls = fit_wcls)

ee_names <- "wcls"
alpha_names <- c("Intercept", control_vars)
beta_names <- c("Intercept", moderator_vars)
num_estimator <- length(ee_names)


#########################

alpha <- matrix(sapply(result, function(l) l$alpha_hat), byrow = TRUE,nrow =nsim)
alpha_se <- matrix(sapply(result, function(l) l$alpha_se),byrow = TRUE,nrow =nsim)
alpha_se_adjusted <- matrix(sapply(result, function(l) l$alpha_se_adjusted),byrow = TRUE, nrow = nsim)

colnames(alpha)= colnames(alpha_se)= colnames(alpha_se_adjusted) = alpha_names

beta <- matrix(sapply(result, function(l) l$beta_hat))
beta_se <- matrix(sapply(result, function(l) l$beta_se))
beta_se_adjusted <- matrix(sapply(result, function(l) l$beta_se_adjusted))

colnames(beta)= colnames(beta_se) = colnames(beta_se_adjusted)= beta_names

#################################

result <- compute_result_beta(beta_true_marginal, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level = 0.05)
result_df <- data.frame(ss = rep(sample_size, num_estimator),
                        est = ee_names,
                        bias = result$bias,
                        se=result$se,
                        se_adjusted = result$se_adjusted,
                        sd = result$sd,
                        rmse = result$rmse,
                        cp.unadj = result$coverage_prob,
                        cp.adj = result$coverage_prob_adjusted)
names(result_df) <- c("ss", "est", "bias","se.unadj","se.adj", "sd", "rmse", "cp.unadj", "cp.adj")
rownames(result_df) <- NULL

arrayid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
save(result_df, file=paste0("fity_", arrayid, "_tmax_", i_ss, ".RData"))

# stopCluster(cl)

