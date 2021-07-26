# for binary outcome

# Tianchen Qian, 2018.07.24

# In this code, I constructed examples where:
# p_t is constant (0.5)
# S_t (covariate) takes 3 values {0,1,2}
# and the log linear GEE is biased for beta.

source("c-functions.R")
library(dplyr)

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}


state_transition = function(z){
  Z = c(0,1,2)
  a = rep(NA,length(z))
  a[which(z==0)]= sample(x=Z,length(which(z==0)),prob = c(2,1,1),replace = T)
  a[which(z==1)]= sample(x=Z,length(which(z==1)),prob = c(1,2,1),replace = T)
  a[which(z==2)]= sample(x=Z,length(which(z==2)),prob = c(1,1,2),replace = T)
  return(a)
}

dgm_binary_categorical_covariate <- function(sample_size, total_T,group_ls) {
  # same DGM as dgm_binary above, but faster
  
  baseline_Y_S0 <- 0.1
  baseline_Y_S1 <- 0.25
  baseline_Y_S2 <- 0.2
  
  beta_0 <- 0.1
  beta_1 <- 0.3
  
  
  prob_a <- 0.2
  
  df_names <- c("userid", "day", "Y", "A", "S", "S2", "prob_Y", "prob_Y_A0", "prob_A")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  group = group_str(group_ls)
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  dta$baseline_group_err = rep(group[["group err"]], each = total_T)
  dta$treatment_group_err = rep(group[["group b_g"]], each = total_T)
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    
    # dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
    dta$S[row_index] = if(t>1){
      state_transition(dta$S[row_index-1])
    }else{
      sample(c(0,1,2), sample_size, replace = TRUE)
    }
    
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0) 
    dta$prob_A[row_index] <- rep(prob_a, sample_size)
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
    #baseline (=when A=0)
    dta$prob_Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                       ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
    
    # E(Y|H,A)
    # +dta$treatment_group_err[row_index]
    #+ dta$baseline_group_err[row_index]
    dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * 
                                                              (beta_0 + beta_1 * dta$S[row_index]+dta$treatment_group_err[row_index])+ 
                                                              dta$baseline_group_err[row_index])
    dta$prob_Y[row_index] = pmin(1,dta$prob_Y[row_index])
    # sample
    dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
  }
  
  return(dta)
}

## true beta for (Intercept, S)
beta_true <- c(0.1, 0.3)

## true beta for Intercept
## log((0.1*exp(0.1)+0.25*exp(0.4)+0.2*exp(0.7))/(0.1+0.25+0.2))
beta_true_marginal <- 0.477

## true alpha for (Intercept, S, S2)
alpha_true <- c(-1.6094379,  0.9162907, -1.1394343)

# try out the range of Y
if (0) {
  set.seed(123)
  dta <- dgm_binary_categorical_covariate(100, 30)
  summary(dta$prob_Y)
  summary(dta$prob_A)
}


# compute marginal beta_true
if (0) {
  
  ### numerically, we have ###
  set.seed(123)
  sample_size <- 500000
  total_T <- 30
  dta <- dgm_binary_categorical_covariate(sample_size, total_T)
  
  tmp <- aggregate(Y ~ day + A, dta, mean)
  beta_true_marginal <- log(sum(tmp$Y[(total_T+1):(2*total_T)]) / sum(tmp$Y[1:total_T])) 
  print(beta_true_marginal)
  # p_a = 0.5, beta_true_marginal = 0.4777817
  # p_a = 0.2, beta_true_marginal = 0.4780562
  
  ### analytically, we have ###
  
  baseline_Y_S0 <- 0.1
  baseline_Y_S1 <- 0.25
  baseline_Y_S2 <- 0.2
  
  beta_0 <- 0.1
  beta_1 <- 0.3
  
  prob_a <- 0.2
  
  numerator <- baseline_Y_S0 * exp(beta_0) + baseline_Y_S1 * exp(beta_0 + beta_1) + baseline_Y_S2 * exp(beta_0 + 2 * beta_1)
  denominator <- baseline_Y_S0 + baseline_Y_S1 + baseline_Y_S2
  
  log(numerator / denominator)
}


# compute true alpha for (Intercept, S, S2)
if (0) {
  baseline_Y_S0 <- 0.1
  baseline_Y_S1 <- 0.25
  baseline_Y_S2 <- 0.2
  
  alpha_0 <- log(baseline_Y_S0)
  alpha_1 <- log(baseline_Y_S1 / baseline_Y_S0)
  alpha_2 <- log(baseline_Y_S2) - log(baseline_Y_S0) - 2 * log(baseline_Y_S1 / baseline_Y_S0)
  print(c(alpha_0, alpha_1, alpha_2))
}