source("c-functions.R")

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
  beta_2 <- -0.1
  
  
  prob_a <- 0.2
  
  df_names <- c("userid", "day", "Y", "A", "S", "S2", "prob_Y", "prob_Y_A0", "prob_A")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  group = group_str(group_ls)
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  dta$baseline_group_err = rep(group[["group err"]], each = total_T)
  dta$treatment_group_err = rep(group[["group b_g"]], each = total_T)
  dta$total_effect = rep(NA,sample_size * total_T)
  
  
  groupid <- group[["group_id"]]
  dta$groupid = rep(groupid,each = total_T)
  groupsize <- group[["group size"]]
  
  A_tj = matrix(NA, nrow = sample_size * total_T, ncol = unique(group[["group size"]]-1))
  userid_prime = matrix(NA, nrow = sample_size * total_T, ncol = unique(group[["group size"]]-1))
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    
    # dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
    dta$S[row_index] = if(t>1){
      state_transition(dta$S[row_index-1])
    }else{
      sample(c(0,1,2), sample_size, replace = TRUE)
    }
    
    #### group level moderator
    S_bar = aggregate(x = dta$S[row_index],                # Specify data column
                      by = list(groupid),              # Specify group indicator
                      FUN = mean)[,2]
    
    S_bar_t = rep(S_bar,groupsize)
    
    
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0) 
    dta$prob_A[row_index] <- rep(prob_a, sample_size)
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
    #baseline (=when A=0)
    dta$prob_Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                       ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
    ### indirect effect
    
    totaleff = rep(NA,sample_size)
    
    for (i in 1:sample_size) {
      g = groupid[i]
      g_index = which(groupid==g)
      other = setdiff(g_index, i)
      g_row_index = row_index[other]
      
      totaleff[i] = sum(dta$A[g_row_index])*(beta_2)
      
      A_tj[(i-1)*total_T+t,] = dta$A[g_row_index]
      userid_prime[(i-1)*total_T+t,] = other
    }
    dta$total_effect[row_index] <- totaleff
    # E(Y|H,A)
    dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * 
                                                              (beta_0+ beta_1*S_bar_t +dta$treatment_group_err[row_index])+ 
                                                              dta$baseline_group_err[row_index]+dta$total_effect[row_index])
    
    dta$prob_Y[row_index] = dta$prob_Y[row_index] /(prob_a*exp(beta_2)+(1-prob_a))^(unique(group[["group size"]]-2))
    
    dta$prob_Y[row_index] = pmin(1,dta$prob_Y[row_index])
    # sample
    dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
  }
  dta =dta[rep(seq_len(nrow(dta)),each = unique(group[["group size"]]-1)), ]
  dta$IE = (as.vector(t(A_tj))-prob_a)*(1-dta$A)
  dta$userid_prime = as.vector(t(userid_prime))
  dta$pairid = paste(dta$userid,dta$userid_prime, sep="_")
  
  return(dta= dta[order(dta$pairid),])
}


## true beta for (Intercept, S)
beta_true <- c(0.1, 0.3)

## true beta for Intercept
beta_true_marginal <- -0.1

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