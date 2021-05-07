library(truncnorm)

## specify group structure
group_all = list()
group_all[["250"]] = list()
group_all[["625"]] = list()
group_all[["1000"]] = list()

# maximum b_g sigma2 can go 0.75^2

group_all[["250"]][["n"]] = 250
group_all[["250"]][["group_id"]] = rep(1:25,each=10)
group_all[["250"]][["e_g sigma2"]] = rep(0,25)
group_all[["250"]][["b_g sigma2"]] = rep(0.5,25)
group_all[["250"]][["max_prob"]] = 0.2*exp(0.1+0.3*2)

group_all[["625"]][["n"]] = 625
group_all[["625"]][["group_id"]] = rep(1:25,each=25)
group_all[["625"]][["e_g sigma2"]] = rep(0,25)
group_all[["625"]][["b_g sigma2"]] = rep(0.5,25)
group_all[["625"]][["max_prob"]] = 0.2*exp(0.1+0.3*2)

group_all[["1000"]][["n"]] = 1000
group_all[["1000"]][["group_id"]] = rep(1:40,each=25)
group_all[["1000"]][["e_g sigma2"]] = rep(0,40)
group_all[["1000"]][["b_g sigma2"]] = rep(0.5,40)
group_all[["1000"]][["max_prob"]] = 0.2*exp(0.1+0.3*2)





group_str = function(group){
  # cluster size (assuming all clusters have the same size)
  group[["group size"]] = unname(table(group[["group_id"]]))
  # the number of clusters
  group[["#groups"]] = length(unique(group[["group_id"]]))
  
  max_prob = group[["max_prob"]]
  
  # correction for adding terms on exp{}
  a = -Inf
  b = log(1/max_prob)
  
  # baseline error term (cluster-level intercept term that does not interact with treatment.)
  err = c()
  
  if(all(group[["e_g sigma2"]]!=0)){
    for (i in 1:group[["#groups"]]){
      sig= sqrt(group[["e_g sigma2"]][i])
      e = rtruncnorm(n=1, a, b, mean = 0, sd = sig)- sig^2/2 - log(pnorm(b/sig-sig)/pnorm(b/sig))
      err = c(err,e)
    }
  }else{
    err = rep(0,group[["#groups"]])
  }
  group[["err"]] = err
  group[["group err"]] = rep(group[["err"]],group[["group size"]])
  
  # random cluster-level intercept term that interacts with treatment 
  b_g = c()
  
  if(all(group[["b_g sigma2"]]!=0)){
    for (i in 1:group[["#groups"]]){
      sig =  sqrt(group[["b_g sigma2"]][i])
      e = rtruncnorm(n=1, a, b, mean = 0, sd =sig)- sig^2/2 - log(pnorm(b/sig-sig)/pnorm(b/sig))
      b_g = c(b_g,e)
    }
  }else{
    b_g = rep(0,group[["#groups"]])
  }
  group[["b_g"]] = b_g
  group[["group b_g"]] = rep(group[["b_g"]],group[["group size"]])
  
  return(group)
} 