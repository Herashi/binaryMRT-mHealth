library(truncnorm)

## specify group structure
group_all = list()
group_all[[1]] = list()
group_all[[2]] = list()
group_all[[3]] = list()
group_all[[4]] = list()
group_all[[5]] = list()
group_all[[6]] = list()
group_all[[7]] = list()
group_all[[8]] = list()
group_all[[9]] = list()
group_all[[10]] = list()
group_all[[11]] = list()


group_all[[11]][["n"]] = 200
group_all[[11]][["group_id"]] = rep(1:40,each=5)
group_all[[11]][["e_g sigma2"]] = rep(0,40)
group_all[[11]][["b_g sigma2"]] = rep(0,40)

group_all[[1]][["n"]] = 200
group_all[[1]][["group_id"]] = rep(1:40,each=5)
group_all[[1]][["e_g sigma2"]] = rep(0,40)
group_all[[1]][["b_g sigma2"]] = rep(0.1^2,40)

group_all[[2]][["n"]] = 200
group_all[[2]][["group_id"]] = rep(1:40,each=5)
group_all[[2]][["e_g sigma2"]] = rep(0,40)
group_all[[2]][["b_g sigma2"]] = rep(0.2^2,40)

group_all[[3]][["n"]] = 200
group_all[[3]][["group_id"]] = rep(1:40,each=5)
group_all[[3]][["e_g sigma2"]] = rep(0,40)
group_all[[3]][["b_g sigma2"]] = rep(0.3^2,40)

group_all[[4]][["n"]] = 200
group_all[[4]][["group_id"]] = rep(1:40,each=5)
group_all[[4]][["e_g sigma2"]] = rep(0,40)
group_all[[4]][["b_g sigma2"]] = rep(0.4^2,40)

group_all[[5]][["n"]] = 200
group_all[[5]][["group_id"]] = rep(1:40,each=5)
group_all[[5]][["e_g sigma2"]] = rep(0,40)
group_all[[5]][["b_g sigma2"]] = rep(0.5^2,40)

group_all[[6]][["n"]] = 200
group_all[[6]][["group_id"]] = rep(1:40,each=5)
group_all[[6]][["e_g sigma2"]] = rep(0,40)
group_all[[6]][["b_g sigma2"]] = rep(0.6^2,40)

group_all[[7]][["n"]] = 200
group_all[[7]][["group_id"]] = rep(1:40,each=5)
group_all[[7]][["e_g sigma2"]] = rep(0,40)
group_all[[7]][["b_g sigma2"]] = rep(0.7^2,40)

group_all[[8]][["n"]] = 200
group_all[[8]][["group_id"]] = rep(1:40,each=5)
group_all[[8]][["e_g sigma2"]] = rep(0,40)
group_all[[8]][["b_g sigma2"]] = rep(0.8^2,40)

group_all[[9]][["n"]] = 200
group_all[[9]][["group_id"]] = rep(1:40,each=5)
group_all[[9]][["e_g sigma2"]] = rep(0,40)
group_all[[9]][["b_g sigma2"]] = rep(0.9^2,40)

group_all[[10]][["n"]] = 200
group_all[[10]][["group_id"]] = rep(1:40,each=5)
group_all[[10]][["e_g sigma2"]] = rep(0,40)
group_all[[10]][["b_g sigma2"]] = rep(1^2,40)

group_str = function(group){
  # cluster size (assuming all clusters have the same size)
  group[["group size"]] = unname(table(group[["group_id"]]))
  # the number of clusters
  group[["#groups"]] = length(unique(group[["group_id"]]))
  
  
  # correction for adding terms on exp{}
  a = -1
  b = 1
  
  # baseline error term (cluster-level intercept term that does not interact with treatment.)
  err = c()
  
  if(all(group[["e_g sigma2"]]!=0)){
    for (i in 1:group[["#groups"]]){
      sig= sqrt(group[["e_g sigma2"]][i])
      e = rtruncnorm(n=1, a, b, mean = 0, sd = sig)- sig^2/2 - 
        log((pnorm(b/sig-sig)-pnorm(a/sig-sig))/(pnorm(b/sig)-pnorm(a/sig)))
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
      e = rtruncnorm(n=1, a, b, mean = 0, sd =sig)- sig^2/2 - 
        log((pnorm(b/sig-sig)-pnorm(a/sig-sig))/(pnorm(b/sig)-pnorm(a/sig)))
      b_g = c(b_g,e)
    }
  }else{
    b_g = rep(0,group[["#groups"]])
  }
  group[["b_g"]] = b_g
  group[["group b_g"]] = rep(group[["b_g"]],group[["group size"]])
  
  return(group)
} 


