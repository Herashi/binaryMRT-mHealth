library(truncnorm)

## specify group structure
group_all = list()
group_all[[1]] = list()
group_all[[2]] = list()
group_all[[3]] = list()
group_all[[4]] = list()
group_all[[5]] = list()
group_all[[6]] = list()





group_all[[1]][["n"]] = 125
group_all[[1]][["group_id"]] = rep(1:25,each=5)
group_all[[1]][["e_g sigma2"]] = rep(0,25)
group_all[[1]][["b_g sigma2"]] = rep(0.5^2,25)

group_all[[2]][["n"]] = 250
group_all[[2]][["group_id"]] = rep(1:25,each=10)
group_all[[2]][["e_g sigma2"]] = rep(0,25)
group_all[[2]][["b_g sigma2"]] = rep(0.5^2,25)

group_all[[3]][["n"]] = 500
group_all[[3]][["group_id"]] = rep(1:50,each=10)
group_all[[3]][["e_g sigma2"]] = rep(0,50)
group_all[[3]][["b_g sigma2"]] = rep(0.5^2,50)

group_all[[4]][["n"]] = 1000
group_all[[4]][["group_id"]] = rep(1:50,each=20)
group_all[[4]][["e_g sigma2"]] = rep(0,50)
group_all[[4]][["b_g sigma2"]] = rep(0.5^2,50)

group_all[[5]][["n"]] = 2000
group_all[[5]][["group_id"]] = rep(1:100,each=20)
group_all[[5]][["e_g sigma2"]] = rep(0,100)
group_all[[5]][["b_g sigma2"]] = rep(0.5^2,100)

group_all[[6]][["n"]] = 2500
group_all[[6]][["group_id"]] = rep(1:100,each=25)
group_all[[6]][["e_g sigma2"]] = rep(0,100)
group_all[[6]][["b_g sigma2"]] = rep(0.5^2,100)



group_str = function(group){
  # cluster size (assuming all clusters have the same size)
  group[["group size"]] = unname(table(group[["group_id"]]))
  # the number of clusters
  group[["#groups"]] = length(unique(group[["group_id"]]))
  
  
  # correction for adding terms on exp{}
  a = -0.8
  b = 0.8
  
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


