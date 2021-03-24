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




group_all[[1]][["n"]] = 100
group_all[[1]][["group_id"]] = rep(1:10,each=10)
group_all[[1]][["baseline sigma2"]] = rep(0,10)
group_all[[1]][["slope sigma2"]] = rep(0.01^2,10)

group_all[[2]][["n"]] = 100
group_all[[2]][["group_id"]] = rep(1:10,each=10)
group_all[[2]][["baseline sigma2"]] = rep(0,10)
group_all[[2]][["slope sigma2"]] = rep(0.02^2,10)

group_all[[3]][["n"]] = 100
group_all[[3]][["group_id"]] = rep(1:10,each=10)
group_all[[3]][["baseline sigma2"]] = rep(0,10)
group_all[[3]][["slope sigma2"]] = rep(0.03^2,10)

group_all[[4]][["n"]] = 100
group_all[[4]][["group_id"]] = rep(1:10,each=10)
group_all[[4]][["baseline sigma2"]] = rep(0,10)
group_all[[4]][["slope sigma2"]] = rep(0.04^2,10)

group_all[[5]][["n"]] = 100
group_all[[5]][["group_id"]] = rep(1:10,each=10)
group_all[[5]][["baseline sigma2"]] = rep(0,10)
group_all[[5]][["slope sigma2"]] = rep(0.05^2,10)

group_all[[6]][["n"]] = 100
group_all[[6]][["group_id"]] = rep(1:10,each=10)
group_all[[6]][["baseline sigma2"]] = rep(0,10)
group_all[[6]][["slope sigma2"]] = rep(0.06^2,10)

group_all[[7]][["n"]] = 100
group_all[[7]][["group_id"]] = rep(1:10,each=10)
group_all[[7]][["baseline sigma2"]] = rep(0,10)
group_all[[7]][["slope sigma2"]] = rep(0.07^2,10)

group_all[[8]][["n"]] = 100
group_all[[8]][["group_id"]] = rep(1:10,each=10)
group_all[[8]][["baseline sigma2"]] = rep(0,10)
group_all[[8]][["slope sigma2"]] = rep(0.08^2,10)

group_all[[9]][["n"]] = 100
group_all[[9]][["group_id"]] = rep(1:10,each=10)
group_all[[9]][["baseline sigma2"]] = rep(0,10)
group_all[[9]][["slope sigma2"]] = rep(0.09^2,10)

group_all[[10]][["n"]] = 100
group_all[[10]][["group_id"]] = rep(1:10,each=10)
group_all[[10]][["baseline sigma2"]] = rep(0,10)
group_all[[10]][["slope sigma2"]] = rep(0.1^2,10)





group_str = function(group){
  # cluster size (assuming all clusters have the same size)
  group[["group size"]] = unname(table(group[["group_id"]]))
  # the number of clusters
  group[["#groups"]] = length(unique(group[["group_id"]]))
  # correction for adding terms on exp{}
  group[["baseline correction"]] = unique(group[["baseline sigma2"]])/2
  group[["slope correction"]] = unique(group[["slope sigma2"]])/2
  
  # baseline error term (cluster-level intercept term that does not interact with treatment.)
  err = c()
  for (i in 1:group[["#groups"]]){
    e = rnorm(1,mean = 0,sd = sqrt(group[["baseline sigma2"]][i]))
    err = c(err,e)
  }
  group[["err"]] = err
  group[["group err"]] = rep(group[["err"]],group[["group size"]])- group[["baseline correction"]]
  
  # random cluster-level intercept term that interacts with treatment 
  slope = c()
  for (i in 1:group[["#groups"]]){
    e = rnorm(1,mean = 0,sd = sqrt(group[["slope sigma2"]][i]))
    slope = c(slope,e)
  }
  group[["slope"]] = slope
  group[["random slope"]] = rep(group[["slope"]],group[["group size"]])-group[["slope correction"]]
  
  
  return(group)
} 