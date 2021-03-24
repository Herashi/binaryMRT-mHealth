## specify group structure
group_all = list()
group_all[["250"]] = list()
group_all[["625"]] = list()
group_all[["1000"]] = list()



group_all[["250"]][["n"]] = 250
group_all[["250"]][["group_id"]] = rep(1:25,each=10)
group_all[["250"]][["baseline sigma2"]] = rep(0.05,25)
group_all[["250"]][["slope sigma2"]] = rep(0.01,25)

group_all[["625"]][["n"]] = 625
group_all[["625"]][["group_id"]] = rep(1:25,each=25)
group_all[["625"]][["baseline sigma2"]] = rep(0.05,25)
group_all[["625"]][["slope sigma2"]] = rep(0.01,25)

group_all[["1000"]][["n"]] = 1000
group_all[["1000"]][["group_id"]] = rep(1:40,each=25)
group_all[["1000"]][["baseline sigma2"]] = rep(0.05,40)
group_all[["1000"]][["slope sigma2"]] = rep(0.01,40)





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

