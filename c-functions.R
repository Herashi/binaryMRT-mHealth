## specify group structure
group_all = list()
group_all[["250"]] = list()



group_all[["250"]][["n"]] = 250
group_all[["250"]][["group_id"]] = rep(1:25,each=10)
group_all[["250"]][["baseline sigma2"]] = rep(0.5,25)
group_all[["250"]][["slope sigma2"]] = rep(0.1,25)




group_str = function(group){
  # cluster size (assuming all clusters have the same size)
  group[["group size"]] = unname(table(group[["group_id"]]))
  # the number of clusters
  group[["#groups"]] = length(unique(group[["group_id"]]))

  # baseline error term (cluster-level intercept term that does not interact with treatment.)
  err = c()
  for (i in 1:group[["#groups"]]){
    e = rnorm(1,mean = 0,sd = sqrt(group[["baseline sigma2"]][i]))
    err = c(err,e)
  }
  group[["err"]] = err
  group[["group err"]] = rep(group[["err"]],group[["group size"]])
  
  # random cluster-level intercept term that interacts with treatment 
  slope = c()
  for (i in 1:group[["#groups"]]){
    e = rnorm(1,mean = 0,sd = sqrt(group[["slope sigma2"]][i]))
    slope = c(slope,e)
  }
  group[["slope"]] = slope
  group[["random slope"]] = rep(group[["slope"]],group[["group size"]])
  
  
  return(group)
} 

