
#Below this should be another job (part2) after the array jobs finish running
#set empty matrix
# SET IT HERE
setwd("~/binary_mrt/Scenario 4")

MASTERLIST =list.files(pattern="^fity")

out <- NULL

for (i in 1:length(MASTERLIST)){
  load(file= MASTERLIST[i])
  out <-rbind(out, result_df)
}

print(length(MASTERLIST))

out$rmse <- with(out, bias^2)


## mean and SD estimate, number of replicates
out <- cbind(aggregate(cbind(bias,se.unadj,se.adj, cp.unadj, cp.adj, rmse) ~
                         est+ss,
                       data = out, FUN = mean),
             sd = aggregate(bias ~est+ss,
                            data = out, FUN=var)$bias)

out$sd = sqrt(out$sd)
out$rmse = sqrt(out$rmse)

save(out,file = "test.RData")