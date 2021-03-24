load("~/binary_mrt/Scenario 2/figure/test.RData")

df = data.frame(matrix(NA, nrow =11 , ncol = 2))
colnames(df) = c("SD","Coverage")

df[,"SD"] = 0:10/100
df[,"Coverage"] = result_df_collected$cp.adj[c(11,1:10)]

ggplot(data = df)+
  geom_line(aes(x=SD,y=Coverage))+
  xlab("Standard Deviation")+
  ylab("Coverage Probability")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0.8,1,0.05),limits = c(0.8,1))+
  geom_hline(yintercept = 0.95,linetype="dashed")




# library(ggplot2)
# library(ggpubr)
# 
# # Two plots
# 
# df = data.frame(matrix(NA, nrow =20, ncol = 3))
# colnames(df) = c("Method","size","Coverage")
# 
# df[,"size"] = rep(c(1,2,4,5,10,6,8,11,13,15),2)
# 
# 
# setwd("C:/Users/herashi/Desktop/figure/group_size")
# 
# for (k in 1:10){
#   load(paste0("test-",k,".RData"))
#   df[k,"Method"] = "WCLS"
#   df[k,"Coverage"] =  omit[1,"cpc"]
# }
# 
# 
# setwd("C:/Users/herashi/Desktop/figure/group_size_Gwcls")
# 
# for (k in 1:10){
#   load(paste0("test-",k,".RData"))
#   df[10+ k,"Method"] = "C-WCLS"
#   df[10+ k,"Coverage"] =  omit[1,"cpc"]
# }
# 
# 
# 
# df$Method = as.factor(df$Method)
# 
# p1 <- ggplot(data = df)+
#   geom_line(aes(x=size,y=Coverage, colour =Method))+
#   xlab("Group Size")+
#   ylab("Coverage Probability")+
#   theme_bw()+
#   scale_y_continuous(breaks = seq(0.6,1,0.05),limits = c(0.6,1))+
#   geom_hline(yintercept = 0.95,linetype="dashed")
