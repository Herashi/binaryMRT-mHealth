load("~/binary_mrt/Scenario 2/figure/test_cwcls_var.RData")
library(ggplot2)
df = data.frame(matrix(NA, nrow =22 , ncol = 2))
colnames(df) = c("SD","Coverage")

df[,"SD"] = rep(0:10/10,2)
df[1:11,"Coverage"] = result_df_collected_1$cp.adj[c(11,1:10)]

load("~/binary_mrt/Scenario 2/figure/test-var.RData")
df[12:22,"Coverage"] = result_df_collected_2$cp.adj[c(11,1:10)]

df[,"Method"] = rep(c("C-WCLS","EMEE"),each = 11)


p1 <- ggplot(data = df)+
     geom_line(aes(x=SD,y=Coverage, colour =Method))+
     xlab("Standard Deviation")+
     ylab("Coverage Probability")+
     theme_bw()+
     scale_y_continuous(breaks = seq(0.55,1,0.05),limits = c(0.55,1))+
     geom_hline(yintercept = 0.95,linetype="dashed")



library(ggpubr)

# Two plots

df = data.frame(matrix(NA, nrow =18, ncol = 3))
colnames(df) = c("Method","size","Coverage")
df[,"Method"] = rep(c("C-WCLS","EMEE"),each = 9)

load("~/binary_mrt/Scenario 2/figure/group_size_cwcls.RData")
df[1:9,"Coverage"] = result_df_collected_1$cp.adj[1:9]

load("~/binary_mrt/Scenario 2/figure/test_groupsize.RData")
df[10:18,"Coverage"] = result_df_collected_2$cp.adj[1:9]
df[,"size"] = rep(c(1,2,4,5,8,10,12,15,18),2)

df$Method = as.factor(df$Method)

p2 <- ggplot(data = df)+
  geom_line(aes(x=size,y=Coverage, colour =Method))+
  xlab("Group Size")+
  ylab("Coverage Probability")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0.55,1,0.05),limits = c(0.55,1))+
  geom_hline(yintercept = 0.95,linetype="dashed")

ggarrange(p2,p1,legend = "top",common.legend = T)
