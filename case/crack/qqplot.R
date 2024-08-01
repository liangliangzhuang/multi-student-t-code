# === Fatigue crack-size data 的可视化 ==== 
library(MASS)      # use 'fitdistr' function to fit a normal distribution
library(QRM)       # use 'fit.st' function to fit a Student-t distribution
library(metRology) # use 'qt.scaled' function to find the quantile of a Student-t distribution
library(numDeriv)
library(robustbase)
library(ggpubr)
library(grid)
library(latex2exp)

cal_dat3 = read.csv("case/crack/result/tidy_crack_dat.csv")
time2 = seq(3,m,3)
newdata = newdata2 = newdata3 = matrix(NA,n,length(time2))

for(i in 1:length(time2)){
  uuu = cal_dat3[cal_dat3$PC=="PC1" & cal_dat3$Time == time2[i],]
  vvv = cal_dat3[cal_dat3$PC=="PC2" & cal_dat3$Time == time2[i],]
  zzz = cal_dat3[cal_dat3$PC=="PC3" & cal_dat3$Time == time2[i],]
  newdata[,i] = uuu[,4] #as.numeric(unlist(uuu[,4]))
  newdata2[,i] = vvv[,4]
  newdata3[,i] = zzz[,4]
}
newdata = data.frame(newdata); newdata2 = data.frame(newdata2); newdata3 = data.frame(newdata3)


# 箱线图（图2）
colnames(newdata) <- time2; colnames(newdata2) <- time2; colnames(newdata3) <- time2
comdata = cbind(rbind(newdata,newdata2,newdata3),rep(paste("PC",1:p,sep=""),each = n))
colnames(comdata) = c(time2,"PC")
box_plot = tidyr::pivot_longer(comdata,cols = !PC,values_to = "value",names_to = "Unit") %>% 
  ggplot(aes(factor(Unit,levels = time2),value,fill=factor(Unit,levels = time2))) + 
  geom_boxplot() +
  facet_wrap(vars(PC),scale="free") +
  scale_fill_viridis(discrete = TRUE,alpha = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") +
  xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches))'))
# ggsave(paste("case/crack/box_plot_new.pdf",sep=""), height = 5, width = 8)
