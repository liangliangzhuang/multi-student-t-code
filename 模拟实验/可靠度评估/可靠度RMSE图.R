# （最终版本 图7）==============================================================
# 汇总可靠度的结果，并画图
# n20 ==========
# III
n1 = c(10,20,30)
mttf_sum = list()
for(hh in 1:length(n1)){
  #III
  load(paste("模拟实验/可靠度评估/result/III-n",n1[hh],"-m10.RData",sep=""))
  III_10 = data.frame(rmse)
  III_10_mean = c(as.numeric(apply(III_10,2,mean)),10)
  III_10$m = 10
  III_10$Scen = "III"
  load(paste("模拟实验/可靠度评估/result/III-n",n1[hh],"-m20.RData",sep=""))
  III_20 = data.frame(rmse)
  III_20_mean = c(as.numeric(apply(III_20,2,mean,na.rm=TRUE)),20)
  III_20$m = 20
  III_20$Scen = "III"
  load(paste("模拟实验/可靠度评估/result/III-n",n1[hh],"-m30.RData",sep=""))
  III_30 = data.frame(rmse)
  III_30_mean = c(as.numeric(apply(III_30,2,mean)),30)
  III_30$m = 30
  III_30$Scen = "III"
  
  # IV
  load(paste("模拟实验/可靠度评估/result/IV-n",n1[hh],"-m10.RData",sep=""))
  IV_10 = data.frame(rmse)
  IV_10_mean = c(as.numeric(apply(IV_10,2,mean)),10)
  IV_10$m = 10
  IV_10$Scen = "IV"
  load(paste("模拟实验/可靠度评估/result/IV-n",n1[hh],"-m20.RData",sep=""))
  IV_20 = data.frame(rmse)
  IV_20_mean = c(as.numeric(apply(IV_20,2,mean)),20)
  IV_20$m = 20
  IV_20$Scen = "IV"
  load(paste("模拟实验/可靠度评估/result/IV-n",n1[hh],"-m30.RData",sep=""))
  IV_30 = data.frame(rmse)
  IV_30_mean = c(as.numeric(apply(IV_30,2,mean)),30)
  IV_30$m = 30
  IV_30$Scen = "IV"
  
  mttf_sum[[hh]] = cbind(rbind(III_10,III_20,III_30,IV_10,IV_20,IV_30),n1[hh])
}


# 结果汇总
mttf_sum_all = bind_rows(mttf_sum)
colnames(mttf_sum_all) = c("Wiener","Student.t","m","Scen","n")
mttf_sum_all %>% pivot_longer(1:2,values_to = "value",names_to = "Model") -> mttf_sum_all2
mttf_sum_all2$m = factor(mttf_sum_all2$m,levels = c(10,20,30))
mttf_sum_all2$Scen = factor(mttf_sum_all2$Scen,levels = c("III","IV"))
mttf_sum_all2$Model = factor(mttf_sum_all2$Model,levels = c("Wiener","Student.t"))


f_names <- list('10' = TeX(c("$m=10$")), '20' = TeX(c("$m=20$")), '30' = TeX(c("$m=30$"))) #,
#"rho" = TeX(c("$\\sigma_{12}$")))
f_labeller <- function(variable, value){return(f_names[value])}

mttf_sum_all2 %>% filter(Scen=="IV") %>% 
  ggplot(aes(factor(n), value, fill = Model)) +
  geom_boxplot(alpha=0.8) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "#008580", position = position_dodge(0.75)) +
  scale_fill_viridis(name="Model", discrete = TRUE, labels=c("Wiener","Student-t")) +
  facet_wrap(vars(factor(m)),nrow=1,scales = "fixed",labeller = f_labeller) +
  theme_bw() + 
  xlab("n") +
  ylab(TeX(r'(RMSE($\times 10^{-2}))')) +
  theme(panel.grid = element_blank(),legend.position = "right",
        legend.key.size = unit(0.6, "cm"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10))

# ggsave(paste("模拟实验/result/MTTF2.pdf",sep=''), width =8, height = 4)
# save.image(file = paste("模拟实验/可靠度评估/可靠度评估-0619.RData", sep = ""))










