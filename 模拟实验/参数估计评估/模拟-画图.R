# 模拟实验中各个参数的RMSE结果（图5 - 6）
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(viridis)
library(latex2exp)

sum_dat = list()
for(i in 1:5){
  sum_dat[[i]] = read.xlsx("模拟实验/参数估计评估/数据汇总.xlsx",sheet = i)
}
 

# 图1 线性 p=2 ===============
data.frame(sum_dat[[1]]) %>% filter(X3 == "RMSE") -> t2_dat #线性 p=2
t2_dat[,4:10] = t2_dat[,4:10] *100 #这里乘了100！

t2_dat$m = rep(c(10,20,30),times=3)
t2_dat$n = rep(c(10,20,30),each=3)
colnames(t2_dat) = c("n","m","RMSE","eta1","eta2","delta1","delta2","sigma1","sigma2","rho")

f_names <- list('eta1' = TeX(c("$\\eta_{1}$")), 'eta2' = TeX(c("$\\eta_{2}$")),
                'delta1' = TeX(c("$\\delta^2_{1}$")), 'delta2' = TeX(c("$\\delta^2_{2}$")),
                'sigma1' = TeX(c("$\\sigma^2_{1}$")), 'sigma2' = TeX(c("$\\sigma^2_{2}$")) #"rho" = TeX(c("$\\rho_{(1,2)}$"
                )
f_labeller <- function(variable, value){return(f_names[value])}
y_ranges <- list(c(2, 8),  c(8, 22), c(8, 20),c(15,35),c(8,15),c(13,28)) 
p1_dat = t2_dat[,1:9] %>% pivot_longer(eta1:sigma2,values_to = "value", names_to = "name") 
p1_dat$name = factor(p1_dat$name,levels = c("eta1","eta2","delta1","delta2",
                                            "sigma1", "sigma2"))
p1 = ggplot(p1_dat, aes(factor(n), value)) + 
  # geom_col(aes(fill = factor(m)),position = "dodge",alpha=1) + 
  geom_point(aes(color = factor(m)),size=2) +
  facet_wrap(vars(name),scales = "free_y",labeller = f_labeller, nrow=2,ncol=4) + 
  theme_bw() + scale_color_viridis(name = "n",discrete = TRUE) + 
  theme(panel.grid = element_blank(),legend.position = "bottom") + xlab("m") + ylab(TeX(r'(RMSE($\times 10^{-2}))')) 
# ggsave(paste("模拟实验/参数估计评估/p2_point(9-6).pdf",sep=''), p1,width =9, height = 6) # 10-5
p1

# 图1 非线性 p=2 ===============
data.frame(sum_dat[[3]]) %>% filter(X3 == "RMSE") -> p2_dat #线性 p=2
p2_dat[,4:12] = p2_dat[,4:12] *100

p2_dat$m = rep(c(10,20,30),times=3)
p2_dat$n = rep(c(10,20,30),each=3)
colnames(p2_dat) = c("n","m","RMSE","eta1","eta2","delta1","delta2","sigma1","sigma2","rho","gamma1","gamma2")
f2_names <- list('eta1' = TeX(c("$\\eta_{1}$")), 'eta2' = TeX(c("$\\eta_{2}$")),
                'delta1' = TeX(c("$\\delta^2_{1}$")), 'delta2' = TeX(c("$\\delta^2_{2}$")),
                'sigma1' = TeX(c("$\\sigma^2_{1}$")), 'sigma2' = TeX(c("$\\sigma^2_{2}$")),
                 "gamma1" = TeX(c("$\\gamma_{1}$")),  "gamma2" = TeX(c("$\\gamma_{2}$"))) #"rho" = TeX(c("$\\rho_{(1,2)}$")),
f2_labeller <- function(variable, value){return(f2_names[value])}


p2_dat = p2_dat %>% dplyr::select(-rho) %>% pivot_longer(eta1:gamma2,values_to = "value", names_to = "name")
  
p2_dat$name = factor(p2_dat$name,levels = c("eta1","eta2","delta1","delta2",
                                            "sigma1", "sigma2", "gamma1", "gamma2"))
p2 = ggplot(p2_dat, aes(factor(n), value)) + 
  # geom_col(aes(fill = factor(m)),position = "dodge",alpha=1) + 
  geom_point(aes(color = factor(m)),size=2) +
  facet_wrap(vars(factor(name)),scales = "free_y",labeller = f2_labeller, nrow=2) + 
  scale_color_viridis(name = "n",discrete = TRUE) + theme_bw() + 
  theme(panel.grid = element_blank(),legend.position = "bottom") + xlab("m") + ylab(TeX(r'(RMSE($\times 10^{-2}$))')) 
# ggsave(paste("模拟实验/参数估计评估/f2_point(9-6).pdf",sep=''), p2, width =9, height = 6)

library(ggpubr)
library(grid)
pp1 = ggarrange(p1 + rremove("ylab") + rremove("xlab"), 
          p2 + rremove("ylab") + rremove("xlab"), 
          labels = c("I","III"), label.y = 1,label.x = 0,
          common.legend = TRUE, legend="top", nrow=2, align ="hv")

annotate_figure(pp1, left = textGrob(TeX(r'(RMSE($\times 10^{-2}))'), rot = 90, vjust = 0.5, gp = gpar(cex = 0.9)),
                bottom = textGrob("m", gp = gpar(cex = 0.9)))
# ggsave(paste("模拟实验/参数估计评估/RMSE_all(9-9).pdf",sep=''), width =9, height = 9)
# ggsave(paste("模拟实验/参数估计评估/RMSE_all(14-6).pdf",sep=''), width =14, height = 6)


# 图三 覆盖率======
#### 版本一 ========
cp_dat = data.frame(sum_dat[[5]]) 
cp_dat[,4:12] = abs((cp_dat[,4:12]-95)/95 * 100) #计算RMSE * 百分比
cp_dat$Scen. = rep(c("t1","t2"),each = 9)
cp_dat$n = rep(rep(c(10,20,30),each=3),2)
# 相对损失-绝对值
colnames(cp_dat) = c("Scen","n","m","eta1","eta2","delta1","delta2","sigma1","sigma2","rho","gamma1","gamma2")
f3_names <- list('eta1' = TeX(c("$\\eta_{1}$")), 'eta2' = TeX(c("$\\eta_{2}$")),
                 'delta1' = TeX(c("$\\delta^2_{1}$")), 'delta2' = TeX(c("$\\delta^2_{2}$")),
                 'sigma1' = TeX(c("$\\sigma^2_{1}$")), 'sigma2' = TeX(c("$\\sigma^2_{2}$")),
                 "gamma1" = TeX(c("$\\gamma_{1}$")), 
                 "gamma2" = TeX(c("$\\gamma_{2}$")))
f3_labeller <- function(variable, value){return(f3_names[value])}

p4_dat = cp_dat %>% filter(Scen=="t2") %>% dplyr::select(!`rho`) %>% pivot_longer(eta1:gamma2,values_to = "value",names_to = "name") 

p4_dat$name = factor(p4_dat$name,levels = c("eta1","eta2","delta1","delta2",
                                            "sigma1", "sigma2", "gamma1", "gamma2"))

p3 = ggplot(p4_dat,aes(factor(n), value)) + 
  geom_point(aes(color = factor(m)),size=2) +
  # geom_col(aes(fill = factor(m)), position = "dodge",alpha=1) +
  facet_wrap(vars(name),scales = "free_y",labeller = f3_labeller, nrow=2) +
  scale_color_viridis(name = "n",discrete = TRUE) + theme_bw() + 
  theme(panel.grid = element_blank()) + xlab("m") + ylab(TeX(r'(Absolute value of relative loss ($\%$))')) 
# ggsave(paste("模拟实验/参数估计评估/cp_point.pdf",sep=''), p3, width = 10, height = 5)
p3

#### 版本二 ========

cp_dat = data.frame(sum_dat[[5]]) 
# cp_dat[,4:12] = abs((cp_dat[,4:12]-95)/95 * 100) #计算RMSE * 百分比
cp_dat$Scen. = rep(c("t1","t2"),each = 9)
cp_dat$n = rep(rep(c(10,20,30),each=3),2)
# 相对损失-绝对值
colnames(cp_dat) = c("Scen","n","m","eta1","eta2","delta1","delta2","sigma1","sigma2","rho","gamma1","gamma2")
f3_labeller <- function(variable, value){return(f3_names[value])}

# Scen I ====
cp_dat1 = cp_dat %>% filter(Scen=="t1") %>% dplyr::select(!c(`gamma1`,`gamma2`)) %>% pivot_longer(eta1:sigma2,values_to = "value",names_to = "name")
cp_dat1$name = factor(cp_dat1$name,levels = c("delta1","delta2","eta1","eta2","sigma1","sigma2"))

f3_names <- list('eta1' = TeX(c("$\\eta_{1}$")), 'eta2' = TeX(c("$\\eta_{2}$")),
                 'delta1' = TeX(c("$\\delta^2_{1}$")), 'delta2' = TeX(c("$\\delta^2_{2}$")),
                 'sigma1' = TeX(c("$\\sigma^2_{1}$")), 'sigma2' = TeX(c("$\\sigma^2_{2}$"))
                 )

p3 = cp_dat1 %>% dplyr::select(!`rho`) %>% ggplot(aes(factor(n), value)) + 
  geom_point(aes(color = factor(m)),size=2) +
  scale_y_continuous(limits = c(80,96)) +
  geom_hline(yintercept = 95,linetype = 2) +
  # geom_col(aes(fill = factor(m)), position = "dodge",alpha=1) +
  facet_wrap(vars(factor(name,levels = )),scales = "free_y",labeller = f3_labeller, ncol=4) +
  scale_color_viridis(name = "n",discrete = TRUE) + theme_bw() + 
  theme(panel.grid = element_blank()) + xlab("m") + ylab("Coverage probability") 
p3 
# ggsave(paste("模拟实验/参数估计评估/cp_point-2.pdf",sep=''), p3, width = 10, height = 5)

# Scen III ====
cp_dat2 = cp_dat %>% filter(Scen=="t2") %>% dplyr::select(!`rho`) %>% pivot_longer(eta1:gamma2,values_to = "value",names_to = "name")
cp_dat2$name = factor(cp_dat2$name,levels = c("eta1","eta2","delta1","delta2","sigma1","sigma2","gamma1","gamma2"))

f4_names <- list('eta1' = TeX(c("$\\eta_{1}$")), 'eta2' = TeX(c("$\\eta_{2}$")),
                 'delta1' = TeX(c("$\\delta^2_{1}$")), 'delta2' = TeX(c("$\\delta^2_{2}$")),
                 'sigma1' = TeX(c("$\\sigma^2_{1}$")), 'sigma2' = TeX(c("$\\sigma^2_{2}$")),
                 "gamma1" = TeX(c("$\\gamma_{1}$")), 
                 "gamma2" = TeX(c("$\\gamma_{2}$")))
f4_labeller <- function(variable, value){return(f4_names[value])}

p4 = cp_dat2 %>% ggplot(aes(factor(n), value)) + 
  geom_point(aes(color = factor(m)),size=2) +
  scale_y_continuous(limits = c(80,96)) +
  geom_hline(yintercept = 95,linetype = 2) +
  # geom_col(aes(fill = factor(m)), position = "dodge",alpha=1) +
  facet_wrap(vars(factor(name)),scales = "free_y",labeller = f4_labeller, nrow=2) +
  scale_color_viridis(name = "n",discrete = TRUE) + theme_bw() + 
  theme(panel.grid = element_blank()) + xlab("m") + ylab("Coverage probability") 
p4

library(ggpubr)
library(grid)
pp2 = ggarrange(p3 + rremove("ylab") + rremove("xlab"), 
                p4 + rremove("ylab") + rremove("xlab"), 
                labels = c("I","III"), label.y = 1,label.x = 0,
                common.legend = TRUE, legend="top", nrow=2, align ="hv")

annotate_figure(pp2, left = textGrob("Coverage probability", rot = 90, vjust = 0.5, gp = gpar(cex = 0.9)),
                bottom = textGrob("m", gp = gpar(cex = 0.9)))

# ggsave(paste("模拟实验/参数估计评估/CP_all(9-9).pdf",sep=''), width = 9, height = 9)
# ggsave(paste("模拟实验/参数估计评估/CP_all(2).pdf",sep=''), width =12, height = 6)

# save.image(file = paste("模拟实验/参数估计评估/参数估计评估-0619.RData", sep = ""))



