# ========= 真实数据分析(Fatigue Crack-size) =========
library(ggplot2)
library(tidyverse)
library(viridis)
library(reshape2)
library(latex2exp)
library(MASS)
library(openxlsx)
source("case_fct.R")
source("data_process.R")
options(scipen = 4) # 设置小数点


# 1.数据读取 =======
fc_dat = fc_dat1 = fc_dat_sub = fc_diff_dat = list()
n = 6
for (i in 1:n) {
  fc_dat[[i]] = read.xlsx("crack/crack.xlsx",sheet = i) *100
  fc_dat_sub[[i]] =  fc_dat[[i]][,] - 0.9*100
  fc_dat1[[i]] =  fc_dat[[i]][-1,] - 0.9*100
  fc_diff_dat[[i]] = apply(fc_dat[[i]],2,diff) 
}
m <- dim(fc_diff_dat[[1]])[1] # diff的维数
p = dim(fc_diff_dat[[1]])[2]


time <- 1 # 观测的时间间隔
t <-(seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)


fc_diff_dat

# 数据处理
LED_dat = Real_path(led_cum = fc_dat_sub, led_diff = fc_diff_dat, p=p, m=m, n=n)

y.diff = LED_dat[[1]]; sumys = LED_dat[[2]]; y =  LED_dat[[3]]
sim_diff_dat = LED_dat[[4]] # 便于绘图（差分）
sim_cum_dat = LED_dat[[5]] # 便于绘图（累计）

# 2.模型选择 ========
# 2.1 绘制退化数据路径图
thr = c(0.9,0.5,0.4)*100  # c(0.45,0.35,0.2)
dummy2 <- data.frame(PC = c("PC1","PC2","PC3"), Z = thr)
g1 = degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") + 
  geom_hline(data = dummy2,aes(yintercept = Z),linetype = "dashed") 
g1
degradation.path.plot(data = sim_cum_dat, leg.pos = "bottom",ech = 5, scale1 = "fixed") +
  scale_color_viridis(discrete = TRUE, name = "Unit",guide = guide_legend(nrow = 1)) + 
  xlab(TeX(r'(Time (hours $\times$ 336))'))  + 
  ylab(TeX(r'(Y(t) (inches $\times 10^{-2}$))') ) 
  # theme(legend.box = "horizontal")
# ggsave("crack/result/crack-dat.pdf", height = 5, width = 9)

# 2.2 箱线图（图2）
data2 = map(sim_cum_dat, ~ mutate(.x, Time = 0:(n()-1)))
merged_df2 <- bind_rows(data2, .id = "Unit")
cal_dat2 = merged_df2 %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
time2 = seq(3,m,3)
newdata = matrix(NA,n,length(time2)); newdata2 = matrix(NA,n,length(time2)); newdata3 = matrix(NA,n,length(time2))

for(i in 1:length(time2)){
  uuu = cal_dat2[cal_dat2$PC=="PC1" & cal_dat2$Time == time2[i],]
  vvv = cal_dat2[cal_dat2$PC=="PC2" & cal_dat2$Time == time2[i],]
  zzz = cal_dat2[cal_dat2$PC=="PC3" & cal_dat2$Time == time2[i],]
  newdata[,i] = as.numeric(unlist(uuu[,4])) #as.numeric(unlist(uuu[,4]))
  newdata2[,i] = as.numeric(unlist(vvv[,4]))
  newdata3[,i] = as.numeric(unlist(zzz[,4]))
}
newdata = data.frame(newdata); newdata2 = data.frame(newdata2); newdata3 = data.frame(newdata3)
colnames(newdata) = colnames(newdata2) = colnames(newdata3) = time2 
comdata = cbind(rbind(newdata, newdata2, newdata3),rep(paste("PC",1:p,sep=""),each = n))
colnames(comdata) = c(time2,"PC")
box_plot = tidyr::pivot_longer(comdata,cols = !PC,values_to = "value",names_to = "Unit") %>% 
  ggplot(aes(factor(Unit,levels = time2),value,fill=factor(Unit,levels = time2))) + 
  geom_boxplot() +
  facet_wrap(vars(PC),scale="free") +
  scale_fill_viridis(discrete = TRUE,alpha = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") +
  xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches $\times 10^{-2}$))') )
box_plot 
# ggsave(paste("crack/result/box_plot.pdf",sep=""), height = 5, width = 8)

# 2.3 PC之间的相关性
cor_re = as.numeric()
for(i in 1:n){
  cor_cal = data.frame("PC1" = sim_diff_dat[[i]][,1], "PC2" = sim_diff_dat[[i]][,2])
  cor_re[i] = cor(cor_cal$PC2, cor_cal$PC1)
  ggplot(cor_cal,aes(x=PC1,y=PC2)) + 
    geom_point(size=1.3,shape=19,color = "#4BA6A3")+
    geom_smooth(method=lm,color ="#978ED6",alpha=0.3) + 
    theme_bw() + 
    theme(panel.grid = element_blank())
}

cor_dat = data.frame("Unit" = 1:n, "value" = abs(cor_re))
ggplot(cor_dat,aes(Unit,value)) + geom_path(alpha=0.8,color = "#4BA6A3") +
  geom_point(color="#978ED6") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0.1,1,0.2)) + 
  theme_bw() + ylab("Correlation") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  theme(panel.grid = element_blank())
# ggsave(paste("crack/result/cor_all.pdf",sep=""), height = 4, width = 5)


# 3. EM算法 ===================================================================

## 模型选择-参数估计
rt_seq = seq(0,100,0.1)

time111 = Sys.time()
item = 100; B_item = 500
source("em/nolinear_T.r")
Sys.time() - time111 
source("em/nolinear_Wiener.r")
# source("em/linear_T.r")
# source("em/linear_Wiener.r")

quan_nolinear_T
quan_nolinear_Wiener
# quan_linear_T
# quan_linear_Wiener

round(em_re_nolinear_T$para,3)
# round(em_re_linear_T$para,3)
em_re_nolinear_T$para
# em_re_linear_T$para
# AIC
c(em_re_nolinear_T$aic, em_re_nolinear_Wiener$aic)
# PARA 区间估计
para_quan = round(cbind(t(quan_nolinear_T),t(quan_nolinear_Wiener)),3)
# write.csv(para_quan,"result/crack_para_quan.csv")


# 拟合路径图(点估计) =======
em_para_final = as.numeric(quan_nolinear_T[2,]) 
Yhat_dat = list()
Yhat_dat[[1]] = matrix(NA,m+1,p)
scen = "exp"

gamma_hat = as.numeric(em_para_final[(length(em_para_final)-p+1):length(em_para_final)])
eta_hat = em_para_final[1:p]
for(j in 1:p){
  if(scen == "linear"){
    mean_t = eta_hat[j] * 0:m # eta * lambda(t) 计算平均退化量
  } else{
    if(scen == "exp"){
      mean_t = eta_hat[j] * (exp(gamma_hat[j] * t[j,])-1) # 非线性——exp
    } else{
      mean_t = eta_hat[j] * t[j,]^gamma_hat[j] # 非线性——power
    }
  }
  Yhat_dat[[1]][,j] = mean_t
  Yhat_dat[[1]] = data.frame(Yhat_dat[[1]])
  colnames(Yhat_dat[[1]]) = paste0("PC",1:p,seq="")
}
mean.path.fit.plot(data=Yhat_dat, true_data = LED_dat[[5]], leg.pos = "none",ech = 10, ci = FALSE) + scale_y_reverse()
# ggsave("result/path-fit-mean.pdf", height = 4, width = 6)

# 可靠度计算 ====================================================================================
# 计算 student-t 和维纳的可靠度
nolinear_T_R_dat = t(round(apply(nolinear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
# linear_T_R_dat = t(round(apply(linear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
nolinear_Wiener_R_dat = t(round(apply(nolinear_Wiener_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))

cdf.MW1 <- nolinear_Wiener_R_dat[,2] # 维纳
# cdf.MW1 = R_cal(para3 = em_re_nolinear_Wiener$para, r_SIG0 = em_re_nolinear_Wiener[[3]], B = B_item, yz = thr, scen ="exp")[[1]][, 2]
cdf.Mt1 <- nolinear_T_R_dat[,2]      # t分布-非线性
# cdf.Mt1 = R_cal(para3 = em_re_nolinear_T$para, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr, scen ="exp")[[1]][, 2]
# cdf.Mt2 <- linear_T_R_dat[,2]      # t分布-线性
# cdf.Mt1 = R_cal(para3 = em_re_nolinear_T$para, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr, scen ="exp")[[1]][, 2]


# cdf.Mt1 = R_cal(para3 = real, r_SIG0 = SIG0 , B = B_item, yz = thr)[[1]][, 2] #真实结果
cdf.90lcl.Mt1 <- nolinear_T_R_dat[,1]  
cdf.90ucl.Mt1 <- nolinear_T_R_dat[,3]

t2 <- rt_seq; t2.t <- rt_seq #时间

# 求出每个产品的失效时间
# y = a * t^b
if(scen=="exp"){
  model_pft <- function(x, a, b) {a * (exp(b*x)-1)}
} else{
  model_pft <- function(x, a, b) {a * x + b}
}

pft_coef_pc =  matrix(NA,n,2)
pft_pp = matrix(NA,n,p)
# pft_span = 100
# PC1 每个产品的拟合系数
for(i in 1:n){
  for(j in 1:p){
    if(scen=="exp"){
      pft_dat = data.frame("x" = t[j,-1], "y" = sim_cum_dat[[i]][-c(1),j]) # exp
    }else{
      pft_dat = data.frame("x" = log(t[j,-1]), "y" = log(sim_cum_dat[[i]][-c(1),j])) # power
      }
    pc1_mod = nls(y ~ model_pft(x, a, b), data = pft_dat, start = list(a = 1, b = 1)) #建立非线性模型
    pft_coef = as.numeric(coef(pc1_mod)) #保存系数
    if(scen=="exp"){
      pft_pp[i,j] = log(thr[j]/pft_coef[1] + 1)/pft_coef[2] # exp # 带入thr，求出伪失效时间
    } else{
      pft_pp[i,j] = exp((log(thr[j])-pft_coef[2]) / pft_coef[1]) # power # 带入thr，求出伪失效时间
    }
  }
}
pft = apply(pft_pp,1,min)
yqq <- 1 - (1:n-0.5) / n

# Anderson–Darling 检验
library(ADGofTest)
# cdf.MW1.PFT <- approx(t2, cdf.MW1, xout = sort(pft))$y 
# ad.test(cdf.MW1.PFT)
cdf.Mt1.PFT <- approx(t2.t, cdf.Mt1, xout = sort(pft))$y
ad.test(cdf.Mt1.PFT)
# cdf.Mt2.PFT <- approx(t2.t, cdf.Mt2, xout = sort(pft))$y
# ad.test(cdf.Mt2.PFT)


# pdf("case/crack/R-com-crack3.pdf",width=8, height=6)
# 画图
plot(t2, cdf.MW1, type = "l", col = "white", xlim = c(4,18),ylim = c(0,1), xlab = "Time", ylab = "Probability", lwd = 2, xaxt = 'n', lty = 2)
polygon(c(t2.t, rev(t2.t)), c(cdf.90lcl.Mt1, rev(cdf.90ucl.Mt1)), col = "skyblue", border = F)
axis(1, at = seq(4,20,by=4)) 
# axis(1, at = c(4e4, 6e4, 7e4, 8e4, 9e4, 2e5, 3e5, 4e5, 6e5, 7e5, 8e5, 9e5), tck = -0.01, label = FALSE)
points(sort(pft), yqq, pch = 19, cex = 0.8)
# lines(t2, cdf.MW1, lty = 2, lwd = 3, col = "red")
lines(t2.t, cdf.Mt1, col = "#440154", lwd = 3)
lines(t2.t, cdf.Mt2, lty = 2, lwd = 3, col = "red")
legend("topright", c(TeX(r'($M_{l}$)'), TeX(r'($M_p$)'), TeX(r'(90% CIs for $M_p$)'), "PFT"),
       col = c("red","#440154","skyblue", 1), lty = c(2, 1, 1, -1), 
       bty = "n", inset = .02, lwd = c(3, 3, 10, -1), pch = c(-1, -1, -1, 19), cex = 1.25) # 7 , 5.5
# legend("topright", c(TeX(r'($M_t$)'), TeX(r'(90% CIs for $M_t$)'), "PFT"), 
#        col = c("#440154","skyblue", 1), lty = c(1, 1, -1), bty = "n", inset = .02, lwd = c(3, 10, -1), pch = c(-1, -1, 19), cex = 1.25) # 7 , 5.5

# dev.off()
# save.image(file = paste("crack/result/0105-crack-new.RData", sep = ""))







