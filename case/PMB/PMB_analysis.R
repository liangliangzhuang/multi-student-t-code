# ======================================================================
# ========================= 真实数据分析 PMB ===========================
# ======================================================================
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggsci)
library(reshape2)
library(latex2exp)
library(MASS)
source("sumqua.R")
source("case/case_fct.R")
options(scipen = 4) # 设置小数点
# 1. 参数设置 =====
n <- 30; m <- 20; p <- 2
v <- 4 # 自由度
time <- 1 ## 观测的时间间隔
eta <- c(5, 6) *1     # eta
delta = c(0.4, 0.6)  # delta sig.w 
scale <- c(1.1,1.2) # gamma=
SIG0 <- matrix(c(1, 0.75,0.75, 1.25), 2)
SIG0cor <- cov2cor(SIG0)
# 时间刻度设置
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)
t.scale <- t^scale
t.diff <- t.scale[, 2:(m+1)]-t.scale[, 1:m]
# 真实参数
real <- c(eta,delta,sqrt(diag(SIG0)),as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)]),v,scale)

# t <- m * time ## 最后一次观测的时间（退化量）
# t.diff <- time ## 观测时间间隔（实际情况中）
# sig <- diag(time, p) ## λ(t)=t,这里已经是差分时间了
# sig_inv <- solve(sig) ## 这里直接全局定义，因为模拟假设同时观测，观测间隔相等
# SIGMA <- diag(t, p) ## 时间尺度矩阵退化量

# 2. 数据模拟 =====
# sim_dat = sim_path(par = real, v = v, SIG0 = SIG0, scen = "non-linear")
# save(sim_dat,file= "case/PMB/sim_dat.RData")

y.diff = sim_dat[[1]]; sumys = sim_dat[[2]]; y = sim_dat[[3]] #EM中常用数据
sim_diff_dat = sim_dat[[4]] # 便于绘图（差分）
sim_cum_dat = sim_dat[[5]] # 便于绘图（累计）


# 2. 探索性分析 =========
# 2.1 绘制退化数据路径图
library(colorspace)
# thr = c(200,300)
thr = c(400,600) 
dummy2 <- data.frame(PC = c("PC1","PC2"), Z = thr)
g1 = degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") + 
  geom_hline(data = dummy2,aes(yintercept = Z),linetype = "dashed") 
g1
degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") +
  # scale_y_reverse() + 
  scale_color_viridis(discrete = TRUE) + 
  theme(legend.position = "none") +
  xlab(TeX(r'(Time(days$\times 3$))'))

# ggsave("case/PMB/PMB-dat-0617.pdf", height = 4, width = 5)
  
# 2.2 箱线图和 QQ 图
data2 = map(sim_cum_dat, ~ mutate(.x, Time = 0:(n()-1)))
merged_df2 <- bind_rows(data2, .id = "Unit")
cal_dat2 = merged_df2 %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
# write.csv(cal_dat2,row.names = FALSE, "case/result/cal_dat2.csv")
source("case/PMB/qqplot.R")

box_plot 
# ggsave(paste("case/PMB/result/box_plot.pdf",sep=""), height = 4, width = 5)
figure[[3]]
# ggsave(paste("case/PMB/result/qqplot/Time-",time2[i],".pdf",sep=""), height = 4, width = 11)

# 2.3 PC之间的相关性
sim3 = do.call(rbind, sim_diff_dat) #数据处理，将所有退化增量汇总
sim3 = data.frame(sim3)
sim3$n = rep(1:n, each = m)
# library(rayshader)
source("case/PMB/PMB-3D.R")

# p_real 
# 转化为3D图形（rayshader包）
# plot_gg(p_real, multicore = TRUE, raytrace = TRUE, width = 7, height = 4,
#         scale = 300, windowsize = c(1400, 866), zoom = 0.55, phi = 30, theta = 45)
# Sys.sleep(0.2)
# render_snapshot(clear = TRUE)
contour_plot
# ggsave("case/result/contour_plot2.pdf", contour_plot, height = 4, width = 5)
# ggsave("case/result/real-2D.pdf", p_real, height = 4, width = 5)
# ggsave("case/result/norm-2D.pdf", p_norm, height = 4, width = 5)


# 3. EM算法 ===================================================================

## 模型选择-参数估计
rt_seq = seq(0,100,0.1)

item = 50; B_item = 5000
source("case/PMB/em/nolinear_T.r")
# item = 100
source("case/PMB/em/linear_T.r")
# 主要比较上面两个
# item = 100
source("case/PMB/em/nolinear_Wiener.r")
# item = 100
source("case/PMB/em/linear_Wiener.r")

# quan_nolinear_T; quan_linear_T; quan_nolinear_Wiener; quan_linear_Wiener

# AIC
c(em_re_linear_T$aic, em_re_nolinear_T$aic, em_re_linear_Wiener$aic, em_re_nolinear_Wiener$aic)
# PARA 区间估计
para_quan = round(cbind(t(quan_linear_T),t(quan_nolinear_T),t(quan_linear_Wiener),t(quan_nolinear_Wiener)),3)
# write.csv(para_quan,"case/PMB/para_quan.csv")
para_quan


# 拟合路径图(点估计) =======
em_para_final = as.numeric(quan_nolinear_T[2,])
Yhat_dat = list()
Yhat_dat[[1]] = matrix(NA,m+1,p)
scen = "non-linear"
for(j in 1:p){
  if(scen == "linear"){
    mean_t = em_para_final[j] * 0:m # eta * lambda(t) 计算平均退化量
  } else{
    gamma_hat = as.numeric(em_para_final[(length(em_para_final)-p+1):length(em_para_final)])
    mean_t = em_para_final[j] * t[j,]^gamma_hat[j] # 非线性
  }
  Yhat_dat[[1]][,j] = mean_t
  Yhat_dat[[1]] = data.frame(Yhat_dat[[1]])
  colnames(Yhat_dat[[1]]) = paste0("PC",1:p,seq="")
}
mean.path.fit.plot(data=Yhat_dat, true_data = sim_dat[[5]], leg.pos = "none",ech = 10, ci = FALSE) +
  xlab(TeX(r'(Time(days$\times 3$))')) + ylab(TeX(r'($Y_{k}(t)$)'))
# ggsave("case/PMB/crack-path-fit-mean2.pdf", height = 3, width = 4)

# 相关系数 ===========================================================================
# 提取sigma_{1,2}
# crack 这部分计算有些问题，应该用自助法的所有结果先计算相关系数，然后求分位数。
rho_hat = matrix(NA,m+1,3)
opt_para = para_quan[,c(4:6)] 
for(i in 1:3){
  sigma_12 = as.numeric(opt_para[7,i] * sqrt(opt_para[5,i]*opt_para[6,i]))
  opt_sigma = opt_para[c(5,6),i]
  opt_delta = opt_para[c(3,4),i]
  opt_gamma = opt_para[c(9,10),i]
  Lambda_new = t^opt_gamma
  rho_hat[,i] = sigma_12/sqrt((opt_sigma[1] + opt_delta[1]/Lambda_new[1,])*(opt_sigma[2] + opt_delta[2]/Lambda_new[2,]))
}
rho_hat = data.frame(cbind(1:m+1, rho_hat))
colnames(rho_hat) = c("Time","5%","50%","95%")
rho_hat %>% pivot_longer(!Time, names_to="class", values_to = "value") %>% 
  ggplot(aes(Time,value,color = fct(class),linetype = fct(class))) + 
  geom_path() + geom_point() +
  scale_color_manual(name="",values = c("#008580","#693476","#008580"), 
                     labels=c("90% CI", "Estimated coeff.", "90% CI")) + 
  scale_linetype_manual(name="",values = c("dashed", "solid", "dashed"), labels=c("90% CI", "Estimated coeff.", "90% CI")) +
  theme_bw() + theme(panel.grid = element_blank(),legend.position = c(0.7,0.25)) +
  ylab(TeX(r'(RMSE($\times 10^{-2}))')) + xlab(TeX(r'(Time(days$\times 3$))'))
# ggsave("case/PMB/rho-hat2.pdf", height = 3, width = 4)


# 可靠度计算 ====================================================================================
thr = c(400,600) 
# em_re_nolinear_T
# R_cal(r_t = seq(0,100,0.1), para3=em_re_nolinear_T$para, B = 500, yz=thr)
  
# 计算t和维纳的可靠度
dim(nolinear_T_R)
nolinear_T_R_dat = t(round(apply(nolinear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
nolinear_Wiener_R_dat = t(round(apply(nolinear_Wiener_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))

# cdf.MW1 <- nolinear_Wiener_R_dat[,2] # 维纳
cdf.MW1 = R_cal(para3 = em_re_nolinear_Wiener$para, r_SIG0 = em_re_nolinear_Wiener[[3]], B = B_item, yz = thr)[[1]][, 2]
cdf.Mt1 <- nolinear_T_R_dat[,2]      # t分布
# cdf.Mt1 = R_cal(para3 = em_re_nolinear_T$para, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr)[[1]][, 2]
# cdf.Mt1 = R_cal(para3 = real, r_SIG0 = SIG0 , B = B_item, yz = thr)[[1]][, 2] #真实结果
cdf.90lcl.Mt1 <- nolinear_T_R_dat[,1]  
cdf.90ucl.Mt1 <- nolinear_T_R_dat[,3]

t2 <- rt_seq; t2.t <- rt_seq #时间

# 求出每个产品的失效时间
# y = a * t^b
model_pft <- function(x, a, b) {a * x + b}

pft_coef_pc =  matrix(NA,n,2)
pft_pp = matrix(NA,n,p)
# pft_span = 100
# PC1 每个产品的拟合系数
for(i in 1:n){
  for(j in 1:p){
    pft_dat = data.frame("x" = log(t[j,-1]), "y" = log(sim_cum_dat[[i]][-1,j])) 
    pc1_mod = nls(y ~ model_pft(x, a, b), data = pft_dat, start = list(a = 1, b = 1)) #建立非线性模型
    pft_coef = as.numeric(coef(pc1_mod)) #保存系数
    pft_pp[i,j] = exp((log(thr[j])-pft_coef[2]) / pft_coef[1]) # 带入thr，求出伪失效时间
  }
}
pft = apply(pft_pp,1,min)
yqq <- 1 - (1:n-0.5) / n


# Anderson–Darling 检验
library(ADGofTest)
cdf.MW1.PFT <- approx(t2, cdf.MW1, xout = sort(pft))$y 
ad.test(cdf.MW1.PFT)
cdf.Mt1.PFT <- approx(t2.t, cdf.Mt1, xout = sort(pft))$y
ad.test(cdf.Mt1.PFT)


# pdf("case/PMB/R-com-pmb.pdf",width=8, height=6)
# 画图
plot(t2, cdf.MW1, type = "l", col = "white", xlim = c(20,60),ylim = c(0,1), xlab = "Time", ylab = "Probability", lwd = 2, xaxt = 'n', lty = 2)
polygon(c(t2.t, rev(t2.t)), c(cdf.90lcl.Mt1, rev(cdf.90ucl.Mt1)), col = "skyblue", border = F)
axis(1, at = seq(20,60,by=10)) 
# axis(1, at = c(4e4, 6e4, 7e4, 8e4, 9e4, 2e5, 3e5, 4e5, 6e5, 7e5, 8e5, 9e5), tck = -0.01, label = FALSE)
points(sort(pft), yqq, pch = 19, cex = 0.8)
lines(t2, cdf.MW1, lty = 2, lwd = 3, col = "red")
lines(t2.t, cdf.Mt1, col = "#440154", lwd = 3)
legend("topright", c(TeX(r'($M^W_{p}$)'), TeX(r'($M_p$)'), TeX(r'(90% CIs for $M_p$)'), "PFT"), 
       col = c("red","#440154", "skyblue", 1), lty = c(2, 1, 1, -1), bty = "n", inset = .02, lwd = c(3, 3, 10, -1), pch = c(-1, -1, -1, 19), cex = 1.25) # 7 , 5.5

# dev.off()
# save.image(file = paste("case/PMB/0619-PMB-final.RData", sep = ""))



