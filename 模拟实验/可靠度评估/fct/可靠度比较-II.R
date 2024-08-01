# ============ 四种方法的可靠度比较（MTTF） ============
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

# itt=100

# 参数设置 ========
# n <- 20; m <- 30; p <- 3
v <- 5 # 自由度
time <- 1 ## 观测的时间间隔
eta <- c(11, 12, 13)    # eta
delta = c(0.4, 0.6, 0.8)  # delta sig.w 
scale <- c(1, 1, 1) # gamma=
SIG0 <- matrix(c(1, 0.75, -1, 0.75, 2.25, 1.2, -1, 1.2, 4), p)
SIG0cor <- cov2cor(SIG0)
# 时间刻度设置
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)
t.scale <- t^scale
t.diff <- t.scale[, 2:(m+1)]-t.scale[, 1:m]
# 真实参数
thr = c(250,250,250) 
real <- c(eta,delta,sqrt(diag(SIG0)),as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)]),v,scale)
# 真实的相关系数（公式4） p=3
real_sigma = diag(SIG0)
real_rho12 = SIG0[1,2]/sqrt((real_sigma[1] + delta[1]/t.scale[1,])*(real_sigma[2] + delta[2]/t.scale[2,])) # 注意第一个为0
real_rho13 = SIG0[1,3]/sqrt((real_sigma[1] + delta[1]/t.scale[1,])*(real_sigma[3] + delta[3]/t.scale[3,])) # 注意第一个为0
real_rho23 = SIG0[2,3]/sqrt((real_sigma[2] + delta[2]/t.scale[2,])*(real_sigma[3] + delta[3]/t.scale[3,])) # 注意第一个为0

NL_Wiener_MTTF = NL_T_MTTF = L_Wiener_MTTF = L_T_MTTF = true_MTTF = numeric()
rho_hat12 = rho_hat13 = rho_hat23 = matrix(NA, itt, m+1)
for(qqq in 1:itt){
  # 1. 数据模拟 =========
  sim_dat = sim_path(par = real, v = v, SIG0 = SIG0, scen = "linear")
  # save(sim_dat,file= "case/PMB/sim_dat.RData")
  y.diff = sim_dat[[1]]; sumys = sim_dat[[2]]; y = sim_dat[[3]] #EM中常用数据
  sim_diff_dat = sim_dat[[4]] # 便于绘图（差分）
  sim_cum_dat = sim_dat[[5]] # 便于绘图（累计）
  
  dummy2 <- data.frame(PC = paste("PC",1:p,sep=""), Z = thr)
  degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") + 
    geom_hline(data = dummy2,aes(yintercept = Z),linetype = "dashed") 
  
  # 3. EM算法 =========
  ## 模型选择-参数估计
  rt_seq = seq(0,100,0.1)
  item = 1; B_item = 1000
  # source("模拟实验/em/nolinear_T.r")
  # source("模拟实验/em/nolinear_Wiener.r")
  source("模拟实验/em/linear_T.r")
  # 估计的相关系数计算（公式4）=======
  opt_para = em_re_linear_T$para  #估计值保存
  sigma_12 = as.numeric(opt_para[10] * sqrt(opt_para[7]*opt_para[8])) # rho_12 / sqrt(simga^2_1 * simga^2_2)
  sigma_13 = as.numeric(opt_para[11] * sqrt(opt_para[7]*opt_para[9])) # rho_13 / sqrt(simga^2_1 * simga^2_3)
  sigma_23 = as.numeric(opt_para[12] * sqrt(opt_para[8]*opt_para[9])) # rho_23 / sqrt(simga^2_2 * simga^2_3)
  opt_sigma = opt_para[7:9]
  opt_delta = opt_para[4:6]
  opt_gamma = opt_para[14:16] # 13是自由度
  Lambda_new = t
  # Lambda_new = t^opt_gamma #因为是线性的，所以gamma为1
  rho_hat12[qqq,] = sigma_12/sqrt((opt_sigma[1] + opt_delta[1]/Lambda_new[1,])*(opt_sigma[2] + opt_delta[2]/Lambda_new[2,]))
  rho_hat13[qqq,] = sigma_13/sqrt((opt_sigma[1] + opt_delta[1]/Lambda_new[1,])*(opt_sigma[3] + opt_delta[3]/Lambda_new[3,]))
  rho_hat23[qqq,] = sigma_23/sqrt((opt_sigma[2] + opt_delta[2]/Lambda_new[2,])*(opt_sigma[3] + opt_delta[3]/Lambda_new[3,]))
  
  source("模拟实验/em/linear_Wiener.r")
  
  # 可靠度计算 =======
  true_R = R_cal(scen = "linear",para3 = real, r_t = rt_seq, r_SIG0 = SIG0, B = B_item, yz = thr)[[1]][, 2]
  # NL_Wiener_R = R_cal(para3 = em_re_nolinear_Wiener[[2]], r_t = rt_seq, r_SIG0 = em_re_nolinear_Wiener[[3]], B = B_item, yz = thr)[[1]][, 2]
  # NL_T_R = R_cal(para3 = em_re_nolinear_T[[2]], r_t = rt_seq, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr)[[1]][, 2]
  L_Wiener_R = R_cal(scen = "linear",para3 =  em_re_linear_Wiener[[2]], r_t = rt_seq, r_SIG0 = em_re_linear_Wiener[[3]], B = B_item, yz = thr)[[1]][, 2]
  L_T_R = R_cal(scen = "linear",para3 = em_re_linear_T[[2]], r_t = rt_seq, r_SIG0 = em_re_linear_T[[3]], B = B_item, yz = thr)[[1]][, 2]
  
  # NL_Wiener_MTTF[qqq] = sum(diff(rt_seq)[1] * NL_Wiener_R)
  # NL_T_MTTF[qqq] = sum(diff(rt_seq)[1] * NL_T_R)
  L_Wiener_MTTF[qqq] = sum(diff(rt_seq)[1] * L_Wiener_R)
  L_T_MTTF[qqq] = sum(diff(rt_seq)[1] * L_T_R)
  true_MTTF[qqq] = sum(diff(rt_seq)[1] * true_R)
  
  print(qqq)
}

# 结果汇总（可靠度）======
all_MTTF = cbind(L_Wiener_MTTF, L_T_MTTF)
#MAE
mae = (all_MTTF  - true_MTTF)/true_MTTF *100 
# RMSE
rmse = sqrt((all_MTTF  - true_MTTF)^2/nrow(all_MTTF)) * 100 
rmse = data.frame(rmse)
order1 = c("Wiener","Student-t")
colnames(rmse) = order1
rmse %>% pivot_longer(everything(),values_to = "value", names_to = "name") %>% 
  ggplot(aes(x = fct_relevel(name, order1), y = value, fill = fct_relevel(name, order1))) + 
  scale_fill_viridis(name ="Model", discrete = TRUE) + 
  xlab("Model") + ylab(TeX(r'(RMSE($\times 10^{-2}))')) +
  geom_boxplot() + theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") 

# 结果汇总（相关系数） =========
# MAE
mae_rho12 = sapply(1:itt, function(i) (rho_hat12[i, ] - real_rho12) / real_rho12)
mae_mean12 = rowMeans(mae_rho12, na.rm = TRUE)[-1] * 100
mae_rho13 = sapply(1:itt, function(i) (rho_hat13[i, ] - real_rho13) / real_rho13)
mae_mean13 = rowMeans(mae_rho13, na.rm = TRUE)[-1] * 100
mae_rho23 = sapply(1:itt, function(i) (rho_hat23[i, ] - real_rho23) / real_rho23)
mae_mean23 = rowMeans(mae_rho23, na.rm = TRUE)[-1] * 100
# RMSE
rmse_rho12 = sapply(1:itt, function(i) (rho_hat12[i, ] - real_rho12)^2)
rmse_mean12 = sqrt(rowMeans(rmse_rho12, na.rm = TRUE))[-1] * 100
rmse_rho13 = sapply(1:itt, function(i) (rho_hat13[i, ] - real_rho13)^2)
rmse_mean13 = sqrt(rowMeans(rmse_rho13, na.rm = TRUE))[-1] * 100
rmse_rho23 = sapply(1:itt, function(i) (rho_hat23[i, ] - real_rho23)^2)
rmse_mean23 = sqrt(rowMeans(rmse_rho23, na.rm = TRUE))[-1] * 100

rmse_mean = rbind(rmse_mean12,rmse_mean13,rmse_mean23)
# save.image(file = paste("模拟实验/可靠度评估/result/II-n",n,"-m",m,".RData", sep = ""))




