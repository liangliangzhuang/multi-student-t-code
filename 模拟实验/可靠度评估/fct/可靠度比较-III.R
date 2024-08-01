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
# n <- 20; m <- 20; p <- 2
v <- 5 # 自由度
time <- 1 ## 观测的时间间隔
eta <- c(11, 12) *1     # eta
delta = c(0.5, 1.5)  # delta sig.w 
scale <- c(1.1,1.2) # gamma=
SIG0 <- matrix(c(1, 0.75,0.75, 2), 2)
SIG0cor <- cov2cor(SIG0)
# 时间刻度设置
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)
t.scale <- t^scale
t.diff <- t.scale[, 2:(m+1)]-t.scale[, 1:m]
# 真实参数
real <- c(eta,delta,sqrt(diag(SIG0)),as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)]),v,scale)
# 真实的相关系数（公式4）
real_sigma = diag(SIG0)
real_rho = SIG0[1,2]/sqrt((real_sigma[1] + delta[1]/t.scale[1,])*(real_sigma[2] + delta[2]/t.scale[2,])) # 注意第一个为0

NL_Wiener_MTTF = NL_T_MTTF = L_Wiener_MTTF = L_T_MTTF = true_MTTF = numeric()
rho_hat = matrix(NA, itt, m+1)
for(qqq in 1:itt){
  # 1. 数据模拟 =========
  sim_dat = sim_path(par = real, v = v, SIG0 = SIG0, scen = "non-linear")
  # save(sim_dat,file= "case/PMB/sim_dat.RData")
  y.diff = sim_dat[[1]]; sumys = sim_dat[[2]]; y = sim_dat[[3]] #EM中常用数据
  sim_diff_dat = sim_dat[[4]] # 便于绘图（差分）
  sim_cum_dat = sim_dat[[5]] # 便于绘图（累计）
  
  
  thr = c(400,500) 
  dummy2 <- data.frame(PC = c("PC1","PC2"), Z = thr)
  degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") + 
    geom_hline(data = dummy2,aes(yintercept = Z),linetype = "dashed") 
  
  # 3. EM算法 =========
  ## 模型选择-参数估计
  rt_seq = seq(0,100,0.1)
  item = 1; B_item = 1000
  source("模拟实验/em/nolinear_T.r")
  para <- list("est_etas" = rep(8, p), "est_delta" = rep(1, p), "est_sig0" = diag(p), "est_v" = 4)
  em_re_nolinear_T <- tryCatch(EM_nonlinear_T(para = para, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys),
                               error = function(e) return(NA))
  
  # 估计的相关系数计算（公式4）=======
  opt_para = em_re_nolinear_T$para  #估计值保存(和I不同)
  sigma_12 = as.numeric(opt_para[7] * sqrt(opt_para[5]*opt_para[6]))
  opt_sigma = opt_para[c(5,6)]
  opt_delta = opt_para[c(3,4)]
  opt_gamma = opt_para[c(9,10)]
  Lambda_new = t^opt_gamma
  # Lambda_new = t^opt_gamma #因为是线性的，所以gamma为1
  rho_hat[qqq,] = sigma_12/sqrt((opt_sigma[1] + opt_delta[1]/Lambda_new[1,])*(opt_sigma[2] + opt_delta[2]/Lambda_new[2,]))
  
  source("模拟实验/em/nolinear_Wiener.r")
  em_re_nolinear_Wiener <- tryCatch(EM_nonlinear_Wiener(para = para, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys),
                                    error = function(e) return(NA))

  # 可靠度计算
  true_R = R_cal(para3 = real, r_t = rt_seq, r_SIG0 = SIG0, B = B_item, yz = thr)[[1]][, 2]
  if(all(is.na(em_re_nolinear_T))){
    NL_T_R = NA
  } else{NL_T_R = R_cal(para3 = em_re_nolinear_T[[2]], r_t = rt_seq, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr)[[1]][, 2]}
  if(all(is.na(em_re_nolinear_Wiener))){
    NL_Wiener_R = NA
  } else{NL_Wiener_R = R_cal(para3 = em_re_nolinear_Wiener[[2]], r_t = rt_seq, r_SIG0 = em_re_nolinear_Wiener[[3]], B = B_item, yz = thr)[[1]][, 2]}
  
  NL_Wiener_MTTF[qqq] = sum(diff(rt_seq)[1] * NL_Wiener_R)
  NL_T_MTTF[qqq] = sum(diff(rt_seq)[1] * NL_T_R)
  true_MTTF[qqq] = sum(diff(rt_seq)[1] * true_R)
   
  print(qqq)
}

# 结果汇总（可靠度）========
# all_MTTF = cbind(NL_Wiener_MTTF, NL_T_MTTF)
# true_MTTF2 = true_MTTF[-which(apply(is.na(all_MTTF), 1, any))]
# all_MTTF2 = matrix(na.omit(all_MTTF),ncol=2)
# #MAE
# mae = (all_MTTF2  - true_MTTF2)/true_MTTF2 *100 
# # RMSE
# rmse = sqrt((all_MTTF2  - true_MTTF2)^2/nrow(all_MTTF2)) * 100 
# rmse = data.frame(rmse)
# order1 = c("Wiener","Student-t")
# colnames(rmse) = order1
# rmse %>% pivot_longer(everything(),values_to = "value", names_to = "name") %>% 
#   ggplot(aes(x = fct_relevel(name, order1), y = value, fill = fct_relevel(name, order1))) + 
#   scale_fill_viridis(name ="Model", discrete = TRUE) + 
#   xlab("Model") + ylab(TeX(r'(RMSE($\times 10^{-2}))')) +
#   geom_boxplot() + theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") 


# 结果汇总（相关系数） =========
# MAE
mae_rho = sapply(1:itt, function(i) (rho_hat[i, ] - real_rho) / real_rho)
mae_mean = rowMeans(mae_rho, na.rm = TRUE)[-1] * 100
# RMSE
rmse_rho = sapply(1:itt, function(i) (rho_hat[i, ] - real_rho)^2)
rmse_mean = sqrt(rowMeans(rmse_rho, na.rm = TRUE))[-1] * 100
# save.image(file = paste("模拟实验/可靠度评估/result/III-n",n,"-m",m,".RData", sep = ""))
