# ========= 真实数据分析(Fatigue Crack-size) =========
library(ggplot2)
library(tidyverse)
library(viridis)
library(reshape2)
library(latex2exp)
library(MASS)
library(openxlsx)
source("EM-alg.R")
source("case/case_fct.R")
options(scipen = 4) # 设置小数点

# 1.数据读取 =======
fc_dat = fc_dat1 = fc_dat_sub = fc_diff_dat = list()
n = 6
for (i in 1:n) {
  fc_dat[[i]] = read.xlsx("case/crack/Fatigue-crack-size.xlsx",sheet = i) 
  fc_dat_sub[[i]] =  fc_dat[[i]][,] - 0.9
  fc_dat1[[i]] =  fc_dat[[i]][-1,] - 0.9
  fc_diff_dat[[i]] = apply(fc_dat[[i]],2,diff) 
}
m <- dim(fc_diff_dat[[1]])[1] 
p = dim(fc_diff_dat[[1]])[2]

time <- 1 ## 观测的时间间隔
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)

# 数据处理
crack_dat = crack_path(led_cum = fc_dat_sub, led_diff = fc_diff_dat, p=p, m=m, n=n)
y.diff = crack_dat[[1]]; 
sumys = crack_dat[[2]]; 
y =  crack_dat[[3]] # EM 中常用数据
sim_diff_dat = crack_dat[[4]] # 便于绘图（差分）
sim_cum_dat = crack_dat[[5]] # 便于绘图（累计）

# 2.探索性分析 ========
# 2.1 绘制退化数据路径图 ====
thr = c(0.9,0.5,0.4)  
dummy2 <- data.frame(PC = c("PC1","PC2","PC3"), Z = thr)
g1 = degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") + 
  geom_hline(data = dummy2,aes(yintercept = Z),linetype = "dashed") 
g1
degradation.path.plot(data = sim_cum_dat, leg.pos = "bottom",ech = 5, scale1 = "fixed") +
  scale_color_viridis(discrete = TRUE, name = "Unit",guide = guide_legend(nrow = 1)) + 
  xlab(TeX(r'(Time (hours $\times$ 336))'))  + 
  ylab(TeX(r'(Y(t) (inches $\times 10^{-2}$))') ) 
  # theme(legend.box = "horizontal")
# ggsave("case/crack/crack-dat.pdf", height = 5, width = 9)

# 2.2 箱线图和 QQ 图 ====
data2 = map(sim_cum_dat, ~ mutate(.x, Time = 0:(n()-1)))
merged_df2 <- bind_rows(data2, .id = "Unit")
tidy_crack_dat = merged_df2 %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
# write.csv(tidy_crack_dat,row.names = FALSE, "case/crack/result/tidy_crack_dat.csv")
source("case/crack/qqplot.R")
box_plot + xlab("Millions of cycles")
# ggsave(paste("case/crack/box_plot_new.pdf",sep=""),  height = 4, width = 8)


# 2.3 PC之间的相关性 ====
cor_re = as.numeric()
for(i in 1:n){
  cor_cal = data.frame("PC1" = sim_diff_dat[[i]][,1], "PC2" = sim_diff_dat[[i]][,2])
  cor_re[i] = cor(cor_cal$PC2, cor_cal$PC1)
  ggplot(cor_cal,aes(x=PC1,y=PC2)) + 
    geom_point(size=1.3,shape=19,color = "#4BA6A3")+
    geom_smooth(method=lm,color ="#978ED6",alpha=0.3) + 
    theme_bw() + 
    theme(panel.grid = element_blank())
  # ggsave(paste("第三章/case/result/cor.pdf",sep=""), height = 4, width = 5)
}
mean(cor_re) #平局相关系数 0.845

cor_dat = data.frame("Unit" = 1:n, "value" = abs(cor_re))
ggplot(cor_dat,aes(Unit,value)) + geom_path(alpha=0.8,color = "#4BA6A3") +
  geom_point(color="#978ED6") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0.1,1,0.2)) + 
  theme_bw() + ylab("Correlation") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  theme(panel.grid = element_blank())
# ggsave(paste("第三章/case/result/cor_all.pdf",sep=""), height = 4, width = 5)


# 3. EM算法 ===================================================================

# R_cal(para3 = em_re_nolinear_T$para, B = 10000, yz = thr); g1
## 模型选择-参数估计
rt_seq = seq(0,100,0.1)

item = 500; B_item = 1000
t111 = Sys.time()
source("case/crack/em/nolinear_T.r") # 8min -> item = 500; B_item = 1000
Sys.time() - t111
source("case/crack/em/nolinear_Wiener.r")

# 其他两种情况结果就不计算了
# item = 10
# source("case/crack/em/linear_T.r")
# source("case/crack/em/linear_Wiener.r")

quan_nolinear_T
# quan_linear_T
quan_nolinear_Wiener
# quan_linear_Wiener


round(em_re_nolinear_T$para,3)
round(em_re_linear_T$para,3)
em_re_nolinear_T$para
em_re_linear_T$para

# AIC
# c(em_re_linear_T$aic, em_re_nolinear_T$aic, em_re_linear_Wiener$aic, em_re_nolinear_Wiener$aic)
c(em_re_nolinear_T$aic, em_re_nolinear_Wiener$aic) #表3
# PARA 区间估计
para_quan = round(cbind(t(quan_nolinear_T),t(quan_nolinear_Wiener)),3)
# para_quan = round(cbind(t(quan_linear_T),t(quan_nolinear_T),t(quan_linear_Wiener),t(quan_nolinear_Wiener)),3)
# write.csv(para_quan,"case/crack/crack_para_quan.csv")


# 4. 拟合路径图(点估计) =======
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
mean.path.fit.plot(data=Yhat_dat, true_data = LED_dat[[5]], leg.pos = "none",ech = 10, ci = FALSE) + 
  scale_y_reverse() + xlab("Millions of cycles")
# ggsave("case/crack/result/path-fit-mean.pdf", height = 4, width = 6)

# 5. 相关系数绘制 ===========================================================================
# 提取sigma_{1,2}
opt_para = para_quan[,c(1:3)]
# 参数整理（自助法所有的结果，为了后续计算相关系数的分位数）
RHO = BT_ALL[, grep("rho", colnames(BT_ALL), value = TRUE)]
SIMGA = BT_ALL[, grep("sigma", colnames(BT_ALL), value = TRUE)]
DELTA = BT_ALL[, grep("delta", colnames(BT_ALL), value = TRUE)]
GAMMA = BT_ALL[, grep("gamma", colnames(BT_ALL), value = TRUE)]
# sigma
sigma_12 = as.numeric(RHO[,1] * sqrt(SIMGA[,1] * SIMGA[,2]))
sigma_13 = as.numeric(RHO[,2] * sqrt(SIMGA[,1] * SIMGA[,3]))
sigma_23 = as.numeric(RHO[,3] * sqrt(SIMGA[,2] * SIMGA[,3]))
# opt
opt_sigma = SIMGA; opt_delta = DELTA; opt_gamma = GAMMA
Lambda_new = list()
for(k in 1:dim(opt_gamma)[1]){Lambda_new[[k]] = exp(opt_gamma[k,] * t) -1 } # exp(gamma * t)-1

rho12_hat = rho13_hat = rho23_hat = matrix(NA, m+1, dim(opt_gamma)[1]) #计算rho(t)的结果，行为时间，列为自助法次数
for(k in 1:dim(opt_gamma)[1]){
  rho12_hat[,k] = sigma_12[k]/sqrt((opt_sigma[k,1] + opt_delta[k,1]/Lambda_new[[k]][1,]) * (opt_sigma[k,2] + opt_delta[k,2]/Lambda_new[[k]][2,]))
  rho13_hat[,k] = sigma_13[k]/sqrt((opt_sigma[k,1] + opt_delta[k,1]/Lambda_new[[k]][1,]) * (opt_sigma[k,3] + opt_delta[k,3]/Lambda_new[[k]][3,]))
  rho23_hat[,k] = sigma_23[k]/sqrt((opt_sigma[k,2] + opt_delta[k,2]/Lambda_new[[k]][2,]) * (opt_sigma[k,3] + opt_delta[k,3]/Lambda_new[[k]][3,]))
}
rho12_hat_qua = apply(rho12_hat, 1, quantile, c(0.25, 0.5, 0.95), na.rm = TRUE)
rho13_hat_qua = apply(rho13_hat, 1, quantile, c(0.25, 0.5, 0.95), na.rm = TRUE)
rho23_hat_qua = apply(rho23_hat, 1, quantile, c(0.25, 0.5, 0.95), na.rm = TRUE)
# rho 整理
rho_hat = data.frame(cbind(1:(m+1), t(rho12_hat_qua), t(rho13_hat_qua), t(rho23_hat_qua)))
colnames(rho_hat) = c("Time", paste('rho',rep(1:3,each=3),"_",c("5%","50%","95%"),sep =""))
# 分面绘制会使用到
f_names <- list('rho1' = TeX(c("$\\rho_{12}(t)$")), 'rho2' = TeX(c("$\\rho_{13}(t)$")),
                'rho3' = TeX(c("$\\rho_{23}(t)$")))
f_labeller <- function(variable, value){return(f_names[value])}
# 绘图
rho_hat %>% pivot_longer(!Time, names_to = c("rho", "percentile"), # 分为 rho 和 percentile 两个新列
                         names_sep = "_", # 列名中 rho 和百分比之间以 _ 分隔
                         values_to = "value") %>% 
  ggplot(aes(Time, value, group = interaction(rho, percentile))) + 
  geom_path(aes(color=percentile, linetype = percentile)) + 
  geom_point(aes(color=percentile)) +
  scale_color_manual(name="",values = c("#008580","#693476","#008580"), 
                     labels=c("90% CI", "Estimated coeff.", "90% CI")) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed"), 
                        name = "", labels = c("90% CI", "Estimated coeff.", "90% CI")) +
  facet_wrap(~ rho, scales = "free_y",labeller = f_labeller) +
  scale_x_continuous(breaks = c(0,5,10)) +
  scale_y_continuous(limits = c(0,0.8)) +
  theme_bw() + theme(panel.grid = element_blank(),legend.position = "top") +
  ylab(TeX(r'(RMSE($\times 10^{-2}))')) + xlab("Millions of cycles")
# ggsave("case/crack/result/rho_crack.pdf", height = 4, width = 8)


# 可靠度计算 ====================================================================================

# em_re_nolinear_T
# R_cal(r_t = seq(0,100,0.1), para3=em_re_nolinear_T$para, B = 500, yz=thr)

# 计算t和维纳的可靠度
dim(nolinear_T_R)

nolinear_T_R_dat = t(round(apply(nolinear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
# linear_T_R_dat = t(round(apply(linear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
nolinear_Wiener_R_dat = t(round(apply(nolinear_Wiener_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))

cdf.MW1 <- nolinear_Wiener_R_dat[,2] # 维纳
# cdf.MW1 = R_cal(para3 = em_re_nolinear_Wiener$para, r_SIG0 = em_re_nolinear_Wiener[[3]], B = B_item, yz = thr, scen ="exp")[[1]][, 2]
cdf.Mt1 <- nolinear_T_R_dat[,2]      # t分布-非线性
# cdf.Mt1 = R_cal(para3 = em_re_nolinear_T$para, r_SIG0 = em_re_nolinear_T[[3]], B = B_item, yz = thr, scen ="exp")[[1]][, 2]
# cdf.Mt2 <- linear_T_R_dat[,2]      # t分布-非线性
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
cdf.MW1.PFT <- approx(t2, cdf.MW1, xout = sort(pft))$y
ad.test(cdf.MW1.PFT)
cdf.Mt1.PFT <- approx(t2.t, cdf.Mt1, xout = sort(pft))$y
ad.test(cdf.Mt1.PFT)

# pdf("case/crack/R-com-crack3.pdf",width=8, height=6)
# 画图
plot(t2, cdf.MW1, type = "l", col = "white", xlim = c(4,18),ylim = c(0,1), xlab = "Time", ylab = "Probability", lwd = 2, xaxt = 'n', lty = 2)
polygon(c(t2.t, rev(t2.t)), c(cdf.90lcl.Mt1, rev(cdf.90ucl.Mt1)), col = "skyblue", border = F)
axis(1, at = seq(4,20,by=4)) 
# axis(1, at = c(4e4, 6e4, 7e4, 8e4, 9e4, 2e5, 3e5, 4e5, 6e5, 7e5, 8e5, 9e5), tck = -0.01, label = FALSE)
points(sort(pft), yqq, pch = 19, cex = 0.8)
lines(t2, cdf.MW1, lty = 2, lwd = 3, col = "red")
lines(t2.t, cdf.Mt1, col = "#440154", lwd = 3)
# lines(t2.t, cdf.Mt2, lty = 2, lwd = 3, col = "red")
legend("topright", c(TeX(r'($M^W_{e}$)'), TeX(r'($M_e$)'), TeX(r'(90% CIs for $M_e$)'), "PFT"),
       col = c("red","#440154","skyblue", 1), lty = c(2, 1, 1, -1), 
       bty = "n", inset = .02, lwd = c(3, 3, 10, -1), pch = c(-1, -1, -1, 19), cex = 1.25) # 7 , 5.5
# legend("topright", c(TeX(r'($M_t$)'), TeX(r'(90% CIs for $M_t$)'), "PFT"), 
#        col = c("#440154","skyblue", 1), lty = c(1, 1, -1), bty = "n", inset = .02, lwd = c(3, 10, -1), pch = c(-1, -1, 19), cex = 1.25) # 7 , 5.5

# dev.off()
# save.image(file = paste("case/crack/0619-crack.RData", sep = ""))






