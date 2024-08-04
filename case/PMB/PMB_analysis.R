# ======================================================================
# ========================= PMB Data Analysis  =========================
# ======================================================================
# Load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, viridis, reshape2, latex2exp, MASS, openxlsx)

source("case/case_fct.R")
source("case/PMB/vis.R")
source("Model_est.R")
options(scipen = 4) # Set decimal point precision

# 1. Data Loading =======
load(file= "case/PMB/pmb_dat.RData")
y.diff = sim_dat[[1]]; sumys = sim_dat[[2]]; y = sim_dat[[3]] 
sim_diff_dat = sim_dat[[4]] 
sim_cum_dat = sim_dat[[5]] 

p <- dim(y.diff)[1]; m <- dim(y.diff)[2]; n <- dim(y.diff)[3];
time <- 1 
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)

# 2. Exploratory Analysis ========
# 2.1 Plot degradation data paths ====
degradation.path.plot(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") +
  scale_color_viridis(discrete = TRUE) + 
  theme(legend.position = "none") +
  xlab(TeX(r'(Time(days$\times 3$))'))
# ggsave("case/result/PMB/path-dat.pdf", height = 4, width = 5)
  
# 2.2 Boxplot of degradation data
data2 = map(sim_cum_dat, ~ mutate(.x, Time = 0:(n()-1)))
merged_df2 <- bind_rows(data2, .id = "Unit")
cal_dat2 = merged_df2 %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
boxplot = boxplot_path(time2 = seq(3,m,3), data = cal_dat2)
boxplot[[1]] + xlab(TeX(r'(Time(days$\times 3$))')) + ylab("Value")
# ggsave("case/result/PMB/boxplot.pdf", height = 4, width = 5)

# 2.3 Contour plot
sim3 = do.call(rbind, sim_diff_dat) 
sim3 = data.frame(sim3)
sim3$n = rep(1:n, each = m)
coutour_plot(sim3 = sim3)

# 2.4 Q-Q plot
qqplot_PC(time2 = seq(3,m,3), newdata = boxplot[[2]], newdata2 = boxplot[[3]])

# 3. Parameter Estimation for Different Models ====================

# 3.1 Point estimate ======
types = "power"
# Initial parameter 
para <- list("est_etas" = rep(8, p), "est_delta" = rep(1, p), "est_sig0" = diag(p), "est_v" = 4)
f1_names <- list('eta1' = TeX(c("$\\hat{\\eta}_{1}$")), 'eta2' = TeX(c("$\\hat{\\eta}_{2}$")),
                 'delta1' = TeX(c("$\\hat{\\delta}^2_{1}$")), 'delta2' = TeX(c("$\\hat{\\delta}^2_{2}$")), 
                 'sigma1' = TeX(c("$\\hat{\\sigma}^2_{1}$")), 'sigma2' = TeX(c("$\\hat{\\sigma}^2_{2}$")),
                 'rho12' = expression(hat(rho)[12]),
                 'v' = TeX(c("$\\hat{\\nu}$")))
index = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
          paste0("sigma",1:p, sep = ""), "rho12", "v")
orders = unique(index)
# M_l
em_re_linear_T = EM_linear_T(max_iter = 5000, eps = 10^-5, y = y, y.diff= y.diff, sumys = sumys)
em_re_linear_T_no_na <- data.frame(em_re_linear_T$para_iter[complete.cases(em_re_linear_T$para_iter), ])
colnames(em_re_linear_T_no_na) = index
EM_iter_plot(para_iter = em_re_linear_T_no_na, f_names = f1_names, orders=orders)
# ggsave("case/result/PMB/PMB-EM-iter-M_l.pdf", height = 3, width = 11)

# M_p
em_re_nolinear_T <- EM_nonlinear_T(para = para, type = types, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
em_re_nolinear_T_no_na <- em_re_nolinear_T$para_iter[complete.cases(em_re_nolinear_T$para_iter), ]
colnames(em_re_nolinear_T_no_na) = index
EM_iter_plot(para_iter = em_re_nolinear_T_no_na, f_names = f1_names, orders=orders)
# ggsave("case/result/PMB/PMB-EM-iter-M_p.pdf", height = 3, width = 11)

# M^W_l
em_re_linear_Wiener = EM_linear_Wiener(max_iter = 5000, eps = 10^-5, y = y, y.diff= y.diff, sumys = sumys)
em_re_linear_Wiener_no_na <- em_re_linear_Wiener$para_iter[complete.cases(em_re_linear_Wiener$para_iter), ]
colnames(em_re_linear_Wiener_no_na) = index
EM_iter_plot(para_iter = em_re_linear_Wiener_no_na[,-8], f_names = f1_names, orders=orders) # omit v
# ggsave("case/result/PMB/PMB-EM-iter-M^W_l.pdf", height = 3, width = 11)

# M^W_p
em_re_nolinear_Wiener <- EM_nonlinear_Wiener(para = para, type=types, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
em_re_nolinear_Wiener_no_na <- em_re_nolinear_Wiener$para_iter[complete.cases(em_re_nolinear_Wiener$para_iter), ]
colnames(em_re_nolinear_Wiener_no_na) = index
EM_iter_plot(para_iter = em_re_nolinear_Wiener_no_na[,-8],f_names = f1_names, orders=orders)
# ggsave("case/result/PMB/PMB-EM-iter-M^W_p.pdf", height = 3, width = 11)

# AIC
c(em_re_linear_T$aic, em_re_nolinear_T$aic, em_re_linear_Wiener$aic, em_re_nolinear_Wiener$aic)
# Point estimates
cbind(em_re_linear_T$para, em_re_nolinear_T$para, em_re_linear_Wiener$para, em_re_nolinear_Wiener$para)

# Interval estimate (Bootstrap) =========
rt_seq = seq(0,100,0.1)
thr = c(400,600) 
item = 500; B_item = 5000
source("case/PMB/Bootstrap/BS_nolinear_T.r")
source("case/PMB/Bootstrap/BS_linear_T.r")
source("case/PMB/Bootstrap/BS_linear_Wiener.r")
source("case/PMB/Bootstrap/BS_nolinear_Wiener.r")

para_quan = round(cbind(t(quan_linear_T),t(quan_nolinear_T),
                        t(quan_linear_Wiener),t(quan_nolinear_Wiener)),3)
# write.csv(para_quan,"case/result/PMB/para_quan.csv")
para_quan

# Fitting path graph (point estimation) =======
em_para_final = as.numeric(quan_nolinear_T[2,])
Yhat_dat = list()
Yhat_dat[[1]] = matrix(NA,m+1,p)
for(j in 1:p){
  gamma_hat = as.numeric(em_para_final[(length(em_para_final)-p+1):length(em_para_final)])
  mean_t = em_para_final[j] * t[j,]^gamma_hat[j] 
  Yhat_dat[[1]][,j] = mean_t
  Yhat_dat[[1]] = data.frame(Yhat_dat[[1]])
  colnames(Yhat_dat[[1]]) = paste0("PC",1:p,seq="")
}
mean.path.fit.plot(data=Yhat_dat, true_data = sim_dat[[5]], leg.pos = "none",ech = 10, ci = FALSE) +
  xlab(TeX(r'(Time(days$\times 3$))')) + ylab(TeX(r'($Y_{k}(t)$)'))
# ggsave("case/result/PMB/crack-path-fit-mean.pdf", height = 3, width = 4)

# Correlation coefficient ===========================================
# sigma_{1,2}
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
rho_hat = data.frame(cbind(1:(m+1), rho_hat))
colnames(rho_hat) = c("Time","5%","50%","95%")
rho_hat %>% pivot_longer(!Time, names_to="class", values_to = "value") %>% 
  ggplot(aes(Time,value,color = fct(class),linetype = fct(class))) + 
  geom_path() + geom_point() +
  scale_color_manual(name="",values = c("#008580","#693476","#008580"), 
                     labels=c("90% CI", "Estimated coeff.", "90% CI")) + 
  scale_linetype_manual(name="",values = c("dashed", "solid", "dashed"), labels=c("90% CI", "Estimated coeff.", "90% CI")) +
  theme_bw() + theme(panel.grid = element_blank(),legend.position = c(0.7,0.25)) +
  ylab(TeX(r'(RMSE($\times 10^{-2}))')) + xlab(TeX(r'(Time(days$\times 3$))'))
# ggsave("case/result/PMB/rho-hat2.pdf", height = 3, width = 4)


# Reliability calculation ===================================================================

# Calculate student-t and Wiener's reliability
dim(nolinear_T_R)
nolinear_T_R_dat = t(round(apply(nolinear_T_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))
nolinear_Wiener_R_dat = t(round(apply(nolinear_Wiener_R,1,quantile,c(0.05,0.5,0.95),na.rm=TRUE),5))

cdf.MW1 <- nolinear_Wiener_R_dat[,2] 
cdf.Mt1 <- nolinear_T_R_dat[,2] 
cdf.90lcl.Mt1 <- nolinear_T_R_dat[,1]  
cdf.90ucl.Mt1 <- nolinear_T_R_dat[,3]

t2 <- rt_seq; t2.t <- rt_seq 

model_pft <- function(x, a, b) {a * x + b}

pft_coef_pc =  matrix(NA,n,2)
pft_pp = matrix(NA,n,p)

for(i in 1:n){
  for(j in 1:p){
    pft_dat = data.frame("x" = log(t[j,-1]), "y" = log(sim_cum_dat[[i]][-1,j])) 
    pc1_mod = nls(y ~ model_pft(x, a, b), data = pft_dat, start = list(a = 1, b = 1)) # nonlinear model
    pft_coef = as.numeric(coef(pc1_mod)) 
    pft_pp[i,j] = exp((log(thr[j])-pft_coef[2]) / pft_coef[1]) 
  }
}
pft = apply(pft_pp,1,min)
yqq <- 1 - (1:n-0.5) / n

# Anderson–Darling test
library(ADGofTest)
cdf.MW1.PFT <- approx(t2, cdf.MW1, xout = sort(pft))$y 
ad.test(cdf.MW1.PFT)
cdf.Mt1.PFT <- approx(t2.t, cdf.Mt1, xout = sort(pft))$y
ad.test(cdf.Mt1.PFT)


pdf("case/result/PMB/R-com-pmb.pdf",width=8, height=6)
plot(t2, cdf.MW1, type = "l", col = "white", xlim = c(20,60),ylim = c(0,1), xlab = TeX(r'(Time(days$\times 3$))'), ylab = "Probability", lwd = 2, xaxt = 'n', lty = 2)
polygon(c(t2.t, rev(t2.t)), c(cdf.90lcl.Mt1, rev(cdf.90ucl.Mt1)), col = "skyblue", border = F)
axis(1, at = seq(20,60,by=10)) 
points(sort(pft), yqq, pch = 19, cex = 0.8)
lines(t2, cdf.MW1, lty = 2, lwd = 3, col = "red")
lines(t2.t, cdf.Mt1, col = "#440154", lwd = 3)
legend("topright", c(TeX(r'($M^W_{p}$)'), TeX(r'($M_p$)'), TeX(r'(90% CIs for $M_p$)'), "PFT"), 
       col = c("red","#440154", "skyblue", 1), lty = c(2, 1, 1, -1), bty = "n", inset = .02, lwd = c(3, 3, 10, -1), pch = c(-1, -1, -1, 19), cex = 1.25) # 7 , 5.5
dev.off()
# save.image(file = paste("case/result/PMB/PMB-final-data.RData", sep = ""))



