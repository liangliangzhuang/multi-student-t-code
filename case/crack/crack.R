# ======== Fatigue Crack Size Data Analysis ========
# The "Fatigue-crack-size.xlsx" dataset, originally introduced in Meeker et al. (2022) and 
# further processed as described in Appendix H of Fang et al. (2022).
# Load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, viridis, reshape2, latex2exp, MASS, openxlsx)
source("case/case_fct.R")
source("Model_est.r")
options(scipen = 4) # Set decimal point precision

# 1. Data Loading =======
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
time <- 1 # Observation interval
t <- (seq(0,m,1)*time)
t <- matrix(rep(t, p), nrow = p, byrow = TRUE)

# Data processing
crack_dat = crack_path(led_cum = fc_dat_sub, led_diff = fc_diff_dat, p=p, m=m, n=n)
y.diff = crack_dat[[1]]; 
sumys = crack_dat[[2]]; 
y =  crack_dat[[3]] 
sim_diff_dat = crack_dat[[4]] # Easy to draw (difference)
sim_cum_dat = crack_dat[[5]] # Easy to draw (cumulative)

# 2. Exploratory Analysis ========
# 2.1 Plot degradation data paths ====
degradation.path.plot(data = sim_cum_dat, leg.pos = "bottom",ech = 5, scale1 = "fixed") +
  scale_color_viridis(discrete = TRUE, name = "Unit",guide = guide_legend(nrow = 1)) + 
  xlab(TeX(r'(Time (hours $\times$ 336))'))  + 
  ylab(TeX(r'(Y(t) (inches))') ) 
# ggsave("case/result/crack/path-dat.pdf", height = 5, width = 9)

# 2.2 Boxplot ====
data2 = map(sim_cum_dat, ~ mutate(.x, Time = 0:(n()-1)))
merged_df2 <- bind_rows(data2, .id = "Unit")
cal_dat2 = merged_df2 %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
time2 = seq(3,m,3)
newdata = newdata2 = newdata3 = matrix(NA,n,length(time2))

for(i in 1:length(time2)){
  uuu = cal_dat2[cal_dat2$PC=="PC1" & cal_dat2$Time == time2[i],]
  vvv = cal_dat2[cal_dat2$PC=="PC2" & cal_dat2$Time == time2[i],]
  zzz = cal_dat2[cal_dat2$PC=="PC3" & cal_dat2$Time == time2[i],]
  newdata[,i] = as.numeric(unlist(uuu[,4]))
  newdata2[,i] = as.numeric(unlist(vvv[,4]))
  newdata3[,i] = as.numeric(unlist(zzz[,4]))
}
newdata = data.frame(newdata); newdata2 = data.frame(newdata2); newdata3 = data.frame(newdata3)
colnames(newdata) <- time2; colnames(newdata2) <- time2; colnames(newdata3) <- time2
comdata = cbind(rbind(newdata,newdata2,newdata3),rep(paste("PC",1:p,sep=""),each = n))
colnames(comdata) = c(time2,"PC")
tidyr::pivot_longer(comdata,cols = !PC,values_to = "value",names_to = "Unit") %>% 
  ggplot(aes(factor(Unit,levels = time2),value,fill=factor(Unit,levels = time2))) + 
  geom_boxplot() +
  facet_wrap(vars(PC),scale="free") +
  scale_fill_viridis(discrete = TRUE,alpha = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") +
  xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches))'))

# ggsave(paste("case/result/crack/boxplot.pdf",sep=""), height = 5, width = 8)


# 3. Parameter Estimation for Different Models =============================================
types = "exp" # exp form
# M_e
para <- list("est_etas" = c(0.01,0.01,0.01)*10, "est_delta" = rep(1, p)*0.2, "est_sig0" = diag(p), "est_v" = 1)
em_re_nolinear_T <- EM_nonlinear_T(para = para, type="exp",max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
em_re_nolinear_T_no_na <- em_re_nolinear_T$para_iter[complete.cases(em_re_nolinear_T$para_iter), ]
em_re_nolinear_T_no_na = em_re_nolinear_T_no_na[,-c(10:12)]
colnames(em_re_nolinear_T_no_na) = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
                                     paste0("sigma",1:p, sep = ""),  "v")
# Plot
f1_names <- list('eta1' = TeX(c("$\\hat{\\eta}_{1}$")), 'eta2' = TeX(c("$\\hat{\\eta}_{2}$")),'eta3' = TeX(c("$\\hat{\\eta}_{3}$")),
                 'delta1' = TeX(c("$\\hat{\\delta}^2_{1}$")), 'delta2' = TeX(c("$\\hat{\\delta}^2_{2}$")), 'delta3' = TeX(c("$\\hat{\\delta}^2_{3}$")),
                 'sigma1' = TeX(c("$\\hat{\\sigma}^2_{1}$")), 'sigma2' = TeX(c("$\\hat{\\sigma}^2_{2}$")),'sigma3' = TeX(c("$\\hat{\\sigma}^2_{3}$")),
                 'v' = TeX(c("$\\hat{\\nu}$")))
orders = unique(colnames(em_re_nolinear_T_no_na))
EM_iter_plot(para_iter = em_re_nolinear_T_no_na, f_names = f1_names, orders=orders) # Plot EM Algorithm Iterations
# ggsave("case/result/crack/FCS-EM-iter-exp_T.pdf", height = 5, width = 9)


# M^W_e
para <- list("est_etas" = rep(0.2, p), "est_delta" = rep(0.01, p), "est_sig0" = diag(p)*0.1, "est_v" = 4)
em_re_nolinear_Wiener <- EM_nonlinear_Wiener(para = para, type="exp", max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
em_re_nolinear_Wiener_no_na <- em_re_nolinear_Wiener$para_iter[complete.cases(em_re_nolinear_Wiener$para_iter), ]
em_re_nolinear_Wiener_no_na = em_re_nolinear_Wiener_no_na[,-c(10:13)] # 删除sigma_{12}等结果和v
colnames(em_re_nolinear_Wiener_no_na) = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
                                          paste0("sigma",1:p, sep = ""))
# Plot
f1_names <- list('eta1' = TeX(c("$\\hat{\\eta}_{1}$")), 'eta2' = TeX(c("$\\hat{\\eta}_{2}$")),'eta3' = TeX(c("$\\hat{\\eta}_{3}$")),
                 'delta1' = TeX(c("$\\hat{\\delta}^2_{1}$")), 'delta2' = TeX(c("$\\hat{\\delta}^2_{2}$")), 'delta3' = TeX(c("$\\hat{\\delta}^2_{3}$")),
                 'sigma1' = TeX(c("$\\hat{\\sigma}^2_{1}$")), 'sigma2' = TeX(c("$\\hat{\\sigma}^2_{2}$")),'sigma3' = TeX(c("$\\hat{\\sigma}^2_{3}$"))
)

EM_iter_plot(para_iter = em_re_nolinear_Wiener_no_na,f_names = f1_names,orders=orders) # 绘制EM迭代图
# ggsave("case/result/crack/FCS-EM-iter-exp_Wiener.pdf", height = 5, width = 9)

# AIC
c(em_re_nolinear_T$aic, em_re_nolinear_Wiener$aic)
# Point estimates
cbind(em_re_nolinear_T$para2, em_re_nolinear_Wiener$para2)

# save.image(file = paste("case/result/crack/crack-final-data.RData", sep = ""))


