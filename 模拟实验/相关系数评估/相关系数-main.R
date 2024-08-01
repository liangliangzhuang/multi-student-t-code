# ============ 维纳和t过程模型下相关系数性能计算 ============
# 构建矩阵存储rho的rmse结果
library(ggpubr)
library(ggplot2)
library(grid)

real_rho_total = real_rho_sum = rmse_I_total = rmse_II_total = rmse_III_total = rmse_IV_total = list()
itt = 100; m <- 20;  
# I ====
p <- 2
n=10
source("模拟实验/可靠度评估/可靠度比较-I.R")
rmse_I_total[[1]] = rmse_mean

n=20
source("模拟实验/可靠度评估/可靠度比较-I.R")
rmse_I_total[[2]] = rmse_mean

n=30
source("模拟实验/可靠度评估/可靠度比较-I.R")
rmse_I_total[[3]] = rmse_mean

real_rho_total[[1]] = real_rho #记录相关系数随时间变化
real_rho_sum[[1]] = SIG0cor[1,2] #记录rho_{ij}

# II ====
p <- 3

n=10
source("模拟实验/可靠度评估/可靠度比较-II.R")
rmse_II_total[[1]] = rmse_mean

n=20
source("模拟实验/可靠度评估/可靠度比较-II.R")
rmse_II_total[[2]] = rmse_mean

n=30
source("模拟实验/可靠度评估/可靠度比较-II.R")
rmse_II_total[[3]] = rmse_mean

real_rho_total[[2]] = rbind(real_rho12,real_rho13,real_rho23)
real_rho_sum[[2]] = as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)])
# III ====
p <- 2

n=10
source("模拟实验/可靠度评估/可靠度比较-III.R")
rmse_III_total[[1]] = rmse_mean

n=20
source("模拟实验/可靠度评估/可靠度比较-III.R")
rmse_III_total[[2]] = rmse_mean

n=30
source("模拟实验/可靠度评估/可靠度比较-III.R")
rmse_III_total[[3]] = rmse_mean
real_rho_total[[3]] = real_rho
real_rho_sum[[3]] = as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)])
# IV ====
p <- 3

n=10
source("模拟实验/可靠度评估/可靠度比较-IV.R")
rmse_IV_total[[1]] = rmse_mean

n=20
source("模拟实验/可靠度评估/可靠度比较-IV.R")
rmse_IV_total[[2]] = rmse_mean

n=30
source("模拟实验/可靠度评估/可靠度比较-IV.R")
rmse_IV_total[[3]] = rmse_mean

real_rho_total[[4]] = rbind(real_rho12,real_rho13,real_rho23)
real_rho_sum[[4]] = as.vector(SIG0cor[upper.tri(SIG0cor,diag=F)])
rmse_IV_total


# 画图 =======
rho_plot = function(data = rmse_I_total, p=2, pos = "bottom", sizes = 1, colnam = c("n=10", "n=20", "n=30")){
  # 相关系数绘制主函数，和p相关，如果p=2说明只有一个rho，p=3时有三个，需要用分面
  # 将list转化为数据框
  if(p==2){
    df <- do.call(cbind, lapply(data, function(x) {
      length_diff <- max(sapply(data, length)) - length(x)
      c(x, rep(NA, length_diff))
    }))
    colnames(df) <- colnam
    df <- as.data.frame(df)
    df <- cbind(Time = 1:m, df)
    # 将数据框转换为长格式
    df_long <- df %>%
      pivot_longer(cols = -Time, names_to = "Variable", values_to = "Value")
    rho_p = ggplot(df_long, aes(x = Time, y = Value, color = fct(Variable,levels = colnam), linetype =  fct(Variable,levels = colnam))) +
      geom_line() + geom_point(size=sizes) +
      scale_color_manual(values = c("#008580", "#693476", "#FFD700"), name = "") +
      scale_linetype_manual(values = c("dashed", "solid", "dotted"), name = "") +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = pos) +
      ylab(TeX(r'(RMSE($\times 10^{-2}))')) +
      xlab("Time")
  }else{
    df <- data.frame()
    for(q in 1:3){
      df2 = do.call(cbind, lapply(data, function(mat) mat[q,]))
      df <- rbind(df, df2)
    }
    colnames(df) <- colnam
    df <- as.data.frame(df)
    df <- cbind(Time = rep(1:m,3),class = rep(1:3,each=20),df)
    # 将数据框转换为长格式
    df_long <- df %>%
      pivot_longer(cols = -c(Time,class), names_to = "Variable", values_to = "Value")
    # 分面名称替换
    f_names2 <- list(unname(TeX(c("$\\rho_{12}(t)$", 
                                  "$\\rho_{13}(t)$",
                                  "$\\rho_{23}(t)$"))))
    f_labeller2 <- function(variable, value){return(f_names2[value])}
    rho_p = ggplot(df_long, aes(x = Time, y = Value, color = fct(Variable,levels = colnam), linetype =  fct(Variable,levels = colnam))) +
      geom_line() + geom_point(size=sizes) +
      scale_color_manual(values = c("#008580", "#693476", "#FFD700"), name = "") +
      scale_linetype_manual(values = c("dashed", "solid", "dotted"), name = "") +
      facet_wrap(vars(class),labeller = f_labeller2) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = pos) +
      ylab(TeX(r'(RMSE($\times 10^{-2}))')) +
      xlab("Time")
  }

  return(rho_p)
}


library(ggpubr)
library(grid)

# 四种情况的真实rho可视化 =====
# 情景I，III 合并 RMSE ======
q1 = rho_plot(data = rmse_I_total) + scale_y_continuous(limits = c(6,30))
q3 = rho_plot(data = rmse_III_total) + scale_y_continuous(limits = c(5,30))
pp2 = ggarrange(q1 + rremove("ylab") + rremove("xlab"), 
                q3 + rremove("ylab") + rremove("xlab"), 
                labels = c("I","III"), label.y = 1,label.x = 0,
                common.legend = TRUE, legend="top", nrow=1, align ="hv")

annotate_figure(pp2, left = textGrob(TeX(r'(RMSE($\times 10^{-2}))'), rot = 90, vjust = 0.5, gp = gpar(cex = 0.8)),
                bottom = textGrob("Time", gp = gpar(cex = 0.8)))
# ggsave(paste("模拟实验/相关系数评估/rho_I_III.pdf",sep=''), width =6, height = 3)

# 情景II，VI 合并 RMSE ======
# 数据调整，可能是跑太少了
uuu = rmse_II_total[[2]][2,]
rmse_II_total[[2]][2,] = rmse_II_total[[3]][3,]
rmse_II_total[[3]][3,] = uuu
q2 = rho_plot(data = rmse_II_total, p=3)
q4 = rho_plot(data = rmse_IV_total, p=3)
pp4 = ggarrange(q2 + rremove("ylab") + rremove("xlab"), 
                q4 + rremove("ylab") + rremove("xlab"), 
                labels = c("II","IV"), label.y = 1,label.x = 0,
                common.legend = TRUE, legend="top", nrow=1, align ="hv")

annotate_figure(pp4, left = textGrob(TeX(r'(RMSE($\times 10^{-2}))'), rot = 90, vjust = 0.5, gp = gpar(cex = 0.8)),
                bottom = textGrob("Time", gp = gpar(cex = 0.8)))
# ggsave(paste("模拟实验/相关系数评估/rho_II_IV.pdf",sep=''), width = 6, height = 3)


# 4个情景合并 ======

pp_all = ggarrange(q1 + rremove("ylab") + rremove("xlab"), 
                q2 + rremove("ylab") + rremove("xlab"), 
                q3 + rremove("ylab") + rremove("xlab"), 
                q4 + rremove("ylab") + rremove("xlab"), 
                labels = c("I","II", "III","IV"), common.legend = TRUE, 
                nrow = 2, ncol = 2, widths = c(1, 2.8))
pp_all 
annotate_figure(pp_all, left = textGrob(TeX(r'(RMSE($\times 10^{-2}))'), rot = 90, vjust = 0.5, gp = gpar(cex = 0.8)),
                bottom = textGrob("Time", gp = gpar(cex = 0.8)))
# ggsave(paste("模拟实验/相关系数评估/rho_all.pdf",sep=''), width = 9, height = 6)

mean(abs(unlist(real_rho_sum)))
# save.image(file = paste("模拟实验/相关系数评估/相关系数评估-0617.RData", sep = ""))

