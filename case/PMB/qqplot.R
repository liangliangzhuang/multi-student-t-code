library(MASS)      # use 'fitdistr' function to fit a normal distribution
library(QRM)       # use 'fit.st' function to fit a Student-t distribution
library(metRology) # use 'qt.scaled' function to find the quantile of a Student-t distribution
library(numDeriv)
library(robustbase)
library(ggpubr)
library(grid)
# library(timeSeries)
# library(DEoptimR)
# library(timeDate)
# library(Rcpp)

cal_dat3 = read.csv("case/PMB/cal_dat2.csv")
time2 = seq(3,m,3)
newdata = newdata2 = matrix(NA,n,length(time2))

for(i in 1:length(time2)){
  uuu = cal_dat3[cal_dat2$PC=="PC1" & cal_dat3$Time == time2[i],]
  vvv = cal_dat3[cal_dat2$PC=="PC2" & cal_dat3$Time == time2[i],]
  newdata[,i] = uuu[,4] #as.numeric(unlist(uuu[,4]))
  newdata2[,i] = vvv[,4]
}
newdata = data.frame(newdata); newdata2 = data.frame(newdata2)


# 箱线图（图2）
colnames(newdata) <- time2; colnames(newdata2) <- time2 
comdata = cbind(rbind(newdata,newdata2),rep(paste("PC",1:2,sep=""),each = n))
colnames(comdata) = c(time2,"PC")
box_plot = pivot_longer(comdata,cols = !PC,values_to = "value",names_to = "Unit") %>% 
  ggplot(aes(factor(Unit,levels = time2),value,fill=factor(Unit,levels = time2))) + 
  geom_boxplot() +
  facet_wrap(vars(PC),scale="free") +
  scale_fill_viridis(discrete = TRUE,alpha = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") +
  xlab("Time") + ylab("Y(t)")
# ggsave(paste("第三章/case/result/box_plot.pdf",sep=""), height = 5, width = 8)

# QQ 图
index2.p.v <- c(.001, .01, .05, .1, .3, .5, .7, .9, .95, .99, .999)
index2.lower.p.v <- c(seq(.001, .01, by = .002), seq(.01, .1, by = .01), seq(.1, .9, by = .1), seq(.9, .99, by = .01), seq(.99, .999, by = .002))

figure = list()
for(i in 1:ncol(newdata)) {
  quan_plot = function(data = newdata){
    
    # 参数估计
    qy <- sort(data[, i])
    par_n <- fitdistr(qy, "normal")
    par_t <- fit.st(qy)
    qline <- 1000
    
    # 理论分位数
    qx_n <- qnorm((1:length(qy) - 0.5) / length(qy), mean = par_n$estimate[1], sd = par_n$estimate[2])
    qx.line <- qnorm((1:qline - 0.5) / qline, mean = par_n$estimate[1], sd = par_n$estimate[2])
    qx.t.line <- qt.scaled(as.matrix((1:qline-0.5)/qline), df = par_t$par.ests[1], mean = par_t$par.ests[2], sd = par_t$par.ests[3])
    # 创建两个数据框
    qx_dat = data.frame("norm" = qx.line, "t" = qx.t.line)
    plot_data = data.frame(Theoretical = qx_n, Sampled = qy)
  
    # 画图
    
    p1 = ggplot() +
      geom_point(data = plot_data, aes(x = Theoretical, y = Sampled),color = "#4BA6A3",alpha = 0.8) +
      geom_line(data = qx_dat, aes(x = norm, y = t), linetype = "dashed", color = "#978ED6", size = 1.2) +
      # facet_wrap(vars(PC),scale="free") + 
      labs(x = "Theoretical", y = "Sample") +
      scale_x_continuous(breaks = qnorm(index2.p.v, mean = par_n$estimate[1], sd = par_n$estimate[2]),
                         labels = index2.p.v) +
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray60") +
      geom_vline(xintercept = qnorm(0.5, mean = par_n$estimate[1], sd = par_n$estimate[2]), linetype = "dashed", color = "gray") +
      # scale_x_continuous(
      #   breaks = qnorm(index2.lower.p.v, mean = par_n$estimate[1], sd = par_n$estimate[2]),
      #   labels = NULL) +
      theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(size = 11))
    
    return(list(qx_dat,plot_data,par_n,"plot" = p1))
  }
  # 单个绘制
  qua_dat1 = quan_plot(data = newdata)
  qua_dat2 = quan_plot(data = newdata2)
  # ggsave(paste("case/result/qqplot/Time-",time2[i],"PC1.pdf",sep=""), qua_dat1$plot, height = 4, width = 5)
  ggsave(paste("case/result/qqplot/Time-",time2[i],"PC2.pdf",sep=""), qua_dat2$plot, height = 4, width = 5)
  # 合并两个图形 并共享x y 轴
  figure[[i]] <- ggarrange(qua_dat1$plot + rremove("ylab") + rremove("xlab"), 
                      qua_dat1$plot + rremove("ylab") + rremove("xlab"),
                      labels = c("PC1","PC2"),
                      hjust = -2,vjust = 5,
                      ncol = 2, nrow = 1,
                      common.legend = TRUE, legend = "bottom",
                      align = "hv", 
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

  annotate_figure(figure[[i]], left = textGrob("Sample", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                  bottom = textGrob("Theoretical", gp = gpar(cex = 1)))
  ggsave(paste("case/result/qqplot/Time-",time2[i],".pdf",sep=""), height = 4, width = 10)
}





