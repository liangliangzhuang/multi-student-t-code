library(ggplot2)
library(rayshader)
library(MASS)
# 二元正态拟合
wrap <- function(parms, dat){
  mymu1  = parms[1]
  mymu2  = parms[2]
  mysig1 = parms[3]
  mysig2 = parms[4]
  myrho  = parms[5]
  myx1 <- dat[,1]
  myx2 <- dat[,2]
  n = length(myx1)
  
  f <- function(x1=myx1, x2=myx2, mu1=mymu1, mu2=mymu2, sig1=mysig1, sig2=mysig2, rho=myrho){
    -n*(log(sig1) + log(sig2) + 0.5*log(1-rho^2)) - 0.5/(1-rho^2)*(
      sum((x1-mu1)^2)/sig1^2 + sum((x2-mu2)^2)/sig2^2 - 2*rho*sum((x1-mu1)*(x2-mu2))/(sig1*sig2)
    )
  }
  f(x1=myx1, x2=myx2, mu1=mymu1, mu2=mymu2, sig1=mysig1, sig2=mysig2, rho=myrho)
  
}
eps <- eps <- .Machine$double.eps  # get a small value for bounding the paramter space to avoid things such as log(0).
numML <- optim(rep(0.5,5), wrap, dat=sim3[,1:2], 
               method="L-BFGS-B", 
               lower = c(-Inf, -Inf, eps, eps, -1+eps), 
               upper = c(Inf, Inf, 100, 100, 1-eps), 
               control = list(fnscale=-1))
numML$par # 参数估计结果（mu，sigma，pho）

# 生成正态分布的数据
bivariate_data <- as.data.frame(mvrnorm(n=10000, #dim(sim3)[1]
                                        mu=numML$par[1:2],
                                        Sigma=matrix(c(numML$par[3], numML$par[5], numML$par[5], numML$par[4]), ncol=2)))
# 绘制数据的2D密度图(真实数据)
p_real <- ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black", bins = 20) +
  xlab("PC1") + ylab("PC2") + theme_bw() + 
  theme(panel.grid = element_blank()) +
  # geom_point() +
  # scale_fill_gradient(low = "#FEFAD1", high = "red", na.value = NA)
  scale_fill_gradientn(colours = terrain.colors(4))
# p_real 
# 转化为3D图形（rayshader包）
# plot_gg(p_real, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
#         scale = 300, windowsize = c(1400, 866), zoom = 0.55, phi = 30, theta = 45)
# Sys.sleep(0.2)
# render_snapshot(clear = TRUE)

# 添加正态分布的2D密度图
p_norm = ggplot(data = bivariate_data,aes(x = V1, y = V2)) + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black", bins = 20) +
  xlab("PC1") + ylab("PC2") + theme_bw() + 
  theme(panel.grid = element_blank()) +
  # scale_fill_gradient(low = "#FEFAD1", high = "red", na.value = NA)
  scale_fill_gradientn(colours = terrain.colors(4))
# p_norm

# plot_gg(p_norm, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
#         scale = 300, windowsize = c(1400, 866), zoom = 0.55, phi = 30, theta = 30) 
# Sys.sleep(0.2)
# render_snapshot(clear = TRUE)


# 将正态结果绘制到真实数据中
library(viridis)
contour_plot = ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
  stat_density_2d(contour = TRUE, color = "black", size=0.7, bins = 20, linetype = 5) + 
  geom_point(alpha=0.4,color="#34AC93") +
  # geom_density_2d(data = bivariate_data, aes(x = V1, y = V2,linetype = "Bivariate Gaussian"),
  #                 color = "#8D82CD", bins = 20,size=0.7, alpha=1) +
  xlab("PC1") + ylab("PC2") + theme_bw() + 
  theme(panel.grid = element_blank(),legend.position = "none")
contour_plot

# contour_plot = ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
#   stat_density_2d(contour = TRUE, color = "black", size=0.7, bins = 20, aes(linetype = "Real")) + 
#   geom_point(alpha=0.4,color="#34AC93") +
#   geom_density_2d(data = bivariate_data, aes(x = V1, y = V2,linetype = "Bivariate Gaussian"),
#                   color = "#8D82CD", bins = 20,size=0.7, alpha=1) +
#   xlab("PC1") + ylab("PC2") + theme_bw() + 
#   theme(panel.grid = element_blank(),legend.position = "none") +
#   annotate(geom = "line",
#            x = c(11.6,13),size=0.7,
#            y = c(12,9), color="#8D82CD", alpha=1,
#            arrow = arrow(angle = 30, length = unit(2, "mm"))) +
#   annotate(geom = "text",
#            x = 13.5, y = 8, size=3,
#            label = "Fitted Bivariate Gaussian",color="#8D82CD", alpha=1) #+
#   # annotate(geom = "line",
#   #          x = c(8,5.5), size=0.7, color = "black", linetype =2,
#   #          y = c(20,23), alpha=0.7,
#   #          arrow = arrow(ends = "first", angle = 30, length = unit(4, "mm"))) +
#   # annotate(geom = "text",
#   #          x = 5, y = 24, size=5,
#   #          label = "Real contour",color="black", alpha=0.7)
#   # scale_color_viridis(discrete = TRUE) 
#   # scale_fill_gradientn(colours = terrain.colors(5)) +
#   # scale_color_manual(values = c("1" = 1, "2" = 2)) +
#   # scale_linetype_manual(values = c("Real" = 1, "Bivariate Gaussian" = 2)) +
#   # guides(linetype = guide_legend(title = "Contour"),
#   #        color = guide_legend(title = "Contour"))

