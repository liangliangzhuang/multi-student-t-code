# ======== Exploratory analysis functions for PMB data ========

# Load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(MASS, QRM, metRology, numDeriv, robustbase, ggpubr, grid)

# Boxplot
boxplot_path = function(time2 = seq(3,m,3), data = cal_dat2){
  
  newdata = newdata2 = matrix(NA,n,length(time2))
  for(i in 1:length(time2)){
    uuu = data[data$PC=="PC1" & data$Time == time2[i],]
    vvv = data[data$PC=="PC2" & data$Time == time2[i],]
    newdata[,i] = as.numeric(unlist(uuu[,4]))
    newdata2[,i] = as.numeric(unlist(vvv[,4]))
  }
  newdata = data.frame(newdata); newdata2 = data.frame(newdata2)
  colnames(newdata) <- time2; colnames(newdata2) <- time2 
  comdata = cbind(rbind(newdata,newdata2),rep(paste("PC",1:2,sep=""),each = n))
  colnames(comdata) = c(time2,"PC")
  box_plot = pivot_longer(comdata,cols = !PC,values_to = "value",names_to = "Unit") %>% 
    ggplot(aes(factor(Unit,levels = time2),value,fill=factor(Unit,levels = time2))) + 
    geom_boxplot() +
    facet_wrap(vars(PC),scale="free") +
    scale_fill_viridis(discrete = TRUE,alpha = 0.8) + 
    theme_bw() + theme(panel.grid = element_blank(),legend.position = "none") #+
    # xlab("Time") + ylab("Y(t)")
  # ggsave(paste("case/PMB/result/box_plot.pdf",sep=""), height = 5, width = 8)
  return(list(box_plot, newdata, newdata2))
}

# Contour
coutour_plot = function(sim3=sim3){
  # Binary normal fitting
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
  numML$par # Parameter estimation result (mu，sigma，pho)
  
  # Generate normally distributed data
  bivariate_data <- as.data.frame(mvrnorm(n=10000, #dim(sim3)[1]
                                          mu=numML$par[1:2],
                                          Sigma=matrix(c(numML$par[3], numML$par[5], numML$par[5], numML$par[4]), ncol=2)))
  # Plot 2D density of data
  p_real <- ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)), colour = "black", bins = 20) +
    xlab("PC1") + ylab("PC2") + theme_bw() + 
    theme(panel.grid = element_blank()) +
    scale_fill_gradientn(colours = terrain.colors(4))
  
  p_norm = ggplot(data = bivariate_data,aes(x = V1, y = V2)) + 
    stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)), colour = "black", bins = 20) +
    xlab("PC1") + ylab("PC2") + theme_bw() + 
    theme(panel.grid = element_blank()) +
    scale_fill_gradientn(colours = terrain.colors(4))
  
  # Plot the normal results to the real data
  library(viridis)
  contour_plot = ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
    stat_density_2d(contour = TRUE, color = "black", size=0.7, bins = 20, linetype = 5) + 
    geom_point(alpha=0.4,color="#34AC93") +
    xlab("PC1") + ylab("PC2") + theme_bw() + 
    theme(panel.grid = element_blank(),legend.position = "none")
  
  return(contour_plot)
}

# Q-Q plot
qqplot_PC = function(time2 = seq(3,m,3), newdata = newdata, newdata2 = newdata2){
  index2.p.v <- c(.001, .01, .05, .1, .3, .5, .7, .9, .95, .99, .999)
  index2.lower.p.v <- c(seq(.001, .01, by = .002), seq(.01, .1, by = .01), seq(.1, .9, by = .1), seq(.9, .99, by = .01), seq(.99, .999, by = .002))
  
  figure = list()
  for(i in 1:ncol(newdata)) {
    quan_plot = function(data = newdata){
      
      # Parameter estimation
      qy <- sort(data[, i])
      par_n <- fitdistr(qy, "normal")
      par_t <- fit.st(qy)
      qline <- 1000
      # Theoretical quantile
      qx_n <- qnorm((1:length(qy) - 0.5) / length(qy), mean = par_n$estimate[1], sd = par_n$estimate[2])
      qx.line <- qnorm((1:qline - 0.5) / qline, mean = par_n$estimate[1], sd = par_n$estimate[2])
      qx.t.line <- qt.scaled(as.matrix((1:qline-0.5)/qline), df = par_t$par.ests[1], mean = par_t$par.ests[2], sd = par_t$par.ests[3])
      
      qx_dat = data.frame("norm" = qx.line, "t" = qx.t.line)
      plot_data = data.frame(Theoretical = qx_n, Sampled = qy)
      
      # Plot
      p1 = ggplot() +
        geom_point(data = plot_data, aes(x = Theoretical, y = Sampled),color = "#4BA6A3",alpha = 0.8) +
        geom_line(data = qx_dat, aes(x = norm, y = t), linetype = "dashed", color = "#978ED6", size = 1.2) +
        labs(x = "Theoretical", y = "Sample") +
        scale_x_continuous(breaks = qnorm(index2.p.v, mean = par_n$estimate[1], sd = par_n$estimate[2]),
                           labels = index2.p.v) +
        geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray60") +
        geom_vline(xintercept = qnorm(0.5, mean = par_n$estimate[1], sd = par_n$estimate[2]), linetype = "dashed", color = "gray") +
        theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(size = 11))
      
      return(list(qx_dat,plot_data,par_n,"plot" = p1))
    }
    qua_dat1 = quan_plot(data = newdata)
    qua_dat2 = quan_plot(data = newdata2)
    
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
    # ggsave(paste("case/PMB/result/Time-",time2[i],".pdf",sep=""), height = 4, width = 10)
  }
  return(figure)
}



