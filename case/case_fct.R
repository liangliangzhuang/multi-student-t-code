
# 数据预处理，将原始数据转化为后续常用的数据形式 =====
crack_path = function(led_cum, led_diff, p=p, m=m, n=n){
  # led_diff 为增量
  # Y增量
  y.diff <-array(0,dim=c(p,m,n))   
  for(i in 1:n){
    for(k in 1:(m)){
      for(j in 1:p){
        y.diff[j,k,i] <- led_diff[[i]][k,j]
      }
    }
    
  }
  # Y累积量
  y <- array(0,dim=c(p,m+1,n))
  y[,1,] <- 0
  for(i in 1:n){
    for(k in 2:(m+1)){
      y[,k,i]<-as.numeric(led_cum[[i]][k,])
    }
  }
  
  sumys = matrix(0,p,n)
  for(i in 1:n) sumys[,i]=apply(y.diff[,,i],1,sum)
  
  # 整洁的数据
  tidy_dat = melt(y.diff)
  colnames(tidy_dat) <- c("p", "m", "n", "y")
  tidy_dat$p <- factor(tidy_dat$p)
  
  
  sim_cum_dat = sim_diff_dat = list()
  for(i in 1:n){
    sim_cum_dat[[i]] = data.frame(led_cum[[i]])
    sim_diff_dat[[i]] = data.frame(led_diff[[i]])
    colnames(sim_cum_dat[[i]]) = colnames(sim_diff_dat[[i]]) = c("PC1","PC2","PC3")
  }
  
  return(list("y.diff" = y.diff, 
              "sumys" = sumys, 
              "y" = y, #前三个为原始数据集，用于后续分析。后面用于方便绘图
              "sim_diff_dat" = sim_diff_dat, 
              "sim_cum_dat" = sim_cum_dat,
              "tidy_dat" = tidy_dat))
}

# 生成伪样本数据 =====
sim_path = function(par = real, v = v, SIG0 = SIG0, scen = "linear"){
  # 设置时间（线性和非线性）
  if(scen == "linear"){
    t.diff = matrix(rep(time,p*m), nrow= p) # 全是1
  } else{
    gamma = par[(length(par)-p+1):length(par)] #gamma
    t <-(seq(0,m,1)*time)
    t <- matrix(rep(t, p), nrow = p, byrow = TRUE)
    if(scen == "exp"){
      t.scale <- exp(t * gamma) - 1
    }else{ t.scale <- t^gamma }
    t.diff <- t.scale[, 2:(m+1)]-t.scale[, 1:m] 
  }
  # 参数设置
  eta = par[1:p]; delta =par[1:p + p]
  
  # 产生模拟数据
  tau<-rgamma(n,shape=v/2,scale=2/v)     #生成n个随机tau
  theta<-matrix(0,p,n) 
  for(i in 1:n){
    theta[,i]<-mvrnorm(1,eta,SIG0/tau[i]) #eta
  }                                     #生成n个对随机向量Θ
  y.diff<-array(0,dim=c(p,m,n))   
  for(i in 1:n){
    for(k in 1:m){
      for(j in 1:p){
        y.diff[j,k,i] <- rnorm(1,theta[j,i]*t.diff[j,k],delta[j]*sqrt(t.diff[j,k]/tau[i]))
      }
    }
  }
  # 计算Y累积量数据集
  y <- array(0,dim=c(p,(m+1),n))
  y[,2,]<-y.diff[,1,]
  for(i in 1:n){
    for(k in 2:m){
      y[,(k+1),i]<-apply(y.diff[,1:k,i],1,sum)
    }
  }
  
  sumys = matrix(0,p,n)
  for(i in 1:n) sumys[,i]=apply(y.diff[,,i],1,sum)
  
  # 整洁的数据
  tidy_dat = melt(y.diff)
  colnames(tidy_dat) <- c("p", "m", "n", "y")
  tidy_dat$p <- factor(tidy_dat$p)
  
  
  sim_cum_dat = sim_diff_dat = list()
  for(i in 1:n){
    sim_cum_dat[[i]] = matrix(NA,m+1,p) 
    sim_diff_dat[[i]] = matrix(NA,m,p)
    for(h in 1:p){
      re = tidy_dat[tidy_dat$n==i & tidy_dat$p==h,]
      sim_cum_dat[[i]][,h] = c(0,cumsum(re[,4]))
      sim_diff_dat[[i]][,h] = re[,4]
      sim_cum_dat[[i]] = data.frame(sim_cum_dat[[i]])
      colnames(sim_cum_dat[[i]]) = colnames(sim_diff_dat[[i]]) = paste0("PC",1:p,sep="")
    }
  }
  
  return(list("y.diff" = y.diff, 
              "sumys" = sumys, 
              "y" = y, #前三个为原始数据集，用于后续分析。后面用于方便绘图
              "sim_diff_dat" = sim_diff_dat, 
              "sim_cum_dat" = sim_cum_dat,
              "tidy_dat" = tidy_dat,
              "tau" = tau,
              "theta" = theta))
}


# 绘制退化数据路径图 ====
degradation.path.plot = function(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free"){
  data1 = map(data, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df <- bind_rows(data1, .id = "Unit")
  cal_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
  p1 = cal_dat %>% 
    ggplot(aes(Time,Value,color = factor(Unit))) + 
    geom_line(alpha=0.8) + geom_point(size=0.8) +
    facet_wrap(vars(factor(PC)),nrow = 1, scales = scale1) + 
    theme_bw() +
    # scale_color_viridis(discrete = T) + 
    scale_x_continuous(breaks = seq(0, m, by = ech),limits = c(0, m)) + 
    theme(legend.position = leg.pos) #panel.grid = element_blank()
  return(p1)
}


# EM迭代图 ====
f1_names = list(expression(hat(eta)[1]), expression(hat(eta)[2]), 
                expression(hat(delta)[1]), expression(hat(delta)[2]),
                expression(hat(sigma)[1]), expression(hat(sigma)[2]), 
                expression(hat(rho)[12]), expression(hat(nu)))


EM_iter_plot = function(para_iter, f_names = f1_names){
  # 添加数学公式
  f_labeller <- function(variable, value){return(f_names[value])}
  d1 = para_iter %>% data.frame() %>% 
    mutate("index" = 1:dim(.)[1]) %>% 
    pivot_longer(cols= !index, names_to = "para", values_to = "value") 
  d1$para = factor(d1$para, ordered = TRUE) #levels = orders, 
  
  p1 = d1 %>% ggplot(aes(index,value)) + #,color = para
    geom_line() + 
    facet_wrap(vars(para), ncol = 4, scales = "free",labeller = f_labeller) +
    # scale_color_aaas(name = "Parameters") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
    xlab("Iteration")
  return(p1)
}


# 拟合路径图 ====
mean.path.fit.plot = function(data = Yhat_dat, true_data = sim_dat[[4]], leg.pos = "none",ech = 5,ci = FALSE){
  
  true_data = map(true_data, as.data.frame)
  data1 = map(true_data, ~ mutate(.x, Time = 0:(n()-1)))
  data_mean = map(data, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df1 <- bind_rows(data1, .id = "Unit")
  merged_df2 <- bind_rows(data_mean, .id = "Unit") # 和上面Unit相冲突
  if(ci == TRUE){
    merged_df2$Unit = rep(c("Low","Mean","Up"),each = length(0:m))
    # 两个数据集合并
    merged_df = rbind(merged_df1,merged_df2)
    # 绘图
    mer_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") 
    p1 = mer_dat %>% ggplot(aes(Time,Value,color = factor(Unit), linetype= factor(Unit), size= factor(Unit))) + 
      geom_line(alpha=0.8) + 
      # geom_ribbon(aes(ymin = Value, ymax = Value), fill = "grey70") + 
      #geom_point(size=0.8) +
      facet_wrap(vars(PC),nrow = 1) + 
      theme_bw() +
      # scale_color_manual(name= "", values = c(rep("gray60",n),"#21908C","#440154","#21908C"))+
      scale_color_manual(name= "", values = c(rep("gray60",n),"#440154"))+ #"#21908C","#440154","#21908C""blue","red","blue"
      scale_linetype_manual(values = c(rep(1,n),5)) +
      scale_size_manual(values = c(rep(0.5,n),2)*0.8) +
      # scale_color_viridis(discrete = T) + 
      scale_x_continuous(breaks = seq(0, m, by = ech), limits = c(0, m)) +
      # scale_y_reverse() +
      theme(legend.position = 'none') +#panel.grid = element_blank()
      xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches))'))
    
  }else{
    merged_df2$Unit = rep("Mean",each = length(0:m))
    # 两个数据集合并
    merged_df = rbind(merged_df1,merged_df2)
    # 绘图
    mer_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") 
    p1 = mer_dat %>% ggplot(aes(Time,Value,color = factor(Unit), linetype= factor(Unit), size= factor(Unit))) + 
      geom_line(alpha=0.8) + 
      # geom_ribbon(aes(ymin = Value, ymax = Value), fill = "grey70") + 
      #geom_point(size=0.8) +
      facet_wrap(vars(PC),nrow = 1) + 
      theme_bw() +
      # scale_color_manual(name= "", values = c(rep("gray60",n),"#21908C","#440154","#21908C"))+
      scale_color_manual(name= "", values = c(rep("gray60",n),"#21908C","#440154","#21908C"))+ #"#21908C","#440154","#21908C""blue","red","blue"
      scale_linetype_manual(values = c(rep(1,n),5,1,2)) +
      scale_size_manual(values = c(rep(0.5,n),2,1.1,1)*0.8) +
      # scale_color_viridis(discrete = T) + 
      scale_x_continuous(breaks = seq(0, m, by = ech), limits = c(0, m)) +
      # scale_y_reverse() +
      theme(legend.position = 'none') +#panel.grid = element_blank() 
      xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches))'))
  }

  return(p1)
}


# 可靠度计算 =====
library(statmod)
R_cal = function(r_t = rt_seq, para3, r_SIG0 = bt_re[[3]], B = 1000, yz, scen = "non-linear"){
  # 可靠度计算
  # 参数值设置
  r_eta <- para3[1:p]
  r_delta <- para3[p+1:p]
  r_v = para3[length(para3)-p]
  r_gamma = para3[(length(para3)-p+1):length(para3)]
  # 根据循环，寻找值
  ft_star = numeric()
  for (b in 1:B) {
    #step1
    r_tau <- rgamma(1, shape = r_v / 2, scale = 2 / r_v)
    #step2
    r_theta <- mvrnorm(1, r_eta, r_SIG0 / r_tau)
    #step3
    mean1 = thr/r_theta; shape1 = thr^2*sqrt(r_tau)/r_delta
    lambda_t = numeric()
    for(j in 1:p){
      lambda_t[j] = rinvgauss(1, mean = mean1[j], shape = shape1[j]) #计算出t函数的随机数
    }
    if(scen=="exp"){
      fail_time = log(lambda_t+1)/r_gamma
    } else if(scen == "linear"){
      fail_time = lambda_t #+1  # 转化到t
    }else{
      fail_time = exp(log(lambda_t)/r_gamma) #+1  # 转化到t
    }
    
    ft_star[b] = min(fail_time)
  }
  ft_star = ft_star[!is.na(ft_star)]

  # 计算可靠度
  R_rate = numeric()
  for (i in 1:length(r_t)) {
    R_rate[i] <- length(which(ft_star >= r_t[i])) / length(ft_star)
  }
  r_dat <- data.frame("Time" = r_t, "value" = R_rate)
  # 绘制可靠度
  r_p1 = ggplot(r_dat, aes(Time, value)) +
    geom_line() +
    scale_x_continuous(limits = c(0, 60)) +
    theme_bw() +
    ylab("Reliability")

  return(list(r_dat,r_p1))
}



