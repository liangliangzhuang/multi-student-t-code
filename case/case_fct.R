
# Generate pseudo sample data =====
sim_path = function(par = real, v = v, SIG0 = SIG0, scen = "linear"){
  # Setting time (linear and non-linear)
  if(scen == "linear"){
    t.diff = matrix(rep(time,p*m), nrow= p) # all one
  } else{
    gamma = par[(length(par)-p+1):length(par)] # gamma
    t <-(seq(0,m,1)*time)
    t <- matrix(rep(t, p), nrow = p, byrow = TRUE)
    if(scen == "exp"){
      t.scale <- exp(t * gamma) - 1
    }else if(scen == "power"){ 
      t.scale <- t^gamma 
      }
    t.diff <- t.scale[, 2:(m+1)]-t.scale[, 1:m] 
  }
  # parameter settings
  eta = par[1:p]; delta =par[1:p + p]
  
  # Generate simulated data
  tau<-rgamma(n,shape=v/2,scale=2/v) # Generate n random tau
  theta<-matrix(0,p,n) 
  for(i in 1:n){
    theta[,i]<-mvrnorm(1,eta,SIG0/tau[i]) 
  }                                    
  y.diff<-array(0,dim=c(p,m,n))   
  for(i in 1:n){
    for(k in 1:m){
      for(j in 1:p){
        y.diff[j,k,i] <- rnorm(1,theta[j,i]*t.diff[j,k],delta[j]*sqrt(t.diff[j,k]/tau[i]))
      }
    }
  }
  # Calculate the Y Cumulative Data Set
  y <- array(0,dim=c(p,(m+1),n))
  y[,2,]<-y.diff[,1,]
  for(i in 1:n){
    for(k in 2:m){
      y[,(k+1),i]<-apply(y.diff[,1:k,i],1,sum)
    }
  }
  
  sumys = matrix(0,p,n)
  for(i in 1:n) sumys[,i]=apply(y.diff[,,i],1,sum)
  
  # Tidy data
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
              "y" = y, # The first three are original data sets for subsequent analysis. The latter are used to facilitate drawing
              "sim_diff_dat" = sim_diff_dat, 
              "sim_cum_dat" = sim_cum_dat,
              "tidy_dat" = tidy_dat,
              "tau" = tau,
              "theta" = theta))
}

# Path plot ====
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

# EM Iteration ====
EM_iter_plot = function(para_iter, f_names = f1_names, orders=orders){
  # Adding mathematical formulas
  f_labeller <- function(variable, value){return(f_names[value])}
  d1 = para_iter %>% data.frame() %>% 
    mutate("index" = 1:dim(.)[1]) %>% 
    pivot_longer(cols= !index, names_to = "para", values_to = "value") 
  d1$para = factor(d1$para, ordered = TRUE, levels = orders) #levels = orders, 
  
  p1 = d1 %>% ggplot(aes(index,value)) + #,color = para
    geom_line() + 
    facet_wrap(vars(para), ncol = 4, scales = "free",labeller = f_labeller) +
    # scale_color_aaas(name = "Parameters") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
    xlab("Iteration")
  return(p1)
}

# Fitted Path Plot ====
mean.path.fit.plot = function(data = Yhat_dat, true_data = sim_dat[[4]], leg.pos = "none",ech = 5,ci = FALSE){
  
  true_data = map(true_data, as.data.frame)
  data1 = map(true_data, ~ mutate(.x, Time = 0:(n()-1)))
  data_mean = map(data, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df1 <- bind_rows(data1, .id = "Unit")
  merged_df2 <- bind_rows(data_mean, .id = "Unit") # Conflicts with the above Unit
  if(ci == TRUE){
    merged_df2$Unit = rep(c("Low","Mean","Up"),each = length(0:m))
    # Merge two datasets
    merged_df = rbind(merged_df1,merged_df2)
    # PLot
    mer_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") 
    p1 = mer_dat %>% ggplot(aes(Time,Value,color = factor(Unit), linetype= factor(Unit), size= factor(Unit))) + 
      geom_line(alpha=0.8) + 
      facet_wrap(vars(PC),nrow = 1) + 
      theme_bw() +
      scale_color_manual(name= "", values = c(rep("gray60",n),"#440154"))+ 
      scale_linetype_manual(values = c(rep(1,n),5)) +
      scale_size_manual(values = c(rep(0.5,n),2)*0.8) +
      # scale_color_viridis(discrete = T) + 
      scale_x_continuous(breaks = seq(0, m, by = ech), limits = c(0, m)) +
      # scale_y_reverse() +
      theme(legend.position = 'none') +#panel.grid = element_blank()
      xlab(TeX(r'(Time (hours $\times$ 336))')) + ylab(TeX(r'(Y(t) (inches))'))
    
  }else{
    merged_df2$Unit = rep("Mean",each = length(0:m))
    merged_df = rbind(merged_df1,merged_df2)

    mer_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") 
    p1 = mer_dat %>% ggplot(aes(Time,Value,color = factor(Unit), linetype= factor(Unit), size= factor(Unit))) + 
      geom_line(alpha=0.8) + 
      facet_wrap(vars(PC),nrow = 1) + 
      theme_bw() +
      scale_color_manual(name= "", values = c(rep("gray60",n),"#21908C","#440154","#21908C"))+ 
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

# Reliability calculation =====
library(statmod)
R_cal = function(r_t = rt_seq, para3, r_SIG0 = bt_re[[3]], B = 1000, yz, scen = "power"){
  # Parameter value setting
  r_eta <- para3[1:p]
  r_delta <- para3[p+1:p]
  r_v = para3[length(para3)-p]
  r_gamma = para3[(length(para3)-p+1):length(para3)]

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
      lambda_t[j] = rinvgauss(1, mean = mean1[j], shape = shape1[j]) 
    }
    if(scen=="exp"){
      fail_time = log(lambda_t+1)/r_gamma
    } else if(scen == "linear"){
      fail_time = lambda_t 
    }else if(scen == "power"){
      fail_time = exp(log(lambda_t)/r_gamma)
    }
    
    ft_star[b] = min(fail_time)
  }
  ft_star = ft_star[!is.na(ft_star)]

  # Calculation reliability
  R_rate = numeric()
  for (i in 1:length(r_t)) {
    R_rate[i] <- length(which(ft_star >= r_t[i])) / length(ft_star)
  }
  r_dat <- data.frame("Time" = r_t, "value" = R_rate)
  # PLot reliability
  r_p1 = ggplot(r_dat, aes(Time, value)) +
    geom_line() +
    scale_x_continuous(limits = c(0, 60)) +
    theme_bw() +
    ylab("Reliability")

  return(list(r_dat,r_p1))
}

# Data preprocessing of Fatigue Crack Size Data 
crack_path = function(led_cum, led_diff, p=p, m=m, n=n){
  # Data preprocessing
  y.diff <-array(0,dim=c(p,m,n))   
  for(i in 1:n){
    for(k in 1:(m)){
      for(j in 1:p){
        y.diff[j,k,i] <- led_diff[[i]][k,j]
      }
    }
    
  }
  # Y Cumulative value
  y <- array(0,dim=c(p,m+1,n))
  y[,1,] <- 0
  for(i in 1:n){
    for(k in 2:(m+1)){
      y[,k,i]<-as.numeric(led_cum[[i]][k,])
    }
  }
  
  sumys = matrix(0,p,n)
  for(i in 1:n) sumys[,i]=apply(y.diff[,,i],1,sum)
  
  # Tidy data
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
              "y" = y, 
              "sim_diff_dat" = sim_diff_dat, 
              "sim_cum_dat" = sim_cum_dat,
              "tidy_dat" = tidy_dat))
}

# Functions required in EM =====
sumqua <- function(A, B) {
  nca <- ncol(A)
  sum0 <- t(A[, 1]) %*% B %*% A[, 1]
  for (s in 2:nca) sum0 <- sum0 + t(A[, s]) %*% B %*% A[, s]
  return(sum0)
}

sumqua0 <- function(A, B, C) { 
  nca <- ncol(A)
  sum0 <- t(A[, 1]) %*% B %*% C[, , 1] %*% A[, 1]
  for (s in 2:nca) sum0 <- sum0 + t(A[, s]) %*% C[, , s] %*% B %*% A[, s]
  return(sum0)
}

sumqua1 <- function(a, A) { # Calculate the solution of eta in the M-step
  nca <- ncol(A)
  sum0 <- a[1] * A[, 1]
  for (s in 2:nca) sum0 <- sum0 + a[s] * A[, s]
  return(sum0 / sum(a))
}

sumqua2 <- function(a, A, B) { # Calculate the solution of sigma in the M-step
  nca <- ncol(A)
  sum0 <- a[1] * A[, 1] %*% t(A[, 1]) + B
  for (s in 2:nca) sum0 <- sum0 + a[s] * A[, s] %*% t(A[, s]) + B
  return(sum0)
}



