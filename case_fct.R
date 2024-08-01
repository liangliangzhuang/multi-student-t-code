# ======== 案例分析中使用到的函数 ========
library(statmod)

## 生成伪样本数据
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

# 绘制退化数据路径图
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


# EM 中需要的函数 ======
sumqua <- function(A, B) {
  nca <- ncol(A)
  sum0 <- t(A[, 1]) %*% B %*% A[, 1]
  for (s in 2:nca) sum0 <- sum0 + t(A[, s]) %*% B %*% A[, s]
  return(sum0)
}

sumqua0 <- function(A, B, C) { ## A、B是矩阵，C是张量
  nca <- ncol(A)
  sum0 <- t(A[, 1]) %*% B %*% C[, , 1] %*% A[, 1]
  for (s in 2:nca) sum0 <- sum0 + t(A[, s]) %*% C[, , s] %*% B %*% A[, s]
  return(sum0)
}

sumqua1 <- function(a, A) { # M步中 计算 eta 的解
  nca <- ncol(A)
  sum0 <- a[1] * A[, 1]
  for (s in 2:nca) sum0 <- sum0 + a[s] * A[, s]
  return(sum0 / sum(a))
}

sumqua2 <- function(a, A, B) { # M步中 计算 sigma 的解
  nca <- ncol(A)
  sum0 <- a[1] * A[, 1] %*% t(A[, 1]) + B
  for (s in 2:nca) sum0 <- sum0 + a[s] * A[, s] %*% t(A[, s]) + B
  return(sum0)
}

# EM算法 =======
# EM算法 —— student-t
EM_t <- function(type = "nonlinear", para, max_iter = 5000, eps = 10^-5,
                 y, y.diff, sumys) {
  
  if(type == "nonlinear"){
    # gamma 估计
    gamma_est <- function(data = y, t = t) {
      ## 时间尺度变换估计 gamma
      model <- function(x, a, b) {
        a + b * x
      }
      est_scale <- matrix(NA, p, n)
      for (j in 1:p) {
        x <- t[j, 2:(m + 1)]
        x <- log(x)
        for (i in 1:n) {
          z <- data[j, 2:(m + 1), i]
          z <- log(z)
          fit <- nls(z ~ model(x, a, b),
                     start = list(a = 1, b = 1),
                     control = nls.control(maxiter = 1000)
          )
          est_scale[j, i] <- coef(fit)[2]
        }
      }
      est_scale <- rowMeans(est_scale) # gamma
      est_t <- t^est_scale
      est_t.diff <- est_t[, 2:(m + 1)] - est_t[, 1:m] # Δλ(t)
      
      SIG <- array(0, dim = c(p, p, m))
      SIG_inv <- array(0, dim = c(p, p, m))
      for (k in 1:m) {
        SIG[, , k] <- diag(est_t.diff[, k])
        SIG_inv[, , k] <- solve(SIG[, , k])
      }
      sumSIG <- diag(est_t[, (m + 1)])
      
      return(list(
        "est_scale" = est_scale, "SIG" = SIG,
        "SIG_inv" = SIG_inv, "sumSIG" = sumSIG,
        "est_t.diff" = est_t.diff
      ))
    }
    hat_gamma <- gamma_est(data = y, t = t)
    est_scale <- hat_gamma[[1]] # gamma
    SIG <- hat_gamma[[2]]
    SIG_inv <- hat_gamma[[3]]
    sumSIG <- hat_gamma[[4]]
    est_t.diff <- hat_gamma[[5]]
    
  } else{
    # gamma 估计
    est_scale <- rep(NA,p) # gamma
    SIG <- array(0, dim = c(p, p, m))
    SIG_inv <- array(0, dim = c(p, p, m))
    for (k in 1:m) {
      SIG[, , k] <- diag(time, p)
      SIG_inv[, , k] <- solve(SIG[, , k])
    }
    sumSIG <- diag(t[, (m + 1)]) #est_t 换成t
    est_t.diff <- matrix(rep(1,n*2),nrow=2)
  }
  
  # EM算法 ===============================================
  est_etas <- para[[1]]; est_delta <- para[[2]]; est_sig0 <- para[[3]]; est_v = para[[4]]
  
  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 4 * p) # 记录迭代的结果 (eta[p], delta[p], sig0 [2p-1], v[1])
  
  while (con == FALSE) {
    iter <- iter + 1
    est_sig0_cor <- cov2cor(est_sig0)
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter, ] <- c(est_etas, est_delta, est_sig0_old,est_v)
    
    k1s <- rep(0, n)
    k2s <- rep(0, n)
    SIGdelta <- diag(est_delta^2)
    
    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(SIGdelta)
    inits5 <- inits4 %*% sumSIG
    sigtheta <- solve(inits1 + inits5) ## 已经是逆了
    inits6 <- sumSIG %*% sigtheta
    
    # E(theta)
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i] # mu_i的第二个符号
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua0(y.diff[, , i], inits4, SIG_inv)
    }
    est_theta <- etheta ## 存储迭代中每一次基于前一步参数值的后验期望
    # 公式13
    etau <- (m * p + est_v) / (k1s - k2s + est_v) 
    est_tau <- etau
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
    est_ltau <- eltau
    
    # M steps
    ## eta ======
    est_etas <- sumqua1(etau, etheta) # M步中 计算 eta 的解
    ## Sig_0 ======
    est_sig0 <- sumqua2(etau, etheta - matrix(est_etas, p, n), sigtheta) / n # M步中 计算 sigma0 的解
    ## Omega_delta ====
    for (j in 1:p) { # M步中计算 Omega_delta 的解
      dj <- 0
      for (i in 1:n) {
        dj <- dj + inits6[j, j] + etau[i] * (sum(y.diff[j, , i]^2 / est_t.diff[j, ]) + etheta[j, i]^2 * sumSIG[j, j] - 2 * sumys[j, i] * etheta[j, i])
      }
      est_delta[j] <- sqrt(dj / m / n)
    }
    ## v ========
    est_v_fun = function(x){ 
      re = - 2 * log(gamma(x/2)) + x * log(v/2) + x/n * sum(eltau - etau)
      return(-re)
    }
    est_v = optimize(est_v_fun,c(0,10))$minimum
    
    
    # 收敛停止条件
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)]) # sigma + rho
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # eta, delta, sig0, v
    if (all((abs(new - old[iter, ])) < eps) | iter > max_iter) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }
  
  Sigma_delta <- diag(diag(est_sig0)) # Sigma_delta
  Sigma_diff_t <- list()
  l_ik <- matrix(NA, n, m)
  for (i in 1:n) {
    Sigma_diff_t[[i]] <- matrix(NA, m, p)
    for (k in 1:m) {
      Sigma_diff_t[[i]][k, ] <- est_t.diff[, k] # Sigma_diff_t
      uuu <- diag(est_t.diff[, k]) # Sigma_diff_t的矩阵形式
      l_ik[i, k] <- t(y.diff[, k, i] - uuu %*% t(t(est_theta[, i]))) %*% solve(Sigma_delta * uuu) %*% (y.diff[, k, i] - uuu %*% t(t(est_theta[, i])))
    }
  }
  
  # # 计算对数似然函数
  l_c <- -((m + 1) * p) / 2 * log(2 * pi) - 1/2 * log(det(est_sig0)) - log(gamma(est_v / 2)) + est_v / 2 * log(est_v / 2)
  l_i0 <- as.numeric((est_theta[, i] - est_etas) %*% solve(est_sig0) %*% (est_theta[, i] - est_etas))
  part1 <- numeric()
  for (i in 1:n) {
    part1[i] <- l_c + (((m + 1) * p + est_v) / 2 - 1) * eltau[i] - m * sum(log(est_delta)) - 1 / 2 * sum(log(Sigma_diff_t[[i]]))
  }
  part2 <- -1/2 * sum(etau * (apply(l_ik, 1, sum) + l_i0 + est_v))
  
  logl <- sum(part1) + part2
  if(type == "nonlinear"){
    aic <- -2 * logl + 2 * (dim(old)[2]-1+p) # 非线性多个p
  } else{
    aic <- -2 * logl + 2 * (dim(old)[2]-1) # 非线性多个p
  }
  
  return(list(
    "para_iter" = old, # 迭代的所有结果
    "para" = c(new, est_scale), # 最终结果
    "est_sig0" = est_sig0, # 估计的sigma0，后续有用
    "logl" = logl, # aic
    "aic" = aic
  ))
}

# EM算法 —— Wiener 
EM_Wiener <- function(type = "nonlinear", para, max_iter = 5000, eps = 10^-5,
                      y, y.diff, sumys) {
  # 非线性的情况========
  if(type == "nonlinear"){
    # gamma 估计
    gamma_est <- function(data = y, t = t) {
      ## 时间尺度变换估计 gamma
      model <- function(x, a, b) {
        a + b * x
      }
      est_scale <- matrix(NA, p, n)
      for (j in 1:p) {
        x <- t[j, 2:(m + 1)]
        x <- log(x)
        for (i in 1:n) {
          z <- data[j, 2:(m + 1), i]
          z <- log(z)
          fit <- nls(z ~ model(x, a, b),
                     start = list(a = 1, b = 1),
                     control = nls.control(maxiter = 1000)
          )
          est_scale[j, i] <- coef(fit)[2]
        }
      }
      est_scale <- rowMeans(est_scale) # gamma
      est_t <- t^est_scale
      est_t.diff <- est_t[, 2:(m + 1)] - est_t[, 1:m] # Δλ(t)
      
      SIG <- array(0, dim = c(p, p, m))
      SIG_inv <- array(0, dim = c(p, p, m))
      for (k in 1:m) {
        SIG[, , k] <- diag(est_t.diff[, k])
        SIG_inv[, , k] <- solve(SIG[, , k])
      }
      sumSIG <- diag(est_t[, (m + 1)])
      
      return(list(
        "est_scale" = est_scale, "SIG" = SIG,
        "SIG_inv" = SIG_inv, "sumSIG" = sumSIG,
        "est_t.diff" = est_t.diff
      ))
    }
    hat_gamma <- gamma_est(data = y, t = t)
    est_scale <- hat_gamma[[1]] # gamma
    SIG <- hat_gamma[[2]]
    SIG_inv <- hat_gamma[[3]]
    sumSIG <- hat_gamma[[4]]
    est_t.diff <- hat_gamma[[5]]
    
  } else{
    # gamma 估计
    est_scale <- rep(NA,p) # gamma
    SIG <- array(0, dim = c(p, p, m))
    SIG_inv <- array(0, dim = c(p, p, m))
    for (k in 1:m) {
      SIG[, , k] <- diag(time, p)
      SIG_inv[, , k] <- solve(SIG[, , k])
      
    }
    sumSIG <- diag(t[, (m + 1)]) #est_t 换成t
    est_t.diff <- matrix(rep(1,n*2),nrow=2)
  }
  
  # EM算法 ===============================================
  est_etas <- para[[1]]; est_delta <- para[[2]]; est_sig0 <- para[[3]]; 
  est_v = 100
  
  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 4 * p) # 记录迭代的结果 (eta[p], delta[p], sig0 [2p-1], v[1])
  
  while (con == FALSE) {
    iter <- iter + 1
    est_sig0_cor <- cov2cor(est_sig0)
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter, ] <- c(est_etas, est_delta, est_sig0_old, est_v)
    
    k1s <- rep(0, n)
    k2s <- rep(0, n)
    SIGdelta <- diag(est_delta^2)
    
    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(SIGdelta)
    inits5 <- inits4 %*% sumSIG
    sigtheta <- solve(inits1 + inits5) ## 已经是逆了
    inits6 <- sumSIG %*% sigtheta
    
    # E(theta)
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i] # mu_i的第二个符号
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua0(y.diff[, , i], inits4, SIG_inv)
    }
    est_theta <- etheta ## 存储迭代中每一次基于前一步参数值的后验期望
    # 公式13
    etau <- (m * p + est_v) / (k1s - k2s + est_v) 
    est_tau <- etau
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
    est_ltau <- eltau
    
    # M steps
    ## eta ======
    est_etas <- sumqua1(etau, etheta) # M步中 计算 eta 的解
    ## Sig_0 ======
    est_sig0 <- sumqua2(etau, etheta - matrix(est_etas, p, n), sigtheta) / n # M步中 计算 sigma0 的解
    ## Omega_delta ====
    for (j in 1:p) { # M步中计算 Omega_delta 的解
      dj <- 0
      for (i in 1:n) {
        dj <- dj + inits6[j, j] + etau[i] * (sum(y.diff[j, , i]^2 / est_t.diff[j, ]) + etheta[j, i]^2 * sumSIG[j, j] - 2 * sumys[j, i] * etheta[j, i])
      }
      est_delta[j] <- sqrt(dj / m / n)
    }
    ## v ========
    # est_v_fun = function(v){ 
    #   re = - 2 * log(gamma(v/2)) + v * log(v/2) + v/n * sum(eltau - etau)
    #   return(-re)
    # }
    # est_v = optimize(est_v_fun,c(0,10))$minimum
    est_v = 100
    
    # 收敛停止条件
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)]) # sigma + rho
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # eta, delta, sig0, v
    if (all((abs(new - old[iter, ])) < eps) | iter > max_iter) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }
  
  
  Sigma_delta <- diag(diag(est_sig0)) # Sigma_delta
  Sigma_diff_t <- list()
  l_ik <- matrix(NA, n, m)
  for (i in 1:n) {
    Sigma_diff_t[[i]] <- matrix(NA, m, p)
    for (k in 1:m) {
      Sigma_diff_t[[i]][k, ] <- est_t.diff[, k] # Sigma_diff_t
      uuu <- diag(est_t.diff[, k]) # Sigma_diff_t的矩阵形式
      l_ik[i, k] <- t(y.diff[, k, i] - uuu %*% t(t(est_theta[, i]))) %*% solve(Sigma_delta * uuu) %*% (y.diff[, k, i] - uuu %*% t(t(est_theta[, i])))
    }
  }
  
  # # 计算对数似然函数
  l_c <- -((m + 1) * p) / 2 * log(2 * pi) - 1/2 * log(det(est_sig0)) - log(gamma(est_v / 2)) + est_v / 2 * log(est_v / 2)
  l_i0 <- as.numeric((est_theta[, i] - est_etas) %*% solve(est_sig0) %*% (est_theta[, i] - est_etas))
  part1 <- numeric()
  for (i in 1:n) {
    part1[i] <- l_c + (((m + 1) * p + est_v) / 2 - 1) * eltau[i] - m * sum(log(est_delta)) - 1 / 2 * sum(log(Sigma_diff_t[[i]]))
  }
  part2 <- -1/2 * sum(etau * (apply(l_ik, 1, sum) + l_i0 + est_v))
  
  logl <- sum(part1) + part2
  if(type == "nonlinear"){
    aic <- -2 * logl + 2 * (dim(old)[2]-1+p) # 非线性多个p
  } else{
    aic <- -2 * logl + 2 * (dim(old)[2]-1) # 非线性多个p
  }
  
  return(list(
    "para_iter" = old, # 迭代的所有结果
    "para" = c(new, est_scale), # 最终结果
    "est_sig0" = est_sig0, # 估计的sigma0，后续有用
    "logl" = logl, # aic
    "aic" = aic
  ))
}


# EM 迭代图
f1_names = list(expression(hat(eta)[1]), expression(hat(eta)[2]), 
                expression(hat(delta)[1]), expression(hat(delta)[2]),
                expression(hat(sigma)[1]), expression(hat(sigma)[2]), 
                expression(hat(rho)[12]), expression(hat(nu)))

EM_iter_plot = function(para_iter, f_names = f1_names){
  orders = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
             paste0("sigma",1:p, sep = ""), 
             paste0("rho",1:(p*(p-1)/2), sep = ""), "v")
  colnames(para_iter) = orders
  # 添加数学公式
  f_labeller <- function(variable, value){return(f_names[value])}
  d1 = para_iter %>% data.frame() %>% 
    mutate("index" = 1:dim(.)[1]) %>% 
    pivot_longer(cols= !index, names_to = "para", values_to = "value") 
  d1$para = factor(d1$para, levels = orders, ordered = TRUE)
  
  p1 = d1 %>% ggplot(aes(index,value)) + #,color = para
    geom_line() + 
    facet_wrap(vars(para), ncol = 4, scales = "free",labeller = f_labeller) +
    # scale_color_aaas(name = "Parameters") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
    xlab("Iteration")
  return(p1)
}


# 拟合路径图 
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
      theme(legend.position = 'none') #panel.grid = element_blank()
    
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
      theme(legend.position = 'none') #panel.grid = element_blank()
    
  }

  return(p1)
}


# 可靠度计算
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

