### 线性+相互依赖 (student-T)
EM_linear_T = function(max_iter = 4000, eps = 10^-5, y = y, y.diff= y.diff, sumys = sumys){
  
  t.diff1 <- time
  t1 <- m * time
  est_t.diff <- matrix(rep(1,m*p),nrow=p)
  
  SIGMA <- diag(t1, p)
  sig <- diag(time, p) ## λ(t)=t,这里已经是差分时间了
  sig_inv <- solve(sig) 
  
  ##EM算法估计参数
  est_etas= rep(8,p)   
  est_delta= rep(1,p)
  est_sig0= diag(p)
  est_v = 4
  # 线性情况下
  est_tau = rep(1, n); est_ltau = rep(5, n)
  # est_theta <- matrix(c(0.2971^2, 0.1544 * 0.2971 * 0.2876, 0.1544 * 0.2971 * 0.2876, 0.2876^2), 2)
  con <- FALSE
  iter = 0
  old = matrix(NA, max_iter+6, 2*p + 1 + p*(p+1)/2) #记录迭代的结果
  while (con == FALSE) {
    iter = iter + 1
    est_etas_old <- est_etas
    est_delta_old <- est_delta
    est_sig0_cor <- cov2cor(est_sig0)
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter,] <- c(est_etas_old, est_delta_old, est_sig0_old,est_v)
    
    k1s <- rep(0, n)
    k2s <- rep(0, n)
    sigw <- diag(est_delta^2)
    
    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(sigw)
    inits5 <- inits4 %*% SIGMA
    sigtheta <- solve(inits1 + inits5) ## 已经是逆了
    inits6 <- SIGMA %*% sigtheta
    
    ## expection of thetas
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i]
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua(y.diff[, , i], inits4 %*% sig_inv)
    }
    est_theta <- etheta ## 存储迭代中每一次基于前一步参数值的后验期望
    etau <- (m * p + v) / (k1s - k2s + v)
    est_tau <- etau
    eltau <- digamma((m * p + v) / 2) - log((k1s - k2s + v) / 2)
    est_ltau <- eltau
    
    #### M steps (updating)
    ## eta ======
    est_etas <- sumqua1(etau, etheta)
    ## Sig_0 ======
    est_sig0 <- sumqua2(etau, etheta - matrix(est_etas, p, n), sigtheta) / n
    ## Omega_delta ====
    for (j in 1:p) {
      wj <- 0
      for (i in 1:n) {
        wj <- wj + inits6[j, j] + etau[i] * (sum(y.diff[j, , i]^2) / t.diff1 + etheta[j, i]^2 * t1 - 2 * sumys[j, i] * etheta[j, i])
      }
      est_delta[j] <- sqrt(wj / m / n)
    }
    ## v ========
    est_v_fun <- function(x) {
      re <- -2 * log(gamma(x / 2)) + x * log(v / 2) + x / n * sum(eltau - etau)
      return(-re)
    }
    est_v <- optimize(est_v_fun, c(0, 10))$minimum
    
    ### 收敛停止条件
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)])
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # EM算法的结果
    if (all((abs(new - old[iter,])) < eps)) {
      con <- TRUE
      old[iter+1,] = new
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
  
  # 计算对数似然函数
  l_c <- -((m + 1) * p) / 2 * log(2 * pi) - 1 / 2 * log(det(est_sig0)) - log(gamma(est_v / 2)) + est_v / 2 * log(est_v / 2)
  l_i0 <- as.numeric((est_theta[, i] - est_etas) %*% solve(est_sig0) %*% (est_theta[, i] - est_etas))
  part1 <- numeric()
  for (i in 1:n) {
    part1[i] <- l_c + (((m + 1) * p + est_v) / 2 - 1) * eltau[i] - m * sum(log(est_delta)) - 1 / 2 * sum(log(Sigma_diff_t[[i]]))
  }
  part2 <- -1 / 2 * sum(etau * (apply(l_ik, 1, sum) + l_i0 + est_v))
  
  logl <- sum(part1) + part2
  aic <- -2 * logl + 2 * (dim(old)[2] - 1) # 非线性多个p
  
  return(list("para_iter" = old,  #迭代的所有结果
              "para" = c(new,rep(NA,p)), # 最终结果
              "est_sig0" = est_sig0, #估计的sigma0，后续有用
              "logl" = logl, # aic
              "aic" = aic
  ))
  
}
em_re_linear_T = EM_linear_T(max_iter = 5000, eps = 10^-5, y = y, y.diff= y.diff, sumys = sumys)

