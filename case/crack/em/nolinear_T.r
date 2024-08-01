### 非线性+相互依赖 (v=2)
library(minpack.lm)
EM_nonlinear_T <- function(para, max_iter = 5000, eps = 10^-5,
                         y, y.diff, sumys) {
  # gamma 估计
  gamma_est <- function(data = y, t = t) {
    ## 时间尺度变换估计 gamma
    # exp 形式
    model <- function(x, a, b) {
      a * (exp(b*x)-1) # eta * (exp(at) - 1)
    }
    est_scale = est_eta1 = c()
    for (j in 1:p) {
      x <- t[j, 2:(m + 1)]  #t
      nls_dat = c()
      for (i in 1:n) {
        z <- data[j, 2:(m + 1), i] #y
        z1 = as.numeric(na.omit(z)); x1 = x[!is.na(z)]
        nls_dat = rbind(nls_dat,cbind(x1,z1))
      }
      
      fit <- nlsLM(z1 ~ model(x1, a, b),
                 start = list(a = 1, b = 1),
                 control = nls.control(maxiter = 1000))
      est_scale[j] <- coef(fit)[2]
      est_eta1[j] = coef(fit)[1]
    }
    # est_scale <- rowMeans(est_scale) # gamma
    # est_t <- t^est_scale # power
    est_t <- exp(est_scale*t)-1 # exp
    
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
      "est_t.diff" = est_t.diff,
      "est_eta" = est_eta1
    ))
  }

  hat_gamma <- gamma_est(data = y, t = t)
  est_scale <- hat_gamma[[1]] # gamma
  SIG <- hat_gamma[[2]] # 和 t 相关的Sigma 
  SIG_inv <- hat_gamma[[3]]
  sumSIG <- hat_gamma[[4]]
  est_t.diff <- hat_gamma[[5]]
  
  # EM算法 ===============================================
  est_etas <- hat_gamma[[6]] # 上面eta的估计作为初始值
  est_delta <- para[[2]]
  est_sig0 <- para[[3]] 
  est_v <- para[[4]]

  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2*p + 1 + p*(p+1)/2) # 记录迭代的结果 (eta[p], delta[p], sig0 [p*(p+1)/2], v[1])

  while (con == FALSE) {
    iter <- iter + 1
    est_sig0_cor <- cov2cor(est_sig0) # \Sigam_0
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter, ] <- c(est_etas, est_delta, est_sig0_old, est_v)

    k1s <- rep(0, n)
    k2s <- rep(0, n)
    SIGdelta <- diag(est_delta^2)

    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(SIGdelta) # \Omega_\delta
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
    est_tau <- etau  #E(\tau_i|xxx)
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
    est_ltau <- eltau #E(\ln{\tau_i}|xxx)

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
    est_v_fun <- function(x) {
      re <- -2 * log(gamma(x / 2)) + x * log(x / 2) + x / n * sum(eltau - etau)
      return(-re)
    }
    est_v <- optimize(est_v_fun, c(0, 100))$minimum


    # 收敛停止条件
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)]) # sigma + rho
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # eta, delta, sig0, v
    #输出sigma12
    est_sig0_new2 <- sqrt(c(diag(est_sig0), est_sig0[upper.tri(est_sig0, diag = F)]))
    new2 = c(est_etas, est_delta, est_sig0_new2, est_v)
    if (all((abs(new - old[iter, ])) < eps) | iter > max_iter) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }
  
  # 计算对数似然函数 ===========
  est_E_tau = etau # est_E_tau - 保存每个i的E[tau_i | ΔY_i, Ψ^(s)]
  Sigma_0_inv = solve(est_sig0) # Sigma_0_inv - Σ_0的逆矩阵
  Omega_delta_inv = solve(diag(est_delta)) # Omega_delta_inv - Ω_δ的逆矩阵
  mu_i_s = etheta
  eta = est_etas # eta - η向量
  Sigma_Theta_i_s = sigtheta # Sigma_Theta_i_s - 每个i的Σ_Θ_i^(s)矩阵
  
  # 计算\Delta Sigma(t_{i,k}) ======
  Delta_Sigma = list()
  for (i in 1:n) {
    Delta_Sigma[[i]] <- list()
    for (k in 1:m) {
      Delta_Sigma[[i]][[k]] <- diag(est_t.diff[, k]) # Sigma_diff_t的矩阵形式
    }
  }
  
  # 计算 est_E_tau_ell_0 和 est_E_tau_ell_k =======
  est_E_tau_ell_0 = numeric()
  est_E_tau_ell_k = list()  
  for (i in 1:n) {
    # E[tau_i ell_i,0 | ΔY_i, Ψ^(s)]
    trace_part = sum(diag(Sigma_0_inv %*% Sigma_Theta_i_s))
    quadratic_part = t(mu_i_s[,i] - eta) %*% Sigma_0_inv %*% (mu_i_s[,i] - eta)
    ell_i0 = as.numeric(trace_part + est_E_tau[i] * quadratic_part)
    # 存储E[tau_i ell_i,k | ΔY_i, Ψ^(s)]的向量
    ell_ik_vector <- numeric()
    for (k in 1:m) {
      Delta_Sigma_inv = solve(Delta_Sigma[[i]][[k]])
      trace_ik = sum(diag(Delta_Sigma[[i]][[k]] %*% Omega_delta_inv %*% Sigma_Theta_i_s))
      quadratic_ik = t(mu_i_s[,i] - Delta_Sigma_inv %*% y.diff[, k, i]) %*% 
        solve(diag(est_delta) %*% Delta_Sigma_inv) %*% 
        (mu_i_s[,i] - Delta_Sigma_inv %*% y.diff[, k, i])
      ell_ik_vector[k] = as.numeric(trace_ik + est_E_tau[i] * quadratic_ik)
    }
    
    # 将结果存入列表
    est_E_tau_ell_0[i] = ell_i0
    est_E_tau_ell_k[[i]] = ell_ik_vector
  }
  
  ### Part1
  l_c <- -((m + 1) * p) / 2 * log(2 * pi) - 1 / 2 * log(det(est_sig0)) - log(gamma(est_v / 2)) + est_v / 2 * log(est_v / 2)
  part1 <- numeric()
  for (i in 1:n) {
    part1[i] <- l_c + (((m + 1) * p + est_v) / 2 - 1) * eltau[i] - 
      m * sum(log(est_delta)) - 
      1 / 2 * sum(log(est_t.diff))
  }
  ### Part2
  part2 <- numeric()
  for (i in 1:n) {
    tau_ell_sum <- sum(est_E_tau_ell_k[[i]])
    part2[i] <- tau_ell_sum + est_E_tau_ell_0[i] + est_v * est_E_tau[i]  # 假设est_E_tau是一个向量
  }
  # 将part1和part2结合起来得到最终的Q函数值
  Q_value <- sum(part1 - 1/2 * part2)
  aic <- -2 * Q_value + 2 * (dim(old)[2] - 1 + p) # 非线性多个p
  
  return(list(
    "para_iter" = old, # 迭代的所有结果
    "para" = c(new, est_scale), # 最终结果
    "est_sig0" = est_sig0, # 估计的sigma0，后续有用
    "logl" = Q_value, # aic
    "aic" = aic,
    "para2" = c(new2, est_scale) # 最终结果
  ))
}

para <- list("est_etas" = c(0.01,0.01,0.01)*10, "est_delta" = rep(1, p)*0.2, "est_sig0" = diag(p), "est_v" = 1)
em_re_nolinear_T <- EM_nonlinear_T(para = para, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
round(em_re_nolinear_T$para,6)
round(em_re_nolinear_T$para2,3)
# real
# EM迭代结果 =========
em_re_nolinear_T_no_na <- em_re_nolinear_T$para_iter[complete.cases(em_re_nolinear_T$para_iter), ]
em_re_nolinear_T_no_na = em_re_nolinear_T_no_na[,-c(10:12)] # 删除sigma_{12}等结果
colnames(em_re_nolinear_T_no_na) = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
                        paste0("sigma",1:p, sep = ""),  "v")

f1_names <- list('eta1' = TeX(c("$\\hat{\\eta}_{1}$")), 'eta2' = TeX(c("$\\hat{\\eta}_{2}$")),'eta3' = TeX(c("$\\hat{\\eta}_{3}$")),
                 'delta1' = TeX(c("$\\hat{\\delta}^2_{1}$")), 'delta2' = TeX(c("$\\hat{\\delta}^2_{2}$")), 'delta3' = TeX(c("$\\hat{\\delta}^2_{3}$")),
                 'sigma1' = TeX(c("$\\hat{\\sigma}^2_{1}$")), 'sigma2' = TeX(c("$\\hat{\\sigma}^2_{2}$")),'sigma3' = TeX(c("$\\hat{\\sigma}^2_{3}$")),
                 # "gamma1" = TeX(c("$\\gamma_{1}$")), "gamma2" = TeX(c("$\\gamma_{2}$")), "gamma3" = TeX(c("$\\gamma_{3}$"))
                 'v' = TeX(c("$\\hat{\\nu}$")))
EM_iter_plot(para_iter = em_re_nolinear_T_no_na,f_names = f1_names) # 绘制EM迭代图
# ggsave("case/crack/EM-iter/FCS-non-linear_T-em.pdf", height = 5, width = 9)

# 自助法 ==================
em_para_final <- em_re_nolinear_T[[2]]
# em_para_final[3] = 0.13
est_v = em_para_final[length(em_para_final)-p]
bt_ci <- matrix(NA, item, 3*p + 1 + p*(p+1)/2) #p=2时  10
nolinear_T_R <- matrix(NA, length(rt_seq), item)

## 生成伪样本数据
for (h in 1:item) { 
  # 1. 产生伪样本数据
  bt_dat <- sim_path(par = em_re_nolinear_T[[2]], v = est_v, SIG0 = em_re_nolinear_T[[3]], scen = "exp")
  y.diff_bt <- bt_dat[[1]]
  sumys_bt <- bt_dat[[2]]
  y_bt <- bt_dat[[3]] # EM中常用数据
  # 2. EM 估计
  para_bt <- list("est_etas" = em_re_nolinear_T[[2]][1:p], "est_delta" = em_re_nolinear_T[[2]][1:p + p],
    "est_sig0" = em_re_nolinear_T[[3]], "est_v" = est_v)

  bt_re <- tryCatch(EM_nonlinear_T(para = para_bt, max_iter = 5000, eps = 10^-5, y = y_bt, y.diff = y.diff_bt, sumys = sumys_bt),
    error = function(e) {
      return(NA)
    }
  )

  if (all(is.na(bt_re))) {
    bt_ci[h, ] <- NA
    nolinear_T_R[, h] <- NA
  } else {
    bt_ci[h, ] <- bt_re$para
    nolinear_T_R[, h] <- R_cal(para3 = bt_re$para, r_t = rt_seq, r_SIG0 = bt_re[[3]], B = B_item, yz = thr, scen = "exp")[[1]][, 2]
  }
  print(h)
}

bt_ci_T <- bt_ci

BT_ALL = bt_ci_T[bt_ci_T[,13] <60 & bt_ci_T[,3]<1,]
quan_nolinear_T <- data.frame(round(apply(bt_ci_T[bt_ci_T[,13] <60 & bt_ci_T[,3]<1,], 2, quantile, c(0.1, 0.5, 0.95), na.rm = TRUE), 5))

colnames(BT_ALL) = colnames(quan_nolinear_T) = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
                              paste0("sigma",1:p, sep = ""), 
                              paste0("rho",1:(p*(p-1)/2), sep = ""), "v", 
                               paste0("gamma",1:p, sep = ""))

quan_nolinear_T


