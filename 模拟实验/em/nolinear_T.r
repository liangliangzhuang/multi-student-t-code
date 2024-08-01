### 非线性+相互依赖 (v=2)
EM_nonlinear_T <- function(para, max_iter = 5000, eps = 10^-5,
                         y, y.diff, sumys) {
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

  # EM算法 ===============================================
  est_etas <- para[[1]]
  est_delta <- para[[2]]
  est_sig0 <- para[[3]]
  est_v <- para[[4]]

  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2*p + 1 + p*(p+1)/2) # 记录迭代的结果 (eta[p], delta[p], sig0 [2p-1], v[1])

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
    est_v_fun <- function(x) {
      re <- -2 * log(gamma(x / 2)) + x * log(x / 2) + x / n * sum(eltau - etau)
      return(-re)
    }
    est_v <- optimize(est_v_fun, c(0, 10))$minimum


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
  l_c <- -((m + 1) * p) / 2 * log(2 * pi) - 1 / 2 * log(det(est_sig0)) - log(gamma(est_v / 2)) + est_v / 2 * log(est_v / 2)
  l_i0 <- as.numeric((est_theta[, i] - est_etas) %*% solve(est_sig0) %*% (est_theta[, i] - est_etas))
  part1 <- numeric()
  for (i in 1:n) {
    part1[i] <- l_c + (((m + 1) * p + est_v) / 2 - 1) * eltau[i] - m * sum(log(est_delta)) - 1 / 2 * sum(log(Sigma_diff_t[[i]]))
  }
  part2 <- -1 / 2 * sum(etau * (apply(l_ik, 1, sum) + l_i0 + est_v))

  logl <- sum(part1) + part2
  aic <- -2 * logl + 2 * (dim(old)[2] - 1 + p) # 非线性多个p
  
  return(list(
    "para_iter" = old, # 迭代的所有结果
    "para" = c(new, est_scale), # 最终结果
    "est_sig0" = est_sig0, # 估计的sigma0，后续有用
    "logl" = logl, # aic
    "aic" = aic
  ))
}

# para <- list("est_etas" = rep(8, p), "est_delta" = rep(1, p), "est_sig0" = diag(p), "est_v" = 4)
# em_re_nolinear_T <- EM_nonlinear_T(para = para, max_iter = 5000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys)
#   

