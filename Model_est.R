library(minpack.lm)

EM_nonlinear_T <- function(para, type, max_iter = 5000, eps = 10^-5,
                           y, y.diff, sumys) {
  # Stage 1 (onlinear least squares estimation)======
  gamma_est <- function(data = y, t = t, type) {
    if (type == "power") {
      model <- function(x, a, b) {
        # power form
        a + b * x
      }
      est_scale <- matrix(NA, p, n)
      est_eta1 <- c()
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
        est_eta1[j] <- coef(fit)[1]
      }
      est_scale <- rowMeans(est_scale) # gamma
      est_t <- t^est_scale # power
    } else if (type == "exp") {
      model <- function(x, a, b) {
        # exp form
        a * (exp(b * x) - 1) # eta * (exp(at) - 1)
      }
      est_scale <- est_eta1 <- c()
      for (j in 1:p) {
        x <- t[j, 2:(m + 1)] # t
        nls_dat <- c()
        for (i in 1:n) {
          z <- data[j, 2:(m + 1), i] # y
          z1 <- as.numeric(na.omit(z))
          x1 <- x[!is.na(z)]
          nls_dat <- rbind(nls_dat, cbind(x1, z1))
        }

        fit <- nlsLM(z1 ~ model(x1, a, b),
          start = list(a = 1, b = 1),
          control = nls.control(maxiter = 1000)
        )
        est_scale[j] <- coef(fit)[2]
        est_eta1[j] <- coef(fit)[1]
      }
      est_t <- exp(est_scale * t) - 1 # exp
    }

    est_t.diff <- est_t[, 2:(m + 1)] - est_t[, 1:m] # Δλ(t)
    SIG <- array(0, dim = c(p, p, m))
    SIG_inv <- array(0, dim = c(p, p, m))
    for (k in 1:m) {
      SIG[, , k] <- diag(est_t.diff[, k])
      SIG_inv[, , k] <- solve(SIG[, , k])
    }
    sumSIG <- diag(est_t[, (m + 1)])

    return(list(
      "est_scale" = est_scale,
      "SIG" = SIG,
      "SIG_inv" = SIG_inv,
      "sumSIG" = sumSIG,
      "est_t.diff" = est_t.diff,
      "est_eta" = est_eta1
    ))
  }
  
  hat_gamma <- gamma_est(data = y, t = t, type = type)
  est_scale <- hat_gamma[[1]] # gamma
  SIG <- hat_gamma[[2]]
  SIG_inv <- hat_gamma[[3]]
  sumSIG <- hat_gamma[[4]]
  est_t.diff <- hat_gamma[[5]]

  # EM algorithm ===============================================
  est_etas <- hat_gamma[[6]] # eta estimates from the previous stage
  est_delta <- para[[2]]
  est_sig0 <- para[[3]]
  est_v <- para[[4]]

  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2 * p + 1 + p * (p + 1) / 2) # (eta[p], delta[p], sig0 [2p-1], v[1])

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
    inits4 <- solve(SIGdelta) # \Omega_\delta
    inits5 <- inits4 %*% sumSIG
    sigtheta <- solve(inits1 + inits5)
    inits6 <- sumSIG %*% sigtheta

    # E(theta)
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i] # The second symbol of mu_i
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua0(y.diff[, , i], inits4, SIG_inv)
    }
    est_theta <- etheta # Stores a posteriori expectation for each iteration based on the parameter values of the previous step
    # Eq (14)
    etau <- (m * p + est_v) / (k1s - k2s + est_v)
    est_tau <- etau # E(\tau_i|xxx)
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
    est_ltau <- eltau # E(\ln{\tau_i}|xxx)

    # M steps
    ## eta =====
    est_etas <- sumqua1(etau, etheta) 
    ## Sig_0 =====
    est_sig0 <- sumqua2(etau, etheta - matrix(est_etas, p, n), sigtheta) / n
    ## Omega_delta ====
    for (j in 1:p) {
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

    # Convergence stop condition
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)]) # sigma + rho
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # eta, delta, sig0, v
    # Output sigma_i,j
    est_sig0_new2 <- sqrt(c(diag(est_sig0), est_sig0[upper.tri(est_sig0, diag = F)]))
    new2 <- c(est_etas, est_delta, est_sig0_new2, est_v)
    if (all((abs(new - old[iter, ])) < eps) | iter > max_iter) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }

  # Q-fun ===========
  est_E_tau <- etau # est_E_tau
  Sigma_0_inv <- solve(est_sig0) # Sigma_0_inv - Σ_0 inverse matrix
  Omega_delta_inv <- solve(diag(est_delta)) # Omega_delta_inv - Ω_δ inverse matrix
  mu_i_s <- etheta
  eta <- est_etas # eta - η vector
  Sigma_Theta_i_s <- sigtheta # Sigma_Theta_i_s - The Σ_Θ_i^(s) matrix for each I

  # \Delta Sigma(t_{i,k}) ======
  Delta_Sigma <- list()
  for (i in 1:n) {
    Delta_Sigma[[i]] <- list()
    for (k in 1:m) {
      Delta_Sigma[[i]][[k]] <- diag(est_t.diff[, k]) # Sigma_diff_t
    }
  }

  # Compute est_E_tau_ell_0 和 est_E_tau_ell_k =======
  est_E_tau_ell_0 <- numeric()
  est_E_tau_ell_k <- list()
  for (i in 1:n) {
    # E[tau_i ell_i,0 | ΔY_i, Ψ^(s)]
    trace_part <- sum(diag(Sigma_0_inv %*% Sigma_Theta_i_s))
    quadratic_part <- t(mu_i_s[, i] - eta) %*% Sigma_0_inv %*% (mu_i_s[, i] - eta)
    ell_i0 <- as.numeric(trace_part + est_E_tau[i] * quadratic_part)
    # E[tau_i ell_i,k | ΔY_i, Ψ^(s)]
    ell_ik_vector <- numeric()
    for (k in 1:m) {
      Delta_Sigma_inv <- solve(Delta_Sigma[[i]][[k]])
      trace_ik <- sum(diag(Delta_Sigma[[i]][[k]] %*% Omega_delta_inv %*% Sigma_Theta_i_s))
      quadratic_ik <- t(mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i]) %*%
        solve(diag(est_delta) %*% Delta_Sigma_inv) %*%
        (mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i])
      ell_ik_vector[k] <- as.numeric(trace_ik + est_E_tau[i] * quadratic_ik)
    }

    # save
    est_E_tau_ell_0[i] <- ell_i0
    est_E_tau_ell_k[[i]] <- ell_ik_vector
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
    part2[i] <- tau_ell_sum + est_E_tau_ell_0[i] + est_v * est_E_tau[i] # est_E_tau
  }
  # Combine part1 and part2 to get the final Q function value
  Q_value <- sum(part1 - 1 / 2 * part2)
  aic <- -2 * Q_value + 2 * (dim(old)[2] - 1 + p) 

  return(list(
    "para_iter" = old, # All results of iteration
    "para" = c(new, est_scale), # final para 
    "est_sig0" = est_sig0,
    "logl" = Q_value, 
    "aic" = aic, # aic
    "para2" = c(new2, est_scale) # final para(2)
  ))
}

EM_linear_T <- function(max_iter = 4000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys) {
  t.diff1 <- time
  t1 <- m * time
  est_t.diff <- matrix(rep(1, m * p), nrow = p)

  SIGMA <- diag(t1, p)
  sig <- diag(time, p) ## λ(t)=t
  sig_inv <- solve(sig)

  ## Inital parameters
  est_etas <- rep(8, p)
  est_delta <- rep(1, p)
  est_sig0 <- diag(p)
  est_v <- 4
  # Linear form
  est_tau <- rep(1, n)
  est_ltau <- rep(5, n)
  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2 * p + 1 + p * (p + 1) / 2) # Record the result of the iteration
  while (con == FALSE) {
    iter <- iter + 1
    est_etas_old <- est_etas
    est_delta_old <- est_delta
    est_sig0_cor <- cov2cor(est_sig0)
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter, ] <- c(est_etas_old, est_delta_old, est_sig0_old, est_v)

    k1s <- rep(0, n)
    k2s <- rep(0, n)
    sigw <- diag(est_delta^2)

    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(sigw)
    inits5 <- inits4 %*% SIGMA
    sigtheta <- solve(inits1 + inits5) 
    inits6 <- SIGMA %*% sigtheta

    ## expection of thetas
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i]
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua(y.diff[, , i], inits4 %*% sig_inv)
    }
    est_theta <- etheta 
    etau <- (m * p + est_v) / (k1s - k2s + est_v)
    est_tau <- etau
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
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
      re <- -2 * log(gamma(x / 2)) + x * log(x / 2) + x / n * sum(eltau - etau)
      return(-re)
    }
    est_v <- optimize(est_v_fun, c(0, 30))$minimum

    ### Convergence stop condition
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)])
    new <- c(est_etas, est_delta, est_sig0_new, est_v) 
    est_sig0_new2 <- sqrt(c(diag(est_sig0), est_sig0[upper.tri(est_sig0, diag = F)]))
    new2 <- c(est_etas, est_delta, est_sig0_new2, est_v)
    if (all((abs(new - old[iter, ])) < eps)) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }

  # Q-fun ===========
  est_E_tau <- etau 
  Sigma_0_inv <- solve(est_sig0) 
  Omega_delta_inv <- solve(diag(est_delta)) 
  mu_i_s <- etheta
  eta <- est_etas 
  Sigma_Theta_i_s <- sigtheta 

  # \Delta Sigma(t_{i,k}) ======
  Delta_Sigma <- list()
  for (i in 1:n) {
    Delta_Sigma[[i]] <- list()
    for (k in 1:m) {
      Delta_Sigma[[i]][[k]] <- diag(est_t.diff[, k])
    }
  }

  # est_E_tau_ell_0 and est_E_tau_ell_k =======
  est_E_tau_ell_0 <- numeric()
  est_E_tau_ell_k <- list()
  for (i in 1:n) {
    # E[tau_i ell_i,0 | ΔY_i, Ψ^(s)]
    trace_part <- sum(diag(Sigma_0_inv %*% Sigma_Theta_i_s))
    quadratic_part <- t(mu_i_s[, i] - eta) %*% Sigma_0_inv %*% (mu_i_s[, i] - eta)
    ell_i0 <- as.numeric(trace_part + est_E_tau[i] * quadratic_part)
    # E[tau_i ell_i,k | ΔY_i, Ψ^(s)]
    ell_ik_vector <- numeric()
    for (k in 1:m) {
      Delta_Sigma_inv <- solve(Delta_Sigma[[i]][[k]])
      trace_ik <- sum(diag(Delta_Sigma[[i]][[k]] %*% Omega_delta_inv %*% Sigma_Theta_i_s))
      quadratic_ik <- t(mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i]) %*%
        solve(diag(est_delta) %*% Delta_Sigma_inv) %*%
        (mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i])
      ell_ik_vector[k] <- as.numeric(trace_ik + est_E_tau[i] * quadratic_ik)
    }

    # Save
    est_E_tau_ell_0[i] <- ell_i0
    est_E_tau_ell_k[[i]] <- ell_ik_vector
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
    part2[i] <- tau_ell_sum + est_E_tau_ell_0[i] + est_v * est_E_tau[i] 
  }
  Q_value <- sum(part1 - 1 / 2 * part2)
  aic <- -2 * Q_value + 2 * (dim(old)[2] - 1) 

  return(list(
    "para_iter" = old, 
    "para" = c(new, rep(NA, p)), 
    "est_sig0" = est_sig0, 
    "logl" = Q_value, 
    "aic" = aic,
    "para2" = c(new2, rep(NA, p)) 
  ))
}

EM_linear_Wiener <- function(max_iter = 1000, eps = 10^-5, y = y, y.diff = y.diff, sumys = sumys) {
  t.diff1 <- time
  t1 <- m * time
  est_t.diff <- matrix(rep(1, m * p), nrow = p)

  SIGMA <- diag(t1, p)
  sig <- diag(time, p) 
  sig_inv <- solve(sig)

  est_etas <- rep(8, p)
  est_delta <- rep(1, p)
  est_sig0 <- diag(p)
  est_v <- 100
  est_tau <- rep(1, n)
  est_ltau <- rep(5, n)
  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2 * p + 1 + p * (p + 1) / 2) 
  while (con == FALSE) {
    iter <- iter + 1
    est_etas_old <- est_etas
    est_delta_old <- est_delta
    est_sig0_cor <- cov2cor(est_sig0)
    est_sig0_old <- c(sqrt(diag(est_sig0)), est_sig0_cor[upper.tri(est_sig0_cor, diag = F)])
    old[iter, ] <- c(est_etas_old, est_delta_old, est_sig0_old, est_v)

    k1s <- rep(0, n)
    k2s <- rep(0, n)
    sigw <- diag(est_delta^2)

    inits1 <- solve(est_sig0)
    inits2 <- inits1 %*% est_etas
    inits3 <- est_etas %*% inits2
    inits4 <- solve(sigw)
    inits5 <- inits4 %*% SIGMA
    sigtheta <- solve(inits1 + inits5) 
    inits6 <- SIGMA %*% sigtheta

    # E(theta)
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i]
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua(y.diff[, , i], inits4 %*% sig_inv)
    }
    est_theta <- etheta
    etau <- (m * p + est_v) / (k1s - k2s + est_v)
    est_tau <- etau
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
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
    est_v <- 100

    ### Convergence stop condition
    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)])
    new <- c(est_etas, est_delta, est_sig0_new, est_v) 
    if (all((abs(new - old[iter, ])) < eps)) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }


  # Q-fun ===========
  est_E_tau <- etau 
  Sigma_0_inv <- solve(est_sig0) 
  Omega_delta_inv <- solve(diag(est_delta)) 
  mu_i_s <- etheta
  eta <- est_etas
  Sigma_Theta_i_s <- sigtheta 

  Delta_Sigma <- list()
  for (i in 1:n) {
    Delta_Sigma[[i]] <- list()
    for (k in 1:m) {
      Delta_Sigma[[i]][[k]] <- diag(est_t.diff[, k])
    }
  }

  est_E_tau_ell_0 <- numeric()
  est_E_tau_ell_k <- list()
  for (i in 1:n) {
    trace_part <- sum(diag(Sigma_0_inv %*% Sigma_Theta_i_s))
    quadratic_part <- t(mu_i_s[, i] - eta) %*% Sigma_0_inv %*% (mu_i_s[, i] - eta)
    ell_i0 <- as.numeric(trace_part + est_E_tau[i] * quadratic_part)
    ell_ik_vector <- numeric()
    for (k in 1:m) {
      Delta_Sigma_inv <- solve(Delta_Sigma[[i]][[k]])
      trace_ik <- sum(diag(Delta_Sigma[[i]][[k]] %*% Omega_delta_inv %*% Sigma_Theta_i_s))
      quadratic_ik <- t(mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i]) %*%
        solve(diag(est_delta) %*% Delta_Sigma_inv) %*%
        (mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i])
      ell_ik_vector[k] <- as.numeric(trace_ik + est_E_tau[i] * quadratic_ik)
    }

    est_E_tau_ell_0[i] <- ell_i0
    est_E_tau_ell_k[[i]] <- ell_ik_vector
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
    part2[i] <- tau_ell_sum + est_E_tau_ell_0[i] + est_v * est_E_tau[i] 
  }
  Q_value <- sum(part1 - 1 / 2 * part2)
  aic <- -2 * Q_value + 2 * (dim(old)[2] - 1) 

  return(list(
    "para_iter" = old, 
    "para" = c(new, rep(NA, p)),
    "est_sig0" = est_sig0,
    "logl" = Q_value, 
    "aic" = aic
  ))
}

EM_nonlinear_Wiener <- function(para, type, max_iter = 5000, eps = 10^-5,
                                y, y.diff, sumys) {

  gamma_est <- function(data = y, t = t, type) {
    if (type == "power") {
      model <- function(x, a, b) {
        # power form
        a + b * x
      }
      est_scale <- matrix(NA, p, n)
      est_eta1 <- c()
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
        est_eta1[j] <- coef(fit)[1]
      }
      est_scale <- rowMeans(est_scale) # gamma
      est_t <- t^est_scale # power
    } else if (type == "exp") {
      model <- function(x, a, b) {
        # exp form
        a * (exp(b * x) - 1) # eta * (exp(at) - 1)
      }
      est_scale <- est_eta1 <- c()
      for (j in 1:p) {
        x <- t[j, 2:(m + 1)] # t
        nls_dat <- c()
        for (i in 1:n) {
          z <- data[j, 2:(m + 1), i] # y
          z1 <- as.numeric(na.omit(z))
          x1 <- x[!is.na(z)]
          nls_dat <- rbind(nls_dat, cbind(x1, z1))
        }

        fit <- nlsLM(z1 ~ model(x1, a, b),
          start = list(a = 1, b = 1),
          control = nls.control(maxiter = 1000)
        )
        est_scale[j] <- coef(fit)[2]
        est_eta1[j] <- coef(fit)[1]
      }
      est_t <- exp(est_scale * t) - 1 # exp
    }

    est_t.diff <- est_t[, 2:(m + 1)] - est_t[, 1:m] # Δλ(t)
    SIG <- array(0, dim = c(p, p, m))
    SIG_inv <- array(0, dim = c(p, p, m))
    for (k in 1:m) {
      SIG[, , k] <- diag(est_t.diff[, k])
      SIG_inv[, , k] <- solve(SIG[, , k])
    }
    sumSIG <- diag(est_t[, (m + 1)])

    return(list(
      "est_scale" = est_scale,
      "SIG" = SIG,
      "SIG_inv" = SIG_inv,
      "sumSIG" = sumSIG,
      "est_t.diff" = est_t.diff,
      "est_eta" = est_eta1
    ))
  }
  hat_gamma <- gamma_est(data = y, t = t, type = type)
  est_scale <- hat_gamma[[1]] # gamma
  SIG <- hat_gamma[[2]]
  SIG_inv <- hat_gamma[[3]]
  sumSIG <- hat_gamma[[4]]
  est_t.diff <- hat_gamma[[5]]

  # EM算法 ===============================================
  est_etas <- hat_gamma[[6]]
  est_delta <- para[[2]]
  est_sig0 <- para[[3]]
  est_v <- 100

  con <- FALSE
  iter <- 0
  old <- matrix(NA, max_iter + 6, 2 * p + 1 + p * (p + 1) / 2) #  (eta[p], delta[p], sig0 [2p-1], v[1])

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
    sigtheta <- solve(inits1 + inits5) 
    inits6 <- sumSIG %*% sigtheta

    # E(theta)
    etheta <- matrix(0, p, n)
    for (i in 1:n) {
      s1 <- inits2 + inits4 %*% sumys[, i] 
      etheta[, i] <- sigtheta %*% s1
      k2s[i] <- t(s1) %*% etheta[, i]
      k1s[i] <- inits3 + sumqua0(y.diff[, , i], inits4, SIG_inv)
    }
    est_theta <- etheta 
    etau <- (m * p + est_v) / (k1s - k2s + est_v)
    est_tau <- etau
    eltau <- digamma((m * p + est_v) / 2) - log((k1s - k2s + est_v) / 2)
    est_ltau <- eltau

    # M steps
    ## eta ======
    est_etas <- sumqua1(etau, etheta) 
    ## Sig_0 ======
    est_sig0 <- sumqua2(etau, etheta - matrix(est_etas, p, n), sigtheta) / n
    ## Omega_delta ====
    for (j in 1:p) { 
      dj <- 0
      for (i in 1:n) {
        dj <- dj + inits6[j, j] + etau[i] * (sum(y.diff[j, , i]^2 / est_t.diff[j, ]) + etheta[j, i]^2 * sumSIG[j, j] - 2 * sumys[j, i] * etheta[j, i])
      }
      est_delta[j] <- sqrt(dj / m / n)
    }
    est_v <- 100

    est_sig0_cor2 <- cov2cor(est_sig0)
    est_sig0_new <- c(sqrt(diag(est_sig0)), est_sig0_cor2[upper.tri(est_sig0_cor2, diag = F)]) 
    new <- c(est_etas, est_delta, est_sig0_new, est_v) # eta, delta, sig0, v
    if (all((abs(new - old[iter, ])) < eps) | iter > max_iter) {
      con <- TRUE
      old[iter + 1, ] <- new
    }
  }

  # Q-fun ===========
  est_E_tau <- etau 
  Sigma_0_inv <- solve(est_sig0) 
  Omega_delta_inv <- solve(diag(est_delta)) 
  mu_i_s <- etheta
  eta <- est_etas 
  Sigma_Theta_i_s <- sigtheta 

  Delta_Sigma <- list()
  for (i in 1:n) {
    Delta_Sigma[[i]] <- list()
    for (k in 1:m) {
      Delta_Sigma[[i]][[k]] <- diag(est_t.diff[, k]) 
    }
  }

  # est_E_tau_ell_0 和 est_E_tau_ell_k =======
  est_E_tau_ell_0 <- numeric()
  est_E_tau_ell_k <- list()
  for (i in 1:n) {
    # E[tau_i ell_i,0 | ΔY_i, Ψ^(s)]
    trace_part <- sum(diag(Sigma_0_inv %*% Sigma_Theta_i_s))
    quadratic_part <- t(mu_i_s[, i] - eta) %*% Sigma_0_inv %*% (mu_i_s[, i] - eta)
    ell_i0 <- as.numeric(trace_part + est_E_tau[i] * quadratic_part)
    # E[tau_i ell_i,k | ΔY_i, Ψ^(s)]
    ell_ik_vector <- numeric()
    for (k in 1:m) {
      Delta_Sigma_inv <- solve(Delta_Sigma[[i]][[k]])
      trace_ik <- sum(diag(Delta_Sigma[[i]][[k]] %*% Omega_delta_inv %*% Sigma_Theta_i_s))
      quadratic_ik <- t(mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i]) %*%
        solve(diag(est_delta) %*% Delta_Sigma_inv) %*%
        (mu_i_s[, i] - Delta_Sigma_inv %*% y.diff[, k, i])
      ell_ik_vector[k] <- as.numeric(trace_ik + est_E_tau[i] * quadratic_ik)
    }

    # Save
    est_E_tau_ell_0[i] <- ell_i0
    est_E_tau_ell_k[[i]] <- ell_ik_vector
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
    part2[i] <- tau_ell_sum + est_E_tau_ell_0[i] + est_v * est_E_tau[i]
  }

  Q_value <- sum(part1 - 1 / 2 * part2)
  aic <- -2 * Q_value + 2 * (dim(old)[2] - 1 + p) 

  return(list(
    "para_iter" = old, 
    "para" = c(new, est_scale), 
    "est_sig0" = est_sig0,
    "logl" = Q_value,
    "aic" = aic
  ))
}
