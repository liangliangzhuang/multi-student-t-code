# 自助法 ==================
em_para_final = em_re_linear_T[[2]]
est_v = em_para_final[length(em_para_final)-p]
bt_ci <- matrix(NA, item, 3*p + 1 + p*(p+1)/2) #p=2时  10
for(h in 1:item){
  # 1. 产生伪样本数据
  bt_dat = sim_path(par = em_re_linear_T[[2]], v = est_v, SIG0 = em_re_linear_T[[3]], scen = "linear")
  y.diff_bt = bt_dat[[1]]; sumys_bt = bt_dat[[2]]; y_bt = bt_dat[[3]] #EM中常用数据
  # y
  # sim_diff_dat = sim_dat[[3]] # 便于绘图（差分）
  # bt_cum_dat = bt_dat[[4]] # 便于绘图（累计）
  # degradation.path.plot(data = bt_cum_dat, leg.pos = "none",ech = 5, scale1 = "free")
  # 2. EM 估计
  para_bt <- list("est_etas" = em_re_linear_T[[2]][1:p], "est_delta" = em_re_linear_T[[2]][1:p + p],
                  "est_sig0" = em_re_linear_T[[3]], "est_v" = est_v)
  
  bt_re = tryCatch(EM_linear_T(max_iter = 5000, eps = 10^-5, y = y_bt, y.diff= y.diff_bt, sumys = sumys_bt),
                   error = function(e) return(NA))
  
  if(all(is.na(bt_re))){bt_ci[h,] = NA} else{bt_ci[h,] = bt_re$para}
  print(h)
}

bt_ci_T2 <- bt_ci
quan_linear_T <- data.frame(round(apply(bt_ci_T2, 2, quantile, c(0.05, 0.5, 0.95), na.rm = TRUE), 5))
colnames(quan_linear_T) = c(paste0("eta",1:p, sep = ""), paste0("delta",1:p, sep = ""),
                            paste0("sigma",1:p, sep = ""), 
                            paste0("rho",1:(p*(p-1)/2), sep = ""), "v", 
                            paste0("gamma",1:p, sep = ""))
quan_linear_T



