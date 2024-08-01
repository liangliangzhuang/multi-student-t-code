Real_path = function(led_cum, led_diff, p=p, m=m, n=n){
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
