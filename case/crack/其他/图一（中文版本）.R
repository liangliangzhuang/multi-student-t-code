library(showtext)
library(latex2exp)
showtext_auto()

# ===============
f_names <- list(TeX(r'(制动扭矩 ($N \cdot m$))'), TeX(r'(衔铁释放时间 (ms))'))
f_labeller <- function(variable, value){return(f_names[value])}

degradation.path.plot2 = function(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free"){
  data1 = map(data, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df <- bind_rows(data1, .id = "Unit")
  cal_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value")
  p1 = cal_dat %>% 
    ggplot(aes(Time,Value,color = factor(Unit))) + 
    geom_line(alpha=0.8) + geom_point(size=0.8) +
    facet_wrap(vars(factor(PC)),nrow = 1, scales = scale1,labeller = f_labeller) + 
    theme_bw() +
    # scale_color_viridis(discrete = T) + 
    scale_x_continuous(breaks = seq(0, m, by = ech),limits = c(0, m)) + 
    theme(legend.position = leg.pos) #panel.grid = element_blank()
  return(p1)
}

degradation.path.plot2(data = sim_cum_dat, leg.pos = "none",ech = 5, scale1 = "free") +
  # scale_y_reverse() + 
  scale_color_viridis(discrete = TRUE) + 
  theme(legend.position = "none") +
  xlab("时间 (天数)") + ylab("退化值")
# ggsave("case/result/PMB-dat（中文）.pdf", height = 4, width = 5)


# ===============

ggplot(sim3[,1:2], aes(x = PC1, y = PC2)) +
  stat_density_2d(contour = TRUE, color = "black", size=0.7, bins = 20, aes(linetype = "Real")) + 
  geom_point(alpha=0.4,color="#34AC93") +
  geom_density_2d(data = bivariate_data, aes(x = V1, y = V2,linetype = "Bivariate Gaussian"),
                  color = "#8D82CD", bins = 20,size=0.7, alpha=1) +
  xlab(TeX(r'(制动扭矩 ($N \cdot m$))')) + 
         ylab(TeX(r'(衔铁释放时间 (ms))')) + theme_bw() + 
  theme(panel.grid = element_blank(),legend.position = "none") +
  annotate(geom = "line",
           x = c(11.6,13),size=0.7,
           y = c(12,9), color="#8D82CD", alpha=1,
           arrow = arrow(angle = 30, length = unit(2, "mm"))) +
  annotate(geom = "text",
           x = 13.5, y = 8, size=4,
           label = "双维纳过程",color="#8D82CD", alpha=1) 


# ggsave("case/result/contour_plot（中文）.pdf", height = 4, width = 5)


