# A Multivariate Student-t Process Model for Dependent Tail-weighted Degradation Data

> This work has been submitted to the IISE transactions and is currently under revision.


The repository related to this paper is structured into two primary directories, each designed for specific functions as outlined below:

1. 模拟实验
  - 参数估计评估：Figure 4-5 以及补充目录 Table S1 --> 模拟-画图.R --> 参数估计评估-0619.RData 
  	> 原始代码不见了，估计可以根据 `可靠度-main.R` 进行修改。
  - 可靠度评估：Figure 7 可靠度RMSE图.R --> 可靠度评估-0619.RData  
  	> 运行 `可靠度-main.R` 得到结果存储在 `result` 中。
  - 相关系数评估： Figure 6 --> 相关系数-main --> 相关系数评估-0617.RData
  - em：四种场景下的 EM 算法函数 (将用到上面三种评估中)
  - 其他图形： 保存模拟实验运行的结果

2. case
  - crack --> crack.R --> 0619-crack.RData
  - PMB --> PMB_analysis.R --> 0617-PMB-final.RData
  - `case_fct.R` 案例中所需的函数
  - result 存储结果

3. `simqua.R` EM中使用到的基础函数，主要计算E步中的内容。

