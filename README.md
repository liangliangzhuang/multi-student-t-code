# Codes for the paper "A Multivariate Student-t Process Model for Dependent Tail-weighted Degradation Data"

The repository structure relevant to this article is outlined below:

- case
  - `case_fct.R`:
    - Functions needed for case analysis, including `degradation.path.plot()`, `EM_iter_plot()`, `mean.path.fit.plot()`, `R_cal()`, `crack_path()`, and sub-functions for EM such as `sumqua()`, `sumqua0()`, `sumqua1()`, `sumqua2()`.
  - crack: 
    - `crack.R`: Main file for analyzing Fatigue Crack Size Data.
    - `Fatigue-crack-size.xlsx`: Corresponding data, originally introduced in [Meeker et al. (2022)](https://www.wiley.com/en-us/Statistical+Methods+for+Reliability+Data%2C+2nd+Edition-p-9781118115459) and further processed as described in Appendix H of [Fang et al. (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0377221721008985).
  - PMB:
    - `PMB_analysis.R`: Main file for PMB analysis.
    - `pmb_dat.RData`: Corresponding data.
    - `vis.R`: Graphical visualization, includes: `boxplot_path()`, `coutour_plot()`, `qqplot_PC()`.
    - Bootstrap: Contains bootstrap codes for four models to calculate interval estimates. Includes R files: `BS_nolinear_T.R`, `BS_linear_T.R`, `BS_nolinear_Wiener.R`, `BS_linear_Wiener.R`.
  - result: Stores results of two case analyses.

- `Model_est.R`: Contains two-stage parameter estimation methods for four models. It will be used for model fitting of two datasets.


## Tutorial

For optimal interaction with these codes, it is recommended to open "`student-t.Rproj`" using RStudio, install all necessary packages as initially specified, and proceed to execute the code sequentially, section by section.

You can run the `crack.R` and `PMB_analysis.R` scripts to analyze the two datasets, respectively. These include: data import and processing, exploratory data analysis, model fitting and selection, among other contents.


| Which results to reproduce              | Data File                            | Code File                                                                                                                                                            | Expected output                                                                                            | Run time at the above-specified computer conditions |
|-----------------------------------------|--------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|------------------------------------------------------|
| Figures 1, 2, 8, 9, Figure S1 and Table 2 | case/PMB/pmb_dat.RData               | Model_est.R<br>case/case_fct.R<br>case/PMB/vis.R<br>case/PMB/Bootstrap/<br>- BS_linear_T.R<br>- BS_nolinear_T.R<br>- BS_linear_Wiener.R<br>- BS_nolinear_Wiener.R | The results are in case/result/PMB and include: Figures 1, 2, 8, 9, Figure S1 and Table 2.                  | 10 minutes                                           |
| Figure 10, Figure S2, and Table 3        | case/crack/Fatigue-crack-size.xlsx   | Model_est.R<br>case/case_fct.R<br>crack/crack.R                                                                                                                      | The results are in case/result/crack and include: Figure 10, Figure S2, and Table 3 (AIC and point estimate) | 1 minute                                            |



## Note 

The paper has been accepted for publication in IISE Transactions. If you use the provided code in this project, please remember to cite the paper accordingly. Detailed citation information is

```bibtex
@article{xu2024multivariate,
  title={A Multivariate Student-t Process Model for Dependent Tail-weighted Degradation Data},
  author={Xu, Ancha and Fang, Guanqi and Zhuang, Liangliang and Gu, Cheng},
  journal={IISE Transactions},
  year={2024},
  publisher={Taylor & Francis}
}
```

If you have any questions or need help with the code, please submit them in the [issue](https://github.com/liangliangzhuang/multi-student-t-code/issues).












