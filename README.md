# Replication codes for paper on *"Pruned skewed Kalman filter and smoother with application to DSGE models"*

This repository contains the replication codes for the paper *"Pruned Skewed Kalman Filter and Smoother with Application to DSGE Models"* authored by Gaygysyz Guljanov, Willi Mutschler, and Mark Trede.
The preprint is available at [this link](https://mutschler.eu/files/papers/GuljanovMutschlerTrede_PSKF.pdf).
Previous versions of the paper were circulated under the title *"Pruned Skewed Kalman Filter and Smoother: With Application to the Yield Curve"*, which can be accessed through [CQE Working Paper 101](https://www.wiwi.uni-muenster.de/cqe/sites/cqe/files/CQE_Paper/cqe_wp_101_2022.pdf) and [Dynare Working Paper no. 78](https://www.dynare.org/wp-repo/dynarewp078.pdf).

## Section 2: Closed skew normal distribution

- The R script `figure_csn_univariate.R` in the `R` folder creates Figure 1 of the paper (`figure_csn_univariate.pdf`).

- The R script `figure_csn_bivariate.R` in the `R` folder creates Figure 2 of the paper (`figure_csn_bivariate.pdf`).

## Section 4: Pruning the skewness dimension

The R script `figure_pruned_distribution.R` in the `R` folder creates Figure 3 of the paper (`figure_pruned_distribution.pdf`).
It also prints the correlation matrix before pruning, the pruned parameters, and the Kullback-Leibler divergence between the original and pruned distributions to the console.

## Section 6: Asymmetric shocks in a New Keynesian DSGE model

- The R script `figure_data.R` in the `R` folder creates Figure 4 of the paper (`figure_data.pdf`).
It also computes the summary statistics and prints them to the console.

- The MATLAB script `main_ireland2004.m` in the `Dynare` folder estimates the Gaussian and CSN models with maximum likelihood and Bayesian methods.
The estimation is done with Dynare in several steps and the results are saved in the `results/ireland2004` folder:
    - Table 2 of the paper is in `results/ireland2004/table2.tex`
    - The likelihood ratio test is in `results/ireland2004/LR_statistic.tex` and `results/ireland2004/LR_pvalue.tex`
    - The data for Figure 5 of the paper is in `results/ireland2004/`, the Figure is created with the R script `figure_irfs.R` in the `R` folder (`figure_irfs.pdf`)
    - Table 3 of the paper is in `results/ireland2004/table3.tex`

## Monte Carlo Study
The Online Appendix includes detailed tables from the Monte Carlo study.
You can replicate these tables using the scripts `table1_table2.m`, `table3.m`, and `table4.m`, which are located in the `MATLAB` folder.
The results can also be found in the `results` folder.


## Publication quality figures
The figures reported in the paper have been generated with the scripts in the folder `R`.
These scripts are written in the `R` programming language and require several packages to be installed.
The figures are saved as `eps` files.

## Replicating results from working paper
This repository also contains replication codes for previous versions of the paper.
Namely, in the working paper versions of the paper, we estimate a Dynamic Nelson-Siegel model (Diebold, Rudebusch and Aruoba, 2006) of the yield curve using maximum likelihood in MATLAB.
Additionally, we provide non-Dynare codes for estimating the Ireland (2004) model with maximum likelihood in MATLAB, though not with Bayesian methods.
While the Dynare toolbox is significantly more user-friendly and widely adopted among quantitatively oriented macroeconomists, our alternative approach aims to offer complementary insights.
