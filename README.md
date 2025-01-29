# Replication codes for paper on *"Pruned skewed Kalman filter and smoother with application to DSGE models"*

This repository contains the replication codes for the paper *"Pruned Skewed Kalman Filter and Smoother with Application to DSGE Models"* authored by Gaygysyz Guljanov, Willi Mutschler, and Mark Trede.
The preprint is available at [this link](https://mutschler.eu/files/papers/GuljanovMutschlerTrede_PSKF.pdf).
Previous versions of the paper were circulated under the title *"Pruned Skewed Kalman Filter and Smoother: With Application to the Yield Curve"*, which can be accessed through [CQE Working Paper 101](https://www.wiwi.uni-muenster.de/cqe/sites/cqe/files/CQE_Paper/cqe_wp_101_2022.pdf) and [Dynare Working Paper no. 78](https://www.dynare.org/wp-repo/dynarewp078.pdf).

## Monte Carlo Study
The Online Appendix includes detailed tables from the Monte Carlo study.
You can replicate these tables using the scripts `table1_table2.m`, `table3.m`, and `table4.m`, which are located in the `MATLAB` folder.
The results can also be found in the `results` folder.

## Ireland (2004) with Dynare
All results are organized in folders and log files corresponding to specific mod files.
You can execute any of the following tasks individually or run them all in the specified order.

- The PSKF toolbox is currently under development and is expected to be released with Dynare 7.0.
You can track the development progress on [Dynare's GitLab](https://git.dynare.org/wmutschl/dynare/-/tree/pskf).
Installation files are available in the artifacts of the pipelines (under `pkg`, download the appropriate artifacts for your platform).
For easy replication and personal studies, you can directly download installation packages in the folder `dynare-installation-files` in this repository.

- Ensure that you add the `matlab` folder to your PATH in MATLAB or Octave.
Refer to the instructions in [the manual](https://www.dynare.org/manual/installation-and-configuration.html#configuration).

- Adjust `CPUnbr` and `MatlabOctavePath` in the `__parallelConf.ini` file according to your computer's core count and the path to the MATLAB/Octave binaries.
For more details, see [Dynare's parallel computing capabilities](https://www.dynare.org/manual/the-configuration-file.html#parallel-configuration).
Some parts of the code (e.g., initial value search) utilize MATLAB's dedicated parallel computing toolbox.

- The reported runtimes are approximate and depend on the number of cores available and the platform of your machine.
The reported times are based on:
  - Apple **MacBook** Pro M2 Max (8 performance, 4 efficiency cores), 64 GB RAM, macOS Sequoia 15.2, MATLAB R2024b Update 3 (24.2.0.2806996) 64-bit (maca64) with Dynare 7-unstable
  - Lenovo ThinkSystem SR655 **Linux server** (AMD EPYC 7402P 24C 2.8GHz), 6x16GB TruDDR4 3200MHz, Pop!_OS 20.04, MATLAB R2024b Update 3 (24.2.0.2806996) 64-bit (glnxa64) with Dynare 7 compiled from source
- Note that comparing the Gaussian and CSN model versions may not be meaningful, as the CSN model has more parameters to estimate.

### Maximum likelihood estimation

#### Task ML 1: maximum likelihood estimation of Gaussian model
Estimate the Gaussian model with maximum likelihood using the PSKF to compute the likelihood.
```matlab
dynare ireland2004_ml_1_gaussian
```
Runtime: 2 minutes on MacBook, 1 minute on Linux server.

#### Task ML 2: initial values search for maximum likelihood estimation of CSN model
Search for initial values for the maximum likelihood estimation of the CSN model (*Task ML 3*) using the PSKF to compute the likelihood.
To this end:
1. Fix model and stderr parameters to Gaussian estimates from maximum likelihood from *Task ML 1*.
1. Create large grid for skew parameters.
1. Evaluate likelihood (computed by PSKF) on grid using MATLAB's parallel computing toolbox.
1. Use best five skew parameter combinations as initial guess and optimize for each initial guess the likelihood (computed by PSKF) over both stderr and skew parameters.
```matlab
dynare ireland2004_ml_2_csn_initval_search
```
Runtime: 23 minutes on MacBook, 33 minutes on Linux server.

#### Task ML 3: maximum likelihood estimation of CSN model
Estimate the CSN model with maximum likelihood using the PSKF to compute the likelihood.
The optimizer is initialized with the best optimized values (highest log-likelihood value) found in *Task ML 2*.
```matlab
dynare ireland2004_ml_3_csn
```
Runtime: 20 minutes on MacBook, 36 minutes on Linux server.

### Bayesian estimation
The easiest estimation strategy is to use the Slice sampler for the full estimation (see *Task Bayes 1* and *Task Bayes 2*), but we also provide several different ways to estimate the model with the Random-Walk Metropolis-Hastings (RWMH) sampler (as this is standard practice).

The Slice sampler is usually more efficient (less auto-correlated draws, lower inefficiency factors) and more importantly requires no fine-tuning.
The downside is a longer runtime, because each draw requires many more posterior function evaluations.

The RWMH sampler requires a mode-finding step (with a positive definite inverse Hessian at the mode) and also a tuning of the *mh_jscale* to get a desired acceptance ratio.
On the upside, once it is initialized and fine-tuned, each draw requires only one posterior function evaluation.

Ultimately, the posterior distributions are nearly identical whether estimated using the Slice sampler or any fine-tuned variant of the RWMH sampler.


#### Task Bayes 1: Bayesian estimation with Slice sampler of Gaussian model
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood. 
The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
No fine-tuning is required, approximately 68 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_1_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 37 minutes on MacBook, 1 hour and 25 minutes on Linux server.

#### Task Bayes 2: Bayesian estimation with Slice sampler of CSN model
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
No fine-tuning is required, approximately 86 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_2_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 14 hours and 49 minutes on MacBook, 18 hours and 3 minutes on Linux server.

#### Task Bayes 3: Bayesian estimation with short Slice sampler of Gaussian model (to get close to posterior mode)
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
These draws are used to initialize the covariance matrix for the RWMH sampler (*Task Bayes 5*) and as initial guess for finding the posterior mode using optimization (*Task Bayes 7*).
No fine-tuning is required, approximately 68 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_3_mode_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 3 minutes on MacBook, 6 minutes on Linux server.

#### Task Bayes 4: Bayesian estimation with short Slice sampler of CSN model (to get close to posterior mode)
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
These draws are used to initialize the covariance matrix for the RWMH sampler (*Task Bayes 6*) and as initial guess for finding the posterior mode using optimization (*Task Bayes 8*).
No fine-tuning is required, approximately 86 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_4_mode_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 42 minutes on MacBook, 59 minutes on Linux server.

#### Task Bayes 5: Bayesian estimation with RWMH sampler of Gaussian model, initialized with short Slice from Task Bayes 3
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in *Task Bayes 3*.
No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 30%).
```matlab
dynare ireland2004_bayes_5_rwmh_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 27 minutes on MacBook, 1 hour and 2 minutes on Linux server.

#### Task Bayes 6: Bayesian estimation with RWMH sampler of CSN model, initialized with short Slice from Task Bayes 4
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in *Task Bayes 4*.
No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 33%).
```matlab
dynare ireland2004_bayes_6_rwmh_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 5 hours and 33 minutes on MacBook, 8 hours and 15 minutes on Linux server.

#### Task Bayes 7: Bayesian mode estimation of Gaussian model
Estimate the mode of the Gaussian model with numerical optimization using the PSKF to compute the likelihood.
The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 3*.
First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 9*.
```matlab
dynare ireland2004_bayes_7_mode_gaussian
```
Runtime: 35 minutes on MacBook, 33 minutes on Linux server.

#### Task Bayes 8: Bayesian mode estimation of CSN model
Estimate the mode of the CSN model with numerical optimization using the PSKF to compute the likelihood.
The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 4*.
First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 10*.
```matlab
dynare ireland2004_bayes_8_mode_csn
```
Runtime: 6 hours and 11 minutes on MacBook, 7 hours on Linux server.

#### Task Bayes 9: Bayesian estimation with RWMH sampler of Gaussian model
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*) in *Task Bayes 7*.
Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 7*.
```matlab
dynare ireland2004_bayes_9_rwmh_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 28 minutes on MacBook, 1 hour and 2 minutes on Linux server.

#### Task Bayes 10: Bayesian estimation with RWMH sampler of CSN model
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*) in *Task Bayes 8*.
Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 8*.
```matlab
dynare ireland2004_bayes_10_rwmh_csn parallel conffile=__parallelConf.ini
```
Runtime: 5 hours and 48 minutes on MacBook, 8 hours and 27 minutes on Linux server.

### Impulse response functions
Compute the impulse response functions of the Gaussian and CSN model using the 16th and 84th percentiles of the estimated distributions.
```matlab
dynare ireland2004_irfs
```
Runtime: 1 minute on MacBook, 1 minute on Linux server.

### Simulation of recessions
Simulates large time series with Gaussian and CSN distributed shocks and computes statistics on recessions.
```matlab
dynare ireland2004_recessions
```
Runtime: 1 minute on MacBook, 1 minute on Linux server.

## Publication quality figures
The figures reported in the paper have been generated with the scripts in the folder `R`.
These scripts are written in the `R` programming language and require several packages to be installed.
The figures are saved as `eps` files.

## Replicating results from working paper
This repository also contains replication codes for previous versions of the paper.
Namely, in the working paper versions of the paper, we estimate a Dynamic Nelson-Siegel model (Diebold, Rudebusch and Aruoba, 2006) of the yield curve using maximum likelihood in MATLAB.
Additionally, we provide non-Dynare codes for estimating the Ireland (2004) model with maximum likelihood in MATLAB, though not with Bayesian methods.
While the Dynare toolbox is significantly more user-friendly and widely adopted among quantitatively oriented macroeconomists, our alternative approach aims to offer complementary insights.
