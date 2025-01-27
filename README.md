# Replication codes for paper on *"Pruned skewed Kalman filter and smoother with application to DSGE models"*

This repository contains the replication codes for the paper *"Pruned skewed Kalman filter and smoother with application to DSGE models"* by Gaygysyz Guljanov, Willi Mutschler, and Mark Trede.
Preprint available at [https://mutschler.eu](https://www.mutschler.eu/files/papers/GuljanovMutschlerTrede_PSKF.pdf).
Previous versions of the paper were circulated under the title *"Pruned Skewed Kalman Filter and Smoother:
With Application to the Yield Curve"* ([CQE Working Paper 101](https://www.wiwi.uni-muenster.de/cqe/sites/cqe/files/CQE_Paper/cqe_wp_101_2022.pdf) and [Dynare Working Paper no. 78](https://www.dynare.org/wp-repo/dynarewp078.pdf)).

## Monte-Carlo study
The Online Appendix contains detailed tables of the Monte-Carlo study.
The tables can be replicated with the scripts `table1_table2.m`, `table3.m`, and `table4.m`.
The required scripts are in the folder `MATLAB` and the results can also be accessed in the folder `results`.

## Ireland (2004) with Dynare
All the results are available in the folders and log files corresponding to the specific mod files.
You can run any of the following tasks individually or all of them in the specified order.

- The PSKF toolbox is still in development and planned for release with Dynare 7.0.
You can follow the development on [Dynare's GitLab](https://git.dynare.org/wmutschl/dynare/-/tree/pskf), the installation files are accessible in the artifacts of the pipelines (under `pkg`, download artifacts for your platform).
For easier replications and your own studies, you can also directly download the installation packages for [macOS arm64], [macOS x86_64] and [Windows x86_64].
For Linux please clone the repository and follow the instructions to [build Dynare from source](https://git.dynare.org/Dynare/dynare#building-dynare-from-source), which of course you can also do for macOS and Windows.

- Note that you need to add the `matlab` folder to your PATH in MATLAB or Octave, see these instructions in [the the manual](https://www.dynare.org/manual/installation-and-configuration.html#configuration).

- Make sure to adjust `CPUnbr` and `MatlabOctavePath` in the file `__parallelConf.ini` with respect to your computer's number of cores and path to the MATLAB/Octave binaries, see [Dynare's parallel computing capabilities](https://www.dynare.org/manual/the-configuration-file.html#parallel-configuration) for more details.
Some parts of the code (e.g. initial value search) also make use of MATLAB's dedicated parallel computing toolbox.

- The reported runtimes are approximate and depend on the number of cores available on and the platform of your machine.
The reported times are based on an Apple MacBook Pro M2 Max (8 performance, 4 efficiency cores), 64 GB RAM, MacOS Sequoia 15.2, MATLAB R2024b Update 3 (24.2.0.2806996) 64-bit (maca64) with Dynare 7 (dynare-7-unstable-2025-01-22-1305-b32e7530-arm64.pkg).
Note that a comparison between the Gaussian and CSN model versions is not really meaningful as the CSN model has more parameters to estimate.

### Maximum likelihood estimation

#### Task ML 1: maximum likelihood estimation of Gaussian model
Estimate the Gaussian model with maximum likelihood using the PSKF to compute the likelihood.
```matlab
dynare ireland2004_ml_1_gaussian
```
Runtime: 2 minutes on MacBook.

#### Task ML 2: initial values search for maximum likelihood estimation of CSN model
Search for initial values for the maximum likelihood estimation of the CSN model using the PSKF to compute the likelihood. To this end:
1. Fix model and stderr parameters to Gaussian estimates from maximum likelihood from *Task ML 1*.
1. Create large grid for skew parameters.
1. Evaluate likelihood (computed by PSKF) on grid using MATLAB's parallel computing toolbox.
1. Use best five skew parameter combinations as initial guess.
1. Optimize for each initial guess the likelihood (computed by PSKF) over both stderr and skew parameters.
```matlab
dynare ireland2004_ml_2_csn_initval_search
```
Runtime: 23 minutes on MacBook.

#### Task ML 3: maximum likelihood estimation of CSN model
Estimate the CSN model with maximum likelihood using the PSKF to compute the likelihood.
The optimizer is initialized with the best initial values (highest log-likelihood value) found in *Task ML 2*.
```matlab
dynare ireland2004_ml_3_csn
```
Runtime: 20 minutes on MacBook.

### Bayesian estimation
The easiest estimation strategy is to use the Slice sampler for the full estimation (see *Task Bayes 1* and *Task Bayes 2*), but we also provide several different ways to estimate the model with the Random-Walk Metropolis-Hastings (RWMH) sampler (as this is standard practice).

The Slice sampler is usually more efficient (less auto-correlated draws, lower inefficiency factors) and more importantly requires no fine-tuning.
The downside is a longer runtime, because each draw requires many more posterior function evaluations.

The RWMH sampler requires a mode-finding step (with a positive definite inverse Hessian at the mode) and also a tuning of the *mh_jscale* to get a desired acceptance ratio.
On the upside, once it is fine-tuned, each draw requires only one posterior function evaluation.

In the end, the posterior distributions are virtually the same when estimating with either the Slice or any fine-tuned RWMH sampler variant.


#### Task Bayes 1: Bayesian estimation with Slice sampler of Gaussian model
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood. 
The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
No fine-tuning is required, approximately 68 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_1_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 37 minutes on MacBook.

#### Task Bayes 2: Bayesian estimation with Slice sampler of CSN model
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
No fine-tuning is required, approximately 86 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_2_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 14 hours and 49 minutes on MacBook.

#### Task Bayes 3: Bayesian estimation with short Slice sampler of Gaussian model (to get close to posterior mode)
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
These draws are used to initialize the covariance matrix for the RWMH sampler (*Task Bayes 5*) and as initial guess for finding the posterior mode using optimization (*Task Bayes 7*).
No fine-tuning is required, approximately 68 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_3_mode_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 3 minutes on MacBook.

#### Task Bayes 4: Bayesian estimation with short Slice sampler of CSN model (to get close to posterior mode)
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
These draws are used to initialize the covariance matrix for the RWMH sampler (*Task Bayes 6*) and as initial guess for finding the posterior mode using optimization (*Task Bayes 8*).
No fine-tuning is required, approximately 86 function evaluations per iteration.
```matlab
dynare ireland2004_bayes_4_mode_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 42 minutes on MacBook.

#### Task Bayes 5: Bayesian estimation with RWMH sampler of Gaussian model, initialized with short Slice from Task Bayes 3
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in *Task Bayes 3*.
No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 30%).
```matlab
dynare ireland2004_bayes_5_rwmh_slice_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 27 minutes on MacBook.

#### Task Bayes 6: Bayesian estimation with RWMH sampler of CSN model, initialized with short Slice from Task Bayes 4
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in *Task Bayes 4*.
No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 33%).
```matlab
dynare ireland2004_bayes_6_rwmh_slice_csn parallel conffile=__parallelConf.ini
```
Runtime: 5 hours and 33 minutes on MacBook.

#### Task Bayes 7: Bayesian mode estimation of Gaussian model
Estimate the mode of the Gaussian model with numerical optimization using the PSKF to compute the likelihood.
The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 3*.
First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 9*.
```matlab
dynare ireland2004_bayes_7_mode_gaussian
```
Runtime: 35 minutes on MacBook.

#### Task Bayes 8: Bayesian mode estimation of CSN model
Estimate the mode of the CSN model with numerical optimization using the PSKF to compute the likelihood.
The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 4*.
First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 10*.
```matlab
dynare ireland2004_bayes_8_mode_csn
```
Runtime: 6 hours and 11 minutes on MacBook.

#### Task Bayes 9: Bayesian estimation with RWMH sampler of Gaussian model
Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=8*) in *Task Bayes 7*.
Tuning of *mh_jscale* required to get a desired acceptance ratio (about 33%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 7*.
```matlab
dynare ireland2004_bayes_9_rwmh_gaussian parallel conffile=__parallelConf.ini
```
Runtime: 28 minutes on MacBook.

#### Task Bayes 10: Bayesian estimation with RWMH sampler of CSN model
Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=8*) in *Task Bayes 8*.
Tuning of *mh_jscale* required to get a desired acceptance ratio (about 33%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 8*.
```matlab
dynare ireland2004_bayes_10_rwmh_csn parallel conffile=__parallelConf.ini
```
Runtime: 5 hours and 48 minutes on MacBook.


## Replicating results from working paper
In the previous working paper versions we used the same filter to estimate a Dynamic Nelson Siegel model of the yield curve with Maximum Likelihood in MATLAB. We also provide non-Dynare codes for estimating the Ireland (2004) model with Maximum Likelihood in MATLAB, but not with Bayesian methods as above.
