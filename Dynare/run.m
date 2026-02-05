% All results are organized in folders and log files corresponding to specific mod files.
% You can execute any of the following tasks individually or run them all in the specified order.
% Ensure that you add the `matlab` folder to your PATH in MATLAB or Octave.
% Refer to the instructions in [the manual](https://www.dynare.org/manual/installation-and-configuration.html#configuration).
% The reported runtimes are approximate and depend on the number of cores available and the platform of your machine.
% The reported times are based on:
%   - Apple **MacBook** Pro M2 Max (8 performance, 4 efficiency cores), 64 GB RAM, macOS Sequoia 15.2, MATLAB R2024b Update 3 (24.2.0.2806996) 64-bit (maca64) with Dynare 7-unstable
%   - Lenovo ThinkSystem SR655 **Linux server** (AMD EPYC 7402P 24C 2.8GHz), 6x16GB TruDDR4 3200MHz, Pop!_OS 20.04, MATLAB R2024b Update 3 (24.2.0.2806996) 64-bit (glnxa64) with Dynare 7 compiled from source
% - Note that comparing the Gaussian and CSN model versions may not be meaningful, as the CSN model has more parameters to estimate.
% Bayesian estimation
% The easiest estimation strategy is to use the Slice sampler for the full estimation (see *Task Bayes 1* and *Task Bayes 2*), but we also provide several different ways to estimate the model with the Random-Walk Metropolis-Hastings (RWMH) sampler (as this is standard practice).
% 
% The Slice sampler is usually more efficient (less auto-correlated draws, lower inefficiency factors) and more importantly requires no fine-tuning.
% The downside is a longer runtime, because each draw requires many more posterior function evaluations.
% 
% The RWMH sampler requires a mode-finding step (with a positive definite inverse Hessian at the mode) and also a tuning of the *mh_jscale* to get a desired acceptance ratio.
% On the upside, once it is initialized and fine-tuned, each draw requires only one posterior function evaluation.
% 
% Ultimately, the posterior distributions are nearly identical whether estimated using the Slice sampler or any fine-tuned variant of the RWMH sampler.
% 
% 
clearvars; clc; close all;


%% SETTINGS
% SET VARIABLES TO 1 TO REPLICATE THE RESULTS, 0 TO SKIP THE REPLICATION STEP AND REUSE PREVIOUSLY COMPUTED RESULTS
% -------------------------------------------% | ------------------------ RUNTIME ---------------------- |
% Replication part                           % | maca64_m2max | maca64_m4pro |   glnxa64   |    win64    |
% -------------------------------------------% | ------------ | ------------ | ----------- | ----------- |
REDO_ML_GAUSSIAN                      = 0;   % |   <00h01m    |   <00h01m    |   <00h00m   |   <00h00m   |
REDO_ML_CSN_INITVAL_SEARCH            = 0;   % |    00h17m    |    00h13m    |    00h00m   |    00h00m   |
REDO_ML_CSN                           = 0;   % |    00h13m    |    00h10m    |    00h00m   |    00h00m   |
% -------------------------------------------% | ------------ | ------------ | ----------- | ----------- |
REDO_BAYES_SLICE_LONG_GAUSSIAN        = 0;   % |    00h00m    |    00h20m    |    00h00m   |    00h00m   |
REDO_BAYES_SLICE_LONG_CSN             = 0;   % |    00h00m    |    10h27m    |    00h00m   |    00h00m   |
REDO_BAYES_SLICE_SHORT_GAUSSIAN       = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_SLICE_SHORT_CSN            = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN  = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_RWMH_SLICE_SHORT_CSN       = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_MODE_GAUSSIAN              = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_MODE_CSN                   = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_RWMH_GAUSSIAN              = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
REDO_BAYES_RWMH_CSN                   = 0;   % |    00h00m    |    00h00m    |    00h00m   |    00h00m   |
% -------------------------------------------% | ------------ | ------------ | ----------- | ----------- |


%% COMMON VARIABLES AND PATHS
% please adjust Dynare path and machine identifiers accordingly
MATLAB_VERSION = ['R' version('-release')]; % used for filenames
ARCH           = computer('arch');          % used for filenames
DYNARE_VERSION = '7-unstable';
if strcmp(ARCH,'glnxa64')
    ARCH = 'glnxa64';
    DYNARE_PATH = [getenv('HOME') '/dynare/' DYNARE_VERSION '/matlab'];
elseif strcmp(ARCH,'maca64')
    DYNARE_PATH = ['/Applications/dynare/' DYNARE_VERSION '-arm64/matlab'];
    [~,ARCH] = system('system_profiler SPHardwareDataType | egrep "Chip"'); % distinguish between MacBook Pro M2 Max (maca64_m2max) and Mac mini M4 Pro (maca64_m4pro)
    ARCH = ['maca64_' lower(regexprep(strtrim(erase(ARCH,'Chip: Apple')),'\s+',''))];
elseif strcmp(ARCH,'win64')
    DYNARE_PATH = ['c:/dynare/' DYNARE_VERSION '/matlab'];
end
if isempty(which('dynare')); addpath(DYNARE_PATH); end % add path if not already available


%% MAXIMUM LIKELIHOOD ESTIMATION OF GAUSSIAN MODEL
% Estimate the Gaussian model with maximum likelihood using the PSKF to compute the likelihood
if REDO_ML_GAUSSIAN
    REDO_ML_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_ml_1_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_ML_GAUSSIAN = toc(REDO_ML_GAUSSIAN);
end


%% INITIAL VALUES SEARCH FOR MAXIMUM LIKELIHOOD ESTIMATION OF CSN MODEL
% Search for initial values for the maximum likelihood estimation of the CSN model using the PSKF to compute the likelihood.
% To this end:
% 1. Fix model and stderr parameters to Gaussian estimates from maximum likelihood from *ireland2004_ml_1_gaussian.mod*.
% 2. Create large grid for skew parameters.
% 3. Evaluate likelihood (computed by PSKF) on grid using MATLAB's parallel computing toolbox.
% 4. Use best five skew parameter combinations as initial guess and optimize for each initial guess the likelihood (computed by PSKF) over both stderr and skew parameters.
if REDO_ML_CSN_INITVAL_SEARCH
    REDO_ML_CSN_INITVAL_SEARCH = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_ml_2_csn_initval_search
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s'); rmdir([M_.fname '_*'],'s');
    REDO_ML_CSN_INITVAL_SEARCH = toc(REDO_ML_CSN_INITVAL_SEARCH);
end


%% MAXIMUM LIKELIHOOD ESTIMATION OF CSN MODEL
% Estimate the CSN model with maximum likelihood using the PSKF to compute the likelihood.
% The optimizer is initialized with the best optimized values (highest log-likelihood value) found in *ireland2004_ml_2_csn_initval_search.log*
if REDO_ML_CSN
    REDO_ML_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_ml_3_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_ML_CSN = toc(REDO_ML_CSN);
end


%% BAYESIAN ESTIMATION WITH LONG SLICE SAMPLER OF GAUSSIAN MODEL
% Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
% The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
% No fine-tuning is required, approximately 68 function evaluations per iteration.
if REDO_BAYES_SLICE_LONG_GAUSSIAN
    REDO_BAYES_SLICE_LONG_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_1_slice_long_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_SLICE_LONG_GAUSSIAN = toc(REDO_BAYES_SLICE_LONG_GAUSSIAN);
end


%% BAYESIAN ESTIMATION WITH LONG SLICER SAMPLER OF CSN MODEL
% Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
% The Slice sampler is used to draw from the posterior distribution (8 chains with 5000 draws each, 50% burn-in).
% No fine-tuning is required, approximately 86 function evaluations per iteration.
if REDO_BAYES_SLICE_LONG_CSN
    REDO_BAYES_SLICE_LONG_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_2_slice_long_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_SLICE_LONG_CSN = toc(REDO_BAYES_SLICE_LONG_CSN);
end


%% BAYESIAN ESTIMATION WITH SHORT SLICE SAMPLER OF GAUSSIAN MODEL (TO GET CLOSE TO POSTERIOR MODE)
% Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
% The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
% The posterior mode draw is used to initialize the covariance matrix for the subsequent RWMH sampler and to find the posterior mode below
% No fine-tuning is required, approximately 68 function evaluations per iteration.
if REDO_BAYES_SLICE_SHORT_GAUSSIAN
    REDO_BAYES_SLICE_SHORT_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_3_slice_short_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_SLICE_SHORT_GAUSSIAN = toc(REDO_BAYES_SLICE_SHORT_GAUSSIAN);
end


%% BAYESIAN ESTIMATION WITH SHORT SLICE SAMPLER OF CSN MODEL (TO GET CLOSE TO POSTERIOR MODE)
% Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
% The Slice sampler is used to draw from the posterior distribution (8 chains with 250 draws each, 20% burn-in).
% The posterior mode draw is used to initialize the covariance matrix for the subsequent RWMH sampler and to find the posterior mode below
% No fine-tuning is required, approximately 86 function evaluations per iteration.
if REDO_BAYES_SLICE_SHORT_CSN
    REDO_BAYES_SLICE_SHORT_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_4_slice_short_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_SLICE_SHORT_CSN = toc(REDO_BAYES_SLICE_SHORT_CSN);
end


%% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF GAUSSIAN MODEL, INITIALIZED AT MODE FROM SHORT SLICE
% Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in ireland2004_bayes_3_mode_slice_gaussian.mod.
% No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 30%).
if REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN
    REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_5_rwmh_slice_short_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN = toc(REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN);
end


%% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF CSN MODEL, INITIALIZED AT MODE FROM SHORT SLICE
% Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the short Slice sampler in ireland2004_bayes_4_mode_slice_csn.mod.
% No additional time-consuming mode-finding step is required, but one still needs to tune *mh_jscale* to get a desired acceptance ratio (about 33%).
if REDO_BAYES_RWMH_SLICE_SHORT_CSN
    REDO_BAYES_RWMH_SLICE_SHORT_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_6_rwmh_slice_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_RWMH_SLICE_SHORT_CSN = toc(REDO_BAYES_RWMH_SLICE_SHORT_CSN);
end


% #### Task Bayes 7: Bayesian mode estimation of Gaussian model
% Estimate the mode of the Gaussian model with numerical optimization using the PSKF to compute the likelihood.
% The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 3*.
% First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
% Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
% The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 9*.
% ```matlab
% dynare ireland2004_bayes_7_mode_gaussian
% ```
% Runtime: 35 minutes on MacBook, 33 minutes on Linux server.
% 
% #### Task Bayes 8: Bayesian mode estimation of CSN model
% Estimate the mode of the CSN model with numerical optimization using the PSKF to compute the likelihood.
% The optimization is initialized at the posterior mode found from the previous short Slice sampler in *Task Bayes 4*.
% First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
% Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
% The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler in *Task Bayes 10*.
% ```matlab
% dynare ireland2004_bayes_8_mode_csn
% ```
% Runtime: 6 hours and 11 minutes on MacBook, 7 hours on Linux server.
% 
% #### Task Bayes 9: Bayesian estimation with RWMH sampler of Gaussian model
% Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*) in *Task Bayes 7*.
% Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 7*.
% ```matlab
% dynare ireland2004_bayes_9_rwmh_gaussian parallel conffile=__parallelConf.ini
% ```
% Runtime: 28 minutes on MacBook, 1 hour and 2 minutes on Linux server.
% 
% #### Task Bayes 10: Bayesian estimation with RWMH sampler of CSN model
% Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*) in *Task Bayes 8*.
% Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine in *Task Bayes 8*.
% ```matlab
% dynare ireland2004_bayes_10_rwmh_csn parallel conffile=__parallelConf.ini
% ```
% Runtime: 5 hours and 48 minutes on MacBook, 8 hours and 27 minutes on Linux server.
% 
% ### Impulse response functions
% Compute the impulse response functions of the Gaussian and CSN model using the 16th and 84th percentiles of the estimated distributions.
% ```matlab
% dynare ireland2004_irfs
% ```
% Runtime: 1 minute on MacBook, 1 minute on Linux server.
% 
% ### Simulation of recessions
% Simulates large time series with Gaussian and CSN distributed shocks and computes statistics on recessions.
% ```matlab
% dynare ireland2004_recessions
% ```
% Runtime: 1 minute on MacBook, 1 minute on Linux server.


%% HOUSEKEEPING

fprintf('\n%s\n* RUNTIMES *\n%s\n', repmat('*',1,12), repmat('*',1,12))
if REDO_ML_GAUSSIAN > 0
    fprintf('- Maximum likelihood estimation with Gaussian shocks: %s\n', dynsec2hms(REDO_ML_GAUSSIAN));
end
if REDO_ML_CSN_INITVAL_SEARCH > 0
    fprintf('- Initial values search for maximum likelihood estimation with CSN shocks: %s\n', dynsec2hms(REDO_ML_CSN_INITVAL_SEARCH));
end
if REDO_ML_CSN > 0
    fprintf('- Maximum likelihood estimation with CSN shocks: %s\n', dynsec2hms(REDO_ML_CSN));
end
if REDO_BAYES_SLICE_LONG_GAUSSIAN > 0
    fprintf('- Bayesian estimation using long Slice sampler with Gaussian shocks: %s\n', dynsec2hms(REDO_BAYES_SLICE_LONG_GAUSSIAN));
end
if REDO_BAYES_SLICE_LONG_CSN > 0
    fprintf('- Bayesian estimation using long Slice sampler with CSN shocks: %s\n', dynsec2hms(REDO_BAYES_SLICE_LONG_CSN));
end
if REDO_BAYES_SLICE_SHORT_GAUSSIAN > 0
    fprintf('- Bayesian estimation using short Slice sampler with Gaussian shocks: %s\n', dynsec2hms(REDO_BAYES_SLICE_SHORT_GAUSSIAN));
end
if REDO_BAYES_SLICE_SHORT_CSN > 0
    fprintf('- Bayesian estimation using short Slice sampler with CSN shocks: %s\n', dynsec2hms(REDO_BAYES_SLICE_SHORT_CSN));
end
if REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN > 0
    fprintf('- Bayesian estimation using RWMH sampler initializes at short Slice sampler mode with Gaussian shocks: %s\n', dynsec2hms(REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN));
end
if REDO_BAYES_RWMH_SLICE_SHORT_CSN > 0
    fprintf('- Bayesian estimation using RWMH sampler initializes at short Slice sampler mode with CSN shocks: %s\n', dynsec2hms(REDO_BAYES_RWMH_SLICE_SHORT_CSN));
end
if REDO_BAYES_MODE_GAUSSIAN > 0
    fprintf('- Bayesian mode estimation with Gaussian shocks: %s\n', dynsec2hms(REDO_BAYES_MODE_GAUSSIAN));
end
if REDO_BAYES_MODE_CSN > 0
    fprintf('- Bayesian mode estimation with CSN shocks: %s\n', dynsec2hms(REDO_BAYES_MODE_CSN));
end
if REDO_BAYES_RWMH_GAUSSIAN > 0
    fprintf('- Bayesian mode estimation with Gaussian shocks: %s\n', dynsec2hms(REDO_BAYES_RWMH_GAUSSIAN));
end
if REDO_BAYES_RWMH_CSN > 0
    fprintf('- Bayesian mode estimation with CSN shocks: %s\n', dynsec2hms(REDO_BAYES_RWMH_CSN));
end
rmpath(DYNARE_PATH);
