% =========================================================================
% Replication script for the Dynare-based estimation of the Ireland (2004)
% model with Gaussian and skew-normal distributed shocks.
% This requires Dynare 7.0 (or later); if installed to default locations
% it will be picked up automatically, otherwise adjust DYNARE_PATH below.
% All results (log files, tables, figures) are organized in the results
% folder. Set the corresponding flag to 1 to recompute results. Tasks can
% be run individually or sequentially in the order given below.
%
% MAXIMUM LIKELIHOOD:
% - REDO_ML_GAUSSIAN           - ML estimation of Gaussian model using PSKF
% - REDO_ML_CSN_INITVAL_SEARCH - grid search + optimization over skew and
%                                stderr parameters to find initial values
%                                for CSN ML estimation (uses parpool)
% - REDO_ML_CSN                - ML estimation of CSN model using PSKF,
%                                initialized at best values from grid search;
%                                also computes likelihood ratio test
%
% BAYESIAN SLICE SAMPLER (no tuning required)
% - REDO_BAYES_SLICE_LONG_GAUSSIAN  - long Slice sampler, 8x5000 draws, 50% burn-in (Gaussian)
% - REDO_BAYES_SLICE_LONG_CSN       - long Slice sampler, 8x5000 draws, 50% burn-in (skew-normal)
% - REDO_BAYES_SLICE_SHORT_GAUSSIAN - short Slice sampler, 8x250 draws, 20% burn-in (Gaussian)
% - REDO_BAYES_SLICE_SHORT_CSN      - short Slice sampler, 8x250 draws, 20% burn-in (skew-normal)
% - REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN - RWMH, 8x250000 draws, initalized at mode found by short Slice sampler (Gaussian)
% - REDO_BAYES_RWMH_SLICE_SHORT_CSN      - RWMH, 8x250000 draws, initalized at mode found by short Slice sampler (skew-normal)
%
% BAYESIAN (standard RWMH workflow: mode-finding then RWMH):
% REDO_BAYES_MODE_GAUSSIAN  - numerical mode-finding (mode_compute=8 then 6) (Gaussian)
% REDO_BAYES_MODE_CSN       - numerical mode-finding (mode_compute=8 then 6) (skew-normal)
% REDO_BAYES_RWMH_GAUSSIAN  - RWMH, 8x250000 draws, mode from mode_compute=6 (Gaussian)
% REDO_BAYES_RWMH_CSN       - RWMH, 8x250000 draws, mode from mode_compute=6 (skew-normal)
%
% ADDITIONAL ANALYSES:
% REDO_RECESSIONS                      - simulate 500000 periods, compute recession statistics
% REDO_COMPARISON_WITH_PARTICLE_FILTER - compare PSKF log-likelihood with particle filter
% REDO_IRFS                            - IRFs at 16th/84th percentiles of ML shock distributions
%
% Runtimes are approximate and depend on the number of cores and platform.
% The reported times below are based on the following machines and software versions:
% - [maca64_m2max] Apple MacBook Pro M2 Max (8P+4E cores), 64 GB, macOS 26.2, MATLAB R2025b, Dynare 7.0
% - [maca64_m4pro] Apple Mac mini M4 Pro (10P+4E cores), 64 GB, macOS 26.2, MATLAB R2025b, Dynare 7.0
% - [glnxa64]      Lenovo ThinkSystem SR655 (AMD EPYC 7402P 24C), 96 GB, Ubuntu 25.10, MATLAB R2025b, Dynare 7.0
% - [win64]        HP Elite Tower 800 G9 (Intel Core i7-12700 12C 2.1GHz), 64 GB, Windows 11 25H2, MATLAB R2025b, Dynare 7.0
% The paper uses the [maca64_m4pro] results.

clearvars; clc; close all;


%% SETTINGS
%========================================================================================================|
% Replication part                           % |                         RUNTIMES                        |
% Set to 1 to recompute results              % | maca64_m2max | maca64_m4pro |   glnxa64   |    win64    |
%========================================================================================================|
REDO_ML_GAUSSIAN                      = 0;   % |   <00h01m    |   <00h01m    |   <00h01m   |   <00h01m   |
REDO_ML_CSN_INITVAL_SEARCH            = 0;   % |    00h18m    |    00h14m    |    00h30m   |    00h23m   |
REDO_ML_CSN                           = 0;   % |    00h13m    |    00h10m    |    00h26m   |    00h09m   |
%--------------------------------------------------------------------------------------------------------|
REDO_BAYES_SLICE_LONG_GAUSSIAN        = 0;   % |    00h27m    |    00h21m    |    00h56m   |    00h53m   |
REDO_BAYES_SLICE_LONG_CSN             = 0;   % |    14h33m    |    10h26m    |    17h46m   |    18h29m   |
REDO_BAYES_SLICE_SHORT_GAUSSIAN       = 0;   % |    00h02m    |    00h01m    |    00h03m   |    00h03m   |
REDO_BAYES_SLICE_SHORT_CSN            = 0;   % |    00h46m    |    00h33m    |    00h57m   |    00h58m   |
REDO_BAYES_RWMH_SLICE_SHORT_GAUSSIAN  = 0;   % |    00h22m    |    00h15m    |    00h43m   |    00h41m   |
REDO_BAYES_RWMH_SLICE_SHORT_CSN       = 0;   % |    08h16m    |    04h31m    |    08h08m   |    08h26m   |
REDO_BAYES_MODE_GAUSSIAN              = 0;   % |    00h10m    |    00h06m    |    00h22m   |    00h13m   |
REDO_BAYES_MODE_CSN                   = 0;   % |    03h13m    |    02h08m    |    07h04m   |    03h45m   |
REDO_BAYES_RWMH_GAUSSIAN              = 0;   % |    00h20m    |    00h15m    |    00h44m   |    00h41m   |
REDO_BAYES_RWMH_CSN                   = 0;   % |    07h17m    |    04h26m    |    08h03m   |    08h33m   |
%--------------------------------------------------------------------------------------------------------|
REDO_RECESSIONS                       = 0;   % |   <00h01m    |   <00h01m    |   <00h01m   |   <00h01m   |
REDO_COMPARISON_WITH_PARTICLE_FILTER  = 0;   % |    00h51m    |    00h35m    |    02h17m   |    00h58m   |
REDO_IRFS                             = 0;   % |   <00h01m    |   <00h01m    |   <00h01m   |   <00h01m   |
%--------------------------------------------------------------------------------------------------------|


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
    dynare ireland2004_bayes_6_rwmh_slice_short_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_RWMH_SLICE_SHORT_CSN = toc(REDO_BAYES_RWMH_SLICE_SHORT_CSN);
end


%% BAYESIAN MODE ESTIMATION OF GAUSSIAN MODEL
% Estimate the mode of the Gaussian model with numerical optimization using the PSKF to compute the likelihood.
% The optimization is initialized at the posterior mode found from the previous short Slice sampler.
% First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
% Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
% The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler below.
if REDO_BAYES_MODE_GAUSSIAN
    REDO_BAYES_MODE_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_7_mode_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_MODE_GAUSSIAN = toc(REDO_BAYES_MODE_GAUSSIAN);
end


%% BAYESIAN MODE ESTIMATION OF CSN MODEL
% Estimate the mode of the CSN model with numerical optimization using the PSKF to compute the likelihood.
% The optimization is initialized at the posterior mode found from the previous short Slice sampler.
% First, a Nelder-Mead simplex-based optimization routine (*mode_compute=8*) is used (fast).
% Second, a Monte-Carlo based optimization routine (*mode_compute=6*) is run that guarantees to yield a positive definite inverse Hessian at the mode (very time-consuming).
% The mode with the highest posterior value is reported in the paper and the mode and covariance matrix from *mode_compute=6* is used to initialize the RWMH sampler below.
if REDO_BAYES_MODE_CSN
    REDO_BAYES_MODE_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_8_mode_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_MODE_CSN = toc(REDO_BAYES_MODE_CSN);
end


%% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF GAUSSIAN MODEL
% Estimate the Gaussian model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*).
% Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine.
if REDO_BAYES_RWMH_GAUSSIAN
    REDO_BAYES_RWMH_GAUSSIAN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_9_rwmh_gaussian
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_RWMH_GAUSSIAN = toc(REDO_BAYES_RWMH_GAUSSIAN);
end


%% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF CSN MODEL
% Estimate the CSN model with Bayesian methods using the PSKF to compute the likelihood.
% The RWMH sampler is used to draw from the posterior distribution (8 chains with 250000 draws each, 50% burn-in).
% The sampler is initialized at the posterior mode and covariance matrix found by the Monte-Carlo optimization (*mode_compute=6*).
% Tuning of *mh_jscale* required to get a desired acceptance ratio (about 30%), which is a by-product of running the Monte-Carlo based optimization routine.
if REDO_BAYES_RWMH_CSN
    REDO_BAYES_RWMH_CSN = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_bayes_10_rwmh_csn
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_BAYES_RWMH_CSN = toc(REDO_BAYES_RWMH_CSN);
end


%% STATISTICS ON RECESSIONS
% Simulate 500000 time periods of model with (i) Gaussian and (ii) CSN distributed shocks
% and count number of total, severe and mild recessions.  Recessions are defined by a period
% where output growth falls below -0.5% per quarter and stays negative for at least two quarters.
% Severe and mild recessions are determined by the top and bottom terciles of the implied
% distribution of peak-to-trough output losses.
if REDO_RECESSIONS
    REDO_RECESSIONS = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_recessions
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_RECESSIONS = toc(REDO_RECESSIONS);
end


%% COMPARISON WITH PARTICLE FILTER
if REDO_COMPARISON_WITH_PARTICLE_FILTER
    REDO_COMPARISON_WITH_PARTICLE_FILTER = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_particle_filter_comparison
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_COMPARISON_WITH_PARTICLE_FILTER = toc(REDO_COMPARISON_WITH_PARTICLE_FILTER);
end


%% IMPULSE RESPONSE FUNCTIONS
% Compute impulse response functions of the Gaussian and CSN model variants
% using the 16th and 84th percentiles of the estimated ML shock distributions.
% Results are stored in a tidy-format CSV file for plotting in R.
if REDO_IRFS
    REDO_IRFS = tic;
    clearvars -except DYNARE_PATH ARCH MATLAB_VERSION REDO_*; clc; close all;
    dynare ireland2004_irfs
    % housekeeping
    pause(1); fclose('all'); movefile([M_.fname '.log'], target_logfile);
    rmdir(['+' M_.fname],'s'); rmdir(M_.fname,'s');
    REDO_IRFS = toc(REDO_IRFS);
end


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
if REDO_RECESSIONS > 0
    fprintf('- Simulations and statistics on recessions: %s\n', dynsec2hms(REDO_RECESSIONS));
end
if REDO_COMPARISON_WITH_PARTICLE_FILTER > 0
    fprintf('- Comparison with particle filter: %s\n', dynsec2hms(REDO_COMPARISON_WITH_PARTICLE_FILTER));
end
if REDO_IRFS > 0
    fprintf('- Impulse response functions: %s\n', dynsec2hms(REDO_IRFS));
end
rmpath(DYNARE_PATH);
