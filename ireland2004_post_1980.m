% MAXIMUM LIKELIHOOD ESTIMATION ON POST 1980 SAMPLE
% =========================================================================
% Copyright © 2023 Willi Mutschler
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
%% SETTINGS
% housekeeping
clearvars; clearvars -global; close all; clc;
addpath('MATLAB');
addpath('MATLAB/ireland2004');

options_.dsge = 1;
options_.datafile = "post_1980";
options_.optim_opt = optimset('display','final','MaxFunEvals',1000000,'MaxIter',10000,'TolFun',1e-4,'TolX',1e-4); % optimization options (we'll use fminsearch)
options_.optim_opt.penalize_objective = 0; % 1: checks whether bounds are violated in objective function and penalizes it
options_.mode_compute = ["fminsearch" "fminsearchbnd" "fminunc" "cmaes"];
options_.kalman.csn.prune_tol = 1e-4; % pruning threshold in Pruned Skewed Kalman filter
options_.kalman.csn.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
options_.kalman.lik_init = 1; % 1: Gaussian distribution, where initial matrix of variance of the error of forecast is set equal to the unconditional variance of the state variables; 2: Gaussian with wide prior with an initial matrix of variance of the error of forecast diagonal with 10 on the diagonal
options_.parameters.skewness_bounds = [-0.95 0.95]; % bounds for skewness coefficient in estim_params file
options_.parameters.use_stderr_skew = 1;
options_.parameters.transform.OMEGA = true;
options_.parameters.transform.ALPHA_X = true;
options_.parameters.transform.RHO_PI = true;
options_.parameters.transform.RHO_G = true;
options_.parameters.transform.RHO_X = true;
options_.parameters.transform.RHO_A = true;
options_.parameters.transform.RHO_E = true;
options_.parameters.transform.stderr_eta_a = true;
options_.parameters.transform.stderr_eta_e = true;
options_.parameters.transform.stderr_eta_z = true;
options_.parameters.transform.stderr_eta_r = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_a = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_e = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_z = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_r = true;
options_.parameters.transform.skew_eta_a = true;
options_.parameters.transform.skew_eta_e = true;
options_.parameters.transform.skew_eta_z = true;
options_.parameters.transform.skew_eta_r = true;
options_.parameters.transform.diag_Gamma_eta_a = true;
options_.parameters.transform.diag_Gamma_eta_e = true;
options_.parameters.transform.diag_Gamma_eta_z = true;
options_.parameters.transform.diag_Gamma_eta_r = true;

%% MODEL PREPROCESSING
% create script files with analytic dynamic Jacobian, which is later evaluated numerically to compute the policy function
% inspired by and very similar to what Dynare's preprocessor does, but in MATLAB
M_ = ireland2004_preprocessing;

% common calibration, values taken from Ireland (2004)
M_.params(ismember(M_.param_names,'BETA'))     = 0.99;
M_.params(ismember(M_.param_names,'PSI'))      = 0.1;
M_.params(ismember(M_.param_names,'ALPHA_PI')) = 0.0001; % fix to very small value as estimate is virtually 0 and this parameter makes even Gaussian ML tricky

%% POST 1980 DATA
load data/ireland2004_gpr.dat; % original dataset
ireland2004_gpr = array2timetable(ireland2004_gpr,...
       'RowTimes',datetime('1948-Q2','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'):calquarters(1):datetime('2003-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'),...
       'VariableNames',{'ghat',...  % output growth (needs to be demeaned), based on quarterly changes in seasonally adjusted real GDP, converted to per capita terms by dividing by the civilian noninsitutional population aged 16 and over
                        'pihat',... % inflation (needs to be demeaned), based on quarterly changes in seasonally adjusted GDP deflator
                        'rhat'} ... % nominal interest rate (needs to be demeaned), based on quarterly averages of daily readings on the three-months US Treasury bill rate
       );
ireland2004_gpr = ireland2004_gpr(find(ireland2004_gpr.Time==datetime('1980-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ')):end,:);
datamat = [];
for jvarobs=1:length(M_.varobs)
    datamat = [datamat (ireland2004_gpr.(M_.varobs(jvarobs))-mean(ireland2004_gpr.(M_.varobs(jvarobs))))]; %demean data
end
clear ireland2004_gpr jvarobs


%% MAXIMUM LIKELIHOOD ESTIMATION WITH SOPHISTICATED SEARCH FOR INITIAL VALUES
dsge_maximum_likelihood_estimation_csn;

%% HOUSEKEEPING
rmpath('MATLAB');
rmpath('MATLAB/ireland2004');



%In Stage 1, the first step is to fix the model parameters and the variance of the shocks to the estimates obtained in Stage 0. After this, a grid is created for the skewness of each shock, and for each combination of variance and skewness, the corresponding combinations of covariance and correlation are recovered. Next, the negative log-likelihood of each value on the grid is computed. Then, a chosen number of the best combinations are used as the initial value, and the shock parameters are optimized over using the Pruned Skewed Kalman filter. This process enables the estimation of the best-fitting parameters for the model and the shocks.

