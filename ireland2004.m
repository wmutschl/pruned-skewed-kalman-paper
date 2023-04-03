% Maximum Likelihood estimation of the DSGE model and data used in
% Ireland (2004) - "Technolog Shocks in The New Keynesian Model"
% The model is estimated first with the Gaussian Kalman filter to replicate
% the original paper; and then with the Pruned Skewed Kalman filter to allow
% for skewed innovations to preference, cost-push, productivity and monetary policy shocks.
% The model is estimated both on the full sample as well as the post-1980 sample.
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

%% housekeeping
clearvars; clearvars -global; close all; clc;
addpath('MATLAB');
addpath('MATLAB/ireland2004');

ONLY_LATEX = false;

if ~ONLY_LATEX
%% MODEL PREPROCESSING
% create script files with analytic dynamic Jacobian, which is later evaluated numerically to compute the policy function
% inspired by and very similar to what Dynare's preprocessor does, but in MATLAB
M_ = ireland2004_preprocessing;

% common calibration, values taken from Ireland (2004)
M_.params(ismember(M_.param_names,'BETA'))     = 0.99;
M_.params(ismember(M_.param_names,'PSI'))      = 0.1;

%% DATA
load data/ireland2004_gpr.dat; % original dataset
ireland2004_gpr = array2timetable(ireland2004_gpr,...
       'RowTimes',datetime('1948-Q2','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'):calquarters(1):datetime('2003-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'),...
       'VariableNames',{'ghat',...  % output growth (needs to be demeaned), based on quarterly changes in seasonally adjusted real GDP, converted to per capita terms by dividing by the civilian noninsitutional population aged 16 and over
                        'pihat',... % inflation (needs to be demeaned), based on quarterly changes in seasonally adjusted GDP deflator
                        'rhat'} ... % nominal interest rate (needs to be demeaned), based on quarterly averages of daily readings on the three-months US Treasury bill rate
       );
ireland2004_post1980 = ireland2004_gpr(find(ireland2004_gpr.Time==datetime('1980-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ')):end,:);
datamat = [];
for jvarobs=1:length(M_.varobs)
    datamat = [datamat (ireland2004_post1980.(M_.varobs(jvarobs))-mean(ireland2004_post1980.(M_.varobs(jvarobs))))]; %demean data on actual sample
end
clear ireland2004_gpr ireland2004_post1980 jvarobs

%% OPTIONS
options_.computer_arch = computer('arch');
options_.dsge = 1; % 1: we are interested in the structural model parameters and not the state-space parameters of the linearized solution
options_.optim_opt = optimset('display','final','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-4,'TolX',1e-4); % optimization options
options_.optim_opt.names = ["cmaes" "cmaes_dsge" "fminsearch" "fminsearchbnd" "fmincon" "fminunc" ]; % names of optimizer that will be used in parallel ("simulannealbnd" and "sa_resampling" are very time-consuming, "cmaes" and "cmaes_dsge" are mildly time-consuming, fminsearch and fminsearchbnd can be mildly time-consuming, "fmincon" and "fminunc" are fast (but not as good) )
options_.optim_opt.penalize_objective = 0; % 1: checks whether bounds are violated in objective function and penalizes likelihood (useful for optimizers that don't support parameter bounds)
options_.kalman.csn.prune_tol = 1e-2; % pruning threshold in Pruned Skewed Kalman filter
options_.kalman.csn.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
options_.kalman.csn.prune_algorithm = "correlations"; % "correlations": prunes according to correlation coefficient; "max_dim": prunes according to correlation coefficient+keeps a maximum of 15 dimensions
options_.kalman.lik_init = 1; % 1: stationary distribution, i.e. Gaussian distribution, where initial matrix of variance of the error of forecast is set equal to the unconditional variance of the state variables;
                              % 2: wide prior, i.e. Gaussian with an initial matrix of variance of the error of forecast diagonal with 10 on the diagonal
                              % 3: run Gaussian Kalman filter (without computing the likelihood) and take smoothed Sigma matrix
options_.parameters.skewness_bounds = [-0.95 0.95]; % bounds for skewness coefficient used in estim_params file;
                                                    % note that the absolute value of the theoretical skewness coefficient of the CSN distribution is bounded by 0.995
                                                    % inference at the bound is problematic due to so-called pile ups; in our experience the likelihood becomes extremely flat with respect to Gamma_eta when close to the bound
options_.parameters.use_stderr_skew = 1;            % 1: optimization algorithm optimizes over skewness coefficient and standard error of CSN shocks instead of Sigma_eta and Gamma_eta,
                                                    % 0: optimization algorithm optimizes over diag(Gamma_eta) and sqrt(diag(Sigma_eta)) instead of skewness coefficient and standard error of CSN shocks
                                                    % note that this has implications for the computation of standard errors
% transform bounded parameters into unbounded ones for the optimization algorithm using the log-odds-transformation
options_.parameters.fix.ALPHA_X            = false;
options_.parameters.fix.ALPHA_PI           = false;
options_.parameters.transform.OMEGA        = true;
options_.parameters.transform.RHO_PI       = false;
options_.parameters.transform.RHO_G        = false;
options_.parameters.transform.RHO_X        = false;
options_.parameters.transform.RHO_A        = false;
options_.parameters.transform.RHO_E        = false;
options_.parameters.transform.skew_eta_a   = true;
options_.parameters.transform.skew_eta_e   = true;
options_.parameters.transform.skew_eta_z   = true;
options_.parameters.transform.skew_eta_r   = true;
options_.parameters.transform.stderr_eta_a = true;
options_.parameters.transform.stderr_eta_e = true;
options_.parameters.transform.stderr_eta_z = true;
options_.parameters.transform.stderr_eta_r = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_a = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_e = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_z = true;
options_.parameters.transform.sqrt_diag_Sigma_eta_r = true;
options_.parameters.transform.diag_Gamma_eta_a = true;
options_.parameters.transform.diag_Gamma_eta_e = true;
options_.parameters.transform.diag_Gamma_eta_z = true;
options_.parameters.transform.diag_Gamma_eta_r = true;
if isfield(options_.parameters.fix,'ALPHA_X') && (options_.parameters.fix.ALPHA_X==1)
    M_.params(ismember(M_.param_names,'ALPHA_X')) = 0;
else
    options_.parameters.transform.ALPHA_X  = true;
end
if isfield(options_.parameters.fix,'ALPHA_PI') && (options_.parameters.fix.ALPHA_PI==1)
    M_.params(ismember(M_.param_names,'ALPHA_PI')) = 0;
else    
    options_.parameters.transform.ALPHA_PI = true;
end

%% ESTIMATION ON POST 1980 SAMPLE
options_.filename = sprintf('results_ireland2004_stderrskew%d_KalmanInit%d_FixAx%d_FixAp%d_R%s_%s',options_.parameters.use_stderr_skew,options_.kalman.lik_init,options_.parameters.fix.ALPHA_X,options_.parameters.fix.ALPHA_PI,version('-release'),options_.computer_arch);
%options_.kalman.csn.initval_search = 1; % 1: try to find good initial values first from Gaussian estimation, then grid on csn shock parameters;
%options_.optim_opt.names = ["fminsearch" "fminsearchbnd"]; % names of optimizer that will be used in parallel ("simulannealbnd" and "sa_resampling" are very time-consuming, "cmaes" and "cmaes_dsge" are mildly time-consuming, fminsearch and fminsearchbnd can be mildly time-consuming, "fmincon" and "fminunc" are fast (but not as good) )
options_.kalman.csn.initval_search = 0; % 0: use initial values provided in estim_params file
[oo_, M_] = dsge_maximum_likelihood_estimation_csn(M_,options_,datamat);
end

if ONLY_LATEX
    load('results/results_ireland2004_stderrskew1_KalmanInit1_FixAx0_FixAp0_R2023a_maci64.mat');
else
    load("results/" + options_.filename);
end

%% display all results and create tables for paper
[~,optim_col] = sort(oo_.csn.neg_log_likelihood); optim_col=optim_col(1);
xparams_csn = oo_.csn.xparams(:,optim_col);
xparams_csn_stderr = oo_.csn.xstderr(:,optim_col);
xparams_gauss = [nan(4,1); oo_.gauss.xparams(:,optim_col)];
xparams_gauss_stderr = [nan(4,1); oo_.gauss.xstderr(:,optim_col)];

tbl_latex_model = string.empty; tbl_latex_shock = string.empty;
for jp=1:length(xparams_csn)
    switch erase(estim_params_csn.names(jp),'transformed_')
        case "skew_eta_a"
            latex_name = "$skew(\eta_a)$";
        case "skew_eta_e"
            latex_name = "$skew(\eta_e)$";
        case "skew_eta_z"
            latex_name = "$skew(\eta_z)$";
        case "skew_eta_r"
            latex_name = "$skew(\eta_r)$";
        case "stderr_eta_a"
            latex_name = "$stderr(\eta_a)$";
        case "stderr_eta_e"
            latex_name = "$stderr(\eta_e)$";
        case "stderr_eta_z"
            latex_name = "$stderr(\eta_z)$";
        case "stderr_eta_r"
            latex_name = "$stderr(\eta_r)$";
        case "OMEGA"
            latex_name = "$\omega$";
        case "ALPHA_X"
            latex_name = "$\alpha_x$";
        case "ALPHA_PI"
            latex_name = "$\alpha_\pi$";
        case "RHO_PI"
            latex_name = "$\rho_\pi$";
        case "RHO_G"
            latex_name = "$\rho_g$";
        case "RHO_X"
            latex_name = "$\rho_x$";
        case "RHO_A"
            latex_name = "$\rho_a$";
        case "RHO_E"
            latex_name = "$\rho_e$";
    end
    if contains(latex_name,"stderr") || contains(latex_name,"skew")
        tbl_latex_shock = [tbl_latex_shock;
                           sprintf('%s & $\\underset{(%.4f)}{%.4f}$ & $\\underset{(%.4f)}{%.4f}$'...
                          ,latex_name ...
                          ,xparams_gauss_stderr(jp),xparams_gauss(jp) ...
                          ,xparams_csn_stderr(jp),xparams_csn(jp) ...
                          )];
    else
        tbl_latex_model = [tbl_latex_model;
                           sprintf('%s & $\\underset{(%.4f)}{%.4f}$ & $\\underset{(%.4f)}{%.4f}$'...
                           ,latex_name ...
                           ,xparams_gauss_stderr(jp),xparams_gauss(jp) ...
                           ,xparams_csn_stderr(jp),xparams_csn(jp) ...
                           )];
    end
end
tbl_latex = tbl_latex_model + "         &$\quad$&$\quad$&         " + tbl_latex_shock + " \\";
tbl_latex = [tbl_latex; "\midrule"; sprintf("\\multicolumn{6}{l}{Value of maximized Log-Likelihood function:} & %.2f & %.2f",-1*oo_.gauss.neg_log_likelihood(optim_col),-1*oo_.csn.neg_log_likelihood(optim_col))];
fprintf('\n\nLatex Table Entries:\n\n')
disp(tbl_latex);

%% HOUSEKEEPING
rmpath('MATLAB');
rmpath('MATLAB/ireland2004');