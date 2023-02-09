% Replicates:
% - Table 1: "Expected losses for filtered states"
% - Table 2: "Expected losses for smoothed states"
% of the paper "Pruned Skewed Kalman Filter and Smoother: With Application
% to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% -------------------------------------------------------------------------
% The results are saved both as MAT as well as log files in the "results" folder
% -------------------------------------------------------------------------
% Note that we use the following re-parametrization for the skewness parameters of the CSN distributed innovations eta:
%   - Gamma_eta = lambda_eta*Sigma_eta^(-1/2)
%   - nu_eta = zeros(nx,1)
%   - Delta_eta = (1-lambda^2)*eye(nx)
% where lambda=]-1;1[ (together with Sigma_eta) scales the skewness of the distribution and nu_eta=0 for lack of identification.
% This re-parametrization enables one to compute the theoretical moments easily in closed-form (alternatively you can make use of csn_mean and csn_var,
% which rely on numerical approximations of the Jacobian and Hessian of certain Gaussian CDF functions)
% =========================================================================
% Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
% =========================================================================

%% Settings

% Housekeeping
clearvars; clearvars -global; clc; close all;
addpath('MATLAB');

REDO = 0; % 0: only create tables, 1: redo DGP1, 2: redo DGP2, 3: redo DGP1 and DGP2

% Common options
SAMPLE_SIZES           = [40 80 110];        % loop over different sample lengths
OPT0.mc_replic         = 2400;               % number of Monte Carlo replications
OPT0.prune_tol         = [0 1e-6 1e-4 1e-2]; % pruning thresholds
OPT0.data_burnin       = 100;                % periods to discard in simulations
OPT0.x_burnin          = 20;                 % burn-in phase for state accuracy checks (as we initialize the Kalman filters with a wide prior)
OPT0.Harvey_factor     = 10;                 % factor for wide prior for Kalman filter initialization, as suggested by Harvey and Phillips(1979)
OPT0.loss_fct.params.a = 1;                  % parameter a of asymmetric loss function: a*|xtilde_t - x_t| for x_t > xtilde_t
OPT0.loss_fct.params.b = 4;                  % parameter b of asymmetric loss function: b*|xtilde_t - x_t| for x_t <= xtilde_t
OPT0.cdfmvna_fct       = "logmvncdf_ME";     % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"

% Parallel environment
if REDO > 0
    OPT0.poolobj = gcp('nocreate');
    if isempty(OPT0.poolobj)
        OPT0.physical_cores_nbr=feature('numcores');
        parpool(OPT0.physical_cores_nbr);
    else
        OPT0.physical_cores_nbr = OPT0.poolobj.NumWorkers;
    end
end
OPT0.computer_arch = computer('arch');

if ~exist('results', 'dir')
    mkdir('results');
end

%% Monte-Carlo study for univariate DGP (1)
if REDO == 1 || REDO == 3
    clear OPT PARAMS RES sample_size

    OPT = OPT0;
    OPT.parameter_set    = 'DGP1';
    OPT.seed_nbr         = 4112460544; % randi([0 2^32],1,1); % draw random seed
    OPT.loss_fct.type.L1 = true; % use median for absolute loss
    OPT.loss_fct.type.L2 = true; % use mean for squared loss
    OPT.loss_fct.type.La = true; % use quantile for asymmetric loss

    PARAMS.G          = 0.8;   % matrix multiplying previous states in state transition equation
    PARAMS.R          = 1;     % matrix multiplying innovations in state transition equation
    PARAMS.F          = 10;    % matrix multiplying current states in measurement equation
    PARAMS.mu_eps     = 1;     % expectation of Gaussian noise
    PARAMS.Sigma_eps  = 1e-2;  % variance of Gaussian noise
    PARAMS.mu_eta     = 0.3;   % location parameter of CSN distributed innovation eta(t)
    PARAMS.Sigma_eta  = 0.64;  % scale parameter of CSN distributed innovation eta(t)
    PARAMS.lambda_eta = -0.89; % auxiliary skewness parameter, between -1 and 1, that determines Gamma_eta and Delta_eta below
    % compute skewness parameters Gamma_eta, nu_eta, Delta_eta from re-parametrization values
    [U,S,V] = svd(PARAMS.Sigma_eta); PARAMS.sqrt_Sigma_eta = U*diag(diag(S).^(1/2))*V'; PARAMS.inv_sqrt_Sigma_eta = U*diag(diag(S).^(-1/2))*V'; % compute Sigma_eta^(-1/2)
    PARAMS.Gamma_eta = PARAMS.lambda_eta*PARAMS.inv_sqrt_Sigma_eta; % first skewness matrix of CSN distributed innovations eta(t), making use of re-parametrization
    PARAMS.nu_eta = zeros(size(PARAMS.Gamma_eta,1),1); % second skewness matrix of CSN distributed innovations eta(t), fix at 0 due to re-parametrization and lack of identification
    PARAMS.Delta_eta = (1-PARAMS.lambda_eta^2)*eye(size(PARAMS.Gamma_eta,1));
    if (norm(PARAMS.Delta_eta + PARAMS.Gamma_eta'*PARAMS.Sigma_eta*PARAMS.Gamma_eta - eye(size(PARAMS.nu_eta,1)),'Inf')>1e-7)
        error('Reparametrization failed')
    end
    clear U S V;

    % loop over sample sizes and perform Monte Carlo analysis
    for sample_size = SAMPLE_SIZES
        OPT.sample_size = sample_size;
        fprintf('%s with T=%d and R=%d\n',OPT.parameter_set,OPT.sample_size,OPT.mc_replic);
        tic
        [RES.Loss,RES.SampleStats] = Accuracy_States(PARAMS,OPT); % simulates data and performs both filtering and smoothing
        RES.time_Accuracy_States = toc;
        OPT.filename = sprintf('results_state_estimation_%s_T%d_R%d_%s',OPT.parameter_set,OPT.sample_size,OPT.mc_replic,OPT.computer_arch);
        save(['results/' OPT.filename],'OPT','PARAMS','RES');
    end
end

%% Monte-Carlo study for multivariate DGP (2)
if REDO == 2 || REDO == 3
    clear OPT PARAMS RES sample_size

    OPT = OPT0;
    OPT.parameter_set    = 'DGP2';
    OPT.seed_nbr         = 1936362161; %randi([0 2^32],1,1); % draw random seed
    OPT.loss_fct.type.L1 = false; % don't use median for absolute loss as in multivariate case this is not obvious
    OPT.loss_fct.type.L2 = true;  % use mean for squared loss
    OPT.loss_fct.type.La = false; % don't use quantile for asymmetric loss as in multivariate case this is not obvious

    rng(OPT.seed_nbr);
    y_nbr = 3; x_nbr = 4;
    is_stationary = 0;
    while is_stationary == 0
        PARAMS.G = 80*randn(x_nbr); % matrix multiplying previous states in state transition equation
        [v,d] = eig(PARAMS.G);
        PARAMS.G = real(v*diag(-1+2.*rand(length(PARAMS.G),1))/v);
        PARAMS.G = round(PARAMS.G,4);
        eigG = eig(PARAMS.G);
        if all(imag(eigG)==0) && all(abs(eigG)<1)
            is_stationary=1;
        end
    end
    PARAMS.R = eye(x_nbr);                    % matrix multiplying innovations in state transition equation
    PARAMS.F = round(randn(y_nbr, x_nbr),4);  % matrix multiplying current states in measurement equation
    PARAMS.mu_eps = round(randn(y_nbr,1),4);  % expectation of Gaussian noise
    PARAMS.Sigma_eps = round(tril(0.5*randn(y_nbr)),2); PARAMS.Sigma_eps = PARAMS.Sigma_eps*PARAMS.Sigma_eps'./1e6; % variance of Gaussian noise
    PARAMS.mu_eta = round(randn(x_nbr,1),4);  % location parameter of CSN distributed innovation eta(t)
    PARAMS.Sigma_eta = round(tril(5*randn(x_nbr)),2); PARAMS.Sigma_eta = PARAMS.Sigma_eta*PARAMS.Sigma_eta';        % scale parameter of CSN distributed innovation eta(t)
    while max(abs(PARAMS.Sigma_eta(:))) > 10
        PARAMS.Sigma_eta = PARAMS.Sigma_eta./10;
    end
    PARAMS.lambda_eta = 0.89; % auxiliary skewness parameter, between -1 and 1, that determines Gamma_eta and Delta_eta below
    % compute skewness parameters Gamma_eta, nu_eta, Delta_eta from re-parametrization values
    [U,S,V] = svd(PARAMS.Sigma_eta); PARAMS.sqrt_Sigma_eta = U*diag(diag(S).^(1/2))*V'; PARAMS.inv_sqrt_Sigma_eta = U*diag(diag(S).^(-1/2))*V'; % compute Sigma_eta^(-1/2)
    PARAMS.Gamma_eta = PARAMS.lambda_eta*PARAMS.inv_sqrt_Sigma_eta; % first skewness matrix of CSN distributed innovations eta(t), making use of re-parametrization
    PARAMS.nu_eta = zeros(size(PARAMS.Gamma_eta,1),1); % second skewness matrix of CSN distributed innovations eta(t), fix at 0 due to re-parametrization and lack of identification
    PARAMS.Delta_eta = (1-PARAMS.lambda_eta^2)*eye(size(PARAMS.Gamma_eta,1));
    if (norm(PARAMS.Delta_eta + PARAMS.Gamma_eta'*PARAMS.Sigma_eta*PARAMS.Gamma_eta - eye(size(PARAMS.nu_eta,1)),'Inf')>1e-7)
        error('Reparametrization failed')
    end
    clear v d eigG is_stationary U S V x_nbr y_nbr

    % loop over sample sizes and perform Monte Carlo analysis
    for sample_size = SAMPLE_SIZES
        OPT.prune_tol = OPT0.prune_tol;
        if sample_size > 40 && any(ismember(OPT.prune_tol,0))
            OPT.prune_tol(ismember(OPT.prune_tol,0)) = []; % remove the uncut version in multivariate case
            warning('In the multivariate case, we do not compute the unpruned version for sample sizes larger than 40. I remove the 0 in prune_tol.')
        end
        OPT.sample_size = sample_size;
        fprintf('%s with T=%d and R=%d\n',OPT.parameter_set,OPT.sample_size,OPT.mc_replic);
        tic        
        [RES.Loss,RES.SampleStats] = Accuracy_States(PARAMS,OPT); % simulates data and performs both filtering and smoothing
        RES.time_Accuracy_States = toc;
        OPT.filename = sprintf('results_state_estimation_%s_T%d_R%d_%s',OPT.parameter_set,OPT.sample_size,OPT.mc_replic,OPT.computer_arch);
        save(['results/' OPT.filename],'OPT','PARAMS','RES');
    end
end

%% display all results and create tables for paper
%OPT0.computer_arch = 'glnxa64';
%OPT0.computer_arch = 'maca64';
%OPT0.computer_arch = 'maci64';
%OPT0.computer_arch = 'win64';
DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP1_T40_R2400_'  OPT0.computer_arch]);
DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP1_T80_R2400_'  OPT0.computer_arch]);
DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP1_T110_R2400_' OPT0.computer_arch]);

DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP2_T40_R2400_'  OPT0.computer_arch]);
DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP2_T80_R2400_'  OPT0.computer_arch]);
DisplayResults_MakeLatexTables_ExpectedLosses(['results_state_estimation_DGP2_T110_R2400_' OPT0.computer_arch]);


%% Clean up
rmpath('MATLAB');