% Replicates:
% - Table 4: "Distribution  of Parameter Estimates"
% of the paper "Pruned Skewed Kalman Filter and Smoother: With Application
% to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% -------------------------------------------------------------------------
% The results are saved both as MAT as well as log files in the "results" folder
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

SAMPLE_SIZES = 100:50:200; % loop over different sample lengths
%SAMPLE_SIZES = []; % leave empty to only re-create tables then only tables will be created

% Common options
OPT.mc_replic     = 1000;                       % number of Monte Carlo replications
OPT.prune_tol     = [1e-6 1e-2];                % pruning thresholds
OPT.data_burnin   = 100;                        % periods to discard in simulations
OPT.Harvey_factor = 10;                         % factor for wide prior for Kalman filter initialization, as suggested by Harvey and Phillips(1979)
OPT.seed_nbr      = 653476096;                  % randi([0 2^32],1,1); % draw random seed
OPT.cdfmvna_fct   = "logmvncdf_ME";             % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
OPT.optim_name    = "fminsearch";               % name of optimization function
OPT.optim_options = optimset('Display', 'off'); % options for optimization algorithm, optimset('TolFun',1e-7,'TolX',1e-7,'MaxIter', 200000, 'MaxFunEvals',200000, 'Display', 'off');


OPT.param_names_fixed = ["G"          "$G$";
                         "F"          "$F$";
                         "R"          "$R$";
                         "mu_eps"     "$\mu_\epsilon$";
                         "Sigma_eps"  "$\Sigma_\epsilon$";
                         "nu_eta"     "$\nu_\eta$";
                         "Delta_eta"  "$\Delta_\eta$";
                        ];
OPT.param_names_estim = ["mu_eta"            "$\mu_\eta$";
                         "diaglogSigma_eta"  "$\Sigma_\eta$";
                         "diagGamma_eta"     "$\Gamma_\eta$";                   
                        ];

% DGP (3)
OPT.parameter_set = 'DGP3';
PARAMS.G = [0.9969    0.1256   -0.4803
           -0.8221    0.0386    0.6687
            0.5605    0.6397   -0.4333];
PARAMS.F = eye(3);
PARAMS.R = eye(3);
PARAMS.mu_eps = zeros(3,1);
PARAMS.Sigma_eps = 1e-4*eye(3);
PARAMS.mu_eta = [0.3; -0.1; 0.2];
PARAMS.Sigma_eta = diag([0.8^2;0.6^2;0.7^2]);
PARAMS.nu_eta = zeros(3,1);
PARAMS.Gamma_eta = diag([5;0;-6]);
PARAMS.Delta_eta = eye(3);

% Parallel environment
if ~isempty(SAMPLE_SIZES)
    OPT.poolobj = gcp('nocreate');
    if isempty(OPT.poolobj)
        OPT.physical_cores_nbr=feature('numcores');
        parpool(OPT.physical_cores_nbr);
    else
        OPT.physical_cores_nbr = OPT.poolobj.NumWorkers;
    end
end
OPT.computer_arch = computer('arch');

if ~exist('results', 'dir')
    mkdir('results');
end

%% Monte-Carlo study for univariate DGP (1)

% loop over sample sizes and perform Monte Carlo analysis
for sample_size = SAMPLE_SIZES
    OPT.sample_size = sample_size;
    fprintf('%s with T=%d and R=%d\n',OPT.optim_name,OPT.sample_size,OPT.mc_replic);
    tic
    [RES.EstimParams,RES.Fvals,RES.exitflags,RES.SampleStats] = Parameter_Estimation(PARAMS,OPT); % simulates data and estimates the parameters of the CSN distribution
    RES.time_Computing_Time = toc;
    OPT.filename = sprintf('results_parameter_estimation_%s_T%d_R%d_%s',OPT.optim_name,OPT.sample_size,OPT.mc_replic,OPT.computer_arch);
    save(['results/' OPT.filename],'OPT','PARAMS','RES');
end



%% display all results and create tables for paper
ONLY_LATEX = false;
%OPT0.computer_arch = 'glnxa64';
%OPT0.computer_arch = 'maca64';
%OPT0.computer_arch = 'maci64';
%OPT0.computer_arch = 'win64';
DisplayResults_MakeLatexTables_ParameterEstimation(['results_parameter_estimation_fminsearch_T100_R1000_'  OPT.computer_arch],ONLY_LATEX);
%DisplayResults_MakeLatexTables_ParameterEstimation(['results_parameter_estimation_fminsearch_T150_R10_'  OPT.computer_arch],ONLY_LATEX);
%DisplayResults_MakeLatexTables_ParameterEstimation(['results_parameter_estimation_fminsearch_T200_R10_'  OPT.computer_arch],ONLY_LATEX);

%% Clean up
rmpath('MATLAB');