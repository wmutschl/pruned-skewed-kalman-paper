function [TIME,SAMPLESTATS] = Computing_Time_Log_Likelihood(PARAMS,OPT)
% function [TIME,SAMPLESTATS] = Computing_Time_Log_Likelihood(PARAMS,OPT)
% -------------------------------------------------------------------------
% Monte-Carlo assessment of computing time to evaluate the log-likelihood
% function of both the Gaussian as well as the Pruned Skewed Kalman Filter.
% Replicates results in section 5.2 of the paper
% "Pruned Skewed Kalman Filter and Smoother: With Application to the Yield Curve"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
%
% Notation for linear state-space system with csn distributed innovations and normally distributed noise:
%   x(t) = G*x(t-1) + R*eta(t)   [state transition equation]
%   y(t) = F*x(t)   + eps(t)     [observation equation]
%   eta(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta) [innovations, shocks]
%   eps(t) ~ N(mu_eps,Sigma_eps) [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% -------------------------------------------------------------------------
% INPUTS
% - PARAMS        [structure]   parameters of linear state-space system, see below for fields
% - OPT           [structure]   options for Monte Carlo analysis, see below for fields
% -------------------------------------------------------------------------
% OUTPUTS
% - TIME          [structure]   results of Monte Carlo analysis with measured time in seconds, for field names see bottom of the function
% - SAMPLESTATS   [structure]   descriptive statistics (sample mean, sample standard deviation and sample skewness coefficient) for simulated data used in Monte Carlo
% =========================================================================
% Copyright (C) 2022-2023 Willi Mutschler
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

%% get options from OPT structure
mc_replic     = OPT.mc_replic;     % number of Monte-Carlo replications
sample_size   = OPT.sample_size;   % sample size
data_burnin   = OPT.data_burnin;   % periods to discard in simulations
Harvey_factor = OPT.Harvey_factor; % factor for wide prior for Kalman filter initialization, as suggested by Harvey and Phillips(1979)
prune_tol     = OPT.prune_tol;     % vector of tolerance levels for cutting dimensions in skewed Kalman filter
seed_nbr      = OPT.seed_nbr;      % number to fix the seed for random numbers below
cdfmvna_fct   = OPT.cdfmvna_fct;   % which function to use to evaluate high-dimensional Gaussian log cdf

%% get parameters of linear state-space system from PARAMS structure
G          = PARAMS.G;
F          = PARAMS.F;
R          = PARAMS.R;
mu_eps     = PARAMS.mu_eps;
Sigma_eps  = PARAMS.Sigma_eps;
mu_eta     = PARAMS.mu_eta;
Sigma_eta  = PARAMS.Sigma_eta;
Gamma_eta  = PARAMS.Gamma_eta;
nu_eta     = PARAMS.nu_eta;
Delta_eta  = PARAMS.Delta_eta;
if isfield(PARAMS,'lambda_eta')
    lambda_eta     = PARAMS.lambda_eta;
    sqrt_Sigma_eta = PARAMS.sqrt_Sigma_eta;
end

% compute expectation and variance of eta as we use these in the Gaussian Kalman filter for expectation and covariance
if isfield(PARAMS,'lambda_eta') && (norm(Delta_eta + Gamma_eta'*Sigma_eta*Gamma_eta - eye(size(nu_eta,1)),'Inf')<1e-7)
    % use re-parametrization to compute unconditional first two moments
    E_eta = mu_eta + (2/sqrt(2*pi)*lambda_eta*sqrt_Sigma_eta)*ones(size(nu_eta,1),1);
    V_eta = Sigma_eta*(1-2/pi*lambda_eta^2);
else % use general functions to compute unconditional first two moments
    E_eta = csnMean(mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, cdfmvna_fct);
    V_eta = csnVar(Sigma_eta, Gamma_eta, nu_eta, Delta_eta,cdfmvna_fct);
end

%% Initialization of Kalman filters
% initialize Kalman filter with a wide Normal prior: x_0 ~ CSN(0,Harvey_factor*eye(nx),0,0,eye(nx))
[y_nbr,x_nbr] = size(F);
mu_0      = zeros(x_nbr,1);
Sigma_0   = Harvey_factor*eye(x_nbr);
Gamma_0   = zeros(x_nbr,x_nbr);
nu_0      = zeros(x_nbr,1);
Delta_0   = eye(x_nbr);
skip_lik  = false; % do compute log-likelihood
skip_loss = true;  % do compute loss function (as we are only interested in computing the log-likelihood)

%% Initialize output structures
prune_tol_nbr     = size(prune_tol,2);
SAMPLESTATS       = nan(mc_replic,y_nbr*3);       % stores sample mean, sample std-dev and sample skewness of y
Time_LogLik_gauss = nan(mc_replic,1);             % stores the time to compute the log-likelihood function with the Gaussian Kalman filter
Time_LogLik_csn   = nan(mc_replic,prune_tol_nbr); % stores the time to compute the log-likelihood function with the Pruned Skewed Kalman filter

%% Run Monte-Carlo analysis
% in order to set the seed in the parfor loop, we follow MATLAB's help on "Repeat Random Numbers in parfor-Loops"
sc = parallel.pool.Constant(RandStream('Threefry','Seed',seed_nbr));
fprintf('Seed settings for parpool:\n');
disp(sc.Value);
fprintf('Note: substream will be set inside the parfor loop to the iteration number\n');
disp(repmat('*',1,60));

%parfor_progress(mc_replic);
parfor r = 1:mc_replic
    stream = sc.Value; stream.Substream = r; RandStream.setGlobalStream(stream); % way to set the seed in parfor loops, see above
    y = simulate_csn_state_space(sample_size,data_burnin, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps);
    SAMPLESTATS(r,:) = [mean(y,2)' std(y,0,2)' skewness(y,1,2)']; % compute sample statistics
    
    % evaluate log-likelihood using Gaussian Kalman filter with E_eta and V_eta as inputs
    tic
    log_lik_gauss = kalman_gaussian(y, mu_0,Sigma_0, G,R,F, E_eta,V_eta, mu_eps,Sigma_eps, skip_lik,skip_loss);
    if isinf(log_lik_gauss) || isnan(log_lik_gauss)
        error('kalman_gaussian: log-likelihood is inf or nan for r=%d',r);
    end
    Time_LogLik_gauss(r,1) = toc;

    % evaluate log-likelihood using Pruned Skewed Kalman filter
    Time_Lik_csn_tmp = nan(1,prune_tol_nbr);
    for j_prune = 1:prune_tol_nbr
        tic
        log_lik_csn = kalman_csn(y, mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, cdfmvna_fct,prune_tol(j_prune),skip_lik,skip_loss);
        if isinf(log_lik_csn) || isnan(log_lik_csn)
            error('kalman_csn: log-likelihood is inf or nan for r=%d',r);
        end
        Time_Lik_csn_tmp(1,j_prune) = toc;
    end
    Time_LogLik_csn(r,:) = Time_Lik_csn_tmp;
    %parfor_progress;
end
%parfor_progress(0);

%% Store into output structure
TIME.gauss = Time_LogLik_gauss;
TIME.csn   = Time_LogLik_csn;