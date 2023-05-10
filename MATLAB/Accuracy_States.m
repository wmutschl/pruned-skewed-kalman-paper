function [LOSS,SAMPLESTATS] = Accuracy_States(PARAMS,OPT)
% function [LOSS,SAMPLESTATS] = Accuracy_States(PARAMS,OPT)
% -------------------------------------------------------------------------
% Monte-Carlo assessment of accuracy of filtered and smoothed states.
% Replicates results in section 5.1 of the paper
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
% - LOSS          [structure]   results of Monte Carlo analysis with expected losses for filtered and smoothed states, for field names see bottom of the function
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
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

%% get options from OPT structure
mc_replic     = OPT.mc_replic;     % number of Monte-Carlo replications
sample_size   = OPT.sample_size;   % sample size
data_burnin   = OPT.data_burnin;   % periods to discard in simulations
x_burnin      = OPT.x_burnin;      % burn-in phase for state accuracy checks (as we initialize the Kalman filters with a wide prior)
Harvey_factor = OPT.Harvey_factor; % factor for wide prior for Kalman filter initialization, as suggested by Harvey and Phillips(1979)
prune_tol     = OPT.prune_tol;     % vector of tolerance levels for pruning dimensions in skewed Kalman filter
loss_fct      = OPT.loss_fct;      % options for loss functions
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
skip_lik  = true;  % do not compute log-likelihood (as we are only interested in filtered and smoothed states)
skip_loss = false; % do compute loss function

%% Initialize output structures
prune_tol_nbr  = size(prune_tol,2);
SAMPLESTATS    = nan(mc_replic,y_nbr*3);       % stores sample mean, sample std-dev and sample skewness of y
filtL1_gauss   = nan(mc_replic,1);             % stores L1 loss for states filtered by the Gaussian Kalman filter
filtL2_gauss   = nan(mc_replic,1);             % stores L2 loss for states filtered by the Gaussian Kalman filter
filtLa_gauss   = nan(mc_replic,1);             % stores La loss for states filtered by the Gaussian Kalman filter
filtL1_csn     = nan(mc_replic,prune_tol_nbr); % stores L1 loss for states filtered by the Pruned Skewed Kalman filter
filtL2_csn     = nan(mc_replic,prune_tol_nbr); % stores L2 loss for states filtered by the Pruned Skewed Kalman filter
filtLa_csn     = nan(mc_replic,prune_tol_nbr); % stores La loss for states filtered by the Pruned Skewed Kalman filter
smoothL1_gauss = nan(mc_replic,1);             % stores L1 loss for states smoothed by the Gaussian Kalman smoother
smoothL2_gauss = nan(mc_replic,1);             % stores L2 loss for states smoothed by the Gaussian Kalman smoother
smoothLa_gauss = nan(mc_replic,1);             % stores La loss for states smoothed by the Gaussian Kalman smoother
smoothL1_csn   = nan(mc_replic,prune_tol_nbr); % stores L1 loss for states smoothed by the Pruned Skewed Kalman smoother
smoothL2_csn   = nan(mc_replic,prune_tol_nbr); % stores L2 loss for states smoothed by the Pruned Skewed Kalman smoother
smoothLa_csn   = nan(mc_replic,prune_tol_nbr); % stores La loss for states smoothed by the Pruned Skewed Kalman smoother

%% Run Monte-Carlo analysis
% in order to set the seed in the parfor loop, we follow MATLAB's help on "Repeat Random Numbers in parfor-Loops"
sc = parallel.pool.Constant(RandStream('Threefry','Seed',seed_nbr));
fprintf('Seed settings for parpool:\n');
disp(sc.Value);
fprintf('Note: substream will be set inside the parfor loop to the iteration number\n');
disp(repmat('*',1,60));

parfor_progress(mc_replic);
parfor r = 1:mc_replic
    stream = sc.Value; stream.Substream = r; RandStream.setGlobalStream(stream); % way to set the seed in parfor loops, see above
    [y, x] = simulate_csn_state_space(sample_size,data_burnin, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps);
    SAMPLESTATS(r,:) = [mean(y,2)' std(y,0,2)' skewness(y,1,2)']; % compute sample statistics

    % Gaussian Kalman filter and smoother with E_eta (NOT mu_eta) and V_eta (NOT Sigma_eta) as inputs
    [~, xfilt_gauss, xsmooth_gauss] = kalman_gaussian(y, mu_0,Sigma_0, G,R,F, E_eta,V_eta, mu_eps,Sigma_eps, skip_lik,skip_loss,loss_fct);

    filtL1_gauss(r,1) = sum(sum(abs(xfilt_gauss.L1(:,x_burnin:end)-x(:,x_burnin:end)),2));
    filtL2_gauss(r,1) = sum(sum((xfilt_gauss.L2(:,x_burnin:end)-x(:,x_burnin:end)).^2,2));
    idxLa1 = (x>xfilt_gauss.La); idxLa2 = ~idxLa1;
    idxLa1(:,1:x_burnin-1) = 0; idxLa2(:,1:x_burnin-1) = 0;
    filtLa_gauss(r,1) = sum(loss_fct.params.a*abs(xfilt_gauss.La(idxLa1)-x(idxLa1))) + sum(loss_fct.params.b*abs(xfilt_gauss.La(idxLa2)-x(idxLa2)));

    smoothL1_gauss(r,1) = sum(sum(abs(xsmooth_gauss.L1(:,x_burnin:end)-x(:,x_burnin:end)),2));
    smoothL2_gauss(r,1) = sum(sum((xsmooth_gauss.L2(:,x_burnin:end)-x(:,x_burnin:end)).^2,2));
    idxLa1 = (x>xsmooth_gauss.La); idxLa2 = ~idxLa1;
    idxLa1(:,1:x_burnin-1) = 0; idxLa2(:,1:x_burnin-1) = 0;
    smoothLa_gauss(r,1) = sum(loss_fct.params.a*abs(xsmooth_gauss.La(idxLa1)-x(idxLa1))) + sum(loss_fct.params.b*abs(xsmooth_gauss.La(idxLa2)-x(idxLa2)));

    % Pruned Skewed Kalman filter
    filtL1_csn_tmp = nan(1,size(prune_tol,2));
    filtL2_csn_tmp = nan(1,size(prune_tol,2));
    filtLa_csn_tmp = nan(1,size(prune_tol,2));
    for j_prune = 1:size(prune_tol,2)
        [~, xfilt_csn] = kalman_csn(y, mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, cdfmvna_fct,prune_tol(j_prune),skip_lik,skip_loss,loss_fct);
        filtL1_csn_tmp(j_prune) = sum(sum(abs(xfilt_csn.L1(:,x_burnin:end)-x(:,x_burnin:end)),2));
        filtL2_csn_tmp(j_prune) = sum(sum((xfilt_csn.L2(:,x_burnin:end)-x(:,x_burnin:end)).^2,2));
        idxLa1 = (x>xfilt_csn.La); idxLa2 = ~idxLa1;
        idxLa1(:,1:x_burnin-1) = 0; idxLa2(:,1:x_burnin-1) = 0;
        filtLa_csn_tmp(j_prune) = sum(loss_fct.params.a*abs(xfilt_csn.La(idxLa1)-x(idxLa1))) + sum(loss_fct.params.b*abs(xfilt_csn.La(idxLa2)-x(idxLa2)));
    end
    filtL1_csn(r,:) = filtL1_csn_tmp;
    filtL2_csn(r,:) = filtL2_csn_tmp;
    filtLa_csn(r,:) = filtLa_csn_tmp;

    % Pruned Skewed Kalman smoother
    % note that for smoothing we need CSN parameters of unpruned predicted and filtered states,
    % i.e. we re-run kalman_csn with prune_tol=0, but skip loss function evaluation as this is extremely time-consuming
    [~, ~, pred, filt] = kalman_csn(y, mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, cdfmvna_fct,0,skip_lik,true,loss_fct);
    smoothL1_csn_tmp = nan(1,size(prune_tol,2));
    smoothL2_csn_tmp = nan(1,size(prune_tol,2));
    smoothLa_csn_tmp = nan(1,size(prune_tol,2));
    for j_prune = 1:size(prune_tol,2)
        xsmooth_csn = kalman_csn_smoother(y, pred,filt, G,R,F, Gamma_eta,Delta_eta, cdfmvna_fct,prune_tol(j_prune),loss_fct);
        smoothL1_csn_tmp(j_prune) = sum(sum(abs(xsmooth_csn.L1(:,x_burnin:end)-x(:,x_burnin:end)),2));
        smoothL2_csn_tmp(j_prune) = sum(sum((xsmooth_csn.L2(:,x_burnin:end)-x(:,x_burnin:end)).^2,2));
        idxLa1 = (x>xsmooth_csn.La); idxLa2 = ~idxLa1;
        idxLa1(:,1:x_burnin-1) = 0; idxLa2(:,1:x_burnin-1) = 0;
        smoothLa_csn_tmp(j_prune) = sum(loss_fct.params.a*abs(xsmooth_csn.La(idxLa1)-x(idxLa1))) + sum(loss_fct.params.b*abs(xsmooth_csn.La(idxLa2)-x(idxLa2)));
    end
    smoothL1_csn(r,:) = smoothL1_csn_tmp;
    smoothL2_csn(r,:) = smoothL2_csn_tmp;
    smoothLa_csn(r,:) = smoothLa_csn_tmp;
    parfor_progress;
end
parfor_progress(0);

%% Store into output structure
LOSS.gauss.filt.L1   = filtL1_gauss;
LOSS.gauss.filt.L2   = filtL2_gauss;
LOSS.gauss.filt.La   = filtLa_gauss;
LOSS.gauss.smooth.L1 = smoothL1_gauss;
LOSS.gauss.smooth.L2 = smoothL2_gauss;
LOSS.gauss.smooth.La = smoothLa_gauss;
LOSS.csn.filt.L1     = filtL1_csn;
LOSS.csn.filt.L2     = filtL2_csn;
LOSS.csn.filt.La     = filtLa_csn;
LOSS.csn.smooth.L1   = smoothL1_csn;
LOSS.csn.smooth.L2   = smoothL2_csn;
LOSS.csn.smooth.La   = smoothLa_csn;