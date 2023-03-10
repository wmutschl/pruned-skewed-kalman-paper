function [ESTIMPARAMS,FVALS,EXITFLAGS,SAMPLESTATS] = Parameter_Estimation(PARAMS,OPT)
% function [ESTIMPARAMS,FVALS,EXITFLAGS,SAMPLESTATS] = Parameter_Estimation(PARAMS,OPT)
% -------------------------------------------------------------------------
% Monte-Carlo assessment of accuracy of estimating parameters of the CSN
% distributed innovations in a linear state-space system with Maximum Likelihood.
% Replicates results in section 5.3 of the paper
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
% - PARAMS        [structure]   results of Monte Carlo analysis with estimated parameters, see below for fields
% - FVALS         [structure]   results of Monte Carlo analysis with log-likelihood values for different variants of the Pruned Skewed Kalman filter and Gaussian Kalman filter, see below for fields
% - EXITFLAGS     [structure]   results of Monte Carlo analysis with exitflags for different variants of the Pruned Skewed Kalman filter and Gaussian Kalman filter, see below for fields
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
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================

% get options
mc_replic     = OPT.mc_replic;     % number of Monte-Carlo replications
sample_size   = OPT.sample_size;   % sample size of replications
data_burnin   = OPT.data_burnin;   % periods to discard in simulations
prune_tol     = OPT.prune_tol;     % vector of tolerance levels for pruning dimensions in skewed Kalman filter
seed_nbr      = OPT.seed_nbr;      % number to fix the seed for random numbers below
cdfmvna_fct   = OPT.cdfmvna_fct;   % which function to use to evaluate high-dimensional Gaussian log cdf
optim_name    = OPT.optim_name;    % name of optimization algorithm
optim_options = OPT.optim_options; % options for optimization algorithm

% get "true" parameters
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

% compute expectation and variance of eta as we use these in the Gaussian Kalman filter for expectation and covariance
if isfield(PARAMS,'lambda_eta') && (norm(Delta_eta + Gamma_eta'*Sigma_eta*Gamma_eta - eye(size(nu_eta,1)),'Inf')<1e-7)
    % use re-parametrization to compute unconditional first two moments
    E_eta = mu_eta + (2/sqrt(2*pi)*lambda_eta*sqrt_Sigma_eta)*ones(size(nu_eta,1),1);
    V_eta = Sigma_eta*(1-2/pi*lambda_eta^2);
else % use general functions to compute unconditional first two moments
    E_eta = csnMean(mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, cdfmvna_fct);
    V_eta = csnVar(Sigma_eta, Gamma_eta, nu_eta, Delta_eta,cdfmvna_fct);
end

%% Initialize estimated parameters at true values
xparam0_csn   = [mu_eta(:); log(diag(Sigma_eta)); diag(Gamma_eta)]; % do log transform on variances
xparam0_gauss = [E_eta(:);  log(diag(V_eta))                     ]; % do log transform on variances

%% Initialize output structures
prune_tol_nbr = size(prune_tol,2); [OPT.y_nbr,OPT.x_nbr] = size(F); OPT.eta_nbr = size(E_eta,1); OPT.eps_nbr = size(mu_eps,1); % store dimensions
SAMPLESTATS   = nan(mc_replic,OPT.y_nbr*3);                       % stores sample mean, sample std-dev and sample skewness of y
xparam1_gauss = nan(length(xparam0_gauss),mc_replic);             % stores parameters estimated by the Gaussian Kalman filter
xparam1_csn   = nan(length(xparam0_csn),prune_tol_nbr,mc_replic); % stores parameters estimated by the Pruned Skewed Kalman filter
fval_gauss    = nan(1,mc_replic);                                 % stores function values of log-likelihood estimated by the Gaussian Kalman filter
fval_csn      = nan(prune_tol_nbr,mc_replic);                     % stores function values of log-likelihood estimated by the Pruned Skewed Kalman filter
EXITFLAGS     = true(mc_replic,1);                                % stores indicators if for the simulated sample both filters run without problems (=true) or with problems (=false)

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
    data_not_skewed_enough = true; % make sure data has high skewness
    while data_not_skewed_enough
        y = simulate_csn_state_space(sample_size,data_burnin, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps);
        if sum(abs(skewness(y,1,2)) > 0.6) == sum(diag(Gamma_eta)~=0)
            data_not_skewed_enough = false;
        end
    end
    SAMPLESTATS(r,:) = [mean(y,2)' std(y,0,2)' skewness(y,1,2)']; % compute sample statistics
    
    try
        % Gaussian Kalman filter
        if strcmp(optim_name,'fminsearch')
            f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,"gaussian");
            [xparam1_gauss(:,r),fval_gauss(r),exitflag] = fminsearch(f,xparam0_gauss,optim_options);
            if ~exitflag; EXITFLAGS(r)=0; end
        elseif strcmp(optim_name,'fmincon')
            f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,"gaussian");
            [xparam1_gauss(:,r),fval_gauss(r),exitflag] = fmincon(f,xparam0_gauss,[],[],[],[],[],[],[],optim_options);
            if ~exitflag; EXITFLAGS(r)=0; end
        elseif strcmp(optim_name,'fminunc')
            f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,"gaussian");
            [xparam1_gauss(:,r),fval_gauss(r),exitflag] = fminunc(f,xparam0_gauss,optim_options);
            if ~exitflag; EXITFLAGS(r)=0; end
        elseif strcmp(optim_name,'simulannealbnd')
            f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,"gaussian");
            [xparam1_gauss(:,r),fval_gauss(r),exitflag] = simulannealbnd(f,xparam0_gauss,[],[],optim_options);
            if ~exitflag; EXITFLAGS(r)=0; end
        end
    catch
        warning('Something wrong for run r=%d',r);
        EXITFLAGS(r)=0;
    end

    try
        % CSN Kalman filter
        xparam1_csn_tmp = nan(length(xparam0_csn),size(prune_tol,2));
        fval_csn_tmp = nan(size(prune_tol,2),1);
        for jprune = 1:size(prune_tol,2)
            if strcmp(optim_name,'fminsearch')
                f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,{"pruned_skewed",prune_tol(jprune)});
                [xparam1_csn_tmp(:,jprune),fval_csn_tmp(jprune),exitflag] = fminsearch(f,xparam0_csn,optim_options);
                if ~exitflag; EXITFLAGS(r)=0; end
            elseif strcmp(optim_name,'fmincon')
                f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,{"pruned_skewed",prune_tol(jprune)});
                [xparam1_csn_tmp(:,jprune),fval_csn_tmp(jprune),exitflag] = fmincon(f,xparam0_csn,[],[],[],[],OPT.lb,OPT.ub,[],optim_options);
                if ~exitflag; EXITFLAGS(r)=0; end
            elseif strcmp(optim_name,'fminunc')
                f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,{"pruned_skewed",prune_tol(jprune)});
                [xparam1_csn_tmp(:,jprune),fval_csn_tmp(jprune),exitflag] = fminunc(f,xparam0_csn,optim_options);
                if ~exitflag; EXITFLAGS(r)=0; end
            elseif strcmp(optim_name,'simulannealbnd')
                f = @(xparam) -1*log_likelihood_objective_function(xparam,y,OPT,PARAMS,{"pruned_skewed",prune_tol(jprune)});
                [xparam1_csn_tmp(:,jprune),fval_csn_tmp(jprune),exitflag] = simulannealbnd(f,xparam0_csn,OPT.lb,OPT.ub,optim_options);
                if ~exitflag; EXITFLAGS(r)=0; end            
            end
        end
        xparam1_csn(:,:,r) = xparam1_csn_tmp;
        fval_csn(:,r) = fval_csn_tmp;
    catch
        warning('Something wrong for run r=%d',r);
        EXITFLAGS(r)=0;
    end
    parfor_progress;
end
parfor_progress(0);

%% Store into output structure
idx_gauss = 0; idx_csn = 0;
idx_diag_eps = reshape(find(repmat(logical(eye(OPT.eps_nbr)),[1,1,mc_replic])),OPT.eps_nbr,mc_replic);
idx_diag_eta = reshape(find(repmat(logical(eye(OPT.eta_nbr)),[1,1,mc_replic])),OPT.eta_nbr,mc_replic);
for jvar = OPT.param_names_estim(:,1)'
    if jvar == "G"
        ESTIMPARAMS.gauss.G = reshape(xparam1_gauss(idx_gauss+(1:OPT.x_nbr^2),:),OPT.x_nbr,OPT.x_nbr,mc_replic);
        for jprune = 1:size(prune_tol,2)
            ESTIMPARAMS.csn.G(:,:,:,jprune) = reshape(squeeze(xparam1_csn(idx_csn+(1:OPT.x_nbr^2),jprune,:)),OPT.x_nbr,OPT.x_nbr,mc_replic);
        end
        idx_gauss=idx_gauss+OPT.x_nbr^2;  idx_csn=idx_csn+OPT.x_nbr^2;
    end
    if jvar == "F"
        ESTIMPARAMS.gauss.F = reshape(xparam1_gauss(idx_gauss+(1:OPT.y_nbr*OPT.x_nbr),:),OPT.y_nbr,OPT.x_nbr,mc_replic);
        for jprune = 1:size(prune_tol,2)
            ESTIMPARAMS.csn.F(:,:,:,jprune) = reshape(squeeze(xparam1_csn(idx_csn+(1:OPT.y_nbr*OPT.x_nbr),jprune,:)),OPT.y_nbr,OPT.x_nbr,mc_replic);
        end
        idx_gauss=idx_gauss+OPT.y_nbr*OPT.x_nbr; idx_csn=idx_csn+OPT.y_nbr*OPT.x_nbr;
    end
    if jvar == "R"
        ESTIMPARAMS.gauss.R = reshape(xparam1_gauss(idx_gauss+(1:OPT.x_nbr*OPT.eta_nbr),:),OPT.x_nbr,OPT.eta_nbr,mc_replic);
        for jprune = 1:size(prune_tol,2)
            ESTIMPARAMS.csn.R(:,:,:,jprune) = reshape(squeeze(xparam1_csn(idx_csn+(1:OPT.x_nbr*OPT.eta_nbr),jprune,:)),OPT.x_nbr,OPT.eta_nbr,mc_replic);
        end
        idx_gauss=idx_gauss+OPT.x_nbr*OPT.eta_nbr; idx_csn=idx_csn+OPT.x_nbr*OPT.eta_nbr;
    end
    if jvar == "mu_eps"
        ESTIMPARAMS.gauss.mu_eps = xparam1_gauss(idx_gauss+(1:OPT.eps_nbr),:);
        for jprune = 1:size(prune_tol,2)
            ESTIMPARAMS.csn.mu_eps(:,:,jprune) = squeeze(xparam1_csn(idx_csn+(1:OPT.eps_nbr),jprune,:));
        end
        idx_gauss=idx_gauss+OPT.eps_nbr; idx_csn=idx_csn+OPT.eps_nbr;
    end
    if jvar == "diaglogSigma_eps" % undo log transform and make full matrix
        ESTIMPARAMS.gauss.Sigma_eps = zeros(OPT.eps_nbr,OPT.eps_nbr,mc_replic); ESTIMPARAMS.gauss.Sigma_eps(idx_diag_eps) = exp(xparam1_gauss(idx_gauss+(1:OPT.eps_nbr),:));
        ESTIMPARAMS.csn.Sigma_eps = zeros(OPT.eps_nbr,OPT.eps_nbr,mc_replic,size(prune_tol,2));
        for jprune = 1:size(prune_tol,2)
            tmp_Sigma_eps = ESTIMPARAMS.csn.Sigma_eps(:,:,:,jprune);
            tmp_Sigma_eps(idx_diag_eps) = squeeze(exp(xparam1_csn(idx_csn+(1:OPT.eps_nbr),jprune,:)));
            ESTIMPARAMS.csn.Sigma_eps(:,:,:,jprune) = tmp_Sigma_eps;            
        end
        idx_gauss=idx_gauss+OPT.eps_nbr; idx_csn=idx_csn+OPT.eps_nbr;
    end
    if jvar == "mu_eta"
        ESTIMPARAMS.gauss.mu_eta = xparam1_gauss(idx_gauss+(1:OPT.eta_nbr),:);
        for jprune = 1:size(prune_tol,2)
            ESTIMPARAMS.csn.mu_eta(:,:,jprune) = squeeze(xparam1_csn(idx_csn+(1:OPT.eta_nbr),jprune,:));
        end
        idx_gauss=idx_gauss+OPT.eta_nbr; idx_csn=idx_csn+OPT.eta_nbr;
    end
    if jvar == "diaglogSigma_eta" % undo log transform and make full matrix
        ESTIMPARAMS.gauss.Sigma_eta = zeros(OPT.eta_nbr,OPT.eta_nbr,mc_replic); ESTIMPARAMS.gauss.Sigma_eta(idx_diag_eta) = exp(xparam1_gauss(idx_gauss+(1:OPT.eta_nbr),:));
        ESTIMPARAMS.csn.Sigma_eta = zeros(OPT.eta_nbr,OPT.eta_nbr,mc_replic,size(prune_tol,2));
        for jprune = 1:size(prune_tol,2)
            tmp_Sigma_eta = ESTIMPARAMS.csn.Sigma_eta(:,:,:,jprune);
            tmp_Sigma_eta(idx_diag_eta) = squeeze(exp(xparam1_csn(idx_csn+(1:OPT.eta_nbr),jprune,:)));
            ESTIMPARAMS.csn.Sigma_eta(:,:,:,jprune) = tmp_Sigma_eta;            
        end
        idx_gauss=idx_gauss+OPT.eta_nbr; idx_csn=idx_csn+OPT.eta_nbr;
    end
    if jvar == "diagGamma_eta" % make full matrix        
        ESTIMPARAMS.gauss.Gamma_eta = [];        
        ESTIMPARAMS.csn.Gamma_eta = zeros(OPT.eta_nbr,OPT.eta_nbr,mc_replic,size(prune_tol,2));
        for jprune = 1:size(prune_tol,2)
            tmp_Gamma_eta = ESTIMPARAMS.csn.Gamma_eta(:,:,:,jprune);
            tmp_Gamma_eta(idx_diag_eta) = squeeze(xparam1_csn(idx_csn+(1:OPT.eta_nbr),jprune,:));
            ESTIMPARAMS.csn.Gamma_eta(:,:,:,jprune) = tmp_Gamma_eta;                
        end
        idx_csn=idx_csn+OPT.eta_nbr;
    end
end
FVALS.gauss       = fval_gauss;
FVALS.csn         = fval_csn;