function [negative_log_likelihood, exit_flag, M_] = negative_log_likelihood(xparams, bounds, datamat, estim_params_, options_, M_)
% function [negative_log_likelihood, exit_flag] = negative_log_likelihood(xparam, datamat, bounds, options_, M_)
% -------------------------------------------------------------------------
% computes the negative log-likelihood of a state-space model with csn distributed innovations eta and normally distributed noise eps:
%   x(t) = G*x(t-1) + R*eta(t)   [state transition equation]
%   y(t) = F*x(t)   +   eps(t)   [observation equation]
%   eta(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta) [innovations, shocks]
%   eps(t) ~ N(0,Sigma_eps) [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% Assumptions:
%   - mu_eta is endogenously determined to ensure that E[eta]=0
%   - due to identifiability, we normalize nu_eta=0, Delta_eta=I
% -------------------------------------------------------------------------
% INPUTS
% - xparam          [estim_param_nbr by 1]   values of estimated parameters
% - bounds          [estim_param_nbr by 2]   values with lower and upper bounds for estimated parameters
% - datamat         [varobs_nbr by nobs)     matrix with data
% - params_         [structure]              information on calibrated parameters
% - estim_params_   [structure]              information on estimated parameters, inspired by Dynare's estimated_params block
% - options_        [structure]              options
% - M_              [structure]              model information (also includes M_.Sigma_eps, M_.Sigma_eta and M_.Gamma_eta)
% -------------------------------------------------------------------------
% OUTPUTS
% - negative_log_likelihood   [double]   value of negative log-likelihood function
% - exit_flag                 [boolean]  1: no errors occured
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

%% INITIALIZATIONS
exit_flag = 1;
large_number = Inf;
kf_variant = "gaussian";

if ~isfield(options_,'first_obs')
    options_.first_obs = 1;
end
if ~isfield(options_,'nobs')
    options_.nobs = size(datamat,1);
end
datamat = datamat(options_.first_obs:options_.nobs,:);

%% CHECK BOUNDS
if options_.optim_opt.penalize_objective
    if any(xparams<bounds(:,1)) || any(xparams>bounds(:,2))
        negative_log_likelihood = large_number;
        exit_flag = 0;
        return
    end
end

%% UPDATE PARAMETERS
if options_.dsge
    M_ = set_params_dsge(xparams,estim_params_,options_,M_);
end
% check whether we need Pruned Skewed Kalman filter
if any(diag(M_.Gamma_eta))
    kf_variant = "pruned_skewed";
end

%% CHECK POSITIVE DEFINITENESS OF COVARIANCE MATRICES
if isdiag(M_.Cov_eta)
    if any(diag(M_.Cov_eta)<0)
        negative_log_likelihood = large_number; exit_flag = 0; warning('Cov_eta not positive definite\n');
        return
    end
elseif ~ispd(M_.Cov_eta)
    negative_log_likelihood = large_number; exit_flag = 0; warning('Cov_eta not positive definite\n');
    return
end
if isdiag(M_.Cov_eps)
    if any(diag(M_.Cov_eps)<0)
        negative_log_likelihood = large_number; exit_flag = 0; warning('Cov_eps not positive definite\n');
        return
    end
elseif ~ispd(M_.Cov_eps)
    negative_log_likelihood = large_number; exit_flag = 0; warning('Cov_eps not positive definite\n');
    return
end

%% CHECK THEORETICAL BOUND ON SKEWNESS COEFFICIENT
for jp=1:M_.exo_nbr
    if abs(csnSkewness_univariate(M_.Sigma_eta(jp,jp),M_.Gamma_eta(jp,jp))) > abs((sqrt(2)*(pi-4))/(pi-2)^(3/2))
        negative_log_likelihood = large_number; exit_flag = 0; warning('Skewness coefficient is out of bounds: %.4f\n',csnSkewness_univariate(M_.Sigma_eta(jp,jp),M_.Gamma_eta(jp,jp)));
        return
    end
end

%% FIRST-ORDER PERTURBATION SOLUTION AROUND DETERMINISTIC STEADY STATE
if options_.dsge
    [M_.G,M_.R,error_indicator] = dsge_perturbation_solution_order_1(M_);
    if error_indicator; negative_log_likelihood = large_number; exit_flag = 0; return; end % Return, if no solution is found    
end

%% COMPUTE NEGATIVE LOG-LIKELIHOOD USING KALMAN FILTER

if options_.kalman.lik_init == 1
    % initialize the state vector at the stationary distribution
    mu_0 = zeros(M_.endo_nbr,1); % note that y are the model variables in deviation from steady-state, so the mean is zero by definition
    %Sigma_0 = reshape( inv(eye(MODEL.endo_nbr*MODEL.endo_nbr) - kron(SOL.gx,SOL.gx))*reshape(SOL.gu*MODEL.COV_eta*SOL.gu',MODEL.endo_nbr*MODEL.endo_nbr,1) ,MODEL.endo_nbr,MODEL.endo_nbr); %analytical, but slow
    Sigma_0 = dlyapdoubling(M_.G,M_.R*M_.Cov_eta*M_.R'); % very fast and numerically accurate
end

if kf_variant == "gaussian"
    negative_log_likelihood = -1*kalman_gaussian(datamat', mu_0,Sigma_0, M_.G,M_.R,M_.F, M_.mu_eta,M_.Sigma_eta, M_.mu_eps,M_.Sigma_eps, false, true);
elseif kf_variant == "pruned_skewed"
    Gamma_0 = zeros(M_.endo_nbr,M_.endo_nbr); % initialize at Gaussian distribution
    nu_0    = zeros(M_.endo_nbr,1); % normalization
    Delta_0 = eye(M_.endo_nbr); % normalization
    %log_likelihood1 = kalman_csn(DATAMAT', mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, SOL.gx,SOL.gu,MODEL.F, MODEL.mu_eta,MODEL.Sigma_eta,MODEL.Gamma_eta,MODEL.nu_eta,MODEL.Delta_eta, MODEL.mu_eps,MODEL.Sigma_eps, OPT.cdfmvna_fct,OPT.prune_tol,false,true);
    if options_.dsge
        negative_log_likelihood = -1*kalman_csn_dsge(...
            M_.G,M_.R,M_.F, ...
            mu_0, Sigma_0, Gamma_0, nu_0, Delta_0, ...
            M_.mu_eps, M_.Sigma_eps, ...
            M_.mu_eta, M_.Sigma_eta, M_.Gamma_eta, M_.nu_eta, M_.Delta_eta, ...
            datamat, true, options_.kalman.csn.prune_tol, options_.kalman.csn.cdfmvna_fct);
    end
end

if isnan(negative_log_likelihood) || isinf(negative_log_likelihood) || ~isreal(negative_log_likelihood)
    negative_log_likelihood = large_number;
    exit_flag = 0;
end

end % main function end