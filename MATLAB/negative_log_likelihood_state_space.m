function [negative_log_likelihood,exit_flag] = negative_log_likelihood_state_space(xparam,datamat,bounds, params_, options, M_)
% function [negative_log_likelihood, exit_flag] = negative_log_likelihood_state_space(xparam, datamat, bounds, params_, options_, M_)
% -------------------------------------------------------------------------
% computes the negative log-likelihood of a DSGE model with csn shocks
% solved with perturbation at first order. Note that the solution takes the
% form of a linear state-space system with csn distributed innovations eta 
% and normally distributed noise eps:
%   x(t) = gx*x(t-1) + gu*eta(t)   [state transition equation]
%   y(t) = F*x(t)    + eps(t)      [observation equation]
%   eta(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta) [innovations, shocks]
%   eps(t) ~ N(0,Sigma_eps) [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% Assumptions:
%   - all elements in eta are independent
%   - mu_eta is endogenously determined to ensure that E[eta]=0
%   - due to identifiability, we normalize nu_eta=0, Delta_eta=I
% -------------------------------------------------------------------------
% INPUTS
% - xparam     [estim_param_nbr by 1]   values of estimated parameters
% - datamat    [varobs_nbr by nobs)     matrix with data
% - bounds     [estim_param_nbr by 2]   values with lower and upper bounds for estimated parameters
% - params_    [structure]              information on calibrated parameters
% - options_   [structure]              options
% - M_         [structure]              model information (also includes M_.Sigma_eps, M_.Sigma_eta and M_.Gamma_eta)
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
%% initializations
exit_flag = 1;
verylarge = Inf;

%% update parameters with xparam1
idx = 0;
for jvar = options_.param_names_estim(:,1)'
    if jvar == "G"
        G = reshape(xparam(idx+(1:options_.x_nbr^2)),options_.x_nbr,options_.x_nbr);
        idx = idx + options_.x_nbr^2;
    end
    if jvar == "F"
        F = reshape(xparam(idx+(1:options_.y_nbr*options_.x_nbr)),options_.y_nbr,options_.x_nbr);
        idx = idx + options_.y_nbr*options_.x_nbr;
    end
    if jvar == "R"
        R = reshape(xparam(idx+(1:options_.x_nbr*options_.eta_nbr)),options_.x_nbr,options_.eta_nbr);
        idx = idx + options_.x_nbr*options_.eta_nbr;
    end
    if jvar == "mu_eps"
        mu_eps = xparam(idx+(1:options_.eps_nbr));
        idx = idx + options_.eps_nbr;
    end
    if jvar == "diaglogSigma_eps"
        Sigma_eps = diag(exp(xparam(idx+(1:options_.eps_nbr)))); %undo log transform
        idx = idx + options_.eps_nbr;
    end
    if jvar == "mu_eta"
        mu_eta = xparam(idx+(1:options_.eta_nbr));
        idx = idx + options_.eta_nbr;
    end
    if jvar == "diaglogSigma_eta"
        Sigma_eta = diag(exp(xparam(idx+(1:options_.eta_nbr)))); %undo log transform
        idx = idx + options_.eta_nbr;
    end    
    if strcmp(kf_variant{1},"pruned_skewed")
        if jvar == "diagGamma_eta"
            Gamma_eta = diag(xparam(idx+(1:options_.eta_nbr))); % make full matrix
            idx = idx + options_.eta_nbr;
        end
    end
end

%% set calibrated parameters
for jvar = options_.param_names_fixed(:,1)'
    if any(ismember(fieldnames(params_),jvar))
        eval(sprintf('%s = PARAMS.%s;',jvar,jvar));
    end
end

%% Check stability and positive definitenes
if sum(abs(eig(G)) >= (1-1e-7))
    % disp("Some eigenvalues of G are outside the unit circle")
    negative_log_likelihood = verylarge;
    exit_flag = 0;
    return
end

% Check if covariance matrices are positive definite
if strcmp(kf_variant{1},"pruned_skewed")
    COVeta = csnVar(Sigma_eta,Gamma_eta,nu_eta,Delta_eta,options_.cdfmvna_fct);
elseif strcmp(kf_variant{1},"gaussian")
    COVeta = Sigma_eta;
end
if isdiag(COVeta)
    if any(diag(COVeta)<0)
        negative_log_likelihood = verylarge;
        exit_flag = 0;
        return
    end
elseif ~ispd(COVeta)
    negative_log_likelihood = verylarge;
    exit_flag = 0;
    return
end

if isdiag(Sigma_eps)
    if any(diag(Sigma_eps)<0)
        negative_log_likelihood = verylarge;
        exit_flag = 0;
        return
    end
elseif ~ispd(Sigma_eps)
    negative_log_likelihood = verylarge;
    exit_flag = 0;
    return
end

%% compute log likelihood
if strcmp(kf_variant{1},"gaussian")
    % initialize Kalman filter with a wide Normal prior: x_0 ~ CSN(0,Harvey_factor*eye(nx),0,0,eye(nx)) = N(0,10*eye(nx))
    mu_0     = zeros(options_.x_nbr,1);
    Sigma_0  = options_.Harvey_factor*eye(options_.x_nbr);
    negative_log_likelihood = -1*kalman_gaussian(datamat, mu_0,Sigma_0, G,R,F, mu_eta,Sigma_eta, mu_eps,Sigma_eps, false, true);
elseif strcmp(kf_variant{1},"pruned_skewed")
    % initialize Kalman filter with a wide Normal prior: x_0 ~ CSN(0,Harvey_factor*eye(nx),0,0,eye(nx)) = N(0,10*eye(nx))
    mu_0     = zeros(options_.x_nbr,1);
    Sigma_0  = options_.Harvey_factor*eye(options_.x_nbr);
    Gamma_0  = zeros(options_.x_nbr,options_.x_nbr);
    nu_0     = zeros(options_.x_nbr,1);
    Delta_0  = eye(options_.x_nbr);
    negative_log_likelihood = -1*kalman_csn(datamat, mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, options_.cdfmvna_fct,kf_variant{2},false,true);
end

if isnan(negative_log_likelihood) || isinf(negative_log_likelihood) || ~isreal(negative_log_likelihood)
    negative_log_likelihood = verylarge;
    exit_flag = 0;
end

end % main function end