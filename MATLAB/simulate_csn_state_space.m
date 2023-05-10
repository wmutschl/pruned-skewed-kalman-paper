function [y, x] = simulate_csn_state_space(obs_nbr,burnin, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps)
% [y, x] = simulate_csn_state_space(obs_nbr,burnin, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps)
% -------------------------------------------------------------------------
% Simulate data from linear state-space system with CSN-distributed innovations:
%   x(t) = G*x(t-1) + R*eta(t)  [state transition equation]
%   y(t) = F*x(t)   + eps(t)    [observation equation]
% with csn distributed innovations eta and normally distributed noise eps:
%   eta(t) ~ CSN(mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta) [innovations, shocks]
%   eps(t) ~ N(mu_eps, Sigma_eps)                                 [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% -------------------------------------------------------------------------
% INPUTS
% - obs_nbr     [scalar]               sample size
% - burnin      [scalar]               number of periods to discard after simulation
% - G           [x_nbr by x_nbr]       state transition matrix mapping previous states to current states
% - R           [x_nbr by eta_nbr]     state transition matrix mapping current innovations to current states
% - F           [y_nbr by x_nbr]       observation equation matrix mapping current states into current observables
% - mu_eta      [eta_nbr by 1]         location parameter of CSN distributed innovations eta
% - Sigma_eta   [eta_nbr by eta_nbr]   scale parameter of CSN distributed innovations eta
% - Gamma_eta   [q by eta_nbr]         first skewness parameter of CSN distributed innovations eta
% - nu_eta      [q by 1]               second skewness parameter of CSN distributed innovations eta
% - Delta_eta   [q by q]               third skewness parameter of CSN distributed innovations eta
% - mu_eps      [y_nbr by 1]           mean of Gaussian measurement error
% - Sigma_eps   [y_nbr by y_nbr]       covariance matrix of Gaussian measurement error
% -------------------------------------------------------------------------
% OUTPUTS
% - y           [y_nbr by obs_nbr]     matrix of simulated observable variables
% - x           [x_nbr by obs_nbr]     matrix of simulated state variables
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler
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

% get dimensions
%q = size(Gamma_eta, 1);  % skewness dimension of eta
[y_nbr,x_nbr] = size(F);  % y_nbr number of observables, x_nbr number of states
obs_nbr = obs_nbr+burnin; % add burnin to effective sample size

% initialize states and observables at zero
x = zeros(x_nbr,obs_nbr);
y = zeros(y_nbr,obs_nbr);

for t=2:obs_nbr
    % simulate state variables
    if all(Gamma_eta==0) % Gaussian case
        x(:,t) = G*x(:,t-1) + R*(mu_eta + chol(Sigma_eta,'lower')*randn(x_nbr,1));
    else % CSN case
        x(:,t) = G*x(:,t-1) + R*csnRandDirect(1, mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta);
    end
    % simulate observables
    y(:,t) = F*x(:,t) + mu_eps + chol(Sigma_eps, 'lower')*randn(y_nbr, 1);
end

% remove burn-in
y = y(:,(burnin+1):end);
x = x(:,(burnin+1):end);


end % simulate_csn_state_space