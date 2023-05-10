function x_smooth = kalman_csn_smoother(Y, pred,filt, G,R,F, Gamma_eta,Delta_eta, cdfmvna_fct,prune_tol,loss_fct)
% x_smooth = kalman_csn_smoother(Y, pred,filt, G,R,F, Gamma_eta,Delta_eta, cdfmvna_fct,prune_tol,loss_fct)
% -------------------------------------------------------------------------
% Evaluate smoothed states of linear state space model with
% csn distributed innovations and normally distributed noise:
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
% - Y               [y_nbr by obs_nbr]             matrix with data
% - pred            [structure]                    CSN parameters of predicted states (mu_t_tm1,Sigma_t_tm1,Gamma_t_tm1,nu_t_tm1,Delta_t_tm1)
% - filt            [structure]                    CSN parameters of filtered states (mu_t_t,Sigma_t_t,Gamma_t_t,nu_t_t,Delta_t_t)
% - G               [x_nbr by x_nbr]               state transition matrix mapping previous states to current states
% - R               [x_nbr by eta_nbr]             state transition matrix mapping current innovations to current states
% - F               [y_nbr by x_nbr]               observation equation matrix mapping current states into current observables
% - Gamma_eta       [skeweta_dim by eta_nbr]       3rd parameter of CSN distributed innovations eta (regulates skewness from Gaussian (Gamma=0) to half normal)
% - Delta_eta       [skeweta_dim by skeweta_dim]   5th parameter of CSN distributed innovations eta (enables closure of CSN distribution und marginalization)
% - cdfmvna_fct     [string]                       name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% - prune_tol       [double]                       correlation threshold to prune redundant skewness dimensions, if set to 0 no pruning will be done
% - loss_fct        [structure]                    Underlying loss function for computing point estimate of filtered states x_t_t that minimizes the expected loss.
%                                                  loss_fct.type possible values:
%                                                  - "L1"   [boolean]   xtilde = median(x_t_t), i.e. absolute loss abs(xtilde-x)
%                                                  - "L2"   [boolean]   xtilde = E[x_t_t], i.e. squared loss (xtilde-x)^2
%                                                  - "La"   [boolean]   xtilde = quantile(x_t_t,a/(a+b)), i.e. asymmetric loss function a*abs(xtilde-x) for x>xtilde and b*abs(xtilde-x) for x<=xtilde, , a and b are taken from loss_fct.params.a and loss_fct.params.b
% -------------------------------------------------------------------------
% OUTPUTS
% - x_smooth        [structure]                    smoothed states according to different loss functions given by loss_fct.type
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov
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
% Possible avenues for speedup
% - get rid of inv(), but we do need the inverse so inv() might be still a good choice, alternatively try out inv_chol.mex

% get dimensions
x_nbr   = size(F,2);
obs_nbr = size(Y,2);

% initialize output structure
x_smooth.L1 = nan(x_nbr,obs_nbr);
x_smooth.L2 = nan(x_nbr,obs_nbr);
x_smooth.La = nan(x_nbr,obs_nbr);    

% initialize some matrices
Gamma_eta = Gamma_eta/(R'*R)*R';

for t = obs_nbr:-1:1
    if t == obs_nbr
        % smoothing step for the last time point, i.e. 'T'
        mu_t_T    = filt.mu(:,obs_nbr);
        Sigma_t_T = filt.Sigma(:,:,obs_nbr);
        Gamma_t_T = filt.Gamma{obs_nbr,1};
        nu_t_T    = filt.nu{obs_nbr,1};
        Delta_t_T = filt.Delta{obs_nbr,1};
    else
        Jt = filt.Sigma(:,:,t)*G'/pred.Sigma(:,:,t+1);
        % first two distribution parameters of smoothed value
        mu_t_T    = filt.mu(:,t) + Jt*(mu_tp1_T-pred.mu(:,t+1));
        Sigma_t_T = filt.Sigma(:,:,t) + Jt*(Sigma_tp1_T-pred.Sigma(:,:,t+1))*Jt';
	    % evaluate the auxiliary matrices
        Mt  = Sigma_tp1_T*Jt'/Sigma_t_T;
        Nt  = -Gamma_eta*G + Gamma_eta*Mt;
        O_t = [Nt; O_tp1*Mt];

        temp_mat      = [Gamma_eta; O_tp1];
        Delta_tilde_t = blkdiag(Delta_eta, Delta_tilde_tp1) + temp_mat*(Sigma_tp1_T - Mt*Sigma_t_T*Mt')*temp_mat';

	    % last three distribution parameters of smoothed value
        Gamma_t_T = [filt.Gamma{t}; O_t];
        nu_t_T    = filt.nu{obs_nbr, 1};
        Delta_t_T = blkdiag(filt.Delta{t, 1}, Delta_tilde_t);
    end

    % pruning redundant skewness dimension to speed up computations of loss function for smoothed states
    if prune_tol > 0
        [Sigma_t_T, Gamma_t_T, nu_t_T, Delta_t_T] = csnPruneParams(Sigma_t_T,Gamma_t_T,nu_t_T,Delta_t_T,prune_tol);
    end

    if loss_fct.type.L1
        % absolue loss, i.e. compute median of CSN distributed x_t_T
        x_smooth.L1(:,t) = csnQuantile(0.5, mu_t_T, Sigma_t_T, Gamma_t_T, nu_t_T, Delta_t_T, cdfmvna_fct);
    end
    if loss_fct.type.L2
        % squared loss, i.e. compute mean of CSN distributed x_t_T
        x_smooth.L2(:,t) = csnMean(mu_t_T, Sigma_t_T, Gamma_t_T, nu_t_T, Delta_t_T, cdfmvna_fct);
    end
    if loss_fct.type.La
        % asymmetric loss, i.e. compute a/(a+b) quantile of CSN distributed x_t_T
        x_smooth.La(:,t) = csnQuantile(loss_fct.params.a/(loss_fct.params.a+loss_fct.params.b), mu_t_T, Sigma_t_T, Gamma_t_T, nu_t_T, Delta_t_T, cdfmvna_fct);
    end

    % assign for the next step
    if t == obs_nbr
        % empty auxiliary matrices for the last time point, i.e. 'T'
        O_tp1 = double.empty(0, x_nbr);
        Delta_tilde_tp1 = [];
    else
        O_tp1 = O_t;
        Delta_tilde_tp1 = Delta_tilde_t;
    end
    mu_tp1_T    = mu_t_T;
    Sigma_tp1_T = Sigma_t_T;
end

end % kalman_csn_smoother