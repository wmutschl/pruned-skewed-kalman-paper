function [mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform(G_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar)
% function [mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform(G_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% 1) make proper csn random variable first (in case of singularity)
% 2) do either rank deficient or full rank linear transformation
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - mu_bar = [mu_t_t; mu_eta]
%   - nu_bar = [nu_t_t; nu_eta];
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
% -----------------------------------------------------------------------
% OUTPUT
% transformed parameters of CSN joint distribution of x_t and eta_t:
% mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr
% =========================================================================
% Copyright (C) 2023 Gaygysyz Guljanov
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
tol = 1e-8;

mu_tr = G_bar * mu_bar;
Sigma_tr = G_bar * Sigma_bar * G_bar';
mu_bar = zeros(size(mu_bar)); % make mu_bar a zero vector

%% singular skew-normal to closed skew-normal   
if rcond(Sigma_bar) < 1e-10
    % Sigma_bar is rank deficient
    [S_hat, ~, Sigma_bar, Gamma_bar] = singular_skew_normal_2_closed_skew_normal(mu_bar, Sigma_bar, Gamma_bar, tol, 0);        
    G_bar = G_bar * S_hat;
end

%% Linear Transformation    
% Check if A_mat is rank deficient
dimensions = size(G_bar);
if dimensions(1) > dimensions(2)
    is_singular = rcond(G_bar' * G_bar) < 1e-8;
else
    is_singular = rcond(Sigma_tr) < 1e-8;
end
% If not full rank, use rank deficient linear transformation
if is_singular
    [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_rank_deficient(G_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar, tol);
    return
end
[Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar);
        
end % csn_statespace_linear_transform