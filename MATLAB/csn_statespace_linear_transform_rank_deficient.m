function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_rank_deficient(G_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar, tol)
% function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform(G_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% helper function for rank deficient linear transformation of CSN distributed random variable
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - nu_bar = [nu_t_t; nu_eta];
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
%   - tol: tolerance level showing which values should be assumed to be numerically zero
% -----------------------------------------------------------------------
% OUTPUT
% transformed parameters of rank-deficient CSN joint distribution of x_t and eta_t:
% Gamma_tr, nu_tr, Delta_tr (returns only skewness parameters, as other parameters are easy to obtain)
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
[nrow, ncol] = size(G_bar);
difference = abs(nrow - ncol);
[S_mat, Lambda_mat, T_mat] = svd(G_bar); % apply singular value decomposition
if nrow >= ncol % more rows
    S_mat = S_mat(:, 1:end-difference);
    Lambda_mat = Lambda_mat(1:end-difference, :);
else
    T_mat = T_mat(:, 1:end-difference);
    Lambda_mat = Lambda_mat(:, 1:end-difference);
end    
hold_vec = abs(diag(Lambda_mat)) > tol;      % which dimensions to delete
S_mat = S_mat(:, hold_vec);                  % delete respective columns
Lambda_mat = Lambda_mat(hold_vec, hold_vec); % delete respective rows and columns    
T_mat = T_mat(:, hold_vec);                  % delete respective columns

% Apply the final formulas of theorem
Tmp_mat = T_mat' * Sigma_bar * T_mat;
Tmp_mat2 = Gamma_bar * Sigma_bar;
Tmp_mat3 = Tmp_mat2 * T_mat / Tmp_mat;
Gamma_tr = (Tmp_mat3 / Lambda_mat) / (S_mat' * S_mat) * S_mat';
nu_tr = nu_bar;
Delta_tr = Delta_bar + Tmp_mat2 * Gamma_bar' - Tmp_mat3 * T_mat' * Sigma_bar * Gamma_bar';

end % csn_statespace_linear_transform_rank_deficient