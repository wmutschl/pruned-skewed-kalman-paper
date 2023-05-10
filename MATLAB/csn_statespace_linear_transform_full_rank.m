function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar)
% function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% helper function for full rank linear transformation of CSN distributed random variable
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Sigma_tr : transformed (rank deficient) parameter matrix
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - nu_bar = [nu_t_t; nu_eta];
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
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

%% more rows
if nrow >= ncol
    nu_tr = nu_bar;
    Delta_tr = Delta_bar;
    if nrow == ncol
        Gamma_tr = Gamma_bar / G_bar;
        return
    end
    Gamma_tr = Gamma_bar / (G_bar' * G_bar) * G_bar';
    return
end

%% more columns
Tmp_mat = Sigma_bar * G_bar' / Sigma_tr;
Gamma_tr = Gamma_bar * Tmp_mat;
nu_tr = nu_bar;
Delta_tr = Delta_bar + Gamma_bar * Sigma_bar * Gamma_bar' - Gamma_bar * Tmp_mat * G_bar * Sigma_bar * Gamma_bar';

end