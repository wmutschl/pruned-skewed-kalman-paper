function [S_hat, mu_hat, Lambda_hat, Gamma_hat] = singular_skew_normal_2_closed_skew_normal(mu_bar, Sigma_bar, Gamma_bar, tol, hmdel)
% function [S_hat, mu_hat, Lambda_hat, Gamma_hat] = singular_skew_normal_2_closed_skew_normal(mu_bar, Sigma_bar, Gamma_bar, tol, hmdel)
% -------------------------------------------------------------------------
% write Singular Skew Normal (SSN) distribution as a linear transformation
% of proper Closed Skew Normal (CSN) distribution
% -----------------------------------------------------------------------
% INPUTS
% mu_bar:      first parameter of the SSN distribution
% Sigma_bar:   second parameter of the SSN distribution
% Gamma_bar:   third parameter of the SSN distribution
% tol:         tolerance level showing which values should be assumed to be  numerically zero
% hmdel:       minimum how many dimensions should be deleted, if known already
% -----------------------------------------------------------------------
% OUTPUT
% parameters of transformed SSN distribution that is CSN
% S_hat:        multiplying matrix, which has more rows than columns
% mu_hat    :   first parameter of the proper CSN
% Lambda_hat:   second parameter of the proper CSN
% Gamma_hat :   third parameter of the proper CSN
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
if sum(mu_bar ~= 0)
    error("mu_bar should be a zero vector")
end
[S_mat, Lambda_mat] = schur(Sigma_bar); % get Schur decomposition
% Create a hold_vec, which tells which dimensions to hold
if hmdel > 0
    % Reorder the diagonal elements of Lambda, in an increasing order
    len = length(Lambda_mat);
    diag_el = diag(Lambda_mat);
    [~, min_indices] = mink(diag_el, len);
    % Should we delete more than hmdel
    hmdel_tmp = hmdel;
    ii = 1;
    while diag_el(min_indices(hmdel_tmp + ii)) / diag_el(min_indices(len)) <= tol
        ii = ii + 1;
    end
    hmdel_tmp = hmdel_tmp + ii - 1;
    if hmdel_tmp >= len
        error("Sigma matrix should at least have rank one");
    end
    % Which vector to hold, i.e. to not cut
    hold_vec = true(len, 1);
    hold_vec(min_indices(1:hmdel_tmp)) = false;
else
    hold_vec = diag(Lambda_mat) > tol;
end
% Hold the relevant dimensions, i.e. delete the unnecessary ones
S_hat = S_mat(:, hold_vec);    
mu_hat = zeros(sum(hold_vec), 1);
Lambda_hat = Lambda_mat(hold_vec, hold_vec);
Gamma_hat = Gamma_bar * S_hat;

end % singular_skew_normal_2_closed_skew_normal
