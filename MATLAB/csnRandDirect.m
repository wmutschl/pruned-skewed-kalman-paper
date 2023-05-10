function [csn_draws] = csnRandDirect(n, mu, Sigma, Gamma, nu, Delta)
% function [res] = csnRandDirect(n, mu, Sigma, D, nu, Delta)
% -------------------------------------------------------------------------
% Draws random numbers from CSN distribution using accept-reject method
% adapted from the "rcsn" command in the "csn" package of the R language
% -------------------------------------------------------------------------
% INPUTS
% - n          [scalar]   number of desired random draws
% - mu         [p by 1]   1st parameter of the CSN distribution (location)
% - Sigma      [p by p]   2nd parameter of the CSN distribution (scale)
% - Gamma      [q by p]   3rd parameter of the CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu         [q by 1]   4th parameter of the CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta      [q by q]   5th parameter of the CSN distribution (enables closure of CSN distribution und marginalization)
% -------------------------------------------------------------------------
% OUTPUTS
% - csn_draws  [p by n]   random draws of the CSN(mu,Sigma,Gamma,nu,Delta) distribution
% =========================================================================
% Copyright (C) 2015 GPL v2 Dmitry Pavlyuk, Eugene Girtcius (rcsn.R function in the "csn" R package)
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

q      = length(nu);
v      = nan(q, n);
Sigma2 = Delta + Gamma * Sigma * Gamma';
L      = chol(Sigma2, 'lower');
for i = 1:n
    draws = -1;
    while sum(draws >= 0) ~= q
        draws = -nu + L * randn(q, 1);
    end
    v(:, i) = draws;
end

tmp = Sigma * Gamma' / Sigma2;
Exp = mu + tmp * (v + nu);
Var = Sigma - tmp * Gamma * Sigma;
L   = chol(Var, 'lower');
csn_draws = Exp + L * randn(length(mu), n);

end % csnRandDirect