function [Sigma, Gamma, nu, Delta] = csnPruneParams(Sigma, Gamma, nu, Delta, tol, debug_dim)
% [Sigma, Gamma, nu, Delta] = csnPruneParams(Sigma, Gamma, nu, Delta, tol, debug_dim)
% -------------------------------------------------------------------------
% Reduces the dimension of a CSN distributed random variable according to
% Algorithm 1 of the paper "Pruned Skewed Kalman Filter and Smoother: With Application
% to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
%
% Idea: 
% 1) Representation of CSN distribution:
%    X ~ CSN(mu,Sigma,Gamma,nu,Delta) is equivalent to X = W|Z>=0 where
%   [W;Z] ~ N([mu;-nu], [Sigma, Sigma*Gamma'; Gamma*Sigma, Delta+Gamma*Sigma*Gamma']
% 2) Correlation introduces skewness:
%    Skewness in CSN is based on the correlation between W and Z:
%    - CSN and N are the same distribution if Gamma=0
%    - if correlation is small, then CSN very close to N
% 3) Prune dimensions if correlation in absolute value is below tol.
%    This reduces the overall skewness dimension to qq < q.
% -------------------------------------------------------------------------
% INPUTS
% - Sigma         [p by p]    2nd parameter of the CSN distribution (scale)
% - Gamma         [q by p]    3rd parameter of the CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu            [q by 1]    4th parameter of the CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta         [q by q]    5th parameter of the CSN distribution (enables closure of CSN distribution und marginalization)
% - tol           [double]    threshold value for correlation (in absolute value) below which correlated variables are pruned
% - debug_dim     [boolean]   optional notification if there are more than debug_dim dimensions left unpruned
% where p is the normal dimension and q the skewness dimension
% -------------------------------------------------------------------------
% OUTPUTS
% - Sigma         [p by p]    2nd parameter of the pruned CSN distribution (scale)
% - Gamma         [qq by p]   3rd parameter of the pruned CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu            [qq by 1]   4th parameter of the pruned CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta         [qq by q]   5th parameter of the pruned CSN distribution (enables closure of CSN distribution und marginalization)
% where qq < q is the skewness dimension of the pruned distribution
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
if nargin < 6
    debug_dim = [];
end

p = length(Sigma); % normal dimension

% create correlation matrix of conditional definition of CSN distributed variable
P = [Sigma, Sigma*Gamma'; Gamma*Sigma, Delta+Gamma*Sigma*Gamma']; P = 0.5*(P + P'); % covariance matrix
n = length(P);
try
    % create correlation matrix from covariance matrix
    normalization = diag(1./sqrt(diag(P)));
    R = abs(normalization*P*normalization);
catch
    % error handling
    Sigma = 0;
    Gamma = 0;
    nu    = 0;
    Delta = 0;
    return
end

% prune dimensions in R if they are in absolute value lower than tol
R = R-diag(repelem(Inf, n));
logi2 = max(R(p+1:end, 1:p), [], 2);
logi2 = (logi2 < tol);

% prune dimensions in CSN parameters
Gamma(logi2, :) = [];
nu(logi2)       = [];
Delta(logi2, :) = [];
Delta(:, logi2) = [];

% Notify if there are more than ... dimensions left unpruned
if ~isempty(debug_dim)
    if sum(logi2 > tol) > debug_dim
        disp('----');
        disp(max(logi2)); disp(sum(logi2 > tol));
        disp('----');
    end    
end

end % csnPruneParams