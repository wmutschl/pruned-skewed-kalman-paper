function Skew = csnSkewness_univariate(Sigma, Gamma)
% function Skew = csnSkewness_univariate(Sigma, Gamma)
% -------------------------------------------------------------------------
% computes the skewness coefficient of a univariate CSN distributed random
% variable X ~ CSN(mu,Sigma,Gamma,nu,Delta) with nu=0 and Delta=1
% according to equation 3.6 of Grabek, Klos, and Koloch (2011) - Skew-Normal
% shocks in the linear state space form DSGE model
% -------------------------------------------------------------------------
% INPUTS
% - Sigma      [double]   scale parameter of CSN distribution
% - Gamma      [double]   skewness parameter of CSN distribution
% -------------------------------------------------------------------------
% OUTPUTS
% - Skew       [double]   theoretical skewness coefficient of X
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
term1 = (4 - pi) / 2;
term2 = (sqrt(2/pi) * Gamma * Sigma / sqrt(1 + Gamma^2 * Sigma))^3;
term3 = (Sigma - (2/pi) * (Gamma^2 * Sigma^2) / (1 + Gamma^2 * Sigma))^(3/2);
Skew  = term1 * term2 / term3;

end