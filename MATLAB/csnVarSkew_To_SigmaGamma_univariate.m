function [Sigma, Gamma, error_indicator] = csnVarSkew_To_SigmaGamma_univariate(Var, Skew, check)
% function [Sigma, Gamma, error_indicator] = csnVarSkew_To_SigmaGamma_univariate(Var, Skew, check)
% -------------------------------------------------------------------------
% recovers the Sigma and Gamma parameters from the Variance and Skewness of
% a univariate CSN distributed random variable X ~ CSN(mu,Sigma,Gamma,nu,Delta)
% where Var=V[X], Skew=skewness[X], nu=0, and Delta=1
% this basically inverts the formulas provided in equation 3.6 of Grabek,
% Klos, and Koloch (2011) - Skew-Normal shocks in the linear state space form DSGE model
% -------------------------------------------------------------------------
% INPUTS
% - Var     [double]   unconditional variance of X
% - Skew    [double]   theoretical skewness coefficient of X
% - check   [boolean]  check whether transformation worked
% -------------------------------------------------------------------------
% OUTPUTS
% - Sigma             [double]    scale parameter of CSN distribution
% - Gamma             [double]    skewness parameter of CSN distribution
% - error_indicator   [boolean]   1: if transformation failed (mostly if Gamma is extremely large and/or Sigma is extremely small)
% =========================================================================
% Copyright (C) 2023 Willi Mutschler
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
if nargin < 3
    check = false;
end
error_indicator = 0;
Sigma = (-2^(2/3) * (4-pi)^(1/3) * abs(Skew)^(2/3) + pi - 4) * Var / (pi - 4);
Gamma = sign(Skew) * ( 2^(1/3) * (4-pi)^(2/3)*sqrt(pi)*sqrt( abs(Skew)^(2/3) / ( ( -2 * 2^(1/3) * (4-pi)^(2/3) * (pi-2) * abs(Skew)^(4/3) + 2^(2/3) * (4-pi)^(7/3) * abs(Skew)^(2/3) + 2 * (pi - 4)^2 ) * Var ) ) );

if check 
    if (abs(csnVar(Sigma,Gamma,0,1) - Var)>1e-13) || (abs(csnSkewness_univariate(Sigma,Gamma) - Skew) > 1e-7)
        fprintf('Sigma: %f\n',Sigma);
        fprintf('Gamma: %f\n',Gamma);
        fprintf('Var: %f\n',csnVar(Sigma,Gamma,0,1));
        fprintf('Skew: %f\n',csnSkewness_univariate(Sigma,Gamma));
        fprintf('abs(csnVar(Sigma,Gamma,0,1) - Var): %f\n',abs(csnVar(Sigma,Gamma,0,1) - Var));
        fprintf('abs(csnSkewness_univariate(Sigma,Gamma) - Skew): %f\n',abs(csnSkewness_univariate(Sigma,Gamma) - Skew));        
        warning('could not transform the parameters')
        error_indicator=1;
    end
end