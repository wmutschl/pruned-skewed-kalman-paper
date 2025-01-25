% =========================================================================
% Copyright (C) 2024-2025 Willi Mutschler
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMUM LIKELIHOOD ESTIMATION OF GAUSSIAN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD                     %
% NOTE THAT WE DON'T USE DYNARE'S KALMAN FILTER FUNCTION BUT THE PSKF EVEN IN GAUSSIAN CASE FOR A FAIR COMPARISON %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#define TRANSFORM_PARAMETERS    = 0  // new feature to transform selected bounded parameters to unbounded domain using logit transform (only during optimization)
@#include "_ireland2004_common.inc"

estimated_params;
stderr eta_a, 3.0200,  0,    10;
stderr eta_e, 0.0200,  0,    10;
stderr eta_z, 0.8900,  0,    10;
stderr eta_r, 0.2800,  0,    10;
OMEGA,        0.0581,  0,     1;
%ALPHA_X,      0.00002, 1e-5,  1;
%ALPHA_PI,     0.00002, 1e-5,  1;
RHO_PI,       0.3866,  0,     1;
RHO_G,        0.3960,  0,     1;
RHO_X,        0.1654,  0,     1;
RHO_A,        0.9048,  0,     1;
RHO_E,        0.9907,  0,     1;
end;
% structure estim_params_.skew_exo needs to be created by the preprocessor
estim_params_.skew_exo = zeros(0, 11);
@#include "__transformParameters.inc" // optional new feature to use logit transform of bounded parameters during optimization only

options_.kalman.pskf.use_in_gaussian_case = true; % if shocks are Gaussian, still use PSKF filter instead of Dynare's Gaussian Kalman filter
options_.hessian.use_onesided = true; % new feature: because ALPHAX and ALPHAPI are extremely close to lower bound, we use a modified hessian.m to compute standard errors;
                                      % note that this might trigger "OPTIMIZATION PROBLEM! (minus) the hessian matrix at the "mode" is not positive definite!",
                                      % therefore: check s.d. and t-stat to judge whether there is an OPTIMIZATION PROBLEM or not
estimation(datafile = 'data/ireland2004_data.m'
          , mode_compute = 8
          , lik_init = 1
          , use_univariate_filters_if_singularity_is_detected = 0
          , keep_kalman_algo_if_singularity_is_detected
          );
format long;
disp(oo_.mle_mode.parameters);
disp(oo_.mle_mode.shocks_std);
format short;

csn = M_.csn;
save([M_.dname filesep 'Output' filesep M_.fname '_shock_params'],'csn','-v6');