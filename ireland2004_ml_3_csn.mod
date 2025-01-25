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

@#define TRANSFORM_PARAMETERS    = 0  // new feature to transform selected bounded parameters to unbounded domain using logit transform (only during optimization)
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMUM LIKELIHOOD ESTIMATION OF CSN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values from ireland2004_ml_2_csn_initval_search.mod (see log file or Output folder)
estimated_params;
stderr eta_a, 2.747697233134916,        0,    10;
stderr eta_e, 0.028882164116404,        0,    10;
stderr eta_z, 0.858196786024231,        0,    10;
stderr eta_r, 0.276196004637038,        0,    10;
OMEGA,        0.058086387283355,        0,     1;
%ALPHA_X,      1.003423775824554e-05, 1e-5,     1;
%ALPHA_PI,     1.020447102105185e-05, 1e-5,     1;
RHO_PI,       0.386477141360770,        0,     1;
RHO_G,        0.396013670264708,        0,     1;
RHO_X,        0.165402611432957,        0,     1;
RHO_A,        0.904795669581000,        0,     1;
RHO_E,        0.990673584292426,        0,     1;
% this is not yet possible with preprocessor
% skew eta_a,  -0.111987050675084, -0.995, 0.995;
% skew eta_e,  -0.289566734765780, -0.995, 0.995;
% skew eta_z,  -0.836095597248922, -0.995, 0.995;
% skew eta_r,   0.520016288312854, -0.995, 0.995;
end;
% manually create structure skew_exo
estim_params_.skew_exo(1,:) = [find(ismember(M_.exo_names,'eta_a')), -0.111987050675084, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(2,:) = [find(ismember(M_.exo_names,'eta_e')), -0.289566734765780, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(3,:) = [find(ismember(M_.exo_names,'eta_z')), -0.836095597248922, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(4,:) = [find(ismember(M_.exo_names,'eta_r')),  0.520016288312854, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
@#include "__transformParameters.inc" // optional new feature to use logit transform of bounded parameters during optimization only

%options_.kalman.pskf.prune_tol = 0.01; % 0.01 is default pruning threshold for Pruned Skewed Kalman filter
options_.hessian.use_onesided = true; % new feature: because skewnesss coefficients (and possibly ALPHAX and ALPHAPI) are extremely close to lower bound,
                                      % we use a modified hessian.m to compute standard errors;
                                      % note that this might trigger "OPTIMIZATION PROBLEM! (minus) the hessian matrix at the "mode" is not positive definite!",
                                      % therefore: check s.d. and t-stat to judge whether there is an OPTIMIZATION PROBLEM or not
estimation(datafile = 'data/ireland2004_data.m'
          , mode_compute = 8
          , additional_optimizer_steps = [8]
          , lik_init = 1
          , use_univariate_filters_if_singularity_is_detected = 0
          , keep_kalman_algo_if_singularity_is_detected
          );

% save shock parameters for graphs
csn = M_.csn;
save([M_.dname filesep 'Output' filesep M_.fname '_shock_params'],'csn','-v6');