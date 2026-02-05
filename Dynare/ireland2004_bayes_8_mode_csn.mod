% =========================================================================
% Copyright (C) 2024-2026 Willi Mutschler
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
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAYESIAN MODE ESTIMATION OF CSN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD, MODE IS COMPUTED BY: %
% (1) RUNNING NUMERICAL OPTIMIZATION (MODE_COMPUTE=8) USING SHORT SLICE MODE AS INITIAL GUESS             %
% (2) RUNNING MONTE-CARLO OPTIMIZATION (MODE_COMPUTE=5) TO GET POSITIVE DEFINITE INVERSE HESSIAN AT MODE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
stderr eta_a, ,  0, 10, inv_gamma_pdf, (sqrt(30)),   (sqrt(30)) ;
stderr eta_e, ,  0, 10, inv_gamma_pdf, (sqrt(0.08)), (sqrt(1))  ;
stderr eta_z, ,  0, 10, inv_gamma_pdf, (sqrt(5)),    (sqrt(15)) ;
stderr eta_r, ,  0, 10, inv_gamma_pdf, (sqrt(0.50)), (sqrt(2))  ;
OMEGA,        ,  0,  1, beta_pdf,      0.20,         0.10       ;
RHO_PI,       ,  0,  1, gamma_pdf,     0.30,         0.10       ;
RHO_G,        ,  0,  1, gamma_pdf,     0.30,         0.10       ;
RHO_X,        ,  0,  1, gamma_pdf,     0.25,         0.0625     ;
RHO_A,        ,  0,  1, beta_pdf,      0.85,         0.10       ;
RHO_E,        ,  0,  1, beta_pdf,      0.85,         0.10       ;
skew eta_a,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_e,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_z,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_r,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
end;

estimation(datafile = '../data/ireland2004_data.m'
          , lik_init = 1
          , mh_replic = 0
          , mode_compute = 8
          , mode_file = '../results/ireland2004/bayes/ireland2004_csn_slice_short_mh_mode.mat'
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter
          );
copyfile([M_.dname '/Output/' M_.fname '_mode.mat'], ['../results/ireland2004/bayes/ireland2004_csn_mode8.mat']);

estimation(datafile = '../data/ireland2004_data.m'
          , lik_init = 1
          , mh_replic = 0
          , mode_compute = 6
          , mode_file = '../results/ireland2004/bayes/ireland2004_csn_mode8.mat'
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter
          );
copyfile([M_.dname '/Output/' M_.fname '_mode.mat'], ['../results/ireland2004/bayes/ireland2004_csn_mode6.mat']);
copyfile([M_.dname '/Output/' M_.fname '_optimal_mh_scale_parameter.mat'], ['../results/ireland2004/bayes/ireland2004_csn_optimal_mh_scale_parameter.mat']);

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
