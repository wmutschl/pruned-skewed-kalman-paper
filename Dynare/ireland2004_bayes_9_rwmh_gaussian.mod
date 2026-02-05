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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF GAUSSIAN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD      %
% RWMH SAMPLER IS INITIALIZED AT MODE AND COVARIANCE MATRIX FROM MONTE-CARLO OPTIMIZATION (MODE_COMPUTE=6) %
% BECAUSE THIS GAVE THE HIGHEST MODE AND A POSITIVE DEFINITE INVERSE HESSIAN AT THE MODE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end;

% run estimation in parallel with 8 workers
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', 8);
end
estimation(datafile = '../data/ireland2004_data.m'
          , lik_init = 1
          , posterior_sampling_method = 'random_walk_metropolis_hastings'
          , posterior_sampler_options = ('proposal_distribution','rand_multivariate_normal')
          , mh_replic = 250000
          , mh_nblocks = 8
          , mode_compute = 0
          , mode_file = '../results/ireland2004/bayes/ireland2004_gaussian_mode6.mat'
          , mh_jscale = 0.5731 % determined as by-product of mode_compute 6, see ireland2004_gaussian_optimal_mh_scale_parameter.mat
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability;
                             % note that Dynare's implementation is much faster because it switches to the steady-state Kalman filter which we have not implemented yet for the PSKF
          );

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
