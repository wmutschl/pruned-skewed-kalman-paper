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
% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF CSN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD           %
% RWMH SAMPLER IS INITIALIZED AT MODE AND COVARIANCE MATRIX FROM MONTE-CARLO OPTIMIZATION (MODE_COMPUTE=6) %
% BECAUSE THIS GAVE THE HIGHEST MODE AND A POSITIVE DEFINITE INVERSE HESSIAN AT THE MODE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
          , mode_file = '../results/ireland2004/bayes/ireland2004_csn_mode6.mat'
          , mh_jscale = 0.4419 % determined as by-product of mode_compute 6, see ireland2004_csn_optimal_mh_scale_parameter.mat
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter
          );

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
