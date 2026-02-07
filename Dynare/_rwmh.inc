%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian estimation with RWMH sampler of @{DISTRIB} version of model using PSKF to compute likelihood %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run estimation in parallel with 8 workers
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', 8);
end
options_.mode_file = sprintf('@{MODE_FILE}_%s_%s.mat', ARCH, MATLAB_VERSION);
estimation(datafile = '../data/ireland2004_data.m'
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability
          , lik_init = 1
          , posterior_sampling_method = 'random_walk_metropolis_hastings'
          , posterior_sampler_options = ('proposal_distribution','rand_multivariate_normal')
          , mh_nblocks = 8
          , mh_replic = 250000
          , mode_compute = 0
          , mh_jscale = @{MH_JSCALE}
          , plot_priors = 0
          );

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
