%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian estimation with Slice sampler of @{DISTRIB} version of model using PSKF to compute likelihood %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run estimation in parallel with 8 workers
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', 8);
end
estimation(datafile = '../data/ireland2004_data.m'
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability
          , lik_init = 1
          , posterior_sampling_method = 'slice'
          , posterior_sampler_options = ('save_tmp_file',0, 'save_iter_info_file',0)
          , mh_nblocks = 8
          , mh_replic = @{MH_REPLIC}
          , mh_drop = @{MH_DROP}
          , plot_priors = 0
          , use_pct
          );
copyfile([M_.dname '/Output/' M_.fname '_mh_mode.mat'], sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_slice_@{MH_REPLIC}_mh_mode_%s_%s.mat', ARCH, MATLAB_VERSION));

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
