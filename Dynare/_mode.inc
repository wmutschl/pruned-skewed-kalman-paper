%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAYESIAN MODE ESTIMATION OF @{DISTRIB} VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD, MODE IS COMPUTED BY: %
% (1) RUNNING NUMERICAL OPTIMIZATION (MODE_COMPUTE=8) USING SHORT SLICE MODE AS INITIAL GUESS                    %
% (2) RUNNING MONTE-CARLO OPTIMIZATION (MODE_COMPUTE=5) TO GET POSITIVE DEFINITE INVERSE HESSIAN AT MODE         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_.mode_file = sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_slice_250_mh_mode_%s_%s.mat', ARCH, MATLAB_VERSION);
estimation(datafile = '../data/ireland2004_data.m'
          , lik_init = 1
          , mh_replic = 0
          , mode_compute = 8, silent_optimizer
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability
          , skewed_kalman_smoother_skip
          );
fprintf('Optimization info (mode_compute=8):\n');
disp(oo_.posterior.optimization.optimization_info);
copyfile([M_.dname '/Output/' M_.fname '_mode.mat'], sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_mode8_%s_%s.mat', ARCH, MATLAB_VERSION));

options_.mode_file = sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_mode8_%s_%s.mat', ARCH, MATLAB_VERSION);
options_.silent_optimizer = false;
estimation(datafile = '../data/ireland2004_data.m'
          , lik_init = 1
          , mh_replic = 0
          , mode_compute = 6
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability
          , skewed_kalman_smoother_skip
          );
fprintf('Optimization info (mode_compute=6):\n');
disp(oo_.posterior.optimization.optimization_info);
copyfile([M_.dname '/Output/' M_.fname '_mode.mat'], sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_mode6_%s_%s.mat', ARCH, MATLAB_VERSION));
copyfile([M_.dname '/Output/' M_.fname '_optimal_mh_scale_parameter.mat'], sprintf('../results/ireland2004/bayes/ireland2004_@{DISTRIB}_optimal_mh_scale_parameter_%s_%s.mat', ARCH, MATLAB_VERSION));

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
