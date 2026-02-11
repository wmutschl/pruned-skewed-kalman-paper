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
@#define BEST_OF = 5
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL VALUE SEARCH FOR MAXIMUM LIKELIHOOD ESTIMATION OF CSN VERSION OF MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) fix model and stderr parameters to Gaussian estimates from maximum likelihood
% (2) create large grid for skew parameters
% (3) evaluate likelihood (computed by PSKF) on grid
% (4) use best five skew parameter combinations as initial guess
% (5) optimize likelihood (computed by PSKF) over stderr and skew parameters

% fix model parameters to Gaussian estimates (see ireland2004_ml_1_gaussian.log)
OMEGA    = 0.058086387283355;
RHO_PI   = 0.386477141360770;
RHO_G    = 0.396013670264708;
RHO_X    = 0.165402611432957;
RHO_A    = 0.904795669581000;
RHO_E    = 0.990673584292426;

% fix stderr parameters to Gaussian estimates, see ireland2004_ml_1_gaussian.log
shocks;
var eta_a; stderr 3.016725412332537;
var eta_e; stderr 0.024764214899394;
var eta_z; stderr 0.886476699633342;
var eta_r; stderr 0.279030970775836;
end;

% estimate only skewness parameters (theoretical bounds are automatically set)
estimated_params;
skew eta_a, 0;
skew eta_e, 0;
skew eta_z, 0;
skew eta_r, 0;
end;

% run estimation command without optimization to initialize structures
estimation(datafile = '../data/ireland2004_data.m'
         , dirname = ireland2004_ml_2_csn_initval_search_0 % store results here
         , mode_compute = 0 % no optimization
         , kalman_algo = 5  % use pruned skewed Kalman filter
         , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
         , cova_compute = 0 % skip Hessian computation
         , skewed_kalman_prune_tol = 0.1 % use high pruning threshold (0.01 is default)
         , skewed_kalman_smoother_skip
         );

% create grid for skewness coefficients
grid.nbr      = 16;   % needs to be even number
grid.endpoint = 0.95; % grid is set between +- this value
Skew_eta_grid = linspace(-abs(grid.endpoint),0,grid.nbr/2);
Skew_eta_grid = [Skew_eta_grid -Skew_eta_grid((end-1):-1:1)];
ranges = repmat({1:(grid.nbr-1)}, 1, estim_params_.nsx);
[C{1:estim_params_.nsx}] = ndgrid(ranges{:}); % use ndgrid to produce an n-dimensional grid where each row is a combination (there are (m-1)^n combinations, each is an n-element tuple)
COMBOS = zeros((grid.nbr-1)^estim_params_.nsx, estim_params_.nsx, 'int16');
for i = 1:estim_params_.nsx
    COMBOS(:, i) = C{i}(:);
end

% evaluate likelihood for grid values
neg_log_likelihood_grid = nan(1,(grid.nbr-1)^estim_params_.nsx);
fprintf('\nEvaluate likelihood on %u grid values for skewness parameter combinations using MATLAB''s parallel toolbox, takes less than 10 minutes with 8 cores\n\n', (grid.nbr-1)^estim_params_.nsx);
tStart = tic;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local'); % open parpool
end
parfor jgrid = 1:(grid.nbr-1)^estim_params_.nsx
    [dataset, datasetInfo, xparam1, ~, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options_.varobs, M_.dname, [], M_, options_, oo_, estim_params_, bayestopt_);
    estim_params.skew_exo(:,2) = Skew_eta_grid(COMBOS(jgrid,:))';
    [dataset, datasetInfo, xparam1, ~, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options.varobs, M.dname, [], M, options, oo, estim_params, bayestopt);
    neg_log_likelihood_grid(jgrid) = dsge_likelihood(xparam1,dataset,datasetInfo,options,M,estim_params,bayestopt,BoundsInfo,oo.dr, oo.steady_state,oo.exo_steady_state,oo.exo_det_steady_state);
end
delete(gcp('nocreate')); % close parpool
fprintf('Best @{BEST_OF} values on grid (runtime %s):\n', dynsec2hms(toc(tStart)));
[~,idx_best_grid] = sort(neg_log_likelihood_grid);
disp(array2table([Skew_eta_grid(COMBOS(idx_best_grid(1:@{BEST_OF}),:))'],...
                 'RowNames',["skew " + M_.exo_names],'VariableNames',"best " + transpose(1:@{BEST_OF})));

% optimize stderr and skewness parameters with best guesses from grid search
warning('off')
@#for j in 1:BEST_OF
skew_eta_a = Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_a'))));
skew_eta_e = Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_e'))));
skew_eta_z = Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_z'))));
skew_eta_r = Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_r'))));

estimated_params(overwrite);
% fix stderr parameters to Gaussian estimates, see ireland2004_ml_1_gaussian.log
stderr eta_a, 3.016725412332537, 0, 10;
stderr eta_e, 0.024764214899394, 0, 10;
stderr eta_z, 0.886476699633342, 0, 10;
stderr eta_r, 0.279030970775836, 0, 10;
skew eta_a, (skew_eta_a);
skew eta_e, (skew_eta_e);
skew eta_z, (skew_eta_z);
skew eta_r, (skew_eta_r);
end;

fprintf('\n%s\n* OPTIMIZE SHOCK PARAMETERS: RUN @{j} / @{BEST_OF} *\n%s\n', repmat('*',1,40),repmat('*',1,40));
estimation(datafile = '../data/ireland2004_data.m'
          , dirname = ireland2004_ml_2_csn_initval_search_@{j} % store results here
          , mode_compute = 8 % other optimizers are similar
          , silent_optimizer % below we display optimization_info, so don't show intermediate optimization output
          , kalman_algo = 5  % use pruned skewed Kalman filter
          , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          , cova_compute = 0 % skip Hessian computation
          , skewed_kalman_prune_tol = 0.01 % switch back to default pruning threshold
          , skewed_kalman_smoother_skip
          );
oo_.initval_search.fval(@{j})          = oo_.posterior.optimization.log_density;
oo_.initval_search.shocks_stderr{@{j}} = oo_.mle_mode.shocks_std;
oo_.initval_search.shocks_skew{@{j}}   = oo_.mle_mode.shocks_skew;
fprintf('Optimization info:\n');
disp(oo_.posterior.optimization.optimization_info);
format long;
fprintf('MLE point estimates for stderr parameters:\n');
disp(oo_.mle_mode.shocks_std);
fprintf('MLE point estimates for skew parameters:\n');
disp(oo_.mle_mode.shocks_skew);
format short;
fprintf('Final value of log-likelihood: %.15f\n', oo_.posterior.optimization.log_density);
@#endfor
warning('on');

% display best results
[~, idx_best] = max(oo_.initval_search.fval);
fprintf('\n%s\n* BEST RESULTS WITH RUN=%d *\n%s\n\n', repmat('*',1,27), idx_best, repmat('*',1,27));
format long
fprintf('- point estimates for stderr parameters:\n');
disp(oo_.initval_search.shocks_stderr{idx_best});
fprintf('- point estimates for skew parameters:\n');
disp(oo_.initval_search.shocks_skew{idx_best});
format short
fprintf('\n');

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
