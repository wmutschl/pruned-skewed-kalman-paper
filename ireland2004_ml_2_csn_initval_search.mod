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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL VALUE SEARCH FOR MAXIMUM LIKELIHOOD ESTIMATION OF CSN VERSION OF MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) fix model and stderr parameters to Gaussian estimates from maximum likelihood
% (2) create large grid for skew parameters
% (3) evaluate likelihood (computed by PSKF) on grid
% (4) use best five skew parameter combinations as initial guess
% (5) optimize likelihood (computed by PSKF) over stderr and skew parameters

@#define TRANSFORM_PARAMETERS    = 0  // new feature to transform selected bounded parameters to unbounded domain using logit transform (only during optimization)
@#define BEST_OF = 5

@#include "_ireland2004_common.inc"

% fix model parameters to Gaussian estimates from ireland2004_ml_1_gaussian.mod (see log file or Output folder)
OMEGA    = 0.058086387283355;
%ALPHA_X  = 1.003423775824554e-05;
%ALPHA_PI = 1.020447102105185e-05;
RHO_PI   = 0.386477141360770;
RHO_G    = 0.396013670264708;
RHO_X    = 0.165402611432957;
RHO_A    = 0.904795669581000;
RHO_E    = 0.990673584292426;
% fix stderr parameters to Gaussian estimates from maximum likelihood
shocks; // note that the shocks block calibrates the covariance matrix; with CSN shocks this covariance matrix is NOT equal to Sigma (unless Gamma = 0)
var eta_a; stderr 3.016725412332537;
var eta_e; stderr 0.024764214899394;
var eta_z; stderr 0.886476699633342;
var eta_r; stderr 0.279030970775836;
end;

% create estimated_params block only for skewness parameters (no interface yet for skeweness parameters)
% this is what preprocessor needs to understand:
% estimated_params;
% skew eta_a, 0, -0.995, 0.995;
% skew eta_e, 0, -0.995, 0.995;
% skew eta_z, 0, -0.995, 0.995;
% skew eta_r, 0, -0.995, 0.995;
% end;
estim_params_.var_exo    = zeros(0, 10); % no stderr parameters for shocks to estimate
estim_params_.var_endo   = zeros(0, 10); % no stderr parameters for measurement errors to estimate
estim_params_.corrx      = zeros(0, 11); % no corr parameters for shocks to estimate
estim_params_.corrn      = zeros(0, 11); % no corr parameters for measurement errors to estimate
estim_params_.param_vals = zeros(0, 10); % no model parameters to estimate
estim_params_.skew_exo(1,:) = [find(ismember(M_.exo_names,'eta_a')), 0, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0]; % estimate skewness coefficient of eta_a
estim_params_.skew_exo(2,:) = [find(ismember(M_.exo_names,'eta_e')), 0, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0]; % estimate skewness coefficient of eta_e
estim_params_.skew_exo(3,:) = [find(ismember(M_.exo_names,'eta_z')), 0, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0]; % estimate skewness coefficient of eta_z
estim_params_.skew_exo(4,:) = [find(ismember(M_.exo_names,'eta_r')), 0, -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0]; % estimate skewness coefficient of eta_r
@#include "__transformParameters.inc" // optional new feature to use logit transform of bounded parameters during optimization only

% run estimation command without optimization to initialize structures
options_.kalman.pskf.prune_tol = 0.1; % use high pruning threshold (0.01 is default)
estimation(datafile = 'data/ireland2004_data.m'
         , mode_compute = 0
         , lik_init = 1
         , use_univariate_filters_if_singularity_is_detected = 0
         , keep_kalman_algo_if_singularity_is_detected
         , dirname = ireland2004_ml_2_csn_initval_search_0
         , cova_compute = 0
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
fprintf('\nEvaluate likelihood on %u grid values for skewness parameter combinations using MATLAB''s parallel toolbox, takes roughly 5 minutes with 8 cores\n\n', (grid.nbr-1)^estim_params_.nsx);
tStart = tic;
parfor jgrid = 1:(grid.nbr-1)^estim_params_.nsx
    [dataset, datasetInfo, xparam1, ~, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options_.varobs, M_.dname, [], M_, options_, oo_, estim_params_, bayestopt_);
    estim_params.skew_exo(:,2) = Skew_eta_grid(COMBOS(jgrid,:))';
    [dataset, datasetInfo, xparam1, ~, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options.varobs, M.dname, [], M, options, oo, estim_params, bayestopt);
    neg_log_likelihood_grid(jgrid) = dsge_likelihood(xparam1,dataset,datasetInfo,options,M,estim_params,bayestopt,BoundsInfo,oo.dr, oo.steady_state,oo.exo_steady_state,oo.exo_det_steady_state);
end
tEnd = toc(tStart)
% display best @{BEST_OF} values on grid
[~,idx_best_grid] = sort(neg_log_likelihood_grid);
disp(array2table([Skew_eta_grid(COMBOS(idx_best_grid(1:@{BEST_OF}),:))'],...
                 'RowNames',["skew " + M_.exo_names],'VariableNames',"best " + transpose(1:@{BEST_OF})));

% optimize stderr and skewness parameters with best guesses from grid search
warning('off')
@#for j in 1:BEST_OF
estimated_params(overwrite);
stderr eta_a, 3.016725412332537, 0, 10;
stderr eta_e, 0.024764214899394, 0, 10;
stderr eta_z, 0.886476699633342, 0, 10;
stderr eta_r, 0.279030970775836, 0, 10;
% this is not yet possible with preprocessor
% skew eta_a, (Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_a'))))), -0.995, 0.995;
% skew eta_e, (Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_e'))))), -0.995, 0.995;
% skew eta_z, (Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_z'))))), -0.995, 0.995;
% skew eta_r, (Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_r'))))), -0.995, 0.995;
end;
% manually create structure skew_exo
estim_params_.skew_exo(1,:) = [find(ismember(M_.exo_names,'eta_a')), Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_a')))), -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(2,:) = [find(ismember(M_.exo_names,'eta_e')), Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_e')))), -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(3,:) = [find(ismember(M_.exo_names,'eta_z')), Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_z')))), -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
estim_params_.skew_exo(4,:) = [find(ismember(M_.exo_names,'eta_r')), Skew_eta_grid(COMBOS(idx_best_grid(@{j}),find(ismember(M_.exo_names,'eta_r')))), -0.995, 0.995, 0, NaN, NaN, NaN, NaN, NaN, 0];
@#include "__transformParameters.inc" // optional new feature to use logit transform of bounded parameters during optimization only

fprintf('OPTIMIZE SHOCK PARAMETERS: RUN @{j} / 5\n');
options_.kalman.pskf.prune_tol = 0.01; % switch back to default pruning threshold
estimation(datafile = 'data/ireland2004_data.m'
          , mode_compute = 8
          , silent_optimizer
          , lik_init = 1
          , use_univariate_filters_if_singularity_is_detected = 0
          , keep_kalman_algo_if_singularity_is_detected
          , dirname = ireland2004_ml_2_csn_initval_search_@{j}
          , cova_compute = 0
          );
oo_.initval_search.fval(@{j})          = oo_.posterior.optimization.log_density;
oo_.initval_search.shocks_stderr{@{j}} = oo_.mle_mode.shocks_std;
oo_.initval_search.shocks_skew{@{j}}   = oo_.mle_mode.shocks_skew;
format long;
disp(oo_.posterior.optimization.log_density);
disp(oo_.mle_mode.shocks_std);
disp(oo_.mle_mode.shocks_skew);
format short;
@#endfor
warning('on');