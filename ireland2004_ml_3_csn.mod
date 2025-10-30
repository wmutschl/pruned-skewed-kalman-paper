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
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMUM LIKELIHOOD ESTIMATION OF CSN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
% taken from ireland2004_ml_1_gaussian.log
OMEGA,        0.058086387283355,        0,     1;
RHO_PI,       0.386477141360770,        0,     1;
RHO_G,        0.396013670264708,        0,     1;
RHO_X,        0.165402611432957,        0,     1;
RHO_A,        0.904795669581000,        0,     1;
RHO_E,        0.990673584292426,        0,     1;
% taken from ireland2004_ml_2_csn_initval_search.log
stderr eta_a, 2.747697233134916,        0,    10;
stderr eta_e, 0.028882164116404,        0,    10;
stderr eta_z, 0.858196786024231,        0,    10;
stderr eta_r, 0.276196004637038,        0,    10;
skew eta_a,  -0.111987050675084;
skew eta_e,  -0.289566734765780;
skew eta_z,  -0.836095597248922;
skew eta_r,   0.520016288312854;
end;

options_.kalman.pskf.prune_tol = 0.01;
options_.kalman.pskf.rank_deficiency_transform = false;
options_.kalman.pskf.skip_smoother = false;
estimation(datafile = 'data/ireland2004_data.m'
          , mode_compute = 5 % fval=1215.850522
          , additional_optimizer_steps = [8] % fval=1215.852847
          , silent_optimizer % below we display optimization_info, so don't show intermediate optimization output
          , kalman_algo = 5  % use pruned skewed Kalman filter
          , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          );
fprintf('Optimization info (stage 1, mode_compute=5):\n');
disp(oo_.posterior.optimization.optimization_info.stage_1);
fprintf('Optimization info (stage 2, mode_compute=8):\n');
disp(oo_.posterior.optimization.optimization_info.stage_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE ESTIMATES SHOCK PARAMETERS AND SMOOTHED SHOCK VALUES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lines = findall(figure(1), 'Type', 'line');
eta_r_t_T = lines(1).YData;
eta_z_t_T = lines(3).YData;
eta_e_t_T = lines(5).YData;
eta_a_t_T = lines(7).YData;
csn = M_.csn;
save([M_.dname filesep 'Output' filesep M_.fname '_shock_params'], 'csn', 'eta_*_t_T', '-v6');

%%%%%%%%%%%%%%%%%%%%%
% ONE-SIDED HESSIAN %
%%%%%%%%%%%%%%%%%%%%%
% skew eta_z is at the theoretical bound, so we compute standard errors using one-sided finite differences, similar to what Ireland (2004) did for ALPHA_X and ALPHA_PI
[dataset, datasetInfo, ~, ~, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options_.varobs, M_.dname, [], M_, options_, oo_, estim_params_, bayestopt_);
xparam1 = oo_.posterior.optimization.mode;
neg_log_likelihood = dsge_likelihood(xparam1,dataset,datasetInfo,options,M,estim_params,bayestopt,BoundsInfo,oo.dr, oo.steady_state,oo.exo_steady_state,oo.exo_det_steady_state);
if ~isequal(neg_log_likelihood, -oo_.posterior.optimization.log_density)
    error('one-sided Hessian: something wrong in dsge_likelihood')
end
addpath('MATLAB/external')
hh = hessian_with_bounds('dsge_likelihood', xparam1, options.gstep, BoundsInfo, ...
                         dataset, datasetInfo, options, M, estim_params, bayestopt, BoundsInfo, oo.dr, oo.steady_state, oo.exo_steady_state, oo.exo_det_steady_state);
rmpath('MATLAB/external');
hh = reshape(hh, length(xparam1), length(xparam1));
hsd = sqrt(diag(hh));
invhess = inv(hh./(hsd*hsd'))./(hsd*hsd');
stdh = sqrt(diag(invhess));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY ESTIMATION RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_names = [["" + fieldnames(oo_.mle_mode.parameters)]; ["stderr " + fieldnames(oo_.mle_mode.shocks_std)]; ["skew " + fieldnames(oo_.mle_mode.shocks_skew)]];
tbl_vals = [struct2array(oo_.mle_mode.parameters) struct2array(oo_.mle_mode.shocks_std) struct2array(oo_.mle_mode.shocks_skew)]';
tbl_vals = [tbl_vals stdh];
fprintf('\n%s\n* ESTIMATION RESULTS *\n%s\n', repmat('*',1,23), repmat('*',1,23));
disp(array2table(tbl_vals, 'RowNames', tbl_names, 'VariableNames', ["Mode", "Std-dev"]))

%%%%%%%%%%%%%%%%%%%%%%%%%
% LIKELIHOOD RATIO TEST %
%%%%%%%%%%%%%%%%%%%%%%%%%
oo_gauss = load('ireland2004_ml_1_gaussian/Output/ireland2004_ml_1_gaussian_results.mat','oo_');
loglik_csn = oo_.posterior.optimization.log_density;            % 1215.85
loglik_gauss = oo_gauss.oo_.posterior.optimization.log_density; % 1207.56
lr_stat = 2 * (loglik_csn - loglik_gauss);
lr_pval = chi2cdf(lr_stat, 4, "upper"); % 4 additional parameters with CSN
lr_str = sprintf('test statistic of $%.2f$ and a p-value of $%.4f$.', lr_stat, lr_pval);
fprintf('\n%s\n* LIKELIHOOD RATIO TEST *\n%s\n', repmat('*',1,25), repmat('*',1,25));
fprintf('- unconstrained log-likelihood (CSN): %.4f\n', loglik_csn);
fprintf('- constrained log-likelihood (Gaussian): %.4f\n', loglik_gauss);
fprintf('- %s\n\n',lr_str);