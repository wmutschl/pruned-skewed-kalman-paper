% Monte-Carlo study of the model of Ireland (2004) with CSN distributed shocks
% =========================================================================
% Copyright (C) 2025-2026 Willi Mutschler
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
addpath('_utils');

@#define DATASETS_NBR     = 1000
@#define DATASETS_LENGTHS = [200, 500]

@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%
% TRUE VALUES %
%%%%%%%%%%%%%%%
ALPHA_X  = 0;
ALPHA_PI = 0;
OMEGA    = 0.0581;
RHO_PI   = 0.3865;
RHO_G    = 0.3960;
RHO_X    = 0.1654;
RHO_A    = 0.9048;
RHO_E    = 0.9907;

shocks;
var eta_a; stderr 2.5232;
var eta_e; stderr 0.0212;
var eta_z; stderr 0.7900;
var eta_r; stderr 0.2838;
@#ifndef GAUSSIAN
skew eta_a = 0;
skew eta_e = -0.40;
skew eta_z = -0.95;
skew eta_r = +0.55;
@#endif
end;

%%%%%%%%%%%%%%%%%
% SIMULATE DATA %
%%%%%%%%%%%%%%%%%
% initialize simulation structures and solve linear state space system
set_dynare_seed(0);
stoch_simul(order=1,periods=0,irf=0,nodecomposition,nomoments,nocorr,nofunctions,nomodelsummary) ghat rhat pihat;
M_.csn = csn_update_specification(M_.Sigma_e, M_.Skew_e);
@#for jlength in DATASETS_LENGTHS
fprintf('SIMULATE @{DATASETS_NBR} DATASETS WITH T=@{jlength}...');
for jdat = 1:@{DATASETS_NBR}
    set_dynare_seed(jdat);
    [exo_simul, diagnostics, best_draws, best_diag] = generate_exo_shocks_matched(@{jlength}, M_.Sigma_e, M_.Skew_e, M_.csn);
    y_ = simult_(M_,options_,oo_.dr.ys,oo_.dr,exo_simul',1);
    oo_.endo_simul = y_(:,2:end);
    datatomfile(sprintf('../data/ireland2004_sim_data_T@{jlength}_%u.m',jdat),{'ghat', 'rhat', 'pihat'});
end
fprintf('DONE!\n');
@#endfor

@#ifdef GAUSSIAN
dist_label = 'gaussian';
@#else
dist_label = 'csn';
@#endif

%%%%%%%%%%%%%%
% ESTIMATION %
%%%%%%%%%%%%%%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end
% declare estimated parameters and initial values
estimated_params;
stderr eta_a, 2.5232;
stderr eta_e, 0.0212;
stderr eta_z, 0.7900;
stderr eta_r, 0.2838;
@#ifndef GAUSSIAN
skew eta_a, 0;
skew eta_e, -0.40;
skew eta_z, -0.80;
skew eta_r, +0.55;
@#endif
end;

% run estimation command without optimization to initialize estimation structures
estimation(datafile = '../data/ireland2004_data.m'
         , mode_compute = 0
         , silent_optimizer % below we display optimization_info, so don't show intermediate optimization output
         , kalman_algo = 5  % use pruned skewed Kalman filter
         , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
         , cova_compute = 0
         , frequentist_smoother = false
         );

% optimize likelihood for grid values
fprintf('\n');
@#for jlength in DATASETS_LENGTHS
fprintf('ESTIMATE @{DATASETS_NBR} DATASETS WITH T=@{jlength}...');
t_opt = tic;
xparam1_opt_T@{jlength} = nan(estim_params_.nvx + estim_params_.nsx, @{DATASETS_NBR});
exitflag_T@{jlength} = nan(1, @{DATASETS_NBR});
current_optimizer = 1;
q = parallel.pool.DataQueue;
h_waitbar = waitbar(0, sprintf('Estimating T=@{jlength}: 0/%d', @{DATASETS_NBR}), 'Name', 'Monte Carlo Progress T=@{jlength}');
t_start_parfor = tic;
afterEach(q, @(~) update_parfor_progress(@{DATASETS_NBR}, h_waitbar, t_start_parfor));
parfor jdat = 1:@{DATASETS_NBR}
    % run dynare_estimation_init twice so parfor loop works using local structures
    [dataset, datasetInfo, xparam1, hh, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options_.varobs, M_.dname, [], M_, options_, oo_, estim_params_, bayestopt_);
    options.datafile = sprintf('../data/ireland2004_sim_data_T@{jlength}_%u.m',jdat);
    options.nobs = @{jlength};
    [dataset, datasetInfo, xparam1, hh, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options.varobs, M.dname, [], M, options, oo, estim_params, bayestopt);
    % run actual optimization
    [xparam1_opt_T@{jlength}(:,jdat), ~, exitflag_T@{jlength}(jdat), ~, ~, ~, ~] = dynare_minimize_objective('dsge_likelihood',xparam1,current_optimizer,options,[BoundsInfo.lb BoundsInfo.ub],bayestopt.name,bayestopt,hh, ...
                                                                                               dataset,datasetInfo,options,M,estim_params,bayestopt,BoundsInfo,oo.dr, oo.steady_state,oo.exo_steady_state,oo.exo_det_steady_state);
    % delete dataset so we can watch progress of the estimation in the ../data folder
    delete(options.datafile);
    send(q, jdat);
end
delete(q);
close(h_waitbar);
clear update_parfor_progress;
fclose('all');
fprintf('DONE (%s)!\n',dynsec2hms(toc(t_opt)));

figure('name','Histogram T=@{jlength} stderr'); sgtitle('Histogram T=@{jlength} stderr');
for j = 1:M_.exo_nbr
    subplot(2,2,j);
    histfit(xparam1_opt_T@{jlength}(j,:));
    title(bayestopt_.name{j});
end

@#ifndef GAUSSIAN
figure('name','Histogram T=@{jlength} skew'); sgtitle('Histogram T=@{jlength} skew');
for j = 1:M_.exo_nbr
    subplot(2,2,j);
    histfit(xparam1_opt_T@{jlength}(4+j,:));
    title(bayestopt_.name{4+j});
end
@#endif

@#endfor

% save results in tidy CSV format
nparam = estim_params_.nvx + estim_params_.nsx;
nrows = @{DATASETS_NBR} * nparam * length([@{DATASETS_LENGTHS}]);
tidy_dataset_id = nan(nrows, 1);
tidy_T          = nan(nrows, 1);
tidy_param_type = cell(nrows, 1);
tidy_shock      = cell(nrows, 1);
tidy_value      = nan(nrows, 1);
tidy_exitflag   = nan(nrows, 1);
row = 0;
@#for jlength in DATASETS_LENGTHS
for jdat = 1:@{DATASETS_NBR}
    for jparam = 1:nparam
        row = row + 1;
        tidy_dataset_id(row) = jdat;
        tidy_T(row) = @{jlength};
        if jparam <= estim_params_.nvx
            tidy_param_type{row} = 'stderr';
        else
            tidy_param_type{row} = 'skew';
        end
        tidy_shock{row} = regexprep(bayestopt_.name{jparam}, '^(stderr|skew)\s+', '');
        tidy_value(row) = xparam1_opt_T@{jlength}(jparam, jdat);
        tidy_exitflag(row) = exitflag_T@{jlength}(jdat);
    end
end
@#endfor
T_results = table(tidy_dataset_id, tidy_T, tidy_param_type, tidy_shock, tidy_value, tidy_exitflag, ...
    'VariableNames', {'dataset_id', 'T', 'parameter_type', 'shock', 'value', 'exitflag'});
writetable(T_results, sprintf('../results/ireland2004/montecarlo/xparam_opt_%s_%s_%s.csv', dist_label, ARCH, MATLAB_VERSION));

% true parameter values
true_stderr = [2.5232; 0.0212; 0.7900; 0.2838];
@#ifndef GAUSSIAN
true_skew = [0; -0.40; -0.95; +0.55];
true_values = [true_stderr; true_skew];
@#else
true_values = true_stderr;
@#endif

% print summary table to console
fprintf('\n=== MONTE CARLO SUMMARY (%s) ===\n', upper(dist_label));
fprintf('%-20s %8s', 'Parameter', 'Truth');
@#for jlength in DATASETS_LENGTHS
fprintf(' %10s %15s %10s', 'Mean_@{jlength}', '[p5;p95]_@{jlength}', 'NRMSE_@{jlength}');
@#endfor
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 20 + 8 + length([@{DATASETS_LENGTHS}])*(10+15+10+3)));
for jparam = 1:nparam
    fprintf('%-20s %8.4f', bayestopt_.name{jparam}, true_values(jparam));
    @#for jlength in DATASETS_LENGTHS
    vals = xparam1_opt_T@{jlength}(jparam,:);
    avg = mean(vals); p5 = prctile(vals,5); p95 = prctile(vals,95);
    if true_values(jparam) ~= 0
        nrmse = sqrt(mean((vals - true_values(jparam)).^2)) / abs(true_values(jparam));
        fprintf(' %10.4f [%5.2f;%5.2f] %10.4f', avg, p5, p95, nrmse);
    else
        fprintf(' %10.4f [%5.2f;%5.2f] %10s', avg, p5, p95, 'N/A');
    end
    @#endfor
    fprintf('\n');
end
fprintf('\n');

% generate LaTeX table
param_latex = cell(nparam, 1);
for jparam = 1:nparam
    pname = bayestopt_.name{jparam};
    % convert "stderr eta_a" -> "\(stderr(\eta_a)\)" and "skew eta_a" -> "\(skew(\eta_a)\)"
    tokens = strsplit(pname);
    param_latex{jparam} = sprintf('\\(%s(\\%s)\\)', tokens{1}, tokens{2});
end
xparam_results = cell(1, length([@{DATASETS_LENGTHS}]));
jT = 0;
@#for jlength in DATASETS_LENGTHS
jT = jT + 1;
xparam_results{jT} = xparam1_opt_T@{jlength};
@#endfor
update_montecarlo_table('../results/ireland2004/montecarlo', dist_label, ...
    true_values, param_latex, [@{DATASETS_LENGTHS}], ARCH, MATLAB_VERSION, ...
    xparam_results{:});

rmpath('_utils');

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s_%s.log', M_.fname, dist_label, ARCH, MATLAB_VERSION);
