% Monte-Carlo study of the model of Ireland (2004) with CSN distributed shocks
% =========================================================================
% Copyright (C) 2025 Willi Mutschler
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
@#define DATASETS_NBR     = 1200
@#define DATASETS_LENGTHS = [100, 500]
@#define DATASETS_BURNIN  = 500

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
var eta_a; stderr 3.0167;
var eta_e; stderr 0.0248;
var eta_z; stderr 0.8865;
var eta_r; stderr 0.2790;
@#ifndef GAUSSIAN
skew eta_a = -0.3;
skew eta_e = -0;
skew eta_z = -0.5;
skew eta_r = +0.8;
@#endif
end;

%%%%%%%%%%%%%%%%%
% SIMULATE DATA %
%%%%%%%%%%%%%%%%%
%set_dynare_seed('clock'); % uncomment for random draws
%set_dynare_seed(75);
@#for jlength in DATASETS_LENGTHS
fprintf('\nSIMULATE @{DATASETS_NBR} DATASETS WITH T=@{jlength}...');
t_sim = tic;
  @#for jdat in 1:DATASETS_NBR
stoch_simul(order=1,periods=@{jlength+DATASETS_BURNIN},irf=0,nodecomposition,nomoments,nocorr,nofunctions,nomodelsummary) ghat rhat pihat;
fprintf('\b');
oo_.endo_simul = oo_.endo_simul(:,@{DATASETS_BURNIN+1}:end);
datatomfile('data/ireland2004_sim_data_T@{jlength}_@{jdat}',{'ghat', 'rhat', 'pihat'});
  @#endfor
fprintf('...DONE (%s)!\n',dynsec2hms(toc(t_sim)));
@#endfor

%%%%%%%%%%%%%%
% ESTIMATION %
%%%%%%%%%%%%%%
estimated_params;
stderr eta_a, 3.0167;
stderr eta_e, 0.0248;
stderr eta_z, 0.8865;
stderr eta_r, 0.2790;
@#ifndef GAUSSIAN
skew eta_a, -0.3;
skew eta_e, -0;
skew eta_z, -0.5;
skew eta_r, +0.8;
@#endif
end;

options_.kalman.pskf.prune_tol = 0.01;
options_.kalman.pskf.rank_deficiency_transform = false;
options_.kalman.pskf.skip_smoother = true;
% run estimation command without optimization to initialize structures
estimation(datafile = 'data/ireland2004_data.m'
         , mode_compute = 0
         , silent_optimizer % below we display optimization_info, so don't show intermediate optimization output
         , kalman_algo = 5  % use pruned skewed Kalman filter
         , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
         , cova_compute = 0
         );

% optimize likelihood for grid values
fprintf('Note: You can watch the progress of the estimation in the data folder:\n      At the end there will be no ireland2004_sim_data_T*_*.m files!\n')
@#for jlength in DATASETS_LENGTHS
fprintf('\nESTIMATE @{DATASETS_NBR} DATASETS WITH T=@{jlength}...');
t_opt = tic;
xparam1_opt_T@{jlength} = nan(estim_params_.nvx + estim_params_.nsx, @{DATASETS_NBR});
current_optimizer = 5;
parfor jdat = 1:@{DATASETS_NBR}
    % run dynare_estimation_init twice so parfor loop works
    [dataset, datasetInfo, xparam1, hh, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options_.varobs, M_.dname, [], M_, options_, oo_, estim_params_, bayestopt_);
    options.datafile = sprintf('data/ireland2004_sim_data_T@{jlength}_%u.m',jdat);
    options.nobs = @{jlength};
    [dataset, datasetInfo, xparam1, hh, M, options, oo, estim_params, bayestopt, BoundsInfo] = dynare_estimation_init(options.varobs, M.dname, [], M, options, oo, estim_params, bayestopt);
    % run optimization
    warning('off');
    xparam1_opt_T@{jlength}(:,jdat) = dynare_minimize_objective('dsge_likelihood',xparam1,current_optimizer,options,[BoundsInfo.lb BoundsInfo.ub],bayestopt.name,bayestopt,hh, ...
                                                                dataset,datasetInfo,options,M,estim_params,bayestopt,BoundsInfo,oo.dr, oo.steady_state,oo.exo_steady_state,oo.exo_det_steady_state);
    % delete dataset so we can watch progress of the estimation
    delete(options.datafile);
end
fclose('all');
fprintf('...DONE (%s)!\n',dynsec2hms(toc(t_opt)));

figure('name','Histogram T=@{jlength} stderr'); sgtitle('Histogram T=@{jlength} stderr');
for j = 1:M_.exo_nbr
    subplot(2,2,j);
    histfit(xparam1_opt_T@{jlength}(j,:));
    title(sprintf('stderr %s', bayestopt_.name{j}))
end

@#ifndef GAUSSIAN
figure('name','Histogram T=@{jlength} skew'); sgtitle('Histogram T=@{jlength} skew');
for j = 1:M_.exo_nbr
    subplot(2,2,j);
    histfit(xparam1_opt_T@{jlength}(4+j,:));
    title(sprintf('skew %s', bayestopt_.name{j}))
end
@#endif

@#endfor

save([M_.dname filesep 'Output' filesep M_.fname '_xparam_opt'], ...
    @#for jlength in DATASETS_LENGTHS
    'xparam1_opt_T@{jlength}', ...
    @#endfor
    '-v6');