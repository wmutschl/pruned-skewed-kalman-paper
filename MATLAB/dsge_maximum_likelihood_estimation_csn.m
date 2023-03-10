%function M_ = dsge_maximum_likelihood_estimation_csn(M_,options_,datamat)
% -------------------------------------------------------------------------
% Maximum Likelihood Estimation of a DSGE model with CSN distributed innovations,
% solved with first-order perturbation techniques, and estimated using the
% Pruned Skewed Kalman filter.
% Note that the solution takes the form of a linear state-space system with
% CSN distributed innovations eta and normally distributed noise eps:
%   x(t) = gx*x(t-1) + gu*eta(t)   [state transition equation]
%   y(t) = F*x(t)    + eps(t)      [observation equation]
%   eta(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta) [innovations, shocks]
%   eps(t) ~ N(0,Sigma_eps) [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% Assumptions:
%   - all elements in eta are independent
%   - due to identifiability, we normalize nu_eta=0, Delta_eta=I
%   - mu_eta is endogenously determined to ensure that E[eta]=0
% -------------------------------------------------------------------------
% Note that we run the estimation in several stages in order to get good initial values:
% Stage 0: - use provided initial values in estimated_params_0 file
%          - run ML with Gaussian Kalman filter
% Stage 1: - fix model parameters to estimates from stage 0
%          - fix Var[eta_j] to estimates of stage 0
%          - create grid for Skew[eta_j]_i
%          - for each combination of Var[eta_j] and Skew[eta_j]_i recover corresponding Sigma_eta_j and Gamma_eta_j combinations
%          - compute negative log-likelihood of each value on grid
%          - for a chosen number of best combinations: use these as initial value and optimize shock parameters with Pruned Skewed Kalman filter
% Stage 2: - use best values from stage 1 as initial value for shock parameters
%          - use values from stage 0 as initial value for model parameters
%          - run optimization with Pruned Skewed Kalman filter to estimate model as well as shock parameters
%          - compute standard errors using the inverse hessian (with possibly fine-tuned numerical steps)
% -------------------------------------------------------------------------
% INPUTS
% - M_         [structure]   information on the model (stripped down version of Dynare's M_ structure)
% - options_   [structure]   options (stripped down version of Dynare's options_ structure)
% -------------------------------------------------------------------------
% OUTPUTS
% =========================================================================
% Copyright © 2023 Willi Mutschler
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
if ~isfield(options_.parameters,'transform')
    options_.parameters.transform = [];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 0
% - use provided initial values in estimated_params_ file
% - run ML with Gaussian Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STAGE = 0;
[estim_params0,xparams0Init,bounds0] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
M_ = set_params_dsge(xparams0Init,estim_params0,options_,M_);
tic_id0 = tic;
[neg_log_likelihood0Init,exit_flag0] = negative_log_likelihood(xparams0Init, bounds0, datamat, estim_params0, options_, M_);
elapsed_time0 = toc(tic_id0);
if exit_flag0 ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
else
    fprintf('Stage 0: initial value of the Gaussian log-likelihood function: %6.4f \nTime required to compute log-likelihood function once: %s \n', -1*neg_log_likelihood0Init,dynsec2hms(elapsed_time0));
end
xparams_stage0 = nan(estim_params0.ntot,length(options_.mode_compute)); neg_log_likelihood_stage0 = nan(1,length(options_.mode_compute));
% run maximum likelihood with Gaussian Kalman filter
parfor jopt=1:length(options_.mode_compute)
    if options_.mode_compute(jopt)=="fminsearch"
        [xparams_stage0(:,jopt), neg_log_likelihood_stage0(jopt)] = fminsearch(   @(x) negative_log_likelihood(x,bounds0,datamat,estim_params0,options_,M_),  xparams0Init,                             options_.optim_opt);
    elseif options_.mode_compute(jopt)=="fminsearchbnd"
        [xparams_stage0(:,jopt), neg_log_likelihood_stage0(jopt)] = fminsearchbnd(@(x) negative_log_likelihood(x,bounds0,datamat,estim_params0,options_,M_),  xparams0Init, bounds0(:,1), bounds0(:,2), options_.optim_opt);
    elseif options_.mode_compute(jopt)=="fminunc"
        [xparams_stage0(:,jopt), neg_log_likelihood_stage0(jopt)] = fminunc(      @(x) negative_log_likelihood(x,bounds0,datamat,estim_params0,options_,M_),  xparams0Init,                             options_.optim_opt);
    elseif options_.mode_compute(jopt)=="cmaes"
        cmaesOptions = cmaes('defaults');
        cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='500'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
        cmaesOptions.LBounds = bounds0(:,1); cmaesOptions.UBounds = bounds0(:,2);
        % Set default search volume (SIGMA)
        cmaesSIGMA = (bounds0(:,2)-bounds0(:,1))*0.2; cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
        while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
            cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
        end
        cmaesOptions.TolFun = options_.optim_opt.TolFun; cmaesOptions.MaxIter = options_.optim_opt.MaxIter; cmaesOptions.MaxFunEvals = options_.optim_opt.MaxFunEvals;
        [~, ~, ~, ~, ~, BESTEVER] = cmaes('negative_log_likelihood',xparams0Init,cmaesSIGMA,cmaesOptions,  bounds0,datamat,estim_params0,options_,M_);
        xparams_stage0(:,jopt)=BESTEVER.x;
        neg_log_likelihood_stage0(jopt)=BESTEVER.f;
    end
end
[~,stage0best] = sort(neg_log_likelihood_stage0);
disp(array2table([xparams_stage0(:,stage0best) xparams0Init; -neg_log_likelihood_stage0(stage0best) -neg_log_likelihood0Init],...
                 'RowNames',[estim_params0.names;'log-lik'],...
                 'VariableNames',[options_.mode_compute(stage0best) "Ireland"]))
xparams0Final = xparams_stage0(:,stage0best(1));
neg_log_likelihood0Final = neg_log_likelihood_stage0(stage0best(1));
M_ = set_params_dsge(xparams0Final,estim_params0,options_,M_); % update model and shock parameters for use in stage 1
V_stage0 = inv(reshape(get_hessian('negative_log_likelihood',xparams0Final,[1e-3;1.0],bounds0,datamat,estim_params0,options_,M_),estim_params0.ntot,estim_params0.ntot)); SE_stage0 = sqrt(diag(V_stage0)); % compute standard errors
fprintf('SUMMARY STAGE 0\n\n')
disp( table(xparams0Final,SE_stage0,xparams0Final./SE_stage0,'Var',{'Estimate','s.d.','t-stat'},'Row',estim_params0.names) );
disp(array2table([sqrt(diag(M_.Cov_eta)) M_.Skew_eta],'RowNames',M_.exo_names,'VariableNames',["stderr", "skew"]));
fprintf('Stage 0: final value of the Gaussian log-likelihood function: %6.4f\n', -1*neg_log_likelihood0Final);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 1
% - fix model parameters to estimates from stage 0
% - fix Var[eta_j] to estimates of stage 0
% - create grid for Skew[eta_j]_i
% - for each combination of Var[eta_j] and Skew[eta_j]_i recover corresponding Sigma_eta_j and Gamma_eta_j combinations
% - compute negative log-likelihood of each value on grid
% - for a chosen number of best combinations: use these as initial value and optimize over shock parameters with Pruned Skewed Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STAGE = 1;
[estim_params1,xparams1Init,bounds1] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_); % note that M_ stores the values from stage 0
M_ = set_params_dsge(xparams1Init,estim_params1,options_,M_);
[neg_log_likelihood1Init,exit_flag1] = negative_log_likelihood(xparams1Init, bounds1, datamat, estim_params1, options_, M_);
fprintf('Stage 1: initial value of the Pruned-Skewed log-likelihood function: %6.4f\n', -1*neg_log_likelihood1Init);
% create an evenly spaced grid of skewness coefficients for each estimated skewness or gamma parameter
grid.nbr      = 16;   % needs to be even number
grid.endpoint = 0.95; % grid is set between +- this value
grid.bestof   = 3;    % how many values to keep in stage 1
Skew_eta_grid = linspace(-abs(grid.endpoint),0,grid.nbr/2); Skew_eta_grid = [Skew_eta_grid -Skew_eta_grid((end-1):-1:1)];
if ~options_.parameters.use_stderr_skew
    % for each skewness coefficient compute required Sigma_eta_j and Gamma_eta_j such that V[eta_j] stays equal to previous stage
    diag_Sigma_eta_grid = zeros(estim_params1.nsx,grid.nbr-1);
    diag_Gamma_eta_grid = zeros(estim_params1.nsx,grid.nbr-1);
    for jexo = 1:M_.exo_nbr
        idx_exo = find(ismember(estim_params1.skew_exo(:,1),jexo));
        if ~isempty(idx_exo)
            for jgrid = 1:(grid.nbr-1)
                [diag_Sigma_eta_grid(idx_exo,jgrid),diag_Gamma_eta_grid(idx_exo,jgrid)] = csnVarSkew_To_SigmaGamma_univariate(M_.Cov_eta(jexo,jexo),Skew_eta_grid(jgrid),1);
            end
        end
    end
end
% compute negative loglikelihood for all possible combinations of grid values
neg_log_likelihood_grid = nan(1,(grid.nbr-1)^estim_params1.nsx); % initialize storage
COMBOS = create_combos(estim_params1.nsx,grid.nbr);
fprintf('\n\nRunning grid search for CSN skewness/gamma parameters:\n');
parfor_progress((grid.nbr-1)^estim_params1.nsx);
parfor jgrid = 1:(grid.nbr-1)^estim_params1.nsx
    M_j = M_;
    combo = COMBOS(jgrid,:);
    for jnsx = 1:estim_params1.nsx
        jexo = estim_params1.skew_exo(jnsx,1);
        if options_.parameters.use_stderr_skew
            M_j.Skew_eta(jexo,1) = Skew_eta_grid(combo(jnsx));
        else
            M_j.Sigma_eta(jexo,jexo) = diag_Sigma_eta_grid(jnsx,combo(jnsx));
            M_j.Gamma_eta(jexo,jexo) = diag_Gamma_eta_grid(jnsx,combo(jnsx));
        end
    end
    [estim_params1_j,xparams1_j,bounds1_j] = feval(str2func(M_j.fname + "_estim_params"), STAGE, M_j, options_);
    neg_log_likelihood_grid(jgrid) = negative_log_likelihood(xparams1_j, bounds1_j, datamat, estim_params1_j, options_, M_j);
    parfor_progress;
end
parfor_progress(0);

% get best values on grid
[~,idx_best_grid] = sort(neg_log_likelihood_grid);
tbl_stderr_grid = nan(M_.exo_nbr,grid.bestof); tbl_skew_grid = nan(M_.exo_nbr,grid.bestof); tbl_loglik_grid = nan(1,grid.bestof); % table for all shock parameters (even non-estimated)
for j=1:grid.bestof
    jgrid = idx_best_grid(j);  M_grid{j} = M_;  combo = COMBOS(jgrid,:);
    for jnsx = 1:estim_params1.nsx
        jexo = estim_params1.skew_exo(jnsx,1);
        if options_.parameters.use_stderr_skew
            M_grid{j}.Skew_eta(jexo,1) = Skew_eta_grid(combo(jnsx));
        else
            M_grid{j}.Sigma_eta(jexo,jexo) = diag_Sigma_eta_grid(jnsx,combo(jnsx));
            M_grid{j}.Gamma_eta(jexo,jexo) = diag_Gamma_eta_grid(jnsx,combo(jnsx));
        end
    end
    [estim_params_grid{j},xparams_grid{j},bounds_grid{j}] = feval(str2func(M_grid{j}.fname + "_estim_params"), STAGE, M_grid{j}, options_);
    [negative_log_likelihood_grid{j}, ~, M_grid{j}] = negative_log_likelihood(xparams_grid{j}, bounds_grid{j}, datamat, estim_params_grid{j}, options_, M_grid{j});
    tbl_stderr_grid(:,j) = sqrt(diag(M_grid{j}.Cov_eta));  tbl_skew_grid(:,j) = M_grid{j}.Skew_eta;  tbl_loglik_grid(:,j) = -1*negative_log_likelihood_grid{j};
end
disp(array2table(tbl_stderr_grid,'RowNames',"stderr "+M_.exo_names,'VariableNames',"best " + num2str(transpose(1:grid.bestof))));
disp(array2table(tbl_skew_grid,'RowNames',"skew "+M_.exo_names,'VariableNames',"best " + num2str(transpose(1:grid.bestof))));
disp(array2table(tbl_loglik_grid,'RowNames',"Log-Lik",'VariableNames',"best " + num2str(transpose(1:grid.bestof))));

% optimize over both Sigma_eta/stderr_eta as well as Gamma_eta/Skew_eta using best values on grid as initial point
xparams_stage1 = nan(estim_params1.ntot,grid.bestof,length(options_.mode_compute)) ; neg_log_likelihood_stage1 = nan(grid.bestof,length(options_.mode_compute));
fprintf('\nOptimize over both Sigma_eta/stderr_eta as well as Gamma_eta/Skew_eta using best values from grid search for CSN skewness/gamma parameter:\n');
% run maximum likelihood with Pruned Skewed Kalman filter
for jgrid=1:grid.bestof
    parfor jopt=1:length(options_.mode_compute)
        if options_.mode_compute(jopt)=="fminsearch"
            [xparams_stage1(:,jgrid,jopt), neg_log_likelihood_stage1(jgrid,jopt)] = fminsearch(   @(x) negative_log_likelihood(x,bounds_grid{jgrid},datamat,estim_params_grid{jgrid},options_,M_grid{jgrid}),  xparams_grid{jgrid},                                                   options_.optim_opt);
        elseif options_.mode_compute(jopt)=="fminsearchbnd"
            [xparams_stage1(:,jgrid,jopt), neg_log_likelihood_stage1(jgrid,jopt)] = fminsearchbnd(@(x) negative_log_likelihood(x,bounds_grid{jgrid},datamat,estim_params_grid{jgrid},options_,M_grid{jgrid}),  xparams_grid{jgrid}, bounds_grid{jgrid}(:,1), bounds_grid{jgrid}(:,2), options_.optim_opt);
        elseif options_.mode_compute(jopt)=="fminunc"
            [xparams_stage1(:,jgrid,jopt), neg_log_likelihood_stage1(jgrid,jopt)] = fminunc(      @(x) negative_log_likelihood(x,bounds_grid{jgrid},datamat,estim_params_grid{jgrid},options_,M_grid{jgrid}),  xparams_grid{jgrid},                                                   options_.optim_opt);
        elseif options_.mode_compute(jopt)=="cmaes"
            cmaesOptions = cmaes('defaults');
            cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='500'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
            cmaesOptions.LBounds = bounds_grid{jgrid}(:,1); cmaesOptions.UBounds = bounds_grid{jgrid}(:,2);
            % Set default search volume (SIGMA)
            cmaesSIGMA = (bounds_grid{jgrid}(:,2)-bounds_grid{jgrid}(:,1))*0.2; cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
            while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
                cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
            end
            cmaesOptions.TolFun = options_.optim_opt.TolFun; cmaesOptions.MaxIter = options_.optim_opt.MaxIter; cmaesOptions.MaxFunEvals = options_.optim_opt.MaxFunEvals;
            [~, ~, ~, ~, ~, BESTEVER] = cmaes('negative_log_likelihood',xparams_grid{jgrid},cmaesSIGMA,cmaesOptions,  bounds_grid{jgrid},datamat,estim_params_grid{jgrid},options_,M_grid{jgrid});
            xparams_stage1(:,jgrid,jopt)=BESTEVER.x; neg_log_likelihood_stage1(jgrid,jopt)=BESTEVER.f;
        end
    end
end
strStage1 = [];
for jgrid=1:grid.bestof
    for jopt=1:length(options_.mode_compute)
        strStage1 = [strStage1 (options_.mode_compute(jopt)+"_init"+num2str(jgrid)) ];
    end
end
xparams_stage1 = reshape(xparams_stage1,estim_params1.ntot,grid.bestof*length(options_.mode_compute));
neg_log_likelihood_stage1 = reshape(neg_log_likelihood_stage1,1,grid.bestof*length(options_.mode_compute));
[~,stage1best] = sort(neg_log_likelihood_stage1(:));
disp(array2table([xparams_stage1(:,stage1best); -neg_log_likelihood_stage1(stage1best)],...
                 'RowNames',[estim_params1.names;'log-lik'],...
                 'VariableNames',strStage1(stage1best)));
xparams1Final = xparams_stage1(:,stage1best(1));
neg_log_likelihood1Final = neg_log_likelihood_stage1(stage1best(1));
M_ = set_params_dsge(xparams1Final,estim_params1,options_,M_);
V_stage1 = inv(reshape(get_hessian('negative_log_likelihood',xparams1Final,[1e-3;1.0],bounds1,datamat,estim_params1,options_,M_),estim_params1.ntot,estim_params1.ntot)); SE_stage1 = sqrt(diag(V_stage1)); % compute standard errors
fprintf('SUMMARY STAGE 1\n\n')
disp( table(xparams_stage1(:,stage1best(1)),SE_stage1,xparams_stage1(:,stage1best(1))./SE_stage1,'Var',{'Estimate','s.d.','t-stat'},'Row',estim_params1.names) );
disp(array2table([sqrt(diag(M_.Cov_eta)) M_.Skew_eta],'RowNames',M_.exo_names,'VariableNames',["stderr", "skew"]));
fprintf('Stage 1: final value of the Pruned-Skewed log-likelihood function: %6.4f\n', -1*neg_log_likelihood1Final)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2
% - use best values from stage 1 as initial value for shock parameters
% - use values from stage 0 as initial value for model parameters
% - run optimization with Pruned Skewed Kalman filter to estimate model as well as shock parameters
% - compute standard errors using the inverse hessian (with possibly fine-tuned numerical steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STAGE = 2;
[estim_params2,xparams2Init,bounds2] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
M_ = set_params_dsge(xparams2Init,estim_params2,options_,M_);
[neg_log_likelihood2Init,exit_flag2] = negative_log_likelihood(xparams2Init, bounds2, datamat, estim_params2, options_, M_);
fprintf('\n\nStage 2: initial value of the Pruned-Skewed log-likelihood function: %6.4f\n', -1*neg_log_likelihood2Init);
xparams_stage2 = nan(estim_params2.ntot,length(options_.mode_compute)); neg_log_likelihood_stage2 = nan(1,length(options_.mode_compute));
fprintf('\n\n run maximum likelihood estimation with Pruned Skewed Kalman filter for all parameters using different optimizers in parallel:\n')
parfor jopt=1:length(options_.mode_compute)
    if options_.mode_compute(jopt)=="fminsearch"
        [xparams_stage2(:,jopt), neg_log_likelihood_stage2(jopt)] = fminsearch(   @(x) negative_log_likelihood(x,bounds2,datamat,estim_params2,options_,M_),  xparams2Init,                             options_.optim_opt);
    elseif options_.mode_compute(jopt)=="fminsearchbnd"
        [xparams_stage2(:,jopt), neg_log_likelihood_stage2(jopt)] = fminsearchbnd(@(x) negative_log_likelihood(x,bounds2,datamat,estim_params2,options_,M_),  xparams2Init, bounds2(:,1), bounds2(:,2), options_.optim_opt);
    elseif options_.mode_compute(jopt)=="fminunc"
        [xparams_stage2(:,jopt), neg_log_likelihood_stage2(jopt)] = fminunc(      @(x) negative_log_likelihood(x,bounds2,datamat,estim_params2,options_,M_),  xparams2Init,                             options_.optim_opt);
    elseif options_.mode_compute(jopt)=="cmaes"
        cmaesOptions = cmaes('defaults');
        cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='500'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
        cmaesOptions.LBounds = bounds2(:,1); cmaesOptions.UBounds = bounds2(:,2);
        % Set default search volume (SIGMA)
        cmaesSIGMA = (bounds2(:,2)-bounds2(:,1))*0.2; cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
        while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
            cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
        end
        cmaesOptions.TolFun = options_.optim_opt.TolFun; cmaesOptions.MaxIter = options_.optim_opt.MaxIter; cmaesOptions.MaxFunEvals = options_.optim_opt.MaxFunEvals;
        [~, ~, ~, ~, ~, BESTEVER] = cmaes('negative_log_likelihood',xparams2Init,cmaesSIGMA,cmaesOptions,  bounds2,datamat,estim_params2,options_,M_);
        xparams_stage2(:,jopt)=BESTEVER.x; neg_log_likelihood_stage2(jopt)=BESTEVER.f;
    end
end
[~,stage2best] = sort(neg_log_likelihood_stage2);
disp(array2table([xparams_stage2(:,stage2best); -neg_log_likelihood_stage2(stage2best)],...
                 'RowNames',[estim_params2.names;'log-lik'],...
                 'VariableNames',options_.mode_compute(stage2best)));
xparams2Final = xparams_stage2(:,stage2best(1));
neg_log_likelihood2Final = neg_log_likelihood_stage2(stage2best(1));
M_ = set_params_dsge(xparams2Final,estim_params2,options_,M_); % update parameters
V_stage2 = inv(reshape(get_hessian('negative_log_likelihood',xparams2Final,[1e-3;1.0],bounds2,datamat,estim_params2,options_,M_),estim_params2.ntot,estim_params2.ntot)); SE_stage2 = sqrt(diag(V_stage2)); % compute standard errors
fprintf('SUMMARY STAGE 2\n\n')
disp( table(xparams2Final,SE_stage2,xparams2Final./SE_stage2,'Var',{'Estimate','s.d.','t-stat'},'Row',estim_params2.names) );
disp(array2table(M_.params,'RowNames',M_.param_names,'VariableNames',["final"]));
disp(array2table([sqrt(diag(M_.Cov_eta)) M_.Skew_eta],'RowNames',M_.exo_names,'VariableNames',["stderr", "skew"]));
fprintf('Final value of the Pruned Skewed log-likelihood function: %6.4f\n', -1*neg_log_likelihood2Final)

lr_stat = 2 * ( neg_log_likelihood0Final - neg_log_likelihood2Final);
lr_pval = chi2cdf(lr_stat, 4, "upper");
fprintf('Likelihood Ratio Test Statistic = %.2f with p-val = %.4f\n',lr_stat,lr_pval);