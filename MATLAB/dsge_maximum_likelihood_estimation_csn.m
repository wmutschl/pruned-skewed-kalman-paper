function [oo_,M_] = dsge_maximum_likelihood_estimation_csn(M_,options_,datamat)
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
% - oo_        [structure]   estimation results and structures of different stages
% - M_         [structure]   updated information on the model (stripped down version of Dynare's M_ structure)
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
% set missing options
if ~isfield(options_.parameters,'transform')
    options_.parameters.transform = [];
end
% parallel environment
options_.poolobj = gcp('nocreate');
if isempty(options_.poolobj)
    options_.physical_cores_nbr=feature('numcores');
    parpool(options_.physical_cores_nbr);
else
    options_.physical_cores_nbr = options_.poolobj.NumWorkers;
end
% create results folder
if ~exist('results', 'dir')
    mkdir('results');
end
% open log file
diary off
options_.logfile = ['results/' options_.filename,'.log'];
if exist(options_.logfile, 'file')
    delete(options_.logfile)
end
diary(options_.logfile)
options_
options_.optim_opt
options_.kalman
options_.kalman.csn
options_.parameters
options_.parameters.transform

options_STDERR = options_; options_STDERR.parameters.transform = []; % we compute standard errors on untransformed parameters
objfct = @negative_log_likelihood;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMUM LIKELIHOOD WITH GAUSSIAN KALMAN FILTER
% - use provided initial values in estimated_params_ file
% - run ML with Gaussian Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STAGE = "gaussian_initval_from_ireland2004_paper";
[estim_params_gauss,xparamsInit_gauss_tr,bounds_gauss_tr] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
M_ = dsge_set_params(xparamsInit_gauss_tr, estim_params_gauss, options_, M_);
% run maximum likelihood with Gaussian Kalman filter
diary off
[xparams_gauss_tr, neg_log_likelihood_gauss, best_gauss] = minimize_objective_in_parallel(objfct, xparamsInit_gauss_tr, bounds_gauss_tr(:,1), bounds_gauss_tr(:,2), options_.optim_names, options_.optim_opt,    bounds_gauss_tr,datamat,estim_params_gauss,options_,M_);
diary(options_.logfile)
% compute standard errors
[xparams_gauss, bounds_gauss] = dsge_untransform(xparams_gauss_tr, bounds_gauss_tr, estim_params_gauss);
xstderr_gauss = standard_errors_inverse_hessian(objfct, xparams_gauss, bounds_gauss,    bounds_gauss,datamat,estim_params_gauss,options_STDERR,M_);
% display summary
fprintf('%s\nSUMMARY GAUSSIAN\n\n',repmat('*',1,100))
fprintf('  Point Estimates\n')
disp(array2table([xparams_gauss(:,best_gauss); -neg_log_likelihood_gauss(best_gauss)], 'RowNames',[erase(estim_params_gauss.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_gauss)));
fprintf('  Standard Errors\n')
disp(array2table([xstderr_gauss(:,best_gauss); -neg_log_likelihood_gauss(best_gauss)], 'RowNames',[erase(estim_params_gauss.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_gauss)));
fprintf('Final value of the best Gaussian log-likelihood function: %6.4f\n', -1*neg_log_likelihood_gauss(best_gauss(1)))
fprintf('%s\n\n',repmat('*',1,100))
% store results
oo_.gauss.xparams_tr = xparams_gauss_tr;
oo_.gauss.xparams = xparams_gauss;
oo_.gauss.xstderr = xstderr_gauss;
oo_.gauss.neg_log_likelihood = neg_log_likelihood_gauss;
% update model and shock parameters
M_ = dsge_set_params(xparams_gauss_tr(:,best_gauss(1)),estim_params_gauss,options_,M_);

if isfield(options_.kalman.csn,'initval_search') && (options_.kalman.csn.initval_search==1)
%% INIITAL VALUE SEARCH
    % store old tolerance
    fprintf('INITIAL VALUES: For initial value search, we use a pruning threshold of 10%%. This will be re-set at the final stage.\n')
    old_prune_tol = options_.kalman.csn.prune_tol;
    options_.kalman.csn.prune_tol = 0.1;
    options_STDERR.csn.prune_tol = 0.1;
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIITAL VALUE SEARCH STAGE 1
    % - fix model parameters to previous estimates
    % - fix Var[eta_j] to previous estimates
    % - create grid for Skew[eta_j]_i
    % - for each combination of Var[eta_j] and Skew[eta_j]_i recover corresponding Sigma_eta_j and Gamma_eta_j combinations
    % - compute negative log-likelihood of each value on grid
    % - for a chosen number of best combinations: use these as initial value and optimize over shock parameters with Pruned Skewed Kalman filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAGE = "csn_shock_params";
    [estim_params_1, xparamsInit_stage1_tr, bounds1_tr] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
    M_ = dsge_set_params(xparamsInit_stage1_tr, estim_params_1, options_, M_);
    
    % create an evenly spaced grid of skewness coefficients for each estimated skewness or gamma parameter
    grid.nbr      = 16;   % needs to be even number
    grid.endpoint = 0.95; % grid is set between +- this value
    grid.bestof   = 3;    % how many values to keep in stage 1
    Skew_eta_grid = linspace(-abs(grid.endpoint),0,grid.nbr/2); Skew_eta_grid = [Skew_eta_grid -Skew_eta_grid((end-1):-1:1)];
    if ~options_.parameters.use_stderr_skew
        % for each skewness coefficient compute required Sigma_eta_j and Gamma_eta_j such that V[eta_j] stays equal to previous estimates
        diag_Sigma_eta_grid = zeros(estim_params_1.nsx,grid.nbr-1);
        diag_Gamma_eta_grid = zeros(estim_params_1.nsx,grid.nbr-1);
        for jexo = 1:M_.exo_nbr
            idx_exo = find(ismember(estim_params_1.skew_exo(:,1),jexo));
            if ~isempty(idx_exo)
                for jgrid = 1:(grid.nbr-1)
                    [diag_Sigma_eta_grid(idx_exo,jgrid),diag_Gamma_eta_grid(idx_exo,jgrid)] = csnVarSkew_To_SigmaGamma_univariate(M_.Cov_eta(jexo,jexo),Skew_eta_grid(jgrid),1);
                end
            end
        end
    else
        diag_Sigma_eta_grid = [];
        diag_Gamma_eta_grid = [];
    end
    % compute negative loglikelihood for all possible combinations of grid values
    neg_log_likelihood_grid = nan(1,(grid.nbr-1)^estim_params_1.nsx); % initialize storage
    COMBOS = create_combos(estim_params_1.nsx,grid.nbr);
    fprintf('\n\nINITIAL VALUES: Running grid search for CSN skewness/gamma parameters:\n');
    parfor_progress((grid.nbr-1)^estim_params_1.nsx);
    diary off;
    parfor jgrid = 1:(grid.nbr-1)^estim_params_1.nsx
        M_j = M_;
        combo = COMBOS(jgrid,:);
        for jnsx = 1:estim_params_1.nsx
            jexo = estim_params_1.skew_exo(jnsx,1);
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
    diary(options_.logfile)
    % get best values on grid
    [~,idx_best_grid] = sort(neg_log_likelihood_grid);
    tbl_stderr_grid = nan(M_.exo_nbr,grid.bestof); tbl_skew_grid = nan(M_.exo_nbr,grid.bestof); tbl_loglik_grid = nan(1,grid.bestof); % table for all shock parameters (even non-estimated)
    for j=1:grid.bestof
        jgrid = idx_best_grid(j);  M_grid{j} = M_;  combo = COMBOS(jgrid,:);
        for jnsx = 1:estim_params_1.nsx
            jexo = estim_params_1.skew_exo(jnsx,1);
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
    % display results on grid
    disp(array2table(tbl_stderr_grid,'RowNames',"stderr "+M_.exo_names,'VariableNames',"best " + num2str(transpose(1:grid.bestof))));
    disp(array2table(tbl_skew_grid,'RowNames',"skew "+M_.exo_names,'VariableNames',"best " + num2str(transpose(1:grid.bestof))));
    disp(array2table(tbl_loglik_grid,'RowNames',"Log-Lik",'VariableNames',"best " + num2str(transpose(1:grid.bestof))));
    
    % optimize over both Sigma_eta/stderr_eta as well as Gamma_eta/Skew_eta using best values on grid as initial point
    xparams_stage1_tr = nan(estim_params_1.ntot,grid.bestof*length(options_.optim_names)) ; neg_log_likelihood_stage1 = nan(1,grid.bestof*length(options_.optim_names));
    fprintf('\nINITIAL VALUES: Optimize over both Sigma_eta/stderr_eta as well as Gamma_eta/Skew_eta using best values from grid search for CSN skewness/gamma parameter:\n');
    % run maximum likelihood with Pruned Skewed Kalman filter
    diary off
    for jgrid=1:grid.bestof
        idx = (jgrid-1)*length(options_.optim_names) + (1:length(options_.optim_names));
        [xparams_stage1_tr(:,idx), neg_log_likelihood_stage1(1,idx)] = minimize_objective_in_parallel(objfct, xparams_grid{jgrid}, bounds_grid{jgrid}(:,1), bounds_grid{jgrid}(:,2), options_.optim_names, options_.optim_opt,    bounds_grid{jgrid},datamat,estim_params_grid{jgrid},options_,M_grid{jgrid});
    end
    diary(options_.logfile)
    [~,best_stage1] = sort(neg_log_likelihood_stage1);
    % compute standard errors
    [xparams_stage1, bounds1] = dsge_untransform(xparams_stage1_tr, bounds1_tr, estim_params_1);
    xstderr_stage1 = standard_errors_inverse_hessian(objfct, xparams_stage1, bounds1,    bounds1,datamat,estim_params_1,options_STDERR,M_);
    % display summary
    strOptimNamesInit = repmat(options_.optim_names,1,grid.bestof) + "_init_" + string(kron(1:grid.bestof,ones(1,length(options_.optim_names))));
    fprintf('%s\nSUMMARY INITVAL STAGE 1\n\n',repmat('*',1,100))
    fprintf('  Point Estimates\n')
    disp(array2table([xparams_stage1(:,best_stage1); -neg_log_likelihood_stage1(best_stage1)], 'RowNames',[erase(estim_params_1.names,"transformed_");'Log-Lik'], 'VariableNames', strOptimNamesInit(best_stage1)));
    fprintf('  Standard Errors\n')
    disp(array2table([xstderr_stage1(:,best_stage1); -neg_log_likelihood_stage1(best_stage1)], 'RowNames',[erase(estim_params_1.names,"transformed_");'Log-Lik'], 'VariableNames', strOptimNamesInit(best_stage1)));
    fprintf('Final value of the best log-likelihood function: %6.4f\n', -1*neg_log_likelihood_stage1(best_stage1(1)))
    fprintf('%s\n\n',repmat('*',1,100))
    % store results
    oo_.initval.stage1.xparams_tr = xparams_stage1_tr;
    oo_.initval.stage1.xparams = xparams_stage1;
    oo_.initval.stage1.xstderr = xstderr_stage1;
    oo_.initval.stage1.neg_log_likelihood = neg_log_likelihood_stage1;
    % update model and shock parameters
    M_ = dsge_set_params(xparams_stage1_tr(:,best_stage1(1)),estim_params_1,options_,M_);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL VALUE SEARCH STAGE 2
    % - use best values from stage 1 as initial value for shock parameters
    % - use values from stage 0 as initial value for model parameters
    % - run optimization with Pruned Skewed Kalman filter to estimate model as well as shock parameters
    % - compute standard errors using the inverse hessian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAGE = "all_params";
    [estim_params_2,xparamsInit_stage2_tr,bounds2_tr] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
    M_ = dsge_set_params(xparamsInit_stage2_tr, estim_params_2, options_, M_);
    % run maximum likelihood with Pruned Skewed Kalman filter
    fprintf('\n\nINITIAL VALUES: Run maximum likelihood estimation with Pruned Skewed Kalman filter for all parameters using different optimizers in parallel:\n')
    diary off
    [xparams_stage2_tr, neg_log_likelihood_stage2, best_stage2] = minimize_objective_in_parallel(objfct, xparamsInit_stage2_tr, bounds2_tr(:,1), bounds2_tr(:,2), options_.optim_names, options_.optim_opt,    bounds2_tr,datamat,estim_params_2,options_,M_);
    diary(options_.logfile)
    % compute standard errors
    [xparams_stage2, bounds2] = dsge_untransform(xparams_stage2_tr, bounds2_tr, estim_params_2);
    xstderr_stage2 = standard_errors_inverse_hessian(objfct, xparams_stage2, bounds2,    bounds2,datamat,estim_params_2,options_STDERR,M_);
    lr_stat_stage2 = 2 * ( neg_log_likelihood_gauss(best_gauss(1)) - neg_log_likelihood_stage2(best_stage2(1)));
    lr_pval_stage2 = chi2cdf(lr_stat_stage2, size(xparams_stage2,1)-size(xparams_gauss,1), "upper");
    % display summary
    fprintf('%s\nSUMMARY INITVAL STAGE 2\n\n',repmat('*',1,100))
    fprintf('  Point Estimates\n')
    disp(array2table([xparams_stage2(:,best_stage2); -neg_log_likelihood_stage2(best_stage2)], 'RowNames',[erase(estim_params_2.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_stage2)));
    fprintf('  Standard Errors\n')
    disp(array2table([xstderr_stage2(:,best_stage2); -neg_log_likelihood_stage2(best_stage2)], 'RowNames',[erase(estim_params_2.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_stage2)));
    fprintf('Final value of the best log-likelihood function: %6.4f\n', -1*neg_log_likelihood_stage2(best_stage2(1)))
    fprintf('Likelihood Ratio Test Statistic = %.2f with p-val = %.4f\n',lr_stat_stage2,lr_pval_stage2);
    fprintf('%s\n\n',repmat('*',1,100))
    % store results
    oo_.initval.stage2.xparams_tr = xparams_stage2_tr;
    oo_.initval.stage2.xparams = xparams_stage2;
    oo_.initval.stage2.xstderr = xstderr_stage2;
    oo_.initval.stage2.neg_log_likelihood = neg_log_likelihood_stage2;
    % update model and shock parameters for use in stage 3
    M_ = dsge_set_params(xparams_stage2_tr(:,best_stage2(1)),estim_params_2,options_,M_);

    % reset pruning threshold    
    options_.kalman.csn.prune_tol = old_prune_tol;
    options_STDERR.csn.prune_tol = old_prune_tol;
end

%% MAXIMUM LIKELIHOOD ESTIMATION WITH PRUNED SKEWED KALMAN FILTER
% - use provided values (or the ones computed from above initval search)
% - ML with Pruned Skewed Kalman filter
% - compute standard errors using the inverse hessian
if isfield(options_.kalman.csn,'initval_search') && (options_.kalman.csn.initval_search==1)
    STAGE = "all_params";
else
    STAGE = "csn_initval";
end
if length(options_.optim_names)==1
    options_.optim_opt.Display = 'iter';
end
[estim_params_csn,xparamsInit_csn_tr,bounds_csn_tr] = feval(str2func(M_.fname + "_estim_params"), STAGE, M_, options_);
M_ = dsge_set_params(xparamsInit_csn_tr, estim_params_csn, options_, M_);
% run maximum likelihood with Pruned Skewed Kalman filter
fprintf('\n\nRun maximum likelihood estimation with Pruned Skewed Kalman filter for all parameters using different optimizers in parallel:\n')
diary off
[xparams_csn_tr, neg_log_likelihood_csn, best_csn] = minimize_objective_in_parallel(objfct, xparamsInit_csn_tr, bounds_csn_tr(:,1), bounds_csn_tr(:,2), options_.optim_names, options_.optim_opt,    bounds_csn_tr,datamat,estim_params_csn,options_,M_);
diary(options_.logfile)
% compute standard errors
[xparams_csn, bounds_csn] = dsge_untransform(xparams_csn_tr, bounds_csn_tr, estim_params_csn);
xstderr_csn = standard_errors_inverse_hessian(objfct, xparams_csn, bounds_csn,    bounds_csn,datamat,estim_params_csn,options_STDERR,M_);
lr_stat = 2 * ( neg_log_likelihood_gauss(best_gauss(1)) - neg_log_likelihood_csn(best_csn(1)));
lr_pval = chi2cdf(lr_stat, size(xparams_csn,1)-size(xparams_gauss,1), "upper");
% display summary
fprintf('%s\nSUMMARY CSN\n\n',repmat('*',1,100))
fprintf('  Point Estimates\n')
disp(array2table([xparams_csn(:,best_csn); -neg_log_likelihood_csn(best_csn)], 'RowNames',[erase(estim_params_csn.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_csn)));
fprintf('  Standard Errors\n')
disp(array2table([xstderr_csn(:,best_csn); -neg_log_likelihood_csn(best_csn)], 'RowNames',[erase(estim_params_csn.names,"transformed_");'Log-Lik'], 'VariableNames', options_.optim_names(best_csn)));
fprintf('Final value of the best CSN log-likelihood function: %6.4f\n', -1*neg_log_likelihood_csn(best_csn(1)))
fprintf('Likelihood Ratio Test Statistic = %.2f with p-val = %.4f\n',lr_stat,lr_pval);
fprintf('%s\n\n',repmat('*',1,100))
% store results
oo_.csn.xparams_tr = xparams_csn_tr;
oo_.csn.xparams = xparams_csn;
oo_.csn.xstderr = xstderr_csn;
oo_.csn.neg_log_likelihood = neg_log_likelihood_csn;
% update model and shock parameters
M_ = dsge_set_params(xparams_csn_tr(:,best_csn(1)),estim_params_csn,options_,M_);

%% houeskeeping 
diary off
save(['results/' options_.filename]);
