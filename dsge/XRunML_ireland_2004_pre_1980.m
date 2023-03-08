%% Maximum Likelihood Estimation of Ireland (2004): "Technology Shocks in The New Keynesian Model", The Review of Economics and Statistics.
% based on Dynare replication codes https://github.com/jpfeifer/DSGE_mod
clear variables; clear global; close all; clc;
addpath('../MATLAB'); % this folder contains our core routines
addpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

%% Options
OPT.modelname = "ireland_2004";
%OPT.datafile  = "data_full_sample";
%OPT.datafile = "data_post_1980";
OPT.datafile = "data_pre_1980";
OPT.optimizer.bounds.penalize_objective = 1; % 1: checks whether bounds are violated in objective function and penalize it
OPT.optimizer.bounds.use_for_optimizer  = 0; % 1: if optimizer supports bounds, use these
OPT.use_stderr_skew_transform           = true;
%OPT.optimizer.name = "fminsearch";
%OPT.optimizer.name = ["cmaes","simulannealbnd","fminsearch", "fminunc", "patternsearch", "csminwel"];
OPT.optimizer.name = ["fminsearch"];
OPT.optimizer.optim_options = optimset('display','final','MaxFunEvals',1000000,'MaxIter',10000);%,'TolFun',1e-5,'TolX',1e-6);
OPT.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
OPT.prune_tol   = 1e-4;           % pruning threshold
%% Preprocessing of model, i.e. create script files with dynamic Jacobians, which can be evaluated numerically to compute the policy function
MODEL = feval(str2func(OPT.modelname + "_preprocessing"));

%% Create DATA and selection matrix F
load ../data/ireland_2004_gpr.dat;
DATA.dates = transpose(datetime('1948-Q2','InputFormat','yyyy-QQQ'):calquarters(1):datetime('2003-Q1','InputFormat','yyyy-QQQ')); DATA.dates.Format = 'yyyy-QQQ';
if OPT.datafile == "data_full_sample"
    date_start = find(DATA.dates==datetime('1948-Q2','InputFormat','yyyy-QQQ'));
    date_end   = find(DATA.dates==datetime('2003-Q1','InputFormat','yyyy-QQQ'));
elseif OPT.datafile == "data_pre_1980"
    date_start = find(DATA.dates==datetime('1948-Q2','InputFormat','yyyy-QQQ'));
    date_end   = find(DATA.dates==datetime('1979-Q4','InputFormat','yyyy-QQQ'));
elseif OPT.datafile == "data_post_1980"
    date_start = find(DATA.dates==datetime('1980-Q1','InputFormat','yyyy-QQQ'));
    date_end   = find(DATA.dates==datetime('2003-Q1','InputFormat','yyyy-QQQ'));
end
DATA.dates = DATA.dates(date_start:date_end);
DATA.ghat  = ireland_2004_gpr(date_start:date_end,1) - mean(ireland_2004_gpr(date_start:date_end,1)); % output growth (demeaned), based on quarterly changes in seasonally adjusted real GDP, converted to per capita terms by dividing by the civilian noninsitutional population aged 16 and over
DATA.pihat = ireland_2004_gpr(date_start:date_end,2) - mean(ireland_2004_gpr(date_start:date_end,2)); % inflation (demeaned), based on quarterly changes in seasonally adjusted GDP deflator
DATA.rhat  = ireland_2004_gpr(date_start:date_end,3) - mean(ireland_2004_gpr(date_start:date_end,3)); % nominal interest rate (demeaned), based on quarterly averages of daily readings on the three-months US Treasury bill rate
OPT.nobs   = length(DATA.dates);
DATA.MAT   = nan(OPT.nobs,MODEL.varobs_nbr); % initialize data matrix corresponding to MODEL.varobs_idx_DR
MODEL.F    = zeros(MODEL.varobs_nbr,MODEL.endo_nbr); % initialize selection matrix F corresponding to MODEL.varobs_idx_DR
for j=1:MODEL.varobs_nbr
    DATA.MAT(:,j) = DATA.(MODEL.varobs{j});
    MODEL.F(j,MODEL.varobs_idx_DR(j)) = 1;
end
clear ireland_2004_gpr date_start date_end idx j jexo jvarobs tmp tmp_varobs_names_DR


%% Stage 0: use reported ML estimates of Ireland (2004) for initial values, re-run ML with Gaussian Kalman filter
OPT_0 = OPT;
MODEL_0 = MODEL;
[xparam_0, PARAM_0, ESTIM_PARAM_0, MODEL_0, OPT_0] = feval(str2func(OPT.modelname + "_params"),  0, MODEL_0, OPT_0);
% compute, check and print objective function
tic_id_0 = tic;
[neg_log_likelihood_0,exit_flag_0] = negative_log_likelihood_dsge(xparam_0,DATA.MAT,PARAM_0,MODEL_0,OPT_0);
elapsed_time_0 = toc(tic_id_0);
if exit_flag_0 ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
else
    fprintf('Initial value of the Gaussian log-likelihood function: %6.4f \nTime required to compute log-likelihood function once: %s \n', -1*neg_log_likelihood_0,dynsec2hms(elapsed_time_0));
    if OPT.optimizer.bounds.use_for_optimizer
        [xparam_0,neg_log_likelihood_0] = fminsearchbnd(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_0,MODEL_0,OPT_0),  xparam_0,OPT_0.optimizer.bounds.lb,OPT_0.optimizer.bounds.ub,OPT_0.optimizer.optim_options);
    else
        [xparam_0,neg_log_likelihood_0] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_0,MODEL_0,OPT_0),  xparam_0,OPT_0.optimizer.optim_options);
    end
    fprintf('Final value of the Gaussian log-likelihood function: %6.4f\n', -1*neg_log_likelihood_0)
end
stderr_eta_0 = nan(MODEL_0.exo_nbr,1);
for jexo = 1:MODEL_0.exo_nbr
    exo_name = MODEL_0.exo_names(jexo,1);
    if ~OPT.use_stderr_skew_transform
        stderr_eta_0(jexo) = xparam_0(find(contains(MODEL_0.param_estim_names,"sqrt_diag_Sigma_"+exo_name)));
    else
        stderr_eta_0(jexo) = xparam_0(find(contains(MODEL_0.param_estim_names,"stderr_"+exo_name)));
    end
end

%% Stage 1:
% Note:
% - V[eta] and skewness(eta) are both a function of Sigma_eta and Gamma_eta, note that all the eta's are independent of each other, so we can make use of univariate formulas
% - eta_j ~ CSN(mu_eta_j,Sigma_eta_j,Gamma_eta_j,nu_eta_j,Delta_eta_j) with nu_eta_j=0, Delta_eta_j=1 and mu_eta_j such that E[eta_j]=0; all eta_j are independent of each other
% At this stage we will:
% - fix model parameters to estimated Gaussian ML values of stage 0
% - fix Var[eta_j]=Var[eta_j]_{Gaussian} to Gaussian estimate of stage 0
% - create grid for Skew[eta], but beware of the theoretical bounds: abs(skewness(eta)) < 0.995 < abs(sqrt(2)*(pi-4))/(pi-2)^(3/2))
% - create grid of corresponding Sigma_eta_j and Gamma_eta_j combinations, i.e. recover Sigma_eta_j and Gamma_eta_j from Skew[eta_j] and Var[eta_j]_{Gaussian}
% - compute negative log-likelihood of each value on grid, keep best 5 parameter combinations

OPT_1 = OPT; MODEL_1 = MODEL; MODEL_1.param_estim_names = MODEL_0.param_estim_names;

% create an evenly spaced grid of skewness coefficients between -0.995 and 0.995 for each shock
grid_nbr = 16; % needs to be even number
best_of  = 3;  % best values to keep
Skew_eta_grid = linspace(-0.95,0,grid_nbr/2); Skew_eta_grid = [Skew_eta_grid -Skew_eta_grid((end-1):-1:1)];
if ~OPT.use_stderr_skew_transform
    % for each skewness coefficient compute required Sigma_eta and Gamma_eta such that V[eta] is equal to Gaussian estimate
    diag_Sigma_eta_grid = zeros(MODEL_0.exo_nbr,grid_nbr-1); diag_Gamma_eta_grid = zeros(MODEL_0.exo_nbr,grid_nbr-1);
    for jexo = 1:MODEL_0.exo_nbr
        for jgrid = 1:(grid_nbr-1)
            [diag_Sigma_eta_grid(jexo,jgrid),diag_Gamma_eta_grid(jexo,jgrid)] = csnVarSkew_To_SigmaGamma(stderr_eta_0(jexo)^2,Skew_eta_grid(jgrid),1);        
        end
    end
end
% compute negative loglikelihood for all possible combinations of grid values
COMBOS = zeros((grid_nbr-1)^4,4,'int16'); idx = 1;
for j1=1:(grid_nbr-1)
    for j2=1:(grid_nbr-1)
        for j3=1:(grid_nbr-1)
            for j4=1:(grid_nbr-1)
                COMBOS(idx,1) = j1; COMBOS(idx,2) = j2; COMBOS(idx,3) = j3; COMBOS(idx,4) = j4;
                idx = idx+1;
            end
        end
    end
end
neg_log_likelihood_grid_1 = nan(1,(grid_nbr-1)^4);
fprintf('Running grid search for CSN parameters:\n');
parfor_progress((grid_nbr-1)^4);
parfor jgrid = 1:(grid_nbr-1)^4
    combo = COMBOS(jgrid,:);
    if OPT.use_stderr_skew_transform
        skew_eta_j = [Skew_eta_grid(combo(1)); Skew_eta_grid(combo(2)); Skew_eta_grid(combo(3)); Skew_eta_grid(combo(4))];
        stderr_eta_j = stderr_eta_0;
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, [], [], stderr_eta_j, skew_eta_j);
    else
        sqrt_Sigma_eta_j = sqrt( [diag_Sigma_eta_grid(1,combo(1)) diag_Sigma_eta_grid(2,combo(2)) diag_Sigma_eta_grid(3,combo(3)) diag_Sigma_eta_grid(4,combo(4))] );    
        diag_Gamma_eta_j = [diag_Gamma_eta_grid(1,combo(1)) diag_Gamma_eta_grid(2,combo(2)) diag_Gamma_eta_grid(3,combo(3)) diag_Gamma_eta_grid(4,combo(4))];
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, sqrt_Sigma_eta_j, diag_Gamma_eta_j);
    end    
    neg_log_likelihood_grid_1(jgrid) = negative_log_likelihood_dsge(xparam_1_j,DATA.MAT,PARAM_1_j,MODEL_1_j,OPT_1_j);
    parfor_progress;
end
parfor_progress(0);
[~,idx_best_1] = sort(neg_log_likelihood_grid_1);
best_neg_log_likelihood_grid_1 = neg_log_likelihood_grid_1(idx_best_1(1:best_of))
diag_Sigma_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); diag_Gamma_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); skew_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); stderr_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of);
for jbest=1:best_of
    combo = COMBOS(idx_best_1(jbest),:);
    if OPT.use_stderr_skew_transform
        skew_eta_grid_1(:,jbest) = [Skew_eta_grid(combo(1)); Skew_eta_grid(combo(2)); Skew_eta_grid(combo(3)); Skew_eta_grid(combo(4))];
        stderr_eta_grid_1(:,jbest) = stderr_eta_0;
    else
        diag_Sigma_eta_grid_1(:,jbest) = transpose( [diag_Sigma_eta_grid(1,combo(1)) diag_Sigma_eta_grid(2,combo(2)) diag_Sigma_eta_grid(3,combo(3)) diag_Sigma_eta_grid(4,combo(4))] );
        diag_Gamma_eta_grid_1(:,jbest) = transpose( [diag_Gamma_eta_grid(1,combo(1)) diag_Gamma_eta_grid(2,combo(2)) diag_Gamma_eta_grid(3,combo(3)) diag_Gamma_eta_grid(4,combo(4))] );
        for jexo=1:MODEL_0.exo_nbr
            stderr_eta_grid_1(jexo,jbest) = sqrt(csnVar(diag_Sigma_eta_grid_1(jexo,jbest),diag_Gamma_eta_grid_1(jexo,jbest),0,1));
            skew_eta_grid_1(jexo,jbest) = skewness_coef_theor(diag_Sigma_eta_grid_1(jexo,jbest),diag_Gamma_eta_grid_1(jexo,jbest));
        end
    end
end
diag_Sigma_eta_grid_1
diag_Gamma_eta_grid_1
stderr_eta_grid_1
skew_eta_grid_1

% optimize over sigma_eta and gamma_eta using best values on grid as initial point
xparam_1 = nan(MODEL_1.exo_nbr + MODEL_1.exo_nbr,best_of) ; neg_log_likelihood_1 = nan(1,best_of);
parfor jbest=1:best_of
    if OPT_1.use_stderr_skew_transform
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, [], [], stderr_eta_0, skew_eta_grid_1(:,jbest));
    else
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, sqrt(diag_Sigma_eta_grid_1(:,jbest)), diag_Gamma_eta_grid_1(:,jbest));        
    end
    if OPT_1.optimizer.bounds.use_for_optimizer
        [xparam_1(:,jbest),neg_log_likelihood_1(jbest)] = fminsearchbnd(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_1_j,MODEL_1_j,OPT_1_j),  xparam_1_j,OPT_1_j.optimizer.bounds.lb,OPT_1_j.optimizer.bounds.ub,OPT_1_j.optimizer.optim_options);
    else
        [xparam_1(:,jbest),neg_log_likelihood_1(jbest)] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_1_j,MODEL_1_j,OPT_1_j),  xparam_1_j,OPT_1_j.optimizer.optim_options);
    end
end
[~,best_of_1] = min(neg_log_likelihood_1);
if OPT.use_stderr_skew_transform
    skew_eta_1 = xparam_1(1:4,best_of_1);
    stderr_eta_1 = xparam_1(5:8,best_of_1);
else
    diag_Gamma_eta_1 = xparam_1(1:4,best_of_1);
    diag_Sigma_eta_1 = xparam_1(5:8,best_of_1).^2;
    stderr_eta_1 = nan(4,1); skew_eta_1 = nan(4,1);
    for jexo=1:4
        stderr_eta_1(jexo,1) = sqrt(csnVar(diag_Sigma_eta_1(jexo,1),diag_Gamma_eta_1(jexo,1),0,1));
        skew_eta_1(jexo,1) = skewness_coef_theor(diag_Sigma_eta_1(jexo,1),diag_Gamma_eta_1(jexo,1));
    end
    diag_Sigma_eta_1
    diag_Gamma_eta_1
end
stderr_eta_1
skew_eta_1

%% Stage 2:
MODEL_2 = MODEL; MODEL_2.param_estim_names = MODEL_0.param_estim_names; OPT_2 = OPT; OPT_2.optimizer.optim_options.Display = 'iter';
if OPT.use_stderr_skew_transform
    [xparam_2, PARAM_2, ESTIM_PARAM_2, MODEL_2, OPT_2] = feval(str2func(OPT_2.modelname + "_params"),  2, MODEL_2, OPT_2, xparam_0, [], [], stderr_eta_1, skew_eta_1);
else
    [xparam_2, PARAM_2, ESTIM_PARAM_2, MODEL_2, OPT_2] = feval(str2func(OPT_2.modelname + "_params"),  2, MODEL_2, OPT_2, xparam_0, sqrt(diag_Sigma_eta_1), diag_Gamma_eta_1);
end
if OPT_2.optimizer.bounds.use_for_optimizer
    [xparam_2,neg_log_likelihood_2] = fminsearchbnd(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2,MODEL_2,OPT_2),  xparam_2,OPT_2.optimizer.bounds.lb,OPT_2.optimizer.bounds.ub,OPT_2.optimizer.optim_options);
else
    [xparam_2,neg_log_likelihood_2] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2,MODEL_2,OPT_2),  xparam_2,OPT_2.optimizer.optim_options);
end

if OPT.use_stderr_skew_transform
    skew_eta_2 = xparam_2(1:4,1);
    stderr_eta_2 = xparam_2(5:8,1);
else
    diag_Gamma_eta_2 = xparam_2(1:4);
    diag_Sigma_eta_2 = xparam_2(5:8).^2;
    stderr_eta_2 = nan(4,1); skew_eta_2 = nan(4,1);
    for jexo=1:4
        stderr_eta_2(jexo,1) = sqrt(csnVar(diag_Sigma_eta_2(jexo,1),diag_Gamma_eta_2(jexo,1),0,1));
        skew_eta_2(jexo,1) = skewness_coef_theor(diag_Sigma_eta_2(jexo,1),diag_Gamma_eta_2(jexo,1));
    end
end
xparam_2(9:end)
stderr_eta_2
skew_eta_2
fprintf('Final value of the CSN log-likelihood function: %6.4f\n', -1*neg_log_likelihood_2)
lr_stat = 2 * ( neg_log_likelihood_0 - neg_log_likelihood_2)
lr_pval = chi2cdf(lr_stat, 4, "upper")

% 
% Compute standard errors
OPT_2.optimizer.bounds.lb = OPT_2.optimizer.bounds.lb -0.5;
OPT_2.optimizer.bounds.ub = OPT_2.optimizer.bounds.ub + 0.5;
hess = get_hessian('negative_log_likelihood_dsge',xparam_2,[1e-3;1.0],    DATA.MAT, PARAM_2, MODEL_2, OPT_2);
hess = reshape(hess,MODEL_2.param_estim_nbr,MODEL_2.param_estim_nbr);
V = inv(hess); % estimated covariance matrix of coefficients
SE_values = sqrt(diag(V))

% housekeeping
rmpath('../MATLAB'); % this folder contains our core routines
rmpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

