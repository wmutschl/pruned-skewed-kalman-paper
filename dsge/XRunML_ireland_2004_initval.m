%% Maximum Likelihood Estimation of Ireland (2004): "Technology Shocks in The New Keynesian Model", The Review of Economics and Statistics.
% based on Dynare replication codes https://github.com/jpfeifer/DSGE_mod
clear variables; clear global; close all; clc;
addpath('../MATLAB'); % this folder contains our core routines
addpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

%% Options
OPT.modelname = "ireland_2004";
%OPT.datafile  = "data_full_sample";
OPT.datafile = "data_post_1980";
%OPT.datafile = "data_pre_1980";
OPT.optimizer.randomize_initval         = 0; % 1: randomize initial values
OPT.optimizer.bounds.do_param_transform = 0; % 1: transforms parameters with bounded support into parameters with unbounded support by using either exp() or logistic transformation
OPT.optimizer.bounds.penalize_objective = 0; % 1: checks whether bounds are violated in objective function and penalize it
OPT.optimizer.bounds.use_for_optimizer  = 0; % 1: if optimizer supports bounds, use these
%OPT.optimizer.name = "fminsearch";
%OPT.optimizer.name = ["cmaes","simulannealbnd","fminsearch", "fminunc", "patternsearch", "csminwel"];
OPT.optimizer.name = ["fminsearch"];
OPT.optimizer.optim_options = optimset('display','final','MaxFunEvals',1000000,'MaxIter',10000);%,'TolFun',1e-5,'TolX',1e-6);
OPT.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
OPT.prune_tol   = 1e-4;           % pruning threshold
OPT.use_stderr_skew_transform = false;
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
    fprintf('Initial value of the Gaussian log-likelihood function: %6.4f \nTime required to compute log-likelihood function once: %s \n', -1*neg_log_likelihood_0,dynsec2hms(elapsed_time_0))
    [xparam_0,neg_log_likelihood_0] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_0,MODEL_0,OPT_0),  xparam_0,OPT_0.optimizer.optim_options);
    fprintf('Final value of the Gaussian log-likelihood function: %6.4f\n', -1*neg_log_likelihood_0)
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
best_of  = 4;  % best values to keep
Skew_eta_grid = linspace(-0.995,0,grid_nbr/2); Skew_eta_grid = [Skew_eta_grid -Skew_eta_grid((end-1):-1:1)];
if ~OPT.use_stderr_skew_transform
    % for each skewness coefficient compute required Sigma_eta and Gamma_eta such that V[eta] is equal to Gaussian estimate
    diag_Sigma_eta_grid = zeros(MODEL_0.exo_nbr,grid_nbr-1); diag_Gamma_eta_grid = zeros(MODEL_0.exo_nbr,grid_nbr-1);
    for jexo = 1:MODEL_0.exo_nbr
        exo_name = MODEL_0.exo_names(jexo,1);
        Var_eta = ESTIM_PARAM_0.(sprintf('sqrt_diag_Sigma_%s',exo_name)){1}^2;
        for jgrid = 1:(grid_nbr-1)
            [diag_Sigma_eta_grid(jexo,jgrid),diag_Gamma_eta_grid(jexo,jgrid)] = csnVarSkew_To_SigmaGamma(Var_eta,Skew_eta_grid(jgrid),1);        
        end
    end
else
    stderr_eta_0 = nan(MODEL_0.exo_nbr,1);
    for jexo = 1:MODEL_0.exo_nbr
        exo_name = MODEL_0.exo_names(jexo,1);
        stderr_eta_0(jexo) = xparam_0(find(contains(MODEL_0.param_estim_names,"stderr_"+exo_name)));
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
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, [], [], stderr_eta_0, skew_eta_j);
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
diag_Sigma_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); diag_Gamma_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); skew_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of); var_eta_grid_1 = nan(MODEL_1.exo_nbr,best_of);
for jbest=1:best_of
    combo = COMBOS(idx_best_1(jbest),:);
    if OPT.use_stderr_skew_transform
        skew_eta_grid_1(:,jbest) = [Skew_eta_grid(combo(1)); Skew_eta_grid(combo(2)); Skew_eta_grid(combo(3)); Skew_eta_grid(combo(4))];
    else
        diag_Sigma_eta_grid_1(:,jbest) = transpose( [diag_Sigma_eta_grid(1,combo(1)) diag_Sigma_eta_grid(2,combo(2)) diag_Sigma_eta_grid(3,combo(3)) diag_Sigma_eta_grid(4,combo(4))] );
        diag_Gamma_eta_grid_1(:,jbest) = transpose( [diag_Gamma_eta_grid(1,combo(1)) diag_Gamma_eta_grid(2,combo(2)) diag_Gamma_eta_grid(3,combo(3)) diag_Gamma_eta_grid(4,combo(4))] );
        for jexo=1:MODEL_0.exo_nbr
            var_eta_grid_1(jexo,jbest) = csnVar(diag_Sigma_eta_grid_1(jexo,jbest),diag_Gamma_eta_grid_1(jexo,jbest),0,1);
            skew_eta_grid_1(jexo,jbest) = skewness_coef_theor(diag_Sigma_eta_grid_1(jexo,jbest),diag_Gamma_eta_grid_1(jexo,jbest));
        end        
    end
end
diag_Sigma_eta_grid_1
diag_Gamma_eta_grid_1
var_eta_grid_1
skew_eta_grid_1

% optimize over sigma_eta and gamma_eta using best values on grid as initial point
xparam_1 = nan(MODEL_1.exo_nbr + MODEL_1.exo_nbr,best_of) ; neg_log_likelihood_1 = nan(1,best_of);
parfor jbest=1:best_of
    if OPT.use_stderr_skew_transform
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, [], [], stderr_eta_0, skew_eta_grid_1(:,jbest));
    else
        [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, sqrt(diag_Sigma_eta_grid_1(:,jbest)), diag_Gamma_eta_grid_1(:,jbest));        
    end
    [xparam_1(:,jbest),neg_log_likelihood_1(jbest)] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_1_j,MODEL_1_j,OPT_1_j),  xparam_1_j,OPT_1_j.optimizer.optim_options);
end
if OPT.use_stderr_skew_transform
    skew_eta_1 = xparam_1(1:4,1);
    stderr_eta_1 = xparam_1(5:8,1);
else
    diag_Gamma_eta_1 = xparam_1(1:4,1);
    diag_Sigma_eta_1 = xparam_1(5:8,1).^2;
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
% for each of the best 5 parameter combinations: re-run ML optimization of all model parameters (initialized at Gaussian values) and of Sigma_eta and Gamma_eta (initialized at best Sigma_eta_j and Gamma_eta_j of the grid)
MODEL_2 = MODEL; MODEL_2.param_estim_names = MODEL_0.param_estim_names; OPT_2 = OPT; OPT_2.optimizer.optim_options.Display = 'iter';
if OPT.use_stderr_skew_transform
    [xparam_2, PARAM_2, ESTIM_PARAM_2, MODEL_2, OPT_2] = feval(str2func(OPT_2.modelname + "_params"),  2, MODEL_2, OPT_2, xparam_0, [], [], stderr_eta_1, skew_eta_1);
else
    [xparam_2, PARAM_2, ESTIM_PARAM_2, MODEL_2, OPT_2] = feval(str2func(OPT_2.modelname + "_params"),  2, MODEL_2, OPT_2, xparam_0, sqrt(diag_Sigma_eta_1), diag_Gamma_eta_1);
end
[xparam_2,neg_log_likelihood_2] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2,MODEL_2,OPT_2),  xparam_2,OPT_2.optimizer.optim_options);

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

lr_stat = 2 * ( neg_log_likelihood_0 - neg_log_likelihood_2);
lr_pval = chi2cdf(lr_stat, 4, "upper")

% 
% Compute standard errors
OPT_2.optimizer.bounds.lb = OPT_2.optimizer.bounds.lb -0.5;
OPT_2.optimizer.bounds.ub = OPT_2.optimizer.bounds.ub + 0.5;
hess = get_hessian('negative_log_likelihood_dsge',xparam_2,[1e-3;1.0],    DATA.MAT, PARAM_2, MODEL_2, OPT_2);
hess = reshape(hess,MODEL_2.param_estim_nbr,MODEL_2.param_estim_nbr);
V = inv(hess); % estimated covariance matrix of coefficients
SE_values = sqrt(diag(V));


%     [xparam_1_j, PARAM_1_j, ESTIM_PARAM_1_j, MODEL_1_j, OPT_1_j] = feval(str2func(OPT_1.modelname + "_params"),  1, MODEL_1, OPT_1, xparam_0, sqrt_Sigma_eta_j, diag_Gamma_eta_j);
% neg_log_likelihood_grid_1(idx_best_1(1:5))
% gamma_eta_grid = gamma_eta_grid(:,idx_best(1:24));
% neg_log_likelihood_grid = nan(1,size(gamma_eta_grid,2));
% x_grid = nan(MODEL.param_estim_nbr,size(gamma_eta_grid,2));
% 
% objective_function = @(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM,MODEL,OPT);
% parfor_progress(size(gamma_eta_grid,2));
% parfor jgrid = 1:size(gamma_eta_grid,2)
%     xpar0 = xparam_0;
%     xpar0(2) = gamma_eta_grid(1,jgrid);
%     xpar0(4) = gamma_eta_grid(2,jgrid);
%     xpar0(6) = gamma_eta_grid(3,jgrid);
%     xpar0(8) = gamma_eta_grid(4,jgrid);
%     [x_grid(:,jgrid),neg_log_likelihood_grid(:,jgrid)] = fminsearch(objective_function,xpar0,OPT.optimizer.optim_options);
%     parfor_progress;
% end
% parfor_progress(0);
% 
% for jgrid = 1:size(gamma_eta_grid,2)
%     Skew_eta_grid(1,jgrid) = skewness_coef_theor(x_grid(1,jgrid),x_grid(2,jgrid));
%     Skew_eta_grid(2,jgrid) = skewness_coef_theor(x_grid(3,jgrid),x_grid(4,jgrid));
%     Skew_eta_grid(3,jgrid) = skewness_coef_theor(x_grid(5,jgrid),x_grid(6,jgrid));
%     Skew_eta_grid(4,jgrid) = skewness_coef_theor(x_grid(7,jgrid),x_grid(8,jgrid));
% end
% Skew_eta_grid
% [~,idx_best] = sort(neg_log_likelihood_grid);
% neg_log_likelihood_grid(idx_best)
% Skew_eta_grid(:,idx_best)
% x_grid(1,idx_best)
% x_grid(2,idx_best)
% x_grid(3,idx_best)
% x_grid(4,idx_best)
% x_grid(5,idx_best)
% x_grid(6,idx_best)
% x_grid(7,idx_best)
% x_grid(8,idx_best)
%     %% check objective function at initial parameters
%     if OPT.optimizer.randomize_initval
%     end
%     
%     if OPT.optimizer.bounds.do_param_transform    
%     end
%     
%     % compute, check and print objective function
%     tic_id = tic;
%     [neg_log_likelihood,exit_flag] = negative_log_likelihood_dsge(xpar0,DATA.MAT,PARAM,MODEL,OPT);
%     elapsed_time = toc(tic_id);
%     if exit_flag ~= 1
%         error('Something wrong with log-likelihood function at initial parameters')
%     end
%     fprintf('Initial value of the log-likelihood function: %6.4f \n', -1*neg_log_likelihood);
%     fprintf('Time required to compute log-likelihood function once: %5.4f seconds \n', elapsed_time);
%     
%     
%     %% Maximum Likelihood Estimation with different optimizers
%     if MODEL.param_estim_nbr > 0
%     
%         % initialize output structure, columns are for different optimizers
%         ML.f  = nan(1,size(OPT.optimizer.name,2)); % stores function value
%         ML.x  = nan(MODEL.param_estim_nbr,size(OPT.optimizer.name,2)); % stores optimized parameters
%         ML.se = nan(MODEL.param_estim_nbr,size(OPT.optimizer.name,2)); % stores standard errors of optimized parameters
%         
%         % deal with bounds
%         if OPT.optimizer.bounds.use_for_optimizer
%             LB = OPT.optimizer.bounds.lb;
%             UB = OPT.optimizer.bounds.ub;
%         else
%             LB = -Inf*ones(MODEL.param_estim_nbr,1);
%             UB = Inf*ones(MODEL.param_estim_nbr,1);
%         end
%         
%         % pass additional arguments to objective function using an anonymous function
%         objective_function = @(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM,MODEL,OPT);
%         
%         for jopt = 1:size(OPT.optimizer.name,2)            
%             tic_opt = tic;
%             switch OPT.optimizer.name(jopt)
%                 case 'fminsearch' % does not care about bounds
%                     [xparam1,fval] = fminsearch(objective_function,xpar0,OPT.optimizer.optim_options);
%                 case 'fminunc' % does not care about bounds
%                     [xparam1,fval] = fminunc(objective_function,xpar0,OPT.optimizer.optim_options);
%                 case 'patternsearch' % does care about bounds
%                     [xparam1,fval] = patternsearch(objective_function,xpar0,[],[],[],[],LB,UB,[],OPT.optimizer.optim_options);
%                 case 'csminwel' % does not care about bounds
%                     [fval,xparam1] = csminwel(objective_function,xpar0,1e-4*eye(MODEL.param_estim_nbr),[],OPT.optimizer.optim_options.TolFun,OPT.optimizer.optim_options.MaxIter);
%                 case 'simulannealbnd' % does care about bounds
%                     [xparam1,fval] = simulannealbnd(objective_function,xpar0,LB,UB,OPT.optimizer.optim_options);
%                 case 'cmaes' % does care about bounds
%                     cmaesOptions.SaveVariables='on'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='10'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
%                     cmaesOptions.TolFun = OPT.optimizer.optim_options.TolFun;
%                     cmaesOptions.TolX = OPT.optimizer.optim_options.TolX;
%                     cmaesOptions.MaxIter = OPT.optimizer.optim_options.MaxIter;
%                     cmaesOptions.MaxFunEvals = OPT.optimizer.optim_options.MaxFunEvals;        
%                     cmaesOptions.LBounds = LB;
%                     cmaesOptions.UBounds = UB;
%                     % Set default search volume (SIGMA)            
%                     cmaesSIGMA = (UB-LB)*0.2;
%                     cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
%                     while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
%                         cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
%                     end
%                     %warning('off','CMAES:NonfinitenessRange');
%                     %warning('off','CMAES:InitialSigma');
%                     [~, ~, ~, ~, ~, BESTEVER] = cmaes('negative_log_likelihood_dsge',xpar0,cmaesSIGMA,cmaesOptions,...
%                                                                               DATA.MAT, PARAM, MODEL, OPT);
%                     xparam1=BESTEVER.x;
%                     fval=BESTEVER.f;
%             end
%             elapsed_time_opt = toc(tic_opt);
%             % Retransform parameters if needed
%             if OPT.optimizer.bounds.do_param_transform
%                 xopt = param_transform_bounded(xparam1,OPT.optimizer.bounds.lb,OPT.optimizer.bounds.ub);
%             else
%                 xopt = xparam1;        
%             end
%             % Compute standard errors
%             old_param_transform = OPT.optimizer.bounds.do_param_transform;
%             OPT.optimizer.bounds.do_param_transform = 0; %turn off for computing standard errors
%             OPT.optimizer.bounds.lb = OPT.optimizer.bounds.lb -0.5;
%             OPT.optimizer.bounds.ub = OPT.optimizer.bounds.ub + 0.5;
%             hess = get_hessian('negative_log_likelihood_dsge',xopt,[1e-3;1.0],    DATA.MAT, PARAM, MODEL, OPT);
%             hess = reshape(hess,MODEL.param_estim_nbr,MODEL.param_estim_nbr);
%             V = inv(hess); % estimated covariance matrix of coefficients
%             SE_values = sqrt(diag(V));
%             OPT.optimizer.do_param_transform = old_param_transform; %turn back on
%         
%             % store into output structure
%             ML.f(jopt) = fval;
%             ML.x(:,jopt) = xopt;
%             ML.se(:,jopt) = SE_values;
%             ML.elapsed_time(jopt) = elapsed_time_opt;
%             disp( table(ML.x(:,jopt),ML.se(:,jopt),ML.x(:,jopt)./ML.se(:,jopt),'Var',{'Estimate','s.d.','t-stat'},'Row',MODEL.param_estim_names) );
%             fprintf('Value of maximimized log-likelihood function: %12.10f.\n',-ML.f(jopt));
%             fprintf('Time required: %5.4f seconds \n', ML.elapsed_time(jopt));    
%             %input(sprintf('Finished optimization with %s. Hit enter for next optimizer',OPT.optimizer.name(jopt)));
%         end
%         
%         % Display estimation results
%         clc; % clear screen
%         [~,jbest] = max(-ML.f);
%         fprintf('SUMMARY\n\n')
%         for jopt=1:size(OPT.optimizer.name,2)
%             if jopt == jbest
%                 fprintf('BEST ')
%             end
%             fprintf('Optimizer %s: ',OPT.optimizer.name(jopt));
%             fprintf('value of maximimized log-likelihood function: %12.10f.\n',-ML.f(jopt));
%             fprintf('Time required: %5.4f seconds \n', ML.elapsed_time(jopt));
%             disp( table(ML.x(:,jopt),ML.se(:,jopt),ML.x(:,jopt)./ML.se(:,jopt),'Var',{'Estimate','s.d.','t-stat'},'Row',MODEL.param_estim_names) );
%         end 
%     end
%     neg_log_likelihood_grid(jgrid) = ML.f(jopt);
%     x_grid(:,jgrid) = xopt;
% end % jgrid


% housekeeping
rmpath('../MATLAB'); % this folder contains our core routines
rmpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

