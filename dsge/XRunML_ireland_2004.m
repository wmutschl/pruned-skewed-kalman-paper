%% Maximum Likelihood Estimation of Ireland (2004): "Technology Shocks in The New Keynesian Model", The Review of Economics and Statistics.
% based on Dynare replication codes https://github.com/jpfeifer/DSGE_mod
clear variables; clear global; close all; clc;
addpath('../MATLAB'); % this folder contains our core routines
addpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

%% ------------------------------------------------------------------------
% USER CHOICES
% -------------------------------------------------------------------------
% Options
OPT.modelname = "ireland_2004";
OPT.datafile  = "data_full_sample";
%OPT.datafile = "data_post_1980";
%OPT.datafile = "data_pre_1980";
OPT.optimizer.randomize_initval         = 0; % 1: randomize initial values
OPT.optimizer.bounds.do_param_transform = 0; % 1: transforms parameters with bounded support into parameters with unbounded support by using either exp() or logistic transformation
OPT.optimizer.bounds.penalize_objective = 0; % 1: checks whether bounds are violated in objective function and penalize it
OPT.optimizer.bounds.use_for_optimizer  = 0; % 1: if optimizer supports bounds, use these
%OPT.optimizer.name = "fminsearch";
%OPT.optimizer.name = ["cmaes","simulannealbnd","fminsearch", "fminunc", "patternsearch", "csminwel"];
OPT.optimizer.name = ["fminsearch"];
OPT.optimizer.optim_options = optimset('display','final','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6);
OPT.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
OPT.prune_tol   = 1e-4;           % pruning threshold

OPT.estim.eta_a  = 1;
OPT.estim.eta_e  = 1;
OPT.estim.eta_z  = 1;
OPT.estim.eta_r  = 1;
OPT.estim.params = 0;

OPT.calib.eta_a  = 0;
OPT.calib.eta_e  = 0;
OPT.calib.eta_z  = 0;
OPT.calib.eta_r  = 0;
OPT.calib.params = 1;

%% Preprocessing of model
% create script files with dynamic Jacobians, which can be evaluated numerically to compute the policy function
MODEL = feval(str2func(OPT.modelname + "_preprocessing"));
% get parameter information
[PARAM, ESTIM_PARAM, MODEL]    = feval(str2func(OPT.modelname + "_params"),MODEL,OPT);
if ~isempty(ESTIM_PARAM)
    MODEL.param_estim_names = fieldnames(ESTIM_PARAM);
else
    MODEL.param_estim_names = [];
end
MODEL.param_estim_nbr   = size(MODEL.param_estim_names,1);

if MODEL.varobs_nbr > (MODEL.exo_nbr + nnz(diag(MODEL.Sigma_eps)))
    error('Stochastic Singularity: You have more observables than shocks+measurement errors.')
end

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

gamma_eta_grid = grid_gamma_eta(ESTIM_PARAM.sqrt_Sigma_eta_a{1},ESTIM_PARAM.sqrt_Sigma_eta_e{1},ESTIM_PARAM.sqrt_Sigma_eta_z{1},ESTIM_PARAM.sqrt_Sigma_eta_r{1});

% create vectors from structures
xparam0 = nan(MODEL.param_estim_nbr,1);
for j = 1:MODEL.param_estim_nbr
    OPT.optimizer.bounds.lb(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){2};
    OPT.optimizer.bounds.ub(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){3};
    xparam0(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){1};
    %PARAM.(MODEL.param_estim_names{j}) = xparam0(j,1);
end

neg_log_likelihood_grid = nan(1,size(gamma_eta_grid,2));
dat = DATA.MAT;

parfor_progress(size(gamma_eta_grid,2));
parfor jgrid = 1:size(gamma_eta_grid,2)
    xpar0 = xparam0;
    xpar0(2) = gamma_eta_grid(1,jgrid);
    xpar0(4) = gamma_eta_grid(2,jgrid);
    xpar0(6) = gamma_eta_grid(3,jgrid);
    xpar0(8) = gamma_eta_grid(4,jgrid);    
    neg_log_likelihood_grid(jgrid) = negative_log_likelihood_dsge(xpar0,dat,PARAM,MODEL,OPT);
    parfor_progress;
end
parfor_progress(0);
[~,idx_best] = sort(neg_log_likelihood_grid);

gamma_eta_grid = gamma_eta_grid(:,idx_best(1:8));
neg_log_likelihood_grid = nan(1,size(gamma_eta_grid,2));
x_grid = nan(MODEL.param_estim_nbr,size(gamma_eta_grid,2));

objective_function = @(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM,MODEL,OPT);
parfor_progress(size(gamma_eta_grid,2));
parfor jgrid = 1:size(gamma_eta_grid,2)
    xpar0 = xparam0;
    xpar0(2) = gamma_eta_grid(1,jgrid);
    xpar0(4) = gamma_eta_grid(2,jgrid);
    xpar0(6) = gamma_eta_grid(3,jgrid);
    xpar0(8) = gamma_eta_grid(4,jgrid);
    [x_grid(:,jgrid),neg_log_likelihood_grid(:,jgrid)] = fminsearch(objective_function,xpar0,OPT.optimizer.optim_options);
    parfor_progress;
end
parfor_progress(0);


    %% check objective function at initial parameters
    if OPT.optimizer.randomize_initval
    end
    
    if OPT.optimizer.bounds.do_param_transform    
    end
    
    % compute, check and print objective function
    tic_id = tic;
    [neg_log_likelihood,exit_flag] = negative_log_likelihood_dsge(xpar0,DATA.MAT,PARAM,MODEL,OPT);
    elapsed_time = toc(tic_id);
    if exit_flag ~= 1
        error('Something wrong with log-likelihood function at initial parameters')
    end
    fprintf('Initial value of the log-likelihood function: %6.4f \n', -1*neg_log_likelihood);
    fprintf('Time required to compute log-likelihood function once: %5.4f seconds \n', elapsed_time);
    
    
    %% Maximum Likelihood Estimation with different optimizers
    if MODEL.param_estim_nbr > 0
    
        % initialize output structure, columns are for different optimizers
        ML.f  = nan(1,size(OPT.optimizer.name,2)); % stores function value
        ML.x  = nan(MODEL.param_estim_nbr,size(OPT.optimizer.name,2)); % stores optimized parameters
        ML.se = nan(MODEL.param_estim_nbr,size(OPT.optimizer.name,2)); % stores standard errors of optimized parameters
        
        % deal with bounds
        if OPT.optimizer.bounds.use_for_optimizer
            LB = OPT.optimizer.bounds.lb;
            UB = OPT.optimizer.bounds.ub;
        else
            LB = -Inf*ones(MODEL.param_estim_nbr,1);
            UB = Inf*ones(MODEL.param_estim_nbr,1);
        end
        
        % pass additional arguments to objective function using an anonymous function
        objective_function = @(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM,MODEL,OPT);
        
        for jopt = 1:size(OPT.optimizer.name,2)            
            tic_opt = tic;
            switch OPT.optimizer.name(jopt)
                case 'fminsearch' % does not care about bounds
                    [xparam1,fval] = fminsearch(objective_function,xpar0,OPT.optimizer.optim_options);
                case 'fminunc' % does not care about bounds
                    [xparam1,fval] = fminunc(objective_function,xpar0,OPT.optimizer.optim_options);
                case 'patternsearch' % does care about bounds
                    [xparam1,fval] = patternsearch(objective_function,xpar0,[],[],[],[],LB,UB,[],OPT.optimizer.optim_options);
                case 'csminwel' % does not care about bounds
                    [fval,xparam1] = csminwel(objective_function,xpar0,1e-4*eye(MODEL.param_estim_nbr),[],OPT.optimizer.optim_options.TolFun,OPT.optimizer.optim_options.MaxIter);
                case 'simulannealbnd' % does care about bounds
                    [xparam1,fval] = simulannealbnd(objective_function,xpar0,LB,UB,OPT.optimizer.optim_options);
                case 'cmaes' % does care about bounds
                    cmaesOptions.SaveVariables='on'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='10'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
                    cmaesOptions.TolFun = OPT.optimizer.optim_options.TolFun;
                    cmaesOptions.TolX = OPT.optimizer.optim_options.TolX;
                    cmaesOptions.MaxIter = OPT.optimizer.optim_options.MaxIter;
                    cmaesOptions.MaxFunEvals = OPT.optimizer.optim_options.MaxFunEvals;        
                    cmaesOptions.LBounds = LB;
                    cmaesOptions.UBounds = UB;
                    % Set default search volume (SIGMA)            
                    cmaesSIGMA = (UB-LB)*0.2;
                    cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
                    while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
                        cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
                    end
                    %warning('off','CMAES:NonfinitenessRange');
                    %warning('off','CMAES:InitialSigma');
                    [~, ~, ~, ~, ~, BESTEVER] = cmaes('negative_log_likelihood_dsge',xpar0,cmaesSIGMA,cmaesOptions,...
                                                                              DATA.MAT, PARAM, MODEL, OPT);
                    xparam1=BESTEVER.x;
                    fval=BESTEVER.f;
            end
            elapsed_time_opt = toc(tic_opt);
            % Retransform parameters if needed
            if OPT.optimizer.bounds.do_param_transform
                xopt = param_transform_bounded(xparam1,OPT.optimizer.bounds.lb,OPT.optimizer.bounds.ub);
            else
                xopt = xparam1;        
            end
            % Compute standard errors
            old_param_transform = OPT.optimizer.bounds.do_param_transform;
            OPT.optimizer.bounds.do_param_transform = 0; %turn off for computing standard errors
            OPT.optimizer.bounds.lb = OPT.optimizer.bounds.lb -0.5;
            OPT.optimizer.bounds.ub = OPT.optimizer.bounds.ub + 0.5;
            hess = get_hessian('negative_log_likelihood_dsge',xopt,[1e-3;1.0],    DATA.MAT, PARAM, MODEL, OPT);
            hess = reshape(hess,MODEL.param_estim_nbr,MODEL.param_estim_nbr);
            V = inv(hess); % estimated covariance matrix of coefficients
            SE_values = sqrt(diag(V));
            OPT.optimizer.do_param_transform = old_param_transform; %turn back on
        
            % store into output structure
            ML.f(jopt) = fval;
            ML.x(:,jopt) = xopt;
            ML.se(:,jopt) = SE_values;
            ML.elapsed_time(jopt) = elapsed_time_opt;
            disp( table(ML.x(:,jopt),ML.se(:,jopt),ML.x(:,jopt)./ML.se(:,jopt),'Var',{'Estimate','s.d.','t-stat'},'Row',MODEL.param_estim_names) );
            fprintf('Value of maximimized log-likelihood function: %12.10f.\n',-ML.f(jopt));
            fprintf('Time required: %5.4f seconds \n', ML.elapsed_time(jopt));    
            %input(sprintf('Finished optimization with %s. Hit enter for next optimizer',OPT.optimizer.name(jopt)));
        end
        
        % Display estimation results
        clc; % clear screen
        [~,jbest] = max(-ML.f);
        fprintf('SUMMARY\n\n')
        for jopt=1:size(OPT.optimizer.name,2)
            if jopt == jbest
                fprintf('BEST ')
            end
            fprintf('Optimizer %s: ',OPT.optimizer.name(jopt));
            fprintf('value of maximimized log-likelihood function: %12.10f.\n',-ML.f(jopt));
            fprintf('Time required: %5.4f seconds \n', ML.elapsed_time(jopt));
            disp( table(ML.x(:,jopt),ML.se(:,jopt),ML.x(:,jopt)./ML.se(:,jopt),'Var',{'Estimate','s.d.','t-stat'},'Row',MODEL.param_estim_names) );
        end 
    end
    neg_log_likelihood_grid(jgrid) = ML.f(jopt);
    x_grid(:,jgrid) = xopt;
end % jgrid


% housekeeping
rmpath('../MATLAB'); % this folder contains our core routines
rmpath('../MATLAB/OptimRoutines/CMAES','../MATLAB/OptimRoutines/CSMINWEL'); % this folder has two useful optimizers: csminwel and cmaes

