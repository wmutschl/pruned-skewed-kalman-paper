%% Preprocessing of model
% create script files with dynamic Jacobians, which can be evaluated numerically to compute the policy function
MODEL = feval(str2func([OPT.modelname '_preprocessing']));
% get parameter information
[PARAM, ESTIM_PARAM, MODEL] = feval(str2func([OPT.modelname '_params']),MODEL);

% Create index for estimated params
MODEL.param_estim_names = fieldnames(ESTIM_PARAM);
MODEL.param_estim_nbr = size(MODEL.param_estim_names,1);
MODEL.param_estim_index = [];
xparam0 = nan(MODEL.param_estim_nbr,1);
for j = 1:MODEL.param_estim_nbr
    MODEL.param_estim_index = [MODEL.param_estim_index find(ismember(fieldnames(PARAM),MODEL.param_estim_names(j)))];
    OPT.optimizer.bounds.lb(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){2};
    OPT.optimizer.bounds.ub(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){3};
    xparam0(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){1};    
end


%% Create data matrix and selection matrix H
load(OPT.datafile,'DATA');
% create data matrix corresponding to MODEL.varobs_idx_DR
DATAMAT = nan(OPT.nobs,MODEL.varobs_nbr);
for j=1:MODEL.varobs_nbr
    DATAMAT(:,j) = DATA.(MODEL.varobs{j});
end
DATAMAT = DATAMAT(OPT.first_obs:OPT.nobs,:);
% create selection matrix H
MODEL.varobs_selection_matrix = zeros(MODEL.varobs_nbr,MODEL.endo_nbr);
for j=1:MODEL.varobs_nbr
    MODEL.varobs_selection_matrix(j,MODEL.varobs_idx_DR(j)) = 1;
end
% Check stochastic singularity
if MODEL.varobs_nbr > (MODEL.exo_nbr + nnz(diag(MODEL.Sigma_e)))
    error('Stochastic Singularity: You have more observables than shocks+measurement errors.')
end

%% check objective function at initial parameters
if OPT.optimizer.randomize_initval
    error_indicator=1;
    fprintf('Try drawing a valid initial vector:\n')
    while error_indicator~=0
        xparam_rnd = mvnrnd(xparam0,eye(length(xparam0)));
        if all(xparam_rnd(:)>OPT.optimizer.bounds.lb) && all(xparam_rnd(:)<OPT.optimizer.bounds.ub)
            for jp = 1:MODEL.param_estim_nbr
                PARAM.(MODEL.param_estim_names{jp}) = xparam_rnd(jp);
            end
            [~,error_indicator] = get_first_order_perturbation_solution(MODEL,PARAM);
        end
    end
    xparam0 = xparam_rnd(:);
    disp(xparam0);
end

if OPT.optimizer.bounds.do_param_transform
    xparam0 = param_transform_unbounded(xparam0,OPT.optimizer.bounds.lb,OPT.optimizer.bounds.ub); % this transforms into unbounded parameters, needs to be undone in objective function before we compute the solution
end

% compute, check and print objective function
tic_id = tic;
[neg_log_likelihood,exit_flag] = dsge_negloglikelihood(xparam0,DATAMAT,PARAM,MODEL,OPT);
elapsed_time = toc(tic_id);
if exit_flag ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
end
fprintf('Initial value of the log-likelihood function: %6.4f \n', -1*neg_log_likelihood);
fprintf('Time required to compute log-likelihood function once: %5.4f seconds \n', elapsed_time);


%% Maximum Likelihood Estimation with different optimizers
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
objective_function = @(x) dsge_negloglikelihood(x,DATAMAT,PARAM,MODEL,OPT);

for jopt = 1:size(OPT.optimizer.name,2)
    clc;
    tic_opt = tic;
    switch OPT.optimizer.name(jopt)
        case 'fminsearch' % does not care about bounds
            [xparam1,fval] = fminsearch(objective_function,xparam0,OPT.optimizer.optim_options);
        case 'fminunc' % does not care about bounds
            [xparam1,fval] = fminunc(objective_function,xparam0,OPT.optimizer.optim_options);
        case 'simulannealbnd' % does care about bounds
            [xparam1,fval] = simulannealbnd(objective_function,xparam0,LB,UB,OPT.optimizer.optim_options);
        case 'patternsearch' % does care about bounds
            [xparam1,fval] = patternsearch(objective_function,xparam0,[],[],[],[],LB,UB,[],OPT.optimizer.optim_options);
        case 'csminwel' % does not care about bounds
            [fval,xparam1] = csminwel(objective_function,xparam0,1e-4*eye(MODEL.param_estim_nbr),[],OPT.optimizer.optim_options.TolFun,OPT.optimizer.optim_options.MaxIter);            
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
            [~, ~, ~, ~, ~, BESTEVER] = cmaes('dsge_negloglikelihood',xparam0,cmaesSIGMA,cmaesOptions,...
                                                                      DATAMAT, PARAM, MODEL, OPT);
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
    hess = get_hessian('dsge_negloglikelihood',xopt,[1e-3;1.0],    DATAMAT, PARAM, MODEL, OPT);
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
    input(sprintf('Finished optimization with %s. Hit enter for next optimizer',OPT.optimizer.name(jopt)));
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