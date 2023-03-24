function [XPARAMS, FVALS, bestjopt] = minimize_objective_in_parallel(objfct, xparamsInit, lb, ub, optim_opt, varargin)
objfct_ = @(x) objfct(x,varargin{:});
tic_id = tic;
[fval_init,exitflag_init] = objfct_(xparamsInit);
elapsed_time = toc(tic_id);
if exitflag_init ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
else
    fprintf('Initial value of the objective function: %6.4f \nTime required to compute objective function once: %s \n', fval_init,dynsec2hms(elapsed_time));
end

XPARAMS = nan(length(xparamsInit),length(optim_opt.names));
FVALS = nan(1,length(optim_opt.names));
warning('off','MATLAB:illConditionedMatrix');
parfor jopt=1:length(optim_opt.names)
    if optim_opt.names(jopt)=="fminsearch"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminsearch(objfct_, xparamsInit, optim_opt);
    elseif optim_opt.names(jopt)=="fminsearchbnd"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminsearchbnd(objfct_, xparamsInit, lb, ub, optim_opt);
    elseif optim_opt.names(jopt)=="fminunc"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminunc(objfct_, xparamsInit, optim_opt);
    elseif optim_opt.names(jopt)=="fmincon"
        [XPARAMS(:,jopt), FVALS(jopt)] = fmincon(objfct_, xparamsInit, [], [], [], [], lb, ub, [], optim_opt);
    elseif optim_opt.names(jopt)=="simulannealbnd"
        [XPARAMS(:,jopt), FVALS(jopt)] = simulannealbnd(objfct_, xparamsInit, lb, ub, optim_opt);
    elseif optim_opt.names(jopt)=="cmaes"
        cmaesOptions = cmaes('defaults');
        cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='300'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
        cmaesOptions.LBounds = lb; cmaesOptions.UBounds = ub;
        % Set default search volume (INSIGMA)
        cmaesINSIGMA = (ub-lb)*0.2; cmaesINSIGMA(~isfinite(cmaesINSIGMA)) = 0.01;
        while max(cmaesINSIGMA)/min(cmaesINSIGMA)>1e6 %make sure initial search volume (INSIGMA) is not badly conditioned
            cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA))=0.9*cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA));
        end
        if ~isempty(optim_opt.TolFun); cmaesOptions.TolFun = optim_opt.TolFun; end
        if ~isempty(optim_opt.MaxIter); cmaesOptions.MaxIter = optim_opt.MaxIter; end
        if ~isempty(optim_opt.MaxFunEvals); cmaesOptions.MaxFunEvals = optim_opt.MaxFunEvals; end
        [~, ~, ~, ~, ~, BESTEVER] = cmaes(func2str(objfct), xparamsInit, cmaesINSIGMA, cmaesOptions, varargin{:});
        XPARAMS(:,jopt) = BESTEVER.x; FVALS(jopt) = BESTEVER.f;
    elseif optim_opt.names(jopt)=="cmaes_dsge"
        cmaesOptions = cmaes_dsge('defaults');
        cmaesOptions.Saving='off'; cmaesOptions.Display ='on'; cmaesOptions.Plotting='off'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.VerboseModulo='300'; cmaesOptions.Resume = 'no';
        cmaesOptions.LBounds = lb; cmaesOptions.UBounds = ub;
        % Set default search volume (INSIGMA)
        cmaesSigma = 1;
        cmaesINSIGMA = (ub-lb)*0.2; cmaesINSIGMA(~isfinite(cmaesINSIGMA)) = 0.01;
        while max(cmaesINSIGMA)/min(cmaesINSIGMA)>1e6 %make sure initial search volume (INSIGMA) is not badly conditioned
            cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA))=0.9*cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA));
        end
        if ~isempty(optim_opt.TolFun); cmaesOptions.TolFun = optim_opt.TolFun; end
        if ~isempty(optim_opt.MaxIter); cmaesOptions.MaxIter = optim_opt.MaxIter; end
        if ~isempty(optim_opt.MaxFunEvals); cmaesOptions.MaxFunEvals = optim_opt.MaxFunEvals; end
        [~, ~, ~, ~, ~, BESTEVER] = cmaes_dsge(func2str(objfct), xparamsInit, cmaesSigma, cmaesINSIGMA, cmaesOptions, varargin{:});
        XPARAMS(:,jopt) = BESTEVER.x; FVALS(jopt) = BESTEVER.f;
    elseif optim_opt.names(jopt)=="sa_resampling"
        sa_opt_max = 0;            % we do minimization
        sa_opt_rt = 0.85;          % The reduction rate in temperature
        sa_opt_eps = 1e-6;         % Function tolerance
        sa_opt_ns = 10;            % Number of cycles
        sa_opt_nt = 20;            % Number of random walkers
        sa_opt_neps = 4;           % Stopping criteria
        sa_opt_maxevl = 1e5;       % Maximum number of function evaluations
        sa_opt_maxresample = 15;   % Maximum number of times we resample a new step size
        sa_opt_maxNaNJump = 5;     % Maximum number of times we jump for a given step size
        sa_opt_c = 2*ones(length(xparamsInit),1); % Default setting for adjusting vm - the stepsizes
        sa_opt_iprint = 0;         % Printing property
        sa_opt_t = 100;            % Initial temperature
        sa_opt_vm = (ub-lb)*0.2; sa_opt_vm(~isfinite(sa_opt_vm)) = 0.01; % initial stepsize
        while max(sa_opt_vm)/min(sa_opt_vm)>1e6
            sa_opt_vm(sa_opt_vm==max(sa_opt_vm))=0.9*sa_opt_vm(sa_opt_vm==max(sa_opt_vm));
        end
        sa_opt_Save_To_File = [];
        if ~isempty(optim_opt.TolFun); sa_opt_eps = optim_opt.TolFun; end
        if ~isempty(optim_opt.MaxFunEvals); sa_opt_maxevl = optim_opt.TolFun; end
        [~,~,~,XPARAMS(:,jopt),FVALS(jopt),~,~,~,~] = sa_resampling(objfct_,length(xparamsInit),xparamsInit,sa_opt_max,sa_opt_rt,sa_opt_eps,sa_opt_ns,sa_opt_nt,sa_opt_neps,sa_opt_maxevl,sa_opt_maxresample,sa_opt_maxNaNJump,lb,ub,sa_opt_c,sa_opt_iprint,sa_opt_t,sa_opt_vm,sa_opt_Save_To_File);
    end
end
warning('on','MATLAB:illConditionedMatrix');
[~, bestjopt] = sort(FVALS);