function [XPARAMS, FVALS, bestjopt] = minimize_objective_in_parallel(objfct, xparamsInit, lb, ub, optim_names, optim_opt, varargin)
objfct_ = @(x) objfct(x,varargin{:});
tic_id = tic;
[fval_init,exitflag_init] = objfct_(xparamsInit);
elapsed_time = toc(tic_id);
if exitflag_init ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
else
    fprintf('Initial value of the objective function: %6.4f \nTime required to compute objective function once: %s \n', fval_init,dynsec2hms(elapsed_time));
end

XPARAMS = nan(length(xparamsInit),length(optim_names));
FVALS = nan(1,length(optim_names));
parfor jopt=1:length(optim_names)
    if optim_names(jopt)=="fminsearch"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminsearch(objfct_, xparamsInit, optim_opt);
    elseif optim_names(jopt)=="fminsearchbnd"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminsearchbnd(objfct_, xparamsInit, lb, ub, optim_opt);
    elseif optim_names(jopt)=="fminunc"
        [XPARAMS(:,jopt), FVALS(jopt)] = fminunc(objfct_, xparamsInit, optim_opt);
    elseif optim_names(jopt)=="fmincon"
        [XPARAMS(:,jopt), FVALS(jopt)] = fmincon(objfct_, xparamsInit, [], [], [], [], lb, ub, [], optim_opt);
    elseif optim_names(jopt)=="simulannealbnd"       
        [XPARAMS(:,jopt), FVALS(jopt)] = simulannealbnd(objfct_, xparamsInit, lb, ub, optim_opt);
    elseif optim_names(jopt)=="cmaes"
        cmaesOptions = cmaes('defaults');
        cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='300'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
        cmaesOptions.LBounds = lb; cmaesOptions.UBounds = ub;
        % Set default search volume (INSIGMA)
        cmaesINSIGMA = (ub-lb)*0.2; cmaesINSIGMA(~isfinite(cmaesINSIGMA)) = 0.01;
        while max(cmaesINSIGMA)/min(cmaesINSIGMA)>1e6 %make sure initial search volume (INSIGMA) is not badly conditioned
            cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA))=0.9*cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA));
        end
        cmaesOptions.TolFun = optim_opt.TolFun; cmaesOptions.MaxIter = optim_opt.MaxIter; cmaesOptions.MaxFunEvals = optim_opt.MaxFunEvals;
        [~, ~, ~, ~, ~, BESTEVER] = cmaes(func2str(objfct), xparamsInit, cmaesINSIGMA, cmaesOptions, varargin{:});
        XPARAMS(:,jopt) = BESTEVER.x; FVALS(jopt) = BESTEVER.f;
    elseif optim_names(jopt)=="cmaes_dsge"
        cmaesOptions = cmaes_dsge('defaults');
        cmaesOptions.Saving='off'; cmaesOptions.Display ='on'; cmaesOptions.Plotting='off'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.VerboseModulo='300'; cmaesOptions.Resume = 'no';
        cmaesOptions.LBounds = lb; cmaesOptions.UBounds = ub;
        % Set default search volume (INSIGMA)
        cmaesSigma = 1;
        cmaesINSIGMA = (ub-lb)*0.2; cmaesINSIGMA(~isfinite(cmaesINSIGMA)) = 0.01;
        while max(cmaesINSIGMA)/min(cmaesINSIGMA)>1e6 %make sure initial search volume (INSIGMA) is not badly conditioned
            cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA))=0.9*cmaesINSIGMA(cmaesINSIGMA==max(cmaesINSIGMA));
        end
        cmaesOptions.TolFun = optim_opt.TolFun; cmaesOptions.MaxIter = optim_opt.MaxIter; cmaesOptions.MaxFunEvals = optim_opt.MaxFunEvals;
        [~, ~, ~, ~, ~, BESTEVER] = cmaes_dsge(func2str(objfct), xparamsInit, cmaesSigma, cmaesINSIGMA, cmaesOptions, varargin{:});
        XPARAMS(:,jopt) = BESTEVER.x; FVALS(jopt) = BESTEVER.f;
    elseif optim_names(jopt)=="sa_resampling"
        rt = 0.85;                 % The reduction rate in temperature
        ns = 10;                   % Number of cycles
        nt = 20;                   % Number of random walkers
        neps = 4;                  % Stopping criteria
        maxresample = 15;          % Maximum number of times we resample a new step size
        maxNaNJump = 5;            % Maximum number of times we jump for a given step size
        c = 2*ones(length(xparamsInit),1);             % Default setting for adjusting vm - the stepsizes
        iprint = 0;                % Printing property
        t = 100;                   % Initial temperature
        vm = (ub-lb)*0.2; vm(~isfinite(vm)) = 0.01;
        while max(vm)/min(vm)>1e6
            vm(vm==max(vm))=0.9*vm(vm==max(vm));
        end
        Save_To_File = [];
        [~,~,~,XPARAMS(:,jopt),FVALS(jopt),~,~,~,~] = sa_resampling(objfct_,length(xparamsInit),xparamsInit,0,rt,optim_opt.TolFun,ns,nt,neps,optim_opt.MaxFunEvals,maxresample,maxNaNJump,lb,ub,c,iprint,t,vm,Save_To_File);
    end
end
[~, bestjopt] = sort(FVALS);