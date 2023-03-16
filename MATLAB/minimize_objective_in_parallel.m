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
    elseif optim_names(jopt)=="cmaes"
        cmaesOptions = cmaes('defaults');
        cmaesOptions.SaveVariables='off'; cmaesOptions.DispFinal='on'; cmaesOptions.WarnOnEqualFunctionValues='no'; cmaesOptions.DispModulo='300'; cmaesOptions.LogModulo='0'; cmaesOptions.LogTime='0'; cmaesOptions.Resume = 0;
        cmaesOptions.LBounds = lb; cmaesOptions.UBounds = ub;
        % Set default search volume (SIGMA)
        cmaesSIGMA = (ub-lb)*0.2; cmaesSIGMA(~isfinite(cmaesSIGMA)) = 0.01;
        while max(cmaesSIGMA)/min(cmaesSIGMA)>1e6 %make sure initial search volume (SIGMA) is not badly conditioned
            cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA))=0.9*cmaesSIGMA(cmaesSIGMA==max(cmaesSIGMA));
        end
        cmaesOptions.TolFun = optim_opt.TolFun; cmaesOptions.MaxIter = optim_opt.MaxIter; cmaesOptions.MaxFunEvals = optim_opt.MaxFunEvals;
        [~, ~, ~, ~, ~, BESTEVER] = cmaes(func2str(objfct), xparamsInit, cmaesSIGMA, cmaesOptions, varargin{:});
        XPARAMS(:,jopt) = BESTEVER.x; FVALS(jopt) = BESTEVER.f;
    end
end
[~, bestjopt] = sort(FVALS);