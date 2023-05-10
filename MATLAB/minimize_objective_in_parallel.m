function [XPARAMS, FVALS, bestjopt] = minimize_objective_in_parallel(objfct, xparamsInit, lb, ub, optim_opt, varargin)
% [XPARAMS, FVALS, bestjopt] = minimize_objective_in_parallel(objfct, xparamsInit, lb, ub, optim_opt, varargin)
% -------------------------------------------------------------------------
% minimizes an objective function using different optimizers in parallel
% -------------------------------------------------------------------------
% INPUTS
% - objfct          [function handle]  name of objective function
% - xparamsInit     [vector]           initial values
% - lb              [vector]           lower bounds
% - ub              [vector]           upper bounds
% - optim_opt       [structure]        options passed to optimizers
% - varargin        [cell]             additional arguments passed to objfct
% -------------------------------------------------------------------------
% OUTPUTS
% - XPARAMS         [matrix]           optimized values for each optimizers
% - FVALS           [matrix]           optimized values of objective function for each optimizers
% - bestjopt        [integer]          indicator which optimizer minimized objective best
% =========================================================================
% Copyright (C) 2023 Willi Mutschler
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
objfct_ = @(x) objfct(x,varargin{:}); % pass additional arguments to objective function
tic_id = tic;
[fval_init,exitflag_init] = objfct_(xparamsInit); % evaluate objfct at initial values to check whether everything works
elapsed_time = toc(tic_id);
if exitflag_init ~= 1
    error('Something wrong with log-likelihood function at initial parameters')
else
    fprintf('Initial value of the objective function: %6.4f \nTime required to compute objective function once: %s \n', fval_init,dynsec2hms(elapsed_time));
end

XPARAMS = nan(length(xparamsInit),length(optim_opt.names)); % initialize storage
FVALS = nan(1,length(optim_opt.names));                     % initialize storage
parfor jopt=1:length(optim_opt.names)
    warning('off','MATLAB:illConditionedMatrix');
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
        % create option structure for cmaes as this is an external optimizer
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
        % create option structure for cmaes_dsge as this is an external optimizer
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
 end
    warning('on','MATLAB:illConditionedMatrix');
end
[~, bestjopt] = sort(FVALS);