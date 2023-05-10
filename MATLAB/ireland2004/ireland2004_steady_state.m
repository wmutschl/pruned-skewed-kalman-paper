function [ys_, params_, error_indicator] = ireland2004_steady_state(ys_,exo_,M_,options_)
% function [ys, params_, error_indicator] = ireland2004_steady_state(ys_,exo_,M_,options_)
% -------------------------------------------------------------------------
% computes the steady-state of the small scale New Keynesian model of Ireland (2004): 
% "Technology Shocks in The New Keynesian Model", The Review of Economics and Statistics
% based on Dynare replication codes kindly provided by Johannes Pfeifer at
% https://github.com/JohannesPfeifer/DSGE_mod/tree/master/Ireland_2004
% -------------------------------------------------------------------------
% INPUTS
% - ys         [endo_nbr by 1]   initial values for steady-state
% - exo        [exo_nbr by 1]    initial values for exogenous
% - M_         [structure]       information on the model (inspired by Dynare's M_ structure)
% - options_   [structure]       options (inspired by Dynare's options_ structure)
% -------------------------------------------------------------------------
% OUTPUTS
% - ys         [endo_nbr by 1]    computed values for steady-state
% - params     [param_nbr by 1]   updated parameter values
% - error_indicator [boolean]     0: no error
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

error_indicator = 0; % initialize no error

% read-out parameters
for iter = 1:M_.param_nbr
    paramname = M_.param_names{iter};
    eval([ paramname ' = M_.params(' int2str(iter) ');']);
end

% compute steady-state of exogenous
eta_a = 0; eta_e = 0; eta_z = 0; eta_r = 0;

% compute steady-state of endogenous
ahat = 0; ehat = 0; zhat = 0; xhat = 0; pihat = 0; yhat = 0; ghat = 0; rhat = 0;

% write out parameters
params_=NaN(M_.param_nbr,1);
for iter = 1:length(M_.params)
    eval([ 'params_(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

% write out endogenous
for iter = 1:M_.endo_nbr
    varname = M_.endo_names{iter};
    eval(['ys_(' int2str(iter) ',1) = ' varname ';']);
end

% check residuals
if any( feval(str2func(M_.fname + '_f'),ys_,params_) > 1e-7 )
    error_indicator = 1;
    warning('Something wrong with the steady-state computations as the residuals are nonzero');
    return
end

end% main function end