function q_alph = gaussianQuantile(alph, mu, Sigma, cdfmvna_fct)
% function q_alph = gaussianQuantile(alph, mu, Sigma, cdfmvna_fct)
% -------------------------------------------------------------------------
% Evaluates the quantile of a N(mu,Sigma) distributed random variable
% for the cumulative probability alph in the interval [0,1]
% Idea:
% - q_alph = min{x: F(x) >= alph}
% - We need to solve F_X(q_alpha) = alpha, to find alpha quantile q_alpha
% - univariate: use norminv
% - multivariate: minimization problem: q_alpha = argmin(log(F_X(q_alpha)) - log(alpha))
% -------------------------------------------------------------------------
% INPUTS
% - alph          [0<alph<1] cumulative probability alph in the interval [0,1]
% - mu            [p by 1]   1st parameter of the Gaussian distribution (mean)
% - Sigma         [p by p]   2nd parameter of the CSN distribution (covariance)
% - cdfmvna_fct   [string]   name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% -------------------------------------------------------------------------
% OUTPUTS
% - q_alph        [p by 1]   alpha quantile vector of the CSN(mu,Sigma,Gamma,nu,Delta) distribution
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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

if nargin < 4
    cdfmvna_fct = 'logmvncdf_ME';
end
if strcmp(cdfmvna_fct,'qsilatmvnv')
    warning('csnQuantile: qsilatmvnv not yet working, switching to ''logmvncdf_ME'' to compute the quantile');
    cdfmvna_fct = 'logmvncdf_ME';
end
if strcmp(cdfmvna_fct,'qsimvnv')
    warning('csnQuantile: qsimvnv not yet working, switching to ''logmvncdf_ME'' to compute the quantile');
    cdfmvna_fct = 'logmvncdf_ME';
end

p = size(mu,1);
if p == 1
    q_alph = norminv(alph,mu,sqrt(Sigma));
else
    % run minimization
    optim_opt = optimset('MaxIter', 10000, 'MaxFunEvals', 10000, 'Display', 'off');
    % prepare helper function to get quantile
    if cdfmvna_fct == "logmvncdf_ME"
        % requires zero mean and correlation matrix as inputs
        normalization = diag(1 ./ sqrt(diag(Sigma)));
        Corr_mat      = normalization * Sigma * normalization; 
        Corr_mat      = 0.5 * (Corr_mat + Corr_mat');
        fun = @(x) ( logmvncdf_ME((normalization * (x - mu)), Corr_mat) - log(alph) )^2;
    elseif cdfmvna_fct == "mvncdf"
        fun = @(x) ( log( mvncdf(x,mu,Sigma) ) - log(alph) )^2;
    elseif cdfmvna_fct == "qsilatmvnv"
        error('csnQuantile: qsilatmvnv not yet working')
        fun = @(x) ( log( qsilatmvnv(1000*p, Sigma, repmat(-Inf,p,1)- mu, x-mu) ) - log(alph) )^2;
    elseif cdfmvna_fct == "qsimvnv"
        error('csnQuantile: qsimvnv not yet working')
        fun = @(x) ( log( qsimvnv(1000*p, Sigma, repmat(-Inf,p,1)- mu, x-mu) ) - log(alph) )^2;
    end
    
    q_alph = fminunc(fun, mu, optim_opt);
end

end % gaussianQuantile