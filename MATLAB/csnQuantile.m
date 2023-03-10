function [q_alph] = csnQuantile(alph, mu, Sigma, Gamma, nu, Delta, cdfmvna_fct, optim_fct)
% function [q_alph] = csnQuantile(alph, mu, Sigma, Gamma, nu, Delta, cdfmvna_fct, optim_fct)
% -------------------------------------------------------------------------
% Evaluates the quantile of a CSN(mu,Sigma,Gamma,nu,Delta) distributed random variable
% for the cumulative probability alph in the interval [0,1]
% Idea:
% - q_alph = min{x: F(x) >= alph}
% - Lemma 2.2.1 of chapter 2 of Genton (2004) - "Skew-elliptical distributions
%   and their applications: a journey beyond normality" gives the cdf, F_X(x)
%   of a CSN distributed random variable X
% - We need to solve F_X(q_alpha) = alpha, to find alpha quantile q_alpha
% - Minimization problem: q_alpha = argmin(log(F_X(q_alpha)) - log(alpha))
% -------------------------------------------------------------------------
% INPUTS
% - alph          [0<alph<1] cumulative probability alph in the interval [0,1]
% - mu            [p by 1]   1st parameter of the CSN distribution (location)
% - Sigma         [p by p]   2nd parameter of the CSN distribution (scale)
% - Gamma         [q by p]   3rd parameter of the CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu            [q by 1]   4th parameter of the CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta         [q by q]   5th parameter of the CSN distribution (enables closure of CSN distribution und marginalization)
% - cdfmvna_fct   [string]   name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% - optim_fct     [string]   name of minimization function, possible values: 'fminunc', 'lsqnonlin'
% -------------------------------------------------------------------------
% OUTPUTS
% - q_alph        [p by 1]   alpha quantile vector of the CSN(mu,Sigma,Gamma,nu,Delta) distribution
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================

if nargin < 7 || (nargin >= 7 && isempty(cdfmvna_fct))
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
if nargin < 8
    optim_fct = 'lsqnonlin';
end

q = size(Delta, 1);
Cov_mat = -Sigma * Gamma';

% evaluate the first log-CDF
Var_Cov1 = Delta - Gamma * Cov_mat;        
if cdfmvna_fct == "logmvncdf_ME"
    normalization1 = diag(1 ./ sqrt(diag(Var_Cov1)));
    eval_point1    = normalization1 * (zeros(q, 1) - nu);
    Corr_mat1      = normalization1 * Var_Cov1 * normalization1; 
    Corr_mat1      = 0.5*(Corr_mat1 + Corr_mat1');
    cdf1           = logmvncdf_ME(eval_point1, Corr_mat1);
elseif cdfmvna_fct == "mvncdf"
    cdf1 = log(mvncdf(zeros(1, q), nu', Var_Cov1));
elseif cdfmvna_fct == "qsilatmvnv"
    if q > 1
        cdf1 = log(qsilatmvnv( 1000*q, Var_Cov1, repmat(-Inf,q,1)-nu, zeros(q, 1)-nu ));
    else
        cdf1 = log(normcdf(0, nu, Var_Cov1));
    end
elseif cdfmvna_fct == "qsimvnv"
    if q > 1
        cdf1 = log(qsimvnv( 1000*q, Var_Cov1, repmat(-Inf,q,1)-nu, zeros(q, 1)-nu ));
    else
        cdf1 = log(normcdf(0, nu, Var_Cov1));
    end
end

% prepare helper function to get quantile
Var_Cov2 = [Sigma, Cov_mat; Cov_mat', Var_Cov1];

if cdfmvna_fct == "logmvncdf_ME"
    % requires zero mean and correlation matrix as inputs
    normalization2 = diag(1 ./ sqrt(diag(Var_Cov2)));
    Corr_mat2      = normalization2 * Var_Cov2 * normalization2;
    Corr_mat2      = 0.5 * (Corr_mat2 + Corr_mat2');
    if strcmp(optim_fct,'lsqnonlin')
        fun = @(x) logmvncdf_ME( (normalization2 * ([x; zeros(q, 1)] - [mu; nu])), Corr_mat2 ) - cdf1 - log(alph);
    elseif strcmp(optim_fct,'fminunc')
        fun = @(x) ( logmvncdf_ME( (normalization2 * ([x; zeros(q, 1)] - [mu; nu])), Corr_mat2 ) - cdf1 - log(alph))^2;
    end
elseif cdfmvna_fct == "mvncdf"
    if strcmp(optim_fct,'lsqnonlin')
        fun = @(x) log(mvncdf([x', zeros(1, q)], [mu', nu'], Var_Cov2)) - cdf1 - log(alph);
    elseif strcmp(optim_fct,'fminunc')
        fun = @(x) (log(mvncdf([x', zeros(1, q)], [mu', nu'], Var_Cov2)) - cdf1 - log(alph))^2;
    end
elseif cdfmvna_fct == "qsilatmvnv"
    error('csnQuantile: qsilatmvnv not yet working')
    fun = @(pp_quant) (log(qsilatmvnv(1000*(size(pp_quant,1)+q),Var_Cov2,repmat(-Inf,size(pp_quant,1)+q,1)-[mu', nu']', [pp_quant', zeros(1, q)]' - [mu', nu']')) ...
                       - cdf1 ...
                       - log(alph))^2;
elseif cdfmvna_fct == "qsimvnv"
    error('csnQuantile: qsimvnv not yet working')
    fun = @(pp_quant) (log(qsimvnv(1000*(size(pp_quant,1)+q),Var_Cov2,repmat(-Inf,size(pp_quant,1)+q,1)-[pp_quant', zeros(1, q)],[pp_quant', zeros(1, q)] - [mu', nu'])) ...
                       - cdf1 ...
                       - log(alph))^2;
end

% run minimization
optim_opt = optimset('MaxIter', 10000, 'MaxFunEvals', 10000, 'Display', 'off');
if strcmp(optim_fct,'lsqnonlin')
    q_alph = lsqnonlin(fun, mu, [], [], optim_opt);
elseif strcmp(optim_fct,'fminunc')
    q_alph = fminunc(fun, mu, optim_opt);
end

end % csnQuantile