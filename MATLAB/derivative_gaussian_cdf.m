function [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = derivative_gaussian_cdf(mu, Sigma, which_element, cdfmvna_fct)
% function [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = cdfder(mu, Sigma, which_element, cdfmvna_fct)
% -------------------------------------------------------------------------
% Differentiates the multivariate normal cdf according to steps given in
% Rodenburger (2015, p.28-29) - On Approximations of Normal Integrals in the Context of Probit based Discrete Choice Modeling (2015)
% - First derivative is equal to term1 * term2
% - Output variables "evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr"
%   are used for differentiating the first derivative for the second time,
%   for details look at the double loop in the code of csnVar() function above
% -------------------------------------------------------------------------
% INPUTS
% - mu              [p by 1]   mean of normal distribution
% - Sigma           [p by p]   covariance matrix of normal distribution
% - which_element   [scalar]   indicator (1,...,p) for which variable in multivariate normal cdf the derivative is computed for
% - cdfmvna_fct     [string]   name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% -------------------------------------------------------------------------
% OUTPUTS
% - term1           [double]           evaluated normal PDF, see the first part of the derivative
% - term2           [double]           evaluated normal CDF, see the second part of the derivative
% - evalp_t1        [1 by 1]           evaluation point of normal PDF in term1
% - covar_t1        [1 by 1]           variance term of normal PDF in term1
% - evalp_t2        [(p-1) by 1]       evaluation point of normal CDF in term2
% - covar_t2        [(p-1) by (p-1)]   covariance matrix of normal CDF in term2
% - mult_matr       [(p-1) by 1]       auxiliary expression Sigma_12/Sigma_22, this expression shows up when evaluating conditional mean and conditional Variance out of partitions of mean and Variance of joint distribution
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov
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

if nargin < 4
    cdfmvna_fct = 'logmvncdf_ME';
end

p = length(mu);
logic_vec = true(p, 1);
logic_vec(which_element, 1) = false;

% permute the matrices
mu2        = zeros(p, 1); 
mu2(1:p-1) = mu(logic_vec);
mu2(p)     = mu(which_element);

Sigma2           = zeros(size(Sigma));
Sigma2(1:p-1, :) = Sigma(logic_vec, :);
Sigma2(p, :)     = Sigma(which_element, :);
Sigma2           = [Sigma2(:, logic_vec), Sigma2(:, which_element)];

% Conditional mean and conditional var-cov matrix
eval_point = -mu2;
mult_matr  = Sigma2(1:p-1, p) / Sigma2(p, p);
cond_mean  = mult_matr * eval_point(p);
cond_var   = Sigma2(1:p-1, 1:p-1) - mult_matr * Sigma2(p, 1:p-1);

evalp_t1 = eval_point(p);
covar_t1 = Sigma2(p, p);

evalp_t2 = eval_point(1:p-1) - cond_mean;
covar_t2 = cond_var;

% Evaluate the first term in the derivative
term1 = normpdf(evalp_t1, 0, sqrt(covar_t1));

% Evaluate the second term in the derivative
if isempty(cond_var)
    term2 = 1;
else
    if cdfmvna_fct == "logmvncdf_ME"
        stdnrd   = diag(1./sqrt(diag(cond_var)));
        haspl    = stdnrd * evalp_t2;
        Uytgesme = stdnrd * cond_var * stdnrd; 
        Uytgesme = 0.5 * (Uytgesme + Uytgesme');
        term2    = exp(logmvncdf_ME(haspl, Uytgesme));
    elseif cdfmvna_fct == "mvncdf"
        term2 = mvncdf(eval_point(1:p-1)', cond_mean', cond_var);
    elseif cdfmvna_fct == "qsilatmvnv"
        if (p-1) > 1
            term2 = qsilatmvnv( 1000*(p-1), cond_var, repmat(-Inf,p-1,1)-cond_mean, eval_point(1:p-1)-cond_mean );
        else
            term2 = normcdf(eval_point(1:p-1), cond_mean, cond_var);
        end
    elseif cdfmvna_fct == "qsimvnv"
        if (p-1) > 1
            term2 = qsimvnv( 1000*(p-1), cond_var, repmat(-Inf,p-1,1)-cond_mean, eval_point(1:p-1)-cond_mean );
        else
            term2 = normcdf(eval_point(1:p-1), cond_mean, cond_var);
        end
    end
end

end % derivative_gaussian_cdf