function V_csn = csnVar(Sigma, Gamma, nu, Delta, cdfmvna_fct)
% function V_csn = csnVar(Sigma, Gamma, nu, Delta, cdfmvna_fct)
% -------------------------------------------------------------------------
% Evaluates the unconditional covariance matrix of a CSN(mu,Sigma,Gamma,nu,Delta)
% distributed random variable based on expressions derived in
% - Dominguez-Molina, Gonzalez-Farias, Gupta (2003, sec. 2.4) - The multivariate closed skew normal distribution
% -------------------------------------------------------------------------
% INPUTS
% - Sigma         [p by p]   2nd parameter of the CSN distribution (scale)
% - Gamma         [q by p]   3rd parameter of the CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu            [q by 1]   4th parameter of the CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta         [q by q]   5th parameter of the CSN distribution (enables closure of CSN distribution und marginalization)
% - cdfmvna_fct   [string]   name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% -------------------------------------------------------------------------
% OUTPUTS
% - V_csn         [p by p]   covariance matrix of the CSN(mu,Sigma,Gamma,nu,Delta) distribution
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

if nargin < 5
    cdfmvna_fct = 'logmvncdf_ME';
end
q = size(Delta, 1); % skewness dimension

% Eq. (1) and (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10-11)
Delta2 = Delta + Gamma * Sigma * Gamma';
Delta2 = 0.5 * (Delta2 + Delta2'); %ensure symmetry
% Evaluate first and second derivatives
first_der  = zeros(q, 1); % first derivatives
second_der = zeros(q, q); % second derivatives
for ii = 1:q
    [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = derivative_gaussian_cdf(nu, Delta2, ii, cdfmvna_fct);
    first_der(ii) = term1 * term2;    
    for jj = 1:q
        if ii ~= jj
            which_element = jj - (ii < jj);
            [expr1, expr2] = derivative_gaussian_cdf(-evalp_t2, covar_t2, which_element, cdfmvna_fct);
            second_der(ii, jj) = term1 * expr1 * expr2;
        else
            expr3 = -evalp_t1 / covar_t1 * term1 * term2;            
            expr4 = zeros(1, q-1);
            for kk = 1:q-1
                [first_term, second_term] = derivative_gaussian_cdf(-evalp_t2, covar_t2, kk, cdfmvna_fct);
                expr4(kk) = first_term * second_term;
            end            
            expr5 = term1 * expr4 * (-mult_matr);
            second_der(jj, jj) = expr3 + expr5;
        end
    end
end

% denominator in psi and Lambda of eq. (1) and (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10-11)
if cdfmvna_fct == "logmvncdf_ME"
    % requires zero mean and correlation matrix as inputs
    normalization3 = diag(1 ./ sqrt(diag(Delta2)));
    eval_point3 = normalization3 * (zeros(q, 1) - nu);
    Corr_mat3 = normalization3 * Delta2 * normalization3;
    Corr_mat3 = 0.5 * (Corr_mat3 + Corr_mat3');
    term3 = exp(logmvncdf_ME(eval_point3, Corr_mat3));
elseif cdfmvna_fct == "mvncdf"
    term3 = mvncdf(zeros(1, q), nu', Delta2);
elseif cdfmvna_fct == "qsilatmvnv"
    if q > 1
        term3 = qsilatmvnv( 1000*q, Delta2, repmat(-Inf,q,1)-nu, zeros(q,1)-nu );
    else
        term3 = normcdf(0, nu, Delta2);
    end
elseif cdfmvna_fct == "qsimvnv"
    if q > 1
        term3 = qsimvnv( 1000*q, Delta2, repmat(-Inf,q,1)-nu, zeros(q,1)-nu );
    else
        term3 = normcdf(0, nu, Delta2);
    end
end

% definition of psi in eq (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)
psi = first_der / term3;
% definition of Lambda in eq (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 11)
Lambda = second_der / term3; Lambda = 0.5 * (Lambda + Lambda');
% eq (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 11)
V_csn = Sigma + Sigma * Gamma' * Lambda * Gamma * Sigma - Sigma * Gamma' * (psi * psi') * Gamma * Sigma;

end %csnVar