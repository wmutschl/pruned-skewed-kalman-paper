function [E_csn] = csnMean(mu, Sigma, Gamma, nu, Delta, cdfmvna_fct)
% function [E_csn] = csnMean(mu, Sigma, Gamma, nu, Delta, cdfmvna_fct)
% -------------------------------------------------------------------------
% Evaluates the unconditional expectation vector of a CSN(mu,Sigma,Gamma,nu,Delta)
% distributed random variable based on expressions derived in
% - Dominguez-Molina, Gonzalez-Farias, Gupta (2003, sec. 2.4) - The multivariate closed skew normal distribution
% - Rodenburger (2015) - On Approximations of Normal Integrals in the Context of Probit based Discrete Choice Modeling (2015)
% -------------------------------------------------------------------------
% INPUTS
% - mu            [p by 1]   1st parameter of the CSN distribution (location)
% - Sigma         [p by p]   2nd parameter of the CSN distribution (scale)
% - Gamma         [q by p]   3rd parameter of the CSN distribution (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu            [q by 1]   4th parameter of the CSN distribution (enables closure of CSN distribution under conditioning)
% - Delta         [q by q]   5th parameter of the CSN distribution (enables closure of CSN distribution und marginalization)
% - cdfmvna_fct   [string]   name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% -------------------------------------------------------------------------
% OUTPUTS
% - E_csn         [p by 1]   expectation vector of the CSN(mu,Sigma,Gamma,nu,Delta) distribution
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

if nargin < 6
    cdfmvna_fct = 'logmvncdf_ME';
end
q       = size(Delta, 1); % skewness dimension
derVect = nan(q, 1);      % initialize gradient vector
Delta1  = Delta + Gamma * Sigma * Gamma'; % var-cov matrix inside the expression "psi" of equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)

% differentiate the normal cdf, steps are given in Rodenburger (2015, p.28-29)
for jj = 1:q
    % permutation of var-cov matrix
    logi   = true(q, 1); logi(jj, 1)=false;
    nu2    = nan(q, 1); nu2(1:q-1)=nu(logi); nu2(q)=nu(jj);
    Delta2 = nan(size(Delta1));
    Delta2(1:q-1, :) = Delta1(logi, :);
    Delta2(q, :) = Delta1(jj, :);
    Delta2 = [Delta2(:, logi), Delta2(:, jj)];

    % conditional mean and var-cov matrix
    nu2      = -nu2;
    condMean = (Delta2(1:q-1, q)*nu2(q))/Delta2(q, q);
    condVar  = Delta2(1:q-1, 1:q-1) - (Delta2(1:q-1, q)/Delta2(q, q))*Delta2(q, 1:q-1);

    % evaluate the gradient vector
    term1 = normpdf(nu2(q), 0, sqrt(Delta2(q, q)));
    if isempty(condVar)
        term2 = 1;
    else
        % different functions to evaluate log Gaussian cdf
        if cdfmvna_fct == "logmvncdf_ME"
            % requires zero mean and correlation matrix as inputs
            normalization2 = diag(1./sqrt(diag(condVar)));
            eval_point2 = normalization2*(nu2(1:q-1) - condMean);
            Corr_mat2 = normalization2*condVar*normalization2; 
            Corr_mat2 = 0.5*(Corr_mat2 + Corr_mat2');
            term2 = exp(logmvncdf_ME(eval_point2, Corr_mat2));
        elseif cdfmvna_fct == "mvncdf"
            term2 = mvncdf(nu2(1:q-1)', condMean', condVar);
        elseif cdfmvna_fct == "qsilatmvnv"
            if (q-1) > 1
                term2 = qsilatmvnv( 1000*(q-1), condVar, repmat(-Inf,q-1,1)-condMean, nu2(1:q-1)-condMean );
            else
                term2 = normcdf(nu2, condMean, condVar);
            end
        elseif cdfmvna_fct == "qsimvnv"
            if (q-1) > 1
                term2 = qsimvnv( 1000*(q-1), condVar, repmat(-Inf,q-1,1)-condMean, nu2(1:q-1)-condMean );
            else
                term2 = normcdf(nu2, condMean, condVar);
            end
        end
    end
    derVect(jj) = term1*term2;
end

% evaluate end result using equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)
Delta1 = 0.5*(Delta1+Delta1');
if isempty(Gamma) || isempty(nu) || isempty(Delta1)
    if isempty(Gamma) && isempty(nu) && isempty(Delta1)
        term3 = 1;
    else
        error("Problem with Gamma, nu, Delta being empty / not empty")
    end
else
    if cdfmvna_fct == "logmvncdf_ME"
        % requires zero mean and correlation matrix as inputs
        normalization3 = diag(1./sqrt(diag(Delta1)));
        eval_point3 = normalization3*(zeros(q, 1) - nu);
        Corr_mat3 = normalization3*Delta1*normalization3;
        Corr_mat3 = 0.5*(Corr_mat3 + Corr_mat3');
        term3 = exp(logmvncdf_ME(eval_point3, Corr_mat3));
    elseif cdfmvna_fct == "mvncdf"
        term3 = mvncdf(zeros(1, q), nu', Delta1);
    elseif cdfmvna_fct == "qsilatmvnv"
        if q > 1
            term3 = qsilatmvnv( 1000*q, Delta1, repmat(-Inf,q,1)-nu, zeros(q, 1)-nu );
        else
            term3 = normcdf(0, nu, Delta1);
        end
    elseif cdfmvna_fct == "qsimvnv"
        if q > 1
            term3 = qsimvnv( 1000*q, Delta1, repmat(-Inf,q,1)-nu, zeros(q, 1)-nu );
        else
            term3 = normcdf(0, nu, Delta1);
        end    
    end
end
E_csn = mu + Sigma * Gamma' * (derVect/term3); % equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)

end