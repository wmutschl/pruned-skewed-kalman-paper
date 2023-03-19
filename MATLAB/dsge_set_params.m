function M_ = dsge_set_params(xparams,estim_params_,options_,M_)
% function M_ = set_params_dsge(xparams,estim_params_,options_,M_)
% -------------------------------------------------------------------------
% updates parameters in model structure for a DSGE model
% -------------------------------------------------------------------------
% INPUTS
% - xparams         [double vector]   numerical values for parameters
% - estim_params_   [structure]       information on estimated parameters, inspired by Dynare's estimated_params block
% - options_        [structure]       options, inspired by Dynare's options_ structure
% - M_              [structure]       model information, inspired by Dynare's M_ structure
% -------------------------------------------------------------------------
% OUTPUTS
% - M_              [structure]       model information, inspired by Dynare's M_ structure
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
idxparams=1;

% skewness parameters of shocks (ordered first in xparams)
for jnsx = 1:estim_params_.nsx
    jexo = estim_params_.skew_exo(jnsx,1);
    if ( isfield(options_.parameters.transform,"skew_"+M_.exo_names(jexo)) && (options_.parameters.transform.("skew_"+M_.exo_names(jexo))==1) ) || (isfield(options_.parameters.transform,"diag_Gamma_"+M_.exo_names(jexo)) && (options_.parameters.transform.("diag_Gamma_"+M_.exo_names(jexo))==1))
        parval = transform_to_bounded(xparams(idxparams),estim_params_.skew_exo(jnsx,3),estim_params_.skew_exo(jnsx,4));
    else
        parval = xparams(idxparams);
    end
    if options_.parameters.use_stderr_skew
        M_.Skew_eta(jexo,1) = parval;
    else
        M_.Gamma_eta(jexo,jexo) = parval;
    end
    idxparams=idxparams+1;    
end

% variance parameters of shocks (ordered second in xparams)
for jnvx = 1:estim_params_.nvx
    jexo = estim_params_.var_exo(jnvx,1);
    if ( isfield(options_.parameters.transform,"stderr_"+M_.exo_names(jexo)) && (options_.parameters.transform.("stderr_"+M_.exo_names(jexo))==1) )|| ( isfield(options_.parameters.transform,"sqrt_diag_Sigma_"+M_.exo_names(jexo)) && (options_.parameters.transform.("sqrt_diag_Sigma_"+M_.exo_names(jexo)))==1 )
        parval = transform_to_bounded(xparams(idxparams),estim_params_.var_exo(jnvx,3),estim_params_.var_exo(jnvx,4))^2;
    else
        parval = xparams(idxparams)^2;
    end
    if options_.parameters.use_stderr_skew        
        M_.Cov_eta(jexo,jexo) = parval;
    else
        M_.Sigma_eta(jexo,jexo) = parval;
    end
    idxparams=idxparams+1;
end

% variance parameters of measurement errors (ordered third in xparams)
for jnvn = 1:estim_params_.nvn
    jmeas = estim_params_.var_endo(jnvn,1);
    if isfield(options_.parameters.transform,"stderr_"+M_.varobs(jmeas)) && (options_.parameters.transform.("stderr_"+M_.varobs(jmeas))==1)
        parval = transform_to_bounded(xparams(idxparams),estim_params_.var_endo(jnvn,3),estim_params_.var_endo(jnvn,4))^2;
    else
        parval = xparams(idxparams)^2;
    end
    M_.Sigma_eps(jmeas,jmeas) = parval;
    idxparams=idxparams+1;
end

% model parameters (ordered fourth in xparams)
for jp = 1:estim_params_.np
    jparam = estim_params_.param_vals(jp,1);
    if isfield(options_.parameters.transform,M_.param_names(jparam)) && (options_.parameters.transform.(M_.param_names(jparam))==1)
        parval = transform_to_bounded(xparams(idxparams),estim_params_.param_vals(jp,3),estim_params_.param_vals(jp,4));
    else
        parval = xparams(idxparams);
    end
    M_.params(jparam) = parval;
    idxparams=idxparams+1;
end
if (idxparams-1)~=length(xparams)
    error('something wrong in updating dsge parameters')
end

% update M_ structure, taking into account which parameter transformation is used and that we assume indepent shocks
for jexo = 1:M_.exo_nbr
    if options_.parameters.use_stderr_skew
        if M_.Skew_eta(jexo,1) == 0 % Gaussian
            M_.Sigma_eta(jexo,jexo) = M_.Cov_eta(jexo,jexo);
            M_.Gamma_eta(jexo,jexo) = 0;
        else % CSN
            if abs(M_.Skew_eta(jexo,1)) >= 0.995
                warning('Skewness parameter very close to theoretical bound')
            else
                [M_.Sigma_eta(jexo,jexo),M_.Gamma_eta(jexo,jexo)] = csnVarSkew_To_SigmaGamma_univariate(M_.Cov_eta(jexo,jexo),M_.Skew_eta(jexo,1),1);
            end
        end
    else
        if M_.Gamma_eta(jexo,jexo) == 0 % Gaussian
            M_.Cov_eta(jexo,jexo) = M_.Sigma_eta(jexo,jexo);
            M_.Skew_eta(jexo,1) = 0;
        else % CSN
            M_.Cov_eta(jexo,jexo) = csnVar(M_.Sigma_eta(jexo,jexo),M_.Gamma_eta(jexo,jexo),M_.nu_eta(jexo,1),M_.Delta_eta(jexo,jexo),options_.kalman.csn.cdfmvna_fct);
            M_.Skew_eta(jexo,1) = csnSkewness_univariate(M_.Sigma_eta(jexo,jexo),M_.Gamma_eta(jexo,jexo));
            if abs(M_.Skew_eta(jexo,1)) >= 0.995
                warning('Skewness parameter very close to theoretical bound')
            end
        end
    end
end

if any(diag(M_.Gamma_eta))
    M_.mu_eta  = -csnMean(zeros(M_.exo_nbr,1),M_.Sigma_eta,M_.Gamma_eta,M_.nu_eta,M_.Delta_eta,options_.kalman.csn.cdfmvna_fct); % ensures E[eta]=0
end