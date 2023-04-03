function [estim_params_,xparams,bounds] = ireland2004_estim_params(STAGE, M_, options_)
% function [estim_params_,xparams,bounds] = ireland2004_estim_params(STAGE, M_, options_)
% -------------------------------------------------------------------------
% information on which parameters to estimate; inspired by Dynare's estimated_params block
% -------------------------------------------------------------------------
% INPUTS
% - STAGE      [integer]     switch between different parameter settings for sophisticated search for initial values
% - M_         [structure]   model information, inspired by Dynare's M_ structure
% - options_   [structure]   options, inspired by Dynare's options_ structure
% -------------------------------------------------------------------------
% OUTPUTS
% - estim_params_   [structure]       information on estimated parameters, inspired by Dynare's estimated_params block
% - xparams         [double vector]   numerical initial values passed to objective function
% - bounds          [double matrix]   numerical lower and upper bounds passed to objective function and optimizer
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
% shock indices
idx_eta_a = find(ismember(M_.exo_names,"eta_a"),1);
idx_eta_e = find(ismember(M_.exo_names,"eta_e"),1);
idx_eta_z = find(ismember(M_.exo_names,"eta_z"),1);
idx_eta_r = find(ismember(M_.exo_names,"eta_r"),1);

% initialize structure (inspired by Dynare's estim_params_)
% - column 1: index in M_.exo_names (to set M_.Sigma_eta and M_.Gamma_eta) for shock parameters; or in M_.endo_names (to set M_.Sigma_eps) for measurement erros; or in M_.param_names (to set M_.params)
% - column 2: initial value
% - column 3: lower bound
% - column 4: upper bound
estim_params_.skew_exo   = double.empty(0,4);
estim_params_.var_exo    = double.empty(0,4);
estim_params_.var_endo   = double.empty(0,4);
estim_params_.param_vals = double.empty(0,4);

if STAGE == "gaussian_initval_from_ireland2004_paper"    
    % values taken from table 3 in Ireland (2004)
    estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_a,                                   0.0302,  0,  1 ] ];
    estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_e,                                   0.0002,  0,  1 ] ];
    estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_z,                                   0.0089,  0,  1 ] ];
    estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_r,                                   0.0028,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"OMEGA"),1),    0.0581,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_X"),1),  0.00002, 1e-5,  1 ] ];    
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_PI"),1), 0.00002, 1e-5,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_PI"),1),   0.3866,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_G"),1),    0.3960,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_X"),1),    0.1654,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_A"),1),    0.9048,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_E"),1),    0.9907,  0,  1 ] ];

elseif STAGE == "csn_initval"
    if options_.parameters.use_stderr_skew
        % skew parameters
        estim_params_.skew_exo    = [ estim_params_.skew_exo;   [ idx_eta_a,  -0.217152129209896,  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo    = [ estim_params_.skew_exo;   [ idx_eta_e,  -0.226113231149225,  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo    = [ estim_params_.skew_exo;   [ idx_eta_z,  -0.923849236134372,  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo    = [ estim_params_.skew_exo;   [ idx_eta_r,   0.805559211425317,  options_.parameters.skewness_bounds ] ];
        % stderr parameters
        estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_a,   0.025239387157507,  0,  1 ] ];
        estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_e,   0.000197716076509666,  0,  1 ] ];
        estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_z,   0.00796603498559645,  0,  1 ] ];
        estim_params_.var_exo    = [ estim_params_.var_exo;     [ idx_eta_r,   0.00283620205731904,  0,  1 ] ];
    end
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"OMEGA"),1),    0.150114849014154,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_X"),1),  0.000193236868878676,  1e-5,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_PI"),1), 1.00000002310016e-05,  1e-5,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_PI"),1),   0.267849413749733,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_G"),1),    0.341003324945699,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_X"),1),    0.284612308301521,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_A"),1),    0.916787807033826,  0,  1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_E"),1),    0.981399127628676,  0,  1 ] ];

elseif STAGE == "csn_shock_params"
    % focus only on shock parameters, use calibrated values from M_ structure as initial values
    if options_.parameters.use_stderr_skew
        % skew parameters
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_a,  M_.Skew_eta(idx_eta_a,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_e,  M_.Skew_eta(idx_eta_e,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_z,  M_.Skew_eta(idx_eta_z,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_r,  M_.Skew_eta(idx_eta_r,1),  options_.parameters.skewness_bounds ] ];
        % stderr parameters
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_a,  sqrt(M_.Cov_eta(idx_eta_a,idx_eta_a)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_e,  sqrt(M_.Cov_eta(idx_eta_e,idx_eta_e)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_z,  sqrt(M_.Cov_eta(idx_eta_z,idx_eta_z)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_r,  sqrt(M_.Cov_eta(idx_eta_r,idx_eta_r)),  0, 1 ] ];
    else
        % diag_Gamma_eta parameters
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_a,  M_.Gamma_eta(idx_eta_a,idx_eta_a),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_e,  M_.Gamma_eta(idx_eta_e,idx_eta_e),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_z,  M_.Gamma_eta(idx_eta_z,idx_eta_z),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_r,  M_.Gamma_eta(idx_eta_r,idx_eta_r),  -Inf, Inf ] ];
        % sqrt_diag_Sigma_eta parameters
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_a,  sqrt(M_.Sigma_eta(idx_eta_a,idx_eta_a)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_e,  sqrt(M_.Sigma_eta(idx_eta_e,idx_eta_e)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_z,  sqrt(M_.Sigma_eta(idx_eta_z,idx_eta_z)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_r,  sqrt(M_.Sigma_eta(idx_eta_r,idx_eta_r)),  0, 1 ] ];
    end

elseif STAGE == "all_params"
    % use values provided in M_ structure as initial values
    if options_.parameters.use_stderr_skew
        % skew parameters
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_a,  M_.Skew_eta(idx_eta_a,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_e,  M_.Skew_eta(idx_eta_e,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_z,  M_.Skew_eta(idx_eta_z,1),  options_.parameters.skewness_bounds ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_r,  M_.Skew_eta(idx_eta_r,1),  options_.parameters.skewness_bounds ] ];
        % stderr parameters
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_a,  sqrt(M_.Cov_eta(idx_eta_a,idx_eta_a)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_e,  sqrt(M_.Cov_eta(idx_eta_e,idx_eta_e)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_z,  sqrt(M_.Cov_eta(idx_eta_z,idx_eta_z)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_r,  sqrt(M_.Cov_eta(idx_eta_r,idx_eta_r)),  0, 1 ] ];
    else
        % diag_Gamma_eta parameters
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_a,  M_.Gamma_eta(idx_eta_a,idx_eta_a),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_e,  M_.Gamma_eta(idx_eta_e,idx_eta_e),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_z,  M_.Gamma_eta(idx_eta_z,idx_eta_z),  -Inf, Inf ] ];
        estim_params_.skew_exo = [ estim_params_.skew_exo;  [ idx_eta_r,  M_.Gamma_eta(idx_eta_r,idx_eta_r),  -Inf, Inf ] ];
        % sqrt_diag_Sigma_eta parameters
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_a,  sqrt(M_.Sigma_eta(idx_eta_a,idx_eta_a)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_e,  sqrt(M_.Sigma_eta(idx_eta_e,idx_eta_e)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_z,  sqrt(M_.Sigma_eta(idx_eta_z,idx_eta_z)),  0, 1 ] ];
        estim_params_.var_exo  = [ estim_params_.var_exo;   [ idx_eta_r,  sqrt(M_.Sigma_eta(idx_eta_r,idx_eta_r)),  0, 1 ] ];        
    end
    % model parameters
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"OMEGA"),1),    M_.params(find(ismember(M_.param_names,"OMEGA"),1)),    0, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_X"),1),  M_.params(find(ismember(M_.param_names,"ALPHA_X"),1)),  1e-5, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"ALPHA_PI"),1), M_.params(find(ismember(M_.param_names,"ALPHA_PI"),1)), 1e-5, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_PI"),1),   M_.params(find(ismember(M_.param_names,"RHO_PI"),1)),   0, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_G"),1),    M_.params(find(ismember(M_.param_names,"RHO_G"),1)),    0, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_X"),1),    M_.params(find(ismember(M_.param_names,"RHO_X"),1)),    0, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_A"),1),    M_.params(find(ismember(M_.param_names,"RHO_A"),1)),    0, 1 ] ];
    estim_params_.param_vals = [ estim_params_.param_vals;  [ find(ismember(M_.param_names,"RHO_E"),1),    M_.params(find(ismember(M_.param_names,"RHO_E"),1)),    0, 1 ] ];
end

% remove fixed parameters
if isfield(options_.parameters.fix,'ALPHA_X') && (options_.parameters.fix.ALPHA_X==1)
   estim_params_.param_vals(find(estim_params_.param_vals(:,1)==find(M_.param_names == 'ALPHA_X',1),1),:) = [];
end
if isfield(options_.parameters.fix,'ALPHA_PI') && (options_.parameters.fix.ALPHA_PI==1)
   estim_params_.param_vals(find(estim_params_.param_vals(:,1)==find(M_.param_names == 'ALPHA_PI',1),1),:) = [];
end

% update numbers
estim_params_.nsx = size(estim_params_.skew_exo,1);
estim_params_.nvx = size(estim_params_.var_exo,1);
estim_params_.nvn = size(estim_params_.var_endo,1);
estim_params_.np  = size(estim_params_.param_vals,1);
estim_params_.ntot = estim_params_.nsx + estim_params_.nvx + estim_params_.nvn + estim_params_.np;
estim_params_.transformed = false(estim_params_.ntot,1);

% create names and optionally transform parameters
estim_params_.names = [];
for jnsx = 1:estim_params_.nsx
    if options_.parameters.use_stderr_skew
        tmpname = "skew_" + M_.exo_names(estim_params_.skew_exo(jnsx,1));
    else
        tmpname = "diag_Gamma_" + M_.exo_names(estim_params_.skew_exo(jnsx,1));        
    end
    if isfield(options_.parameters.transform,tmpname) && (options_.parameters.transform.(tmpname)==1)
        if estim_params_.skew_exo(jnsx,2) == estim_params_.skew_exo(jnsx,3)
            estim_params_.skew_exo(jnsx,2) = estim_params_.skew_exo(jnsx,3)+0.0001;
        end
        if estim_params_.skew_exo(jnsx,2) == estim_params_.skew_exo(jnsx,4)
            estim_params_.skew_exo(jnsx,2) = estim_params_.skew_exo(jnsx,4)-0.0001;
        end
        estim_params_.skew_exo(jnsx,2) = transform_to_unbounded(estim_params_.skew_exo(jnsx,2),estim_params_.skew_exo(jnsx,3),estim_params_.skew_exo(jnsx,4));
        estim_params_.names = [estim_params_.names; "transformed_"+tmpname];        
    else
        estim_params_.names = [estim_params_.names; tmpname];
    end
end
for jnvx = 1:estim_params_.nvx
    if options_.parameters.use_stderr_skew
        tmpname = "stderr_" + M_.exo_names(estim_params_.var_exo(jnvx,1));
    else
        tmpname = "sqrt_diag_Sigma_" + M_.exo_names(estim_params_.var_exo(jnvx,1));        
    end
    if isfield(options_.parameters.transform,tmpname) && (options_.parameters.transform.(tmpname)==1)
        if estim_params_.var_exo(jnvx,2) == estim_params_.var_exo(jnvx,3)
            estim_params_.var_exo(jnvx,2) = estim_params_.var_exo(jnvx,3)+0.0001;
        end
        if estim_params_.var_exo(jnvx,2) == estim_params_.var_exo(jnvx,4)
            estim_params_.var_exo(jnvx,2) = estim_params_.var_exo(jnvx,4)-0.0001;
        end
        estim_params_.var_exo(jnvx,2) = transform_to_unbounded(estim_params_.var_exo(jnvx,2),estim_params_.var_exo(jnvx,3),estim_params_.var_exo(jnvx,4));
        estim_params_.names = [estim_params_.names; "transformed_"+tmpname];        
    else
        estim_params_.names = [estim_params_.names; tmpname];
    end
end
for jnvn = 1:estim_params_.nvn
    tmpname = "stderr_" + M_.varobs(estim_params_.var_endo(jnvn,1));
    if isfield(options_.parameters.transform,tmpname) && (options_.parameters.transform.(tmpname)==1)
        if estim_params_.var_endo(jnvn,2) == estim_params_.var_endo(jnvn,3)
            estim_params_.var_endo(jnvn,2) = estim_params_.var_endo(jnvn,3)+0.0001;
        end
        if estim_params_.var_endo(jnvn,2) == estim_params_.var_endo(jnvn,4)
            estim_params_.var_endo(jnvn,2) = estim_params_.var_endo(jnvn,4)-0.0001;
        end
        estim_params_.var_endo(jnvn,2) = transform_to_unbounded(estim_params_.var_endo(jnvn,2),estim_params_.var_endo(jnvn,3),estim_params_.var_endo(jnvn,4));
        estim_params_.names = [estim_params_.names; "transformed_"+tmpname];        
    else
        estim_params_.names = [estim_params_.names; tmpname];
    end
end
for jnp = 1:estim_params_.np
    tmpname = M_.param_names(estim_params_.param_vals(jnp,1));
    if isfield(options_.parameters.transform,tmpname) && (options_.parameters.transform.(tmpname)==1)
        if estim_params_.param_vals(jnp,2) == estim_params_.param_vals(jnp,3)
            estim_params_.param_vals(jnp,2) = estim_params_.param_vals(jnp,3)+0.0001;
        end
        if estim_params_.param_vals(jnp,2) == estim_params_.param_vals(jnp,4)
            estim_params_.param_vals(jnp,2) = estim_params_.param_vals(jnp,4)-0.0001;
        end        
        estim_params_.param_vals(jnp,2) = transform_to_unbounded(estim_params_.param_vals(jnp,2),estim_params_.param_vals(jnp,3),estim_params_.param_vals(jnp,4));
        estim_params_.names = [estim_params_.names; "transformed_"+tmpname];
    else
        estim_params_.names = [estim_params_.names; tmpname];
    end
end
estim_params_.transformed(contains(estim_params_.names,"transformed_")) = true;
% create numerical vector for initial values
xparams = [estim_params_.skew_exo(  :,2);
           estim_params_.var_exo(   :,2);
           estim_params_.var_endo(  :,2);
           estim_params_.param_vals(:,2)];
% create matrix for bounds
bounds = [estim_params_.skew_exo(  :,3)  estim_params_.skew_exo(  :,4);
          estim_params_.var_exo(   :,3)  estim_params_.var_exo(   :,4);
          estim_params_.var_endo(  :,3)  estim_params_.var_endo(  :,4);
          estim_params_.param_vals(:,3)  estim_params_.param_vals(:,4)];
bounds(estim_params_.transformed,1) = -Inf;
bounds(estim_params_.transformed,2) = Inf;
