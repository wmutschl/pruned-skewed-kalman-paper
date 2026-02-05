function [fval,info,exit_flag,DLIK,Hess,ys,trend_coeff,M_,bayestopt_,dr] = non_linear_dsge_likelihood(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [fval,info,exit_flag,DLIK,Hess,ys,trend_coeff,M_,options_,bayestopt_,dr] = non_linear_dsge_likelihood(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Evaluates the posterior kernel of a dsge model using a non linear filter.
%
% INPUTS
% - xparam1                 [double]              n×1 vector, estimated parameters.
% - dataset_                [struct]              MATLAB's structure containing the dataset
% - dataset_info            [struct]              MATLAB's structure describing the dataset
% - options_                [struct]              MATLAB's structure describing the options
% - M_                      [struct]              MATLAB's structure describing the M_
% - estim_params_           [struct]              MATLAB's structure describing the estimated_parameters
% - bayestopt_              [struct]              MATLAB's structure describing the priors
% - BoundsInfo              [struct]              MATLAB's structure specifying the bounds on the parameter values
% - dr                      [structure]           Reduced form model.
% - endo_steady_state       [vector]              steady state value for endogenous variables
% - exo_steady_state        [vector]              steady state value for exogenous variables
% - exo_det_steady_state    [vector]              steady state value for exogenous deterministic variables
%
% OUTPUTS
% - fval                    [double]              scalar, value of the likelihood or posterior kernel.
% - info                    [integer]             4×1 vector, information on whether solution and likelihood could be computed
% - exit_flag               [integer]             scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - DLIK                    [double]              Empty array.
% - Hess                    [double]              Empty array.
% - ys                      [double]              Empty array.
% - trend_coeff             [double]              Empty array.
% - M_                      [struct]              Updated M_ structure described in INPUTS section.
% - options_                [struct]              Updated options_ structure described in INPUTS section.
% - bayestopt_              [struct]              See INPUTS section.
% - dr                      [struct]              decision rule structure described in INPUTS section.

% Copyright © 2010-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% Initialization of the returned arguments.
fval            = [];
ys              = [];
trend_coeff     = [];
exit_flag       = 1;
DLIK            = [];
Hess            = [];

% Ensure that xparam1 is a column vector.
% (Don't do the transformation if xparam1 is empty, otherwise it would become a
%  0×1 matrix, which create issues with older MATLABs when comparing with [] in
%  check_bounds_and_definiteness_estimation)
if ~isempty(xparam1)
    xparam1 = xparam1(:);
end

% Issue an error if loglinear option is used.
if options_.loglinear
    error('non_linear_dsge_likelihood: It is not possible to use a non linear filter with the option loglinear!')
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

M_ = set_all_parameters(xparam1,estim_params_,M_);

[fval,info,exit_flag,Q,H]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, BoundsInfo);
if info(1)
    return
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic steadystate and extract the matrices of the state equation (T and R).
[dr, info, M_.params] = resol(0, M_, options_, dr , endo_steady_state, exo_steady_state, exo_det_steady_state);

if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 || ...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85
        %meaningful second entry of output that can be used
        fval = Inf;
        if ~isfinite(info(2))
            info(4) = 0.1;
        else
            info(4) = info(2);
        end
        exit_flag = 0;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end
end

% Define a vector of indices for the observed variables. Is this really usefull?...
bayestopt_.mf = bayestopt_.mf1;

% Get needed information for Kalman filter routines.
start = options_.presample+1;
Y = transpose(dataset_.data);

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------

mf0 = bayestopt_.mf0;
restrict_variables_idx = dr.restrict_var_list;
state_variables_idx = restrict_variables_idx(bayestopt_.mf0);
number_of_state_variables = length(mf0);

ReducedForm.steadystate = dr.ys(dr.order_var(restrict_variables_idx));
if options_.order == 1
    ReducedForm.constant = ReducedForm.steadystate;
else
    ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
end
ReducedForm.state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));
ReducedForm.Q = Q;
ReducedForm.H = H;
ReducedForm.mf0 = mf0;
ReducedForm.mf1 = bayestopt_.mf1;

if options_.order>3
    ReducedForm.use_k_order_solver = true;
    ReducedForm.dr = dr;
    ReducedForm.udr = folded_to_unfolded_dr(dr, M_, options_);
    if pruning
        error('Pruning is not available for orders > 3');
    end
else
    ReducedForm.use_k_order_solver = false;
    ReducedForm.ghx  = dr.ghx(restrict_variables_idx,:);
    ReducedForm.ghu  = dr.ghu(restrict_variables_idx,:);
    if options_.order > 1
        ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
        ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
        ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
        ReducedForm.ghs2 = dr.ghs2(restrict_variables_idx,:);
    end
    if options_.order==3
        ReducedForm.ghxxx = dr.ghxxx(restrict_variables_idx,:);
        ReducedForm.ghuuu = dr.ghuuu(restrict_variables_idx,:);
        ReducedForm.ghxxu = dr.ghxxu(restrict_variables_idx,:);
        ReducedForm.ghxuu = dr.ghxuu(restrict_variables_idx,:);
        ReducedForm.ghxss = dr.ghxss(restrict_variables_idx,:);
        ReducedForm.ghuss = dr.ghuss(restrict_variables_idx,:);
    end
end

% Set initial condition.
switch options_.particle.initialization
  case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
    StateVectorMean = ReducedForm.constant(mf0);
    [A,B] = kalman_transition_matrix(dr,dr.restrict_var_list,dr.restrict_columns);
    StateVectorVariance = lyapunov_symm(A, B*Q*B', options_.lyapunov_fixed_point_tol, ...
                                        options_.qz_criterium, options_.lyapunov_complex_threshold, [], options_.debug);
    StateVectorVariance = StateVectorVariance(mf0,mf0);
  case 2% Initial state vector covariance is a monte-carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
    StateVectorMean = ReducedForm.constant(mf0);
    old_DynareOptionsperiods = options_.periods;
    options_.periods = 5000;
    old_DynareOptionspruning =  options_.pruning;
    options_.pruning = options_.particle.pruning;
    y_ = simult(endo_steady_state, dr,M_,options_);
    y_ = y_(dr.order_var(state_variables_idx),2001:5000); %state_variables_idx is in dr-order while simult_ is in declaration order
    if any(any(isnan(y_))) ||  any(any(isinf(y_))) && ~ options_.pruning
        fval = Inf;
        info(1) = 202;
        info(4) = 0.1;
        exit_flag = 0;
        return;        
    end
    StateVectorVariance = cov(y_');       
    options_.periods = old_DynareOptionsperiods;
    options_.pruning = old_DynareOptionspruning;
    clear('old_DynareOptionsperiods','y_');
  case 3% Initial state vector covariance is a diagonal matrix (to be used
        % if model has stochastic trends).
    StateVectorMean = ReducedForm.constant(mf0);
    StateVectorVariance = options_.particle.initial_state_prior_std*eye(number_of_state_variables);
  otherwise
    error('Unknown initialization option!')
end
ReducedForm.StateVectorMean = StateVectorMean;
ReducedForm.StateVectorVariance = StateVectorVariance;

if ~options_.particle.use_reduced_rank_cholesky
    [~, flag] = chol(ReducedForm.StateVectorVariance);%reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
    if flag 
        fval = Inf;
        info(1) = 201;
        info(4) = 0.1;
        exit_flag = 0;    
        return;
    end
end

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------
options_.warning_for_steadystate = 0;
[s1,s2,current_stream] = get_dynare_random_generator_state();

startTime = tic;
LIK = feval(options_.particle.algorithm, ReducedForm, Y, start, options_.particle, options_.threads, options_, M_);
fprintf('Time to evaluate likelihood with particle filter with N=%d: %s\n', options_.particle.number_of_particles, dynsec2hms(toc(startTime)));
set_dynare_random_generator_state(s1,s2,current_stream);
if imag(LIK)
    fval = Inf; info(1) = 46; info(4) = 0.1; exit_flag = 0;
    return
elseif isnan(LIK)
    fval = Inf; info(1) = 45; info(4) = 0.1; exit_flag = 0;
    return
elseif isinf(LIK)
    fval = Inf; info(1) = 50; info(4) = 0.1; exit_flag = 0;
    return
else
    likelihood = LIK;
end
options_.warning_for_steadystate = 1;
% ------------------------------------------------------------------------------
% Adds prior if necessary
% ------------------------------------------------------------------------------
lnprior = priordens(xparam1(:),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
if isinf(lnprior)
    fval = Inf; info(1) = 40; info(4) = 0.1; exit_flag = 0;
    return
end
if isnan(lnprior)
    fval = Inf; info(1) = 47; info(4) = 0.1; exit_flag = 0;
    return
end
if imag(lnprior)~=0
    fval = Inf; info(1) = 48; info(4) = 0.1; exit_flag = 0;
    return
end
fval = (likelihood-lnprior);