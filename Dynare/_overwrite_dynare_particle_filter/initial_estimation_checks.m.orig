function oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,dataset_info,M_,estim_params_,options_,bayestopt_,BoundsInfo,oo_)
% function oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,dataset_info,M_,estim_params_,options_,bayestopt_,BoundsInfo,oo_)
% Checks data (complex values, ML evaluation, initial values, BK conditions,..)
%
% INPUTS
%   objective_function  [function handle] of the objective function
%   xparam1             [vector] of parameters to be estimated
%   dataset_            [dseries] object storing the dataset
%   dataset_info        [structure] storing information about the sample.
%   M_                  [structure] describing the model
%   estim_params_       [structure] characterizing parameters to be estimated
%   options_            [structure] describing the options
%   bayestopt_          [structure] describing the priors
%   BoundsInfo          [structure] containing prior bounds
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   oo_                 [structure] storing the results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2025 Dynare Team
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

%get maximum number of simultaneously observed variables for stochastic
%singularity check
maximum_number_non_missing_observations=max(sum(~isnan(dataset_.data(2:end,:)),2));
init_number_non_missing_observations=sum(~isnan(dataset_.data(1,:)),2);

if options_.heteroskedastic_filter
    if options_.order>1 
        error('initial_estimation_checks:: heteroskedastic shocks are only supported with the Kalman filter/smoother')
    end
    observations_by_period=sum(~isnan(dataset_.data),2);
    base_shocks = find(diag(M_.Sigma_e));
    zero_shocks = ~ismember(1:M_.exo_nbr,base_shocks);
    non_zero_shocks_by_period=repmat(length(base_shocks),size(observations_by_period));
    % check periods for which base shocks are scaled to zero
    non_zero_shocks_by_period = non_zero_shocks_by_period-sum(M_.heteroskedastic_shocks.Qscale(base_shocks,1:options_.nobs)==0,1)';
    % check periods for which base shocks are set to zero
    non_zero_shocks_by_period = non_zero_shocks_by_period-sum(M_.heteroskedastic_shocks.Qvalue(base_shocks,1:options_.nobs)==0,1)';
    % check periods for which other shocks are set to a positive number
    non_zero_shocks_by_period = non_zero_shocks_by_period+sum(M_.heteroskedastic_shocks.Qvalue(zero_shocks,1:options_.nobs)>0,1)';
end

if options_.order>1
    if any(any(isnan(dataset_.data)))
        error('initial_estimation_checks:: particle filtering does not support missing observations')
    end
    if options_.prefilter==1
        error('initial_estimation_checks:: particle filtering does not support the prefilter option')
    end
    if bayestopt_.with_trend
        error('initial_estimation_checks:: particle filtering does not support trends')
    end
    if options_.order>3 && options_.particle.pruning==1
        error('initial_estimation_checks:: the particle filter with order>3 does not support pruning')
    end
    if options_.particle.pruning~=options_.pruning
        warning('initial_estimation_checks:: the pruning settings differ between the particle filter and the one used for IRFs/simulations. Make sure this is intended.')
    end
end

if options_.occbin.likelihood.status || options_.occbin.smoother.status
    if options_.prefilter
        error('initial_estimation_checks:: OccBin is incompatible with the prefilter option due to the sample mean generally not corresponding to the steady state with an occasionally binding constraint.')
    end
    if ~options_.occbin.likelihood.inversion_filter && (options_.kalman_algo==2 || options_.kalman_algo==4)
        error('initial_estimation_checks:: OccBin is incompatible with the selected univariate Kalman filter.')        
    end
    if options_.fast_kalman_filter
        error('initial_estimation_checks:: OccBin is incompatible with the fast Kalman filter.')        
    end
    if options_.bayesian_irf
        error('initial_estimation_checks:: OccBin is incompatible with the bayesian_irf option.')                
    end
    if options_.moments_varendo
        error('initial_estimation_checks:: OccBin is incompatible with the moments_varendo option.')                
    end
end

if (options_.occbin.likelihood.status && options_.occbin.likelihood.inversion_filter) || (options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter)
    err_index= find(diag(M_.Sigma_e)~=0);
    if length(err_index)~=length(options_.varobs)
        fprintf('initial_estimation_checks:: The IVF requires exactly as many shocks as observables.')
    end
    var_index=find(any(isnan(dataset_.data)));
    if ~isempty(var_index)
        fprintf('initial_estimation_checks:: The IVF requires exactly as many shocks as observables.\n')
        fprintf('initial_estimation_checks:: The data series %s contains NaN, I am therefore dropping shock %s for these time points.\n',...
            options_.varobs{var_index},M_.exo_names{options_.occbin.likelihood.IVF_shock_observable_mapping(var_index)})
    end
end

if options_.order>1 || (options_.order==1 && ~ischar(options_.mode_compute) && options_.mode_compute==11)
    if options_.order==1 && options_.mode_compute==11
        disp_string='mode_compute=11';
    else
        disp_string='particle filtering';
    end
    if M_.H==0
        error('initial_estimation_checks:: %s requires measurement error on the observables',disp_string)
    else
        if sum(diag(M_.H)>0)<length(options_.varobs)
            error('initial_estimation_checks:: %s requires as many measurement errors as observed variables',disp_string)
        else
            [~,flag]=chol(M_.H);
            if flag
                error('initial_estimation_checks:: the measurement error matrix must be positive definite')
            end
        end
    end
end

non_zero_ME=length(estim_params_.H_entries_to_check_for_positive_definiteness);

print_init_check_warning=false;
if maximum_number_non_missing_observations>M_.exo_nbr+non_zero_ME
    error('initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables!')
end
if init_number_non_missing_observations>M_.exo_nbr+non_zero_ME
    if options_.no_init_estimation_check_first_obs
        print_init_check_warning=true;
    else
        error('initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables in first period!')
    end
end

if options_.heteroskedastic_filter
    if any(observations_by_period>(non_zero_shocks_by_period+non_zero_ME))
        error('initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance: Check heteroskedastic block and shocks calibration!')
    end
else
    if maximum_number_non_missing_observations>length(find(diag(M_.Sigma_e)))+non_zero_ME
        error('initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance!')
    end
end
if init_number_non_missing_observations>length(find(diag(M_.Sigma_e)))+non_zero_ME
    if options_.no_init_estimation_check_first_obs
        print_init_check_warning=true;
    else
        error('initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance in first period!')
    end
end
if print_init_check_warning
    fprintf('ESTIMATION_CHECKS: You decided to ignore test of stochastic singularity in first_obs.\n');
    fprintf('ESTIMATION_CHECKS: If this was not done on purpose (typically when observing a stock variable [capital] in first period, on top of its flow [investment]),\n');
    fprintf('ESTIMATION_CHECKS: it may lead to a crash or provide undesired/wrong results later on!\n');
end

if (any(bayestopt_.pshape  >0 ) && options_.mh_replic) && options_.mh_nblck<1
    error('initial_estimation_checks:: Bayesian estimation cannot be conducted with mh_nblocks=0.')
end

if options_.mh_drop<0 || options_.mh_drop>=1
    error('initial_estimation_checks:: mh_drop must be in [0,1).')    
end

% check and display warnings if steady-state solves static model (except if diffuse_filter == 1) and if steady-state changes estimated parameters
[oo_.steady_state] = check_steady_state_changes_parameters(M_,estim_params_,oo_,options_, [options_.diffuse_filter==0 options_.diffuse_filter==0] );

% check and display warning if prior implies that stderr, corr, or skew params are outside theoretical bounds
if any(bayestopt_.pshape)
    check_prior_stderr_corr_skew(estim_params_,bayestopt_);
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(M_);

[~,~,~,info]= priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
if any(info)
    fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',bayestopt_.name{info,1})
    error('The initial value of the prior is -Inf')
end

if isfield(M_,'filter_initial_state') && ~isempty(M_.filter_initial_state)
    state_indices=oo_.dr.order_var(oo_.dr.restrict_var_list(bayestopt_.mf0));
    for ii=1:size(state_indices,1)
        if ~isempty(M_.filter_initial_state{state_indices(ii),1})
            try
                evaluate_expression(M_.filter_initial_state{state_indices(ii),2},M_,oo_)
            catch
                fprintf('Unable to evaluate the expression\n %s \nfor the filter_initial_state of variable %s\n',M_.filter_initial_state{state_indices(ii),2},M_.endo_names(state_indices(ii),:))
            end
        end
    end
end

if ~isreal(dataset_.data)
    error('initial_estimation_checks: the data contains complex values.')
end

% Evaluate the likelihood.
ana_deriv = options_.analytic_derivation;
options_.analytic_derivation=0;
if ~isequal(options_.mode_compute,11) || ...
        (isequal(options_.mode_compute,11) && isequal(options_.order,1))
    %shut off potentially automatic switch to diffuse filter for the
    %purpose of checking stochastic singularity
    use_univariate_filters_if_singularity_is_detected_old=options_.use_univariate_filters_if_singularity_is_detected;
    options_.use_univariate_filters_if_singularity_is_detected=0;
    [fval,info] = feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    if info(1)==50
        fprintf('\ninitial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.\n')
        fprintf('initial_estimation_checks:: This is often a sign of stochastic singularity, but can also sometimes happen by chance\n')
        fprintf('initial_estimation_checks:: for a particular combination of parameters and data realizations.\n')
        fprintf('initial_estimation_checks:: If you think the latter is the case, you should try with different initial values for the estimated parameters.\n')
        error('initial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.')
    end
    if info(1)==201 || info(1)==202
        message=get_error_message(info,options_);
        error('initial_estimation_checks:: %s.',message)
    end
    %reset options
    options_.use_univariate_filters_if_singularity_is_detected=use_univariate_filters_if_singularity_is_detected_old;
else
    info=0;
    fval = 0;
end
if options_.debug
    oo_.likelihood_at_initial_parameters=fval;
end
options_.analytic_derivation=ana_deriv;

if options_.mode_compute==13
    error('Options mode_compute=13 is only compatible with quadratic objective functions')
end

if isnan(fval)
    error('The initial value of the likelihood is NaN')
elseif imag(fval)
    error('The initial value of the likelihood is complex')
end

if info(1) > 0
    if options_.order>1
        [eigenvalues_] = check(M_,options_, oo_);
        if any(abs(1-abs(eigenvalues_))<abs(options_.qz_criterium-1))
            error('Your model has at least one unit root and you are using a nonlinear filter. Please set nonlinear_filter_initialization=3.')
        end
    else
        disp('Error in computing likelihood for initial parameter values')
        print_info(info, options_.noprint, options_)
    end
end

if options_.prefilter==1
    if (~options_.loglinear && any(abs(oo_.steady_state(bayestopt_.mfys))>1e-9)) || (options_.loglinear && any(abs(log(oo_.steady_state(bayestopt_.mfys)))>1e-9))
        disp('You are trying to estimate a model with a non zero steady state for the observed endogenous')
        disp('variables using demeaned data!')
        error('You should change something in your mod file...')
    end
end

if ~isequal(options_.mode_compute,11)
    disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
end

if options_.mh_tune_jscale.status && (options_.mh_tune_jscale.maxiter<options_.mh_tune_jscale.stepsize)
    warning('You specified mh_tune_jscale, but the maximum number of iterations is smaller than the step size. No update will take place.')
end

if ~isempty(options_.conditional_variance_decomposition) && ~options_.moments_varendo
    disp('The conditional_variance_decomposition-option will be ignored. You need to set moments_varendo');
end

function evaluate_expression(expression,M_,oo_)%#ok<INUSL>
% function evaluate_expression(expression,M_,oo_)
%evaluates expressions relying on M_ and oo_ having their original names
eval(expression);
