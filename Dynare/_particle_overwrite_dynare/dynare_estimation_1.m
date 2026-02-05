function dynare_estimation_1(var_list_,dname)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info

dispString = 'Estimation::mcmc';

if issmc(options_)
    options_.mode_compute = 0;
    options_.mh_replic = 0;
    options_.mh_recover = false;
    options_.load_mh_file = false;
    options_.load_results_after_load_mh = false;
end
if isdime(options_) && options_.prior_trunc
    options_.prior_trunc = 0;
    fprintf('%s: DIME requires no prior truncation. Resetting options_.prior_trunc=0.\n', dispString);
end

if ~isfolder([M_.dname filesep 'Output'])
    mkdir(M_.dname,'Output');
end

if isempty(estim_params_)
    mode_compute_o = options_.mode_compute;
    mh_replic_o = options_.mh_replic;
    options_.mode_compute = 0;
    options_.mh_replic = 0;
    reset_options_related_to_estimation = true;
else
    reset_options_related_to_estimation = false;
end

%store qz_criterium
qz_criterium_old=options_.qz_criterium;
if isnan(options_.first_obs)
    first_obs_nan_indicator=true;
else
    first_obs_nan_indicator=false;
end

% Set particle filter flag.
if options_.order > 1 && ~options_.particle.status
    error('For estimating the model with a higher-order approximation using a nonlinear filter, one should have options_.particle.status=true;')
end
if options_.particle.status
    skipline()
        disp('Estimation using a non linear filter!')
    skipline()
    if ~options_.nointeractive && ismember(options_.mode_compute,[1,3,4]) && ~strcmpi(options_.particle.filter_algorithm,'gf')% Known gradient-based optimizers
        disp('You are using a gradient-based mode-finder. Particle filtering introduces discontinuities in the')
        disp('objective function w.r.t the parameters. Thus, should use a non-gradient based optimizer.')
        fprintf('\nPlease choose a mode-finder:\n')
        fprintf('\t 0 - Continue using gradient-based method (it is most likely that you will no get any sensible result).\n')
        fprintf('\t 6 - Monte Carlo based algorithm\n')
        fprintf('\t 7 - Nelder-Mead simplex based optimization routine (MATLAB optimization toolbox required)\n')
        fprintf('\t 8 - Nelder-Mead simplex based optimization routine (Dynare''s implementation)\n')
        fprintf('\t 9 - CMA-ES (Covariance Matrix Adaptation Evolution Strategy) algorithm\n')
        choice = [];
        while isempty(choice)
            choice = input('Please enter your choice: ');
            if isnumeric(choice) && isint(choice) && ismember(choice,[0 6 7 8 9])
                if choice
                    options_.mode_compute = choice;
                end
            else
                fprintf('\nThis is an invalid choice (you have to choose between 0, 6, 7, 8 and 9).\n')
                choice = [];
            end
        end
    end
end

%
% set objective function
%
if ~options_.dsge_var
    if options_.particle.status
        objective_function = str2func('non_linear_dsge_likelihood');
        [options_.particle] = check_particle_filter_options(options_.particle);
    else
        if options_.occbin.likelihood.status && options_.occbin.likelihood.inversion_filter
            objective_function = str2func('occbin.IVF_posterior');
        elseif options_.conditional_likelihood.status && options_.conditional_likelihood.order==1
            objective_function = str2func('dsge_conditional_likelihood_1');
        else
            objective_function = str2func('dsge_likelihood');
        end
    end
else
    objective_function = str2func('dsge_var_likelihood');
end

[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = ...
    dynare_estimation_init(var_list_, dname, [], M_, options_, oo_, estim_params_, bayestopt_);

if options_.dsge_var
    check_dsge_var_model(M_, estim_params_, bayestopt_);
    if dataset_info.missing.state
        error('Estimation::DsgeVarLikelihood: I cannot estimate a DSGE-VAR model with missing observations!')
    end
    if options_.noconstant
        dataset_info=var_sample_moments(options_.dsge_varlag, -1, dataset_, dataset_info);
    else
        % The steady state is non zero ==> a constant in the VAR is needed!
        dataset_info=var_sample_moments(options_.dsge_varlag, 0, dataset_, dataset_info);
    end
end

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
M_.sigma_e_is_diagonal = true;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = false;
end

data = dataset_.data;
data_index = dataset_info.missing.aindex;
missing_value = dataset_info.missing.state;

% Set number of observations
gend = dataset_.nobs;

% Get the number of parameters to be estimated.
nx = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+estim_params_.np; % Total number of parameters to be estimated.

if ~isempty(estim_params_)
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end

%
% perform initial estimation checks;
%

try
    oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,dataset_info,M_,estim_params_,options_,bayestopt_,bounds,oo_);
catch e % if check fails, provide info on using calibration if present
    if estim_params_.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('ESTIMATION_CHECKS: There was an error in computing the likelihood for initial parameter values.\n')
        fprintf('ESTIMATION_CHECKS: If this is not a problem with the setting of options (check the error message below),\n')
        fprintf('ESTIMATION_CHECKS: you should try using the calibrated version of the model as starting values. To do\n')
        fprintf('ESTIMATION_CHECKS: this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')
        fprintf('ESTIMATION_CHECKS: command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(e);
end

%
% Run smoother if no estimation or mode-finding are requested
%

if isequal(options_.mode_compute,0) && isempty(options_.mode_file) && ~options_.mh_posterior_mode_estimation && ~issmc(options_)
    if options_.order==1 && ~options_.particle.status
        if options_.smoother
            if options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter
                [~, info, ~, ~, ~, ~, ~, ~, oo_.dr, atT, innov, oo_.occbin.smoother.regime_history] = occbin.IVF_posterior(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,prior_bounds(bayestopt_,options_.prior_trunc),oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                if ismember(info(1),[303,304,306])
                    fprintf('\nIVF: smoother did not succeed. No results will be written to oo_.\n')
                else
                    updated_variables = atT*nan;
                    measurement_error=[];
                    ys = oo_.dr.ys;
                    trend_coeff = zeros(length(options_.varobs_id),1);
                    bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
                    options_nk=options_.nk;
                    options_.nk=[]; %unset options_.nk and reset it later
                    [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff);
                    options_.nk=options_nk;
                end
            else
                if options_.occbin.smoother.status
                    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,~,~,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_.mf,a0T,state_uncertainty0] = occbin.DSGE_smoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info);
                    if oo_.occbin.smoother.error_flag(1)==0
                        [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty,a0T,state_uncertainty0);
                    else
                        fprintf('\nOccBin: smoother did not succeed. No results will be written to oo_.\n')
                    end
                else
                    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,~,~,P,PK,decomp,Trend,state_uncertainty,oo_.dr,bayestopt_.mf,a0T,state_uncertainty0] = DsgeSmoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_);
                    [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty,a0T,state_uncertainty0);
                end
            end
            if options_.forecast > 0
                oo_.forecast = forecasts.run(var_list_,M_,options_,oo_,'smoother',dataset_info);
            end
        end
        %reset qz_criterium
        options_.qz_criterium=qz_criterium_old;
        return
    else %allow to continue, e.g. with MCMC_jumping_covariance
        if options_.smoother
            error('Estimation:: Particle Smoothers are not yet implemented.')
        end
    end
end

%
% Estimation of the posterior mode or likelihood mode
%

if ~isequal(options_.mode_compute,0) && ~options_.mh_posterior_mode_estimation && ~issmc(options_)
    optimizer_vec = [options_.mode_compute;num2cell(options_.additional_optimizer_steps)];
    for optim_iter = 1:length(optimizer_vec)
        current_optimizer = optimizer_vec{optim_iter};
        if isnumeric(current_optimizer)
            if current_optimizer==5
                if options_.analytic_derivation
                    old_analytic_derivation = options_.analytic_derivation;
                    options_.analytic_derivation=-1; %force analytic outer product gradient Hessian for each iteration
                end
            end
        end
        [xparam1, fval, ~, hh, Scale, new_rat_hess_info, optimization_info] = ...
            dynare_minimize_objective(objective_function,xparam1,current_optimizer,options_,[bounds.lb bounds.ub],bayestopt_.name,bayestopt_,hh,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        fprintf('\nFinal value of minus the log posterior (or likelihood):%f \n', fval);
        if length(optimizer_vec) > 1
            oo_.posterior.optimization.optimization_info.(sprintf('stage_%d', optim_iter)) = optimization_info;
        else
            oo_.posterior.optimization.optimization_info = optimization_info;
        end
        if isnumeric(current_optimizer)
            if current_optimizer==5
                newratflag = new_rat_hess_info.newratflag;
                new_rat_hess_info = new_rat_hess_info.new_rat_hess_info;
                if options_.analytic_derivation
                    options_.analytic_derivation = old_analytic_derivation;
                end
            elseif current_optimizer==6 %save scaling factor
                save([M_.dname filesep 'Output' filesep M_.fname '_optimal_mh_scale_parameter.mat'],'Scale');
                options_.mh_jscale = Scale;
                bayestopt_.jscale(:) = options_.mh_jscale;
            end
        end
        if ~isnumeric(current_optimizer) || ~isequal(current_optimizer,6) %always already computes covariance matrix
            if options_.cova_compute == 1 %user did not request covariance not to be computed
                if options_.analytic_derivation && strcmp(func2str(objective_function),'dsge_likelihood')
                    ana_deriv_old = options_.analytic_derivation;
                    options_.analytic_derivation = 2;
                    [~,~,~,~,hh] = feval(objective_function,xparam1, ...
                                         dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                    options_.analytic_derivation = ana_deriv_old;
                elseif ~isnumeric(current_optimizer) || ~(isequal(current_optimizer,5) && newratflag~=1 && strcmp(func2str(objective_function),'dsge_likelihood'))
                    % enter here if i) not mode_compute_5, ii) if mode_compute_5 and newratflag==1;
                    % with flag==0 or 2 and dsge_likelihood, we force to use
                    % the Hessian from outer product gradient of optimizer 5 below
                    if options_.hessian.use_penalized_objective
                        penalized_objective_function = str2func('penalty_objective_function');
                        hh = hessian(penalized_objective_function, xparam1, options_.gstep, objective_function, fval, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                    else
                        hh = hessian(objective_function, xparam1, options_.gstep, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                    end 
                    hh = reshape(hh, nx, nx);
                elseif isnumeric(current_optimizer) && isequal(current_optimizer,5)
                    % other numerical Hessian options available with optimizer
                    % 5 and dsge_likelihood
                    %
                    % if newratflag == 0
                    % compute outer product gradient of optimizer 5
                    %
                    % if newratflag == 2
                    % compute 'mixed' outer product gradient of optimizer 5
                    % with diagonal elements computed with numerical second order derivatives
                    %
                    % uses univariate filters, so to get max # of available
                    % densities for outer product gradient
                    kalman_algo0 = options_.kalman_algo;
                    compute_hessian = 1;
                    if ~((options_.kalman_algo == 2) || (options_.kalman_algo == 4))
                        options_.kalman_algo=2;
                        if options_.lik_init == 3
                            options_.kalman_algo=4;
                        end
                    elseif newratflag==0 % hh already contains outer product gradient with univariate filter
                        compute_hessian = 0;
                    end
                    if compute_hessian
                        crit = options_.newrat.tolerance.f;
                        newratflag = newratflag>0;
                        hh = reshape(mr_hessian(xparam1,objective_function,fval,newratflag,crit,new_rat_hess_info,[bounds.lb bounds.ub],bayestopt_.p2,0,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state), nx, nx);
                    end
                    options_.kalman_algo = kalman_algo0;
                end
            end
        end
        parameter_names = bayestopt_.name;
    end
    if options_.cova_compute || current_optimizer==5 || current_optimizer==6
        save([M_.dname filesep 'Output' filesep M_.fname '_mode.mat'],'xparam1','hh','parameter_names','fval');
    else
        save([M_.dname filesep 'Output' filesep M_.fname '_mode.mat'],'xparam1','parameter_names','fval');
    end
end

if ~options_.mh_posterior_mode_estimation && options_.cova_compute && ~issmc(options_)
    check_hessian_at_the_mode(hh, xparam1, M_, estim_params_, options_, bounds);
end

%
% create mode_check_plots
%

if options_.mode_check.status && ~options_.mh_posterior_mode_estimation && ~issmc(options_)
    ana_deriv_old = options_.analytic_derivation;
    options_.analytic_derivation = 0;
    mode_check(objective_function,xparam1,hh,options_,M_,estim_params_,bayestopt_,bounds,false,...
               dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    options_.analytic_derivation = ana_deriv_old;
end

oo_.posterior.optimization.mode = [];
oo_.posterior.optimization.Variance = [];
oo_.posterior.optimization.log_density=[];

invhess = [];

if ~issmc(options_)
    if ~options_.mh_posterior_mode_estimation
        oo_.posterior.optimization.mode = xparam1;
        if exist('fval','var')
            oo_.posterior.optimization.log_density=-fval;
        end
        if options_.cova_compute
            hsd = sqrt(diag(hh));
            invhess = inv(hh./(hsd*hsd'))./(hsd*hsd');
            stdh = sqrt(diag(invhess));
            oo_.posterior.optimization.Variance = invhess;
        end
    else
        variances = bayestopt_.p2.*bayestopt_.p2;
        idInf = isinf(variances);
        variances(idInf) = 1;
        invhess = options_.mh_posterior_mode_estimation*diag(variances);
        xparam1 = bayestopt_.p5;
        idNaN = isnan(xparam1);
        xparam1(idNaN) = bayestopt_.p1(idNaN);
        outside_bound_pars=find(xparam1 < bounds.lb | xparam1 > bounds.ub);
        xparam1(outside_bound_pars) = bayestopt_.p1(outside_bound_pars);
    end
end

if ~options_.cova_compute
    stdh = NaN(length(xparam1),1);
end

if ~issmc(options_) && any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    % display results table and store parameter estimates and standard errors in results
    oo_ = display_estimation_results_table(xparam1, stdh, M_, options_, estim_params_, bayestopt_, oo_, prior_dist_names, 'Posterior', 'posterior');
    % Laplace approximation to the marginal log density:
    if options_.cova_compute
        estim_params_nbr = size(xparam1,1);
        if ispd(invhess)
            log_det_invhess = log(det(invhess./(stdh*stdh')))+2*sum(log(stdh));
            likelihood = feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
            oo_.MarginalDensity.LaplaceApproximation = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess - likelihood;
        else
            oo_.MarginalDensity.LaplaceApproximation = NaN;
        end
        skipline()
        dprintf('Log data density [Laplace approximation] is %f.', oo_.MarginalDensity.LaplaceApproximation)
        skipline()
    end
    if options_.dsge_var
        [~,~,~,~,~,~,~,oo_.dsge_var.posterior_mode.PHI_tilde,oo_.dsge_var.posterior_mode.SIGMA_u_tilde,oo_.dsge_var.posterior_mode.iXX,oo_.dsge_var.posterior_mode.prior] =...
            feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    end
elseif ~issmc(options_) && ~any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    oo_=display_estimation_results_table(xparam1, stdh, M_, options_, estim_params_, bayestopt_, oo_, prior_dist_names, 'Maximum Likelihood', 'mle');
end

if ~issmc(options_)
    invhess = set_mcmc_jumping_covariance(invhess, nx, options_.MCMC_jumping_covariance, bayestopt_, 'dynare_estimation_1');
end

%
% Run SMC sampler.
%

if issmc(options_)
    [posterior_sampler_options, options_, bayestopt_] = check_posterior_sampler_options([], M_.fname, M_.dname, options_, bounds, bayestopt_);
    options_.posterior_sampler_options.current_options = posterior_sampler_options;
    if isonline(options_)
        online_auxiliary_filter(xparam1, dataset_, options_, M_, estim_params_, bayestopt_, oo_);
    elseif ishssmc(options_)
        oo_.MarginalDensity.hssmc = hssmc(objective_function, bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_);
    elseif isdime(options_)
        dime(objective_function, xparam1, bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    elseif isdsmh(options_)
        oo_.MarginalDensity.dsmh = dsmh(objective_function, bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_);
    end 
end

%
% Run MCMC and compute posterior statistics.
%

if issmc(options_) || (any(bayestopt_.pshape>0) && options_.mh_replic) ||  (any(bayestopt_.pshape>0) && options_.load_mh_file) % not ML estimation
    if ~issmc(options_)
        % Reset bounds as lb and ub must only be operational during mode-finding
        bounds = set_mcmc_prior_bounds(xparam1, bayestopt_, options_, 'dynare_estimation_1');
        % Tune the jumping distribution's scale parameter
        if options_.mh_tune_jscale.status
            if strcmp(options_.posterior_sampler_options.posterior_sampling_method, 'random_walk_metropolis_hastings')
                options_.mh_jscale = tune_mcmc_mh_jscale_wrapper(invhess, options_, M_, objective_function, xparam1, bounds,...
                                                                 dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                bayestopt_.jscale(:) = options_.mh_jscale;
                fprintf('mh_jscale has been set equal to %s\n', num2str(options_.mh_jscale));
            else
                warning('mh_tune_jscale is only available with Random Walk Metropolis Hastings!')
            end
        end
        % Run MCMC
        if options_.mh_replic || options_.load_mh_file
            posterior_sampler_options = options_.posterior_sampler_options.current_options;
            posterior_sampler_options.invhess = invhess;
            [posterior_sampler_options, options_, bayestopt_] = check_posterior_sampler_options(posterior_sampler_options, M_.fname, M_.dname, options_, bounds, bayestopt_);
            % store current options in global
            options_.posterior_sampler_options.current_options = posterior_sampler_options;
            if options_.mh_replic
                ana_deriv_old = options_.analytic_derivation;
                options_.analytic_derivation = 0;
                posterior_sampler(objective_function,posterior_sampler_options.proposal_distribution,xparam1,posterior_sampler_options,bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,dispString);
                options_.analytic_derivation = ana_deriv_old;
            end
        end
        % Discard first mh_drop percent of the draws:
        CutSample(M_, options_, dispString);
    end
    if options_.mh_posterior_mode_estimation || (issmc(options_) && options_.smc_posterior_mode_estimation)
        [~, covariance, posterior_mode] = compute_posterior_covariance_matrix(bayestopt_.name, M_.fname, M_.dname, options_);
        oo_ = fill_mh_mode(posterior_mode, sqrt(diag(covariance)), M_, options_, estim_params_, oo_, 'posterior');
        %reset qz_criterium
        options_.qz_criterium = qz_criterium_old;
        return
    else
        % Get stored results if required
        if ~issmc(options_) && options_.load_mh_file && options_.load_results_after_load_mh
            oo_load_mh = load(sprintf('%s/%s/%s_results', M_.dname, 'Output', M_.fname), 'oo_');
        end
        % Compute MCMC convergence diagnostics
        if ~issmc(options_) && ~options_.nodiagnostic
            if (options_.mh_replic>0 || (options_.load_mh_file && ~options_.load_results_after_load_mh))
                oo_= mcmc_diagnostics(options_, estim_params_, M_,oo_);
            elseif options_.load_mh_file && options_.load_results_after_load_mh
                if isfield(oo_load_mh.oo_,'convergence')
                    oo_.convergence=oo_load_mh.oo_.convergence;
                end
            end
        elseif isdime(options_) && ~options_.nodiagnostic
            % provide plot of log densities over iterations
            oo_.lprob = trace_plot_dime(options_, M_);
        end
        % Estimation of the marginal density from the Mh draws:
        if isdsmh(options_) || ishssmc(options_) || isonline(options_) || isdime(options_) || options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            if ~issmc(options_)
                [~, oo_] = marginal_density(M_, options_, estim_params_, oo_, bayestopt_);
            end
            % Store posterior statistics by parameter name
            oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_, prior_dist_names);
            if ~options_.nograph
                oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_);
            end
            % Store posterior mean in a vector and posterior variance in
            % a matrix
            [oo_.posterior.metropolis.mean,oo_.posterior.metropolis.Variance] = GetPosteriorMeanVariance(options_, M_);
        elseif options_.load_mh_file && options_.load_results_after_load_mh
            % load fields from previous MCMC run stored in results-file
            field_names={'posterior_mode','posterior_std_at_mode',...% fields set by marginal_density
                         'posterior_mean','posterior_hpdinf','posterior_hpdsup','posterior_median','posterior_variance','posterior_std','posterior_deciles','posterior_density',...% fields set by GetPosteriorParametersStatistics
                         'prior_density',...%fields set by PlotPosteriorDistributions
                        };
            for field_iter=1:size(field_names,2)
                if isfield(oo_load_mh.oo_,field_names{1,field_iter})
                    oo_.(field_names{1,field_iter})=oo_load_mh.oo_.(field_names{1,field_iter});
                end
            end
            % field set by marginal_density
            if isfield(oo_load_mh.oo_,'MarginalDensity') && isfield(oo_load_mh.oo_.MarginalDensity,'ModifiedHarmonicMean')
                oo_.MarginalDensity.ModifiedHarmonicMean=oo_load_mh.oo_.MarginalDensity.ModifiedHarmonicMean;
            end
            % field set by GetPosteriorMeanVariance
            if isfield(oo_load_mh.oo_,'posterior') && isfield(oo_load_mh.oo_.posterior,'metropolis')
                oo_.posterior.metropolis=oo_load_mh.oo_.posterior.metropolis;
            end
        end
        [options_.sub_draws, error_flag]=set_number_of_subdraws(M_,options_); %check whether number is feasible
        if ~(~isempty(options_.sub_draws) && options_.sub_draws==0)
            if options_.bayesian_irf
                if error_flag
                    error('%s: I cannot compute the posterior IRFs!',dispString)
                end
                if options_.occbin.likelihood.status
                    fprintf('%s: the bayesian_irf option is not compatible with the use of OccBin.',dispString)
                else
                    oo_=PosteriorIRF('posterior',options_,estim_params_,oo_,M_,bayestopt_,dataset_,dataset_info,dispString);
                end                
            end
            if options_.moments_varendo
                if error_flag
                    error('%s: I cannot compute the posterior moments for the endogenous variables!',dispString)
                end
                if options_.load_mh_file
                    if issmc(options_)
                        error('%s: SMC does not yet support the load_mh_file option.\n',dispString);
                    elseif options_.mh_replic==0 %user wants to recompute results for standard MCMC
                        [MetropolisFolder, info] = CheckPath('metropolis',M_.dname);
                        if ~info
                            generic_post_data_file_name={'Posterior2ndOrderMoments','decomposition','PosteriorVarianceDecomposition','correlation','PosteriorCorrelations','conditional decomposition','PosteriorConditionalVarianceDecomposition'};
                            for ii=1:length(generic_post_data_file_name)
                                delete_stale_file([MetropolisFolder filesep M_.fname '_' generic_post_data_file_name{1,ii} '*']);
                            end
                            % restore compatibility for loading pre-4.6.2 runs where estim_params_ was not saved; see 6e06acc7 and !1944
                            NumberOfDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_posterior_draws*' ]));
                            if NumberOfDrawsFiles>0
                                temp=load([M_.dname '/metropolis/' M_.fname '_posterior_draws1']);
                                if ~isfield(temp,'estim_params_')
                                    for file_iter=1:NumberOfDrawsFiles
                                        save([M_.dname '/metropolis/' M_.fname '_posterior_draws' num2str(file_iter)],'estim_params_','-append')
                                    end
                                end
                            end
                        end
                    end
                end
                if options_.occbin.likelihood.status
                    fprintf('%s: the moments_varendo option is not compatible with the use of OccBin.',dispString)
                else
                    oo_ = compute_moments_varendo('posterior',options_,M_,oo_,estim_params_,var_list_);
                end
            end
            if options_.smoother || ~isempty(options_.filter_step_ahead) || options_.forecast
                if error_flag
                    error('%s: I cannot compute the posterior statistics!',dispString)
                end
                if options_.order==1 && ~options_.particle.status
                    oo_=prior_posterior_statistics('posterior',dataset_,dataset_info,M_,oo_,options_,estim_params_,bayestopt_,dispString); %get smoothed and filtered objects and forecasts
                else
                    error('%s: Particle Smoothers are not yet implemented.',dispString)
                end
            end
        else
            fprintf('%s: sub_draws was set to 0. Skipping posterior computations.\n',dispString);
        end
        xparam1 = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
        M_ = set_all_parameters(xparam1,estim_params_,M_);
    end
else
    M_ = set_all_parameters(xparam1,estim_params_,M_); %set to posterior mode
end

if options_.particle.status
    %reset qz_criterium
    options_.qz_criterium=qz_criterium_old;
    return
end

%Run and store classical smoother if needed
if options_.smoother && ... %Bayesian smoother requested before
    (any(bayestopt_.pshape > 0) && options_.mh_replic || ... % Bayesian with MCMC run
    any(bayestopt_.pshape > 0) && options_.load_mh_file) % Bayesian with loaded MCMC
    % nothing to do
elseif options_.partial_information ||...
    options_.order>1 %no particle smoother
    % smoothing not yet supported
else
    %% ML estimation, or posterior mode without Metropolis-Hastings or Metropolis without Bayesian smoothed variables
    oo_=save_display_classical_smoother_results(xparam1,M_,oo_,options_,bayestopt_,dataset_,dataset_info,estim_params_);
end

if options_.forecast == 0 || options_.mh_replic > 0 || options_.load_mh_file
    % nothing to do
elseif options_.order>2 && M_.exo_det_nbr > 0
    %forecasting not yet supported
else
    oo_.forecast = forecasts.run(var_list_,M_,options_,oo_,'smoother',dataset_info);
end

%reset qz_criterium
options_.qz_criterium=qz_criterium_old;

if reset_options_related_to_estimation
    options_.mode_compute = mode_compute_o;
    options_.mh_replic = mh_replic_o;
end
if first_obs_nan_indicator
    options_.first_obs=NaN;
end
