function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,M_,dr] = dsge_likelihood(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
% [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,M_,oo_] = dsge_likelihood(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,oo_,derivatives_info)
% Evaluates the negative of the posterior kernel of a DSGE model using the specified
% kalman_algo; the resulting posterior includes the 2*pi constant of the
% likelihood function
%
% INPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - dataset_            [structure]     dataset after transformations
% - dataset_info        [structure]     storing information about the
%                                       sample; not used but required for interface
% - options_            [structure]     MATLAB's structure describing the current options
% - M_                  [structure]     MATLAB's structure describing the model
% - estim_params_       [structure]     characterizing parameters to be estimated
% - bayestopt_          [structure]     describing the priors
% - BoundsInfo          [structure]     containing prior bounds
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state [vector]       steady state value for exogenous deterministic variables
% - derivatives_info    [structure]     derivative info for identification
%
% OUTPUTS
% - fval                    [double]        scalar, value of minus the likelihood or posterior kernel.
% - info                    [integer]       4×1 vector, information on whether solution and likelihood could be computed
% - exit_flag               [integer]       scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - DLIK                    [double]        Vector with score of the likelihood
% - Hess                    [double]        asymptotic Hessian matrix.
% - SteadyState             [double]        steady state level for the endogenous variables
% - trend_coeff             [double]        Matrix of doubles, coefficients of the deterministic trend in the measurement equation.
% - M_                      [struct]        Updated M_ structure described in INPUTS section.
% - dr                      [structure]     Reduced form model.
%
% This function is called by: dynare_estimation_1, mode_check,
% initial_estimation_checks, hessian, mr_hessian, hssmc, dime, dsmh,
% calibrate_mh_scale_parameter, posterior_sampler_initialization,
% identification.analysis, TaRB_optimizer_wrapper,
% posterior_sampler_iteration, slice_sampler, tempered_likelihood
% This function calls: dynare_resolve, lyapunov_symm, lyapunov_solver, compute_Pinf_Pstar, kalman_filter_d, missing_observations_kalman_filter_d,
% univariate_kalman_filter_d, kalman_steady_state, get_perturbation_params_deriv, kalman_filter, missing_observations_kalman_filter, univariate_kalman_filter, priordens

% Copyright © 2004-2025 Dynare Team
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

% Initial author: stephane DOT adjemian AT univ DASH lemans DOT FR

% Initialization of the returned variables and others...
fval        = [];
SteadyState = [];
trend_coeff = [];
exit_flag   = 1;
info        = zeros(4,1);

if options_.analytic_derivation
    DLIK        = NaN(1,length(xparam1));
else
    DLIK        = [];
end
Hess        = [];

% Ensure that xparam1 is a column vector.
% (Don't do the transformation if xparam1 is empty, otherwise it would become a
%  0×1 matrix, which create issues with older MATLABs when comparing with [] in
%  check_bounds_and_definiteness_estimation)
if ~isempty(xparam1)
    xparam1 = xparam1(:);
end

% Set flag related to analytical derivatives.
analytic_derivation = options_.analytic_derivation;

if analytic_derivation
    if options_.loglinear
        error('The analytic_derivation and loglinear options are not compatible')
    end
    if options_.endogenous_prior
        error('The analytic_derivation and endogenous_prior options are not compatible')
    end
end

if nargout==1
    analytic_derivation=0;
end

if analytic_derivation
    kron_flag=options_.analytic_derivation_mode;
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
is_restrict_state_space = true;
if options_.occbin.likelihood.status
    [occbin_options, options_.occbin.filter.state_covariance] = set_occbin_options(options_);
    if occbin_options.opts_simul.restrict_state_space
        [T,R,SteadyState,info,dr, M_.params,TTx,RRx,CCx, T0, R0] = ...
            occbin.dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,[],'restrict');
    else
        is_restrict_state_space = false;
        oldoo.restrict_var_list = dr.restrict_var_list;
        oldoo.restrict_columns = dr.restrict_columns;
        dr.restrict_var_list = bayestopt_.smoother_var_list;
        dr.restrict_columns = bayestopt_.smoother_restrict_columns;

        % Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
        [T,R,SteadyState,info,M_,dr, M_.params,TTx,RRx,CCx, T0, R0] = ...
            occbin.dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);

        dr.restrict_var_list = oldoo.restrict_var_list;
        dr.restrict_columns = oldoo.restrict_columns;

    end
    occbin_.status = true;
    occbin_.info= {options_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, M_, occbin_options, TTx, RRx, CCx,T0,R0};
else
    % Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
    if options_.kalman_algo == 5 % pruned skewed Kalman filter
        is_restrict_state_space = false;
        [T,R,SteadyState,info,dr, M_.params] = dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    else
        [T,R,SteadyState,info,dr, M_.params] = dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,'restrict');
    end
    occbin_.status = false;
end

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86 || ...
                info(1) == 401 || info(1) == 402 || info(1) == 403 || ... %cycle reduction
                info(1) == 411 || info(1) == 412 || info(1) == 413 % logarithmic reduction
        %meaningful second entry of output that can be used
        fval = Inf;
        if ~isfinite(info(2))
            info(4) = 0.1;
        else
            info(4) = info(2);
        end
        exit_flag = 0;
        if analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    end
end

% check endogenous prior restrictions
info=endogenous_prior_restrictions(T,R,M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
if info(1)
    fval = Inf;
    if ~isfinite(info(2))
        info(4) = 0.1;
    else
        info(4) = info(2);
    end
    exit_flag = 0;
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

if is_restrict_state_space
%% Define a vector of indices for the observed variables
    bayestopt_.mf = bayestopt_.mf1;
else
%get location of observed variables and requested smoothed variables in
%decision rules
    bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
end

% Define the constant vector of the measurement equation.
if options_.noconstant
    constant = zeros(dataset_.vobs,1);
else
    if options_.loglinear
        constant = log(SteadyState(bayestopt_.mfys));
    else
        constant = SteadyState(bayestopt_.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if bayestopt_.with_trend
    [trend_addition, trend_coeff]=compute_trend_coefficients(M_,options_,dataset_.vobs,dataset_.nobs);
    trend = repmat(constant,1,dataset_.nobs)+trend_addition;
else
    trend_coeff = zeros(dataset_.vobs,1);
    trend = repmat(constant,1,dataset_.nobs);
end

% Get needed information for Kalman filter routines.
start = options_.presample+1;
Z = bayestopt_.mf;           %selector for observed variables
no_missing_data_flag = ~dataset_info.missing.state;
mm = length(T);             %number of states
pp = dataset_.vobs;    %number of observables
rr = length(Q);             %number of shocks
kalman_tol = options_.kalman_tol;
diffuse_kalman_tol = options_.diffuse_kalman_tol;
riccati_tol = options_.riccati_tol;
Y = transpose(dataset_.data)-trend;
smpl = size(Y,2);

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------
kalman_algo = options_.kalman_algo;


diffuse_periods = 0;
expanded_state_vector_for_univariate_filter=0;
singular_diffuse_filter = 0;
if options_.heteroskedastic_filter
    Qvec=get_Qvec_heteroskedastic_filter(Q,smpl,M_);
end

switch options_.lik_init
  case 1% Standard initialization with the steady state of the state equation.
    if (kalman_algo~=2) && (kalman_algo~=5)
        % Use standard Kalman filter except if the univariate filter or pruned skewed filter is explicitly chosen.
        kalman_algo = 1;
    end
    Pstar=lyapunov_solver(T,R,Q,options_);
    Pinf  = [];
    a     = zeros(mm,1);
    a=set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
    a_0_given_tm1=T*a; %set state prediction for first Kalman step;

    if options_.occbin.likelihood.status || (kalman_algo == 5)
        Z =zeros(length(bayestopt_.mf),size(T,1));
        for i = 1:length(bayestopt_.mf)
            Z(i,bayestopt_.mf(i))=1;
        end
        Zflag = 1;
    else
        Zflag = 0;
    end
  case 2% Initialization with large numbers on the diagonal of the covariance matrix of the states (for non stationary models).
    if (kalman_algo~=2) && (kalman_algo~=5)
        % Use standard Kalman filter except if the univariate filter or pruned skewed filter is explicitly chosen.
        kalman_algo = 1;
    end
    Pstar = options_.Harvey_scale_factor*eye(mm);
    Pinf  = [];
    a     = zeros(mm,1);
    a = set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
    a_0_given_tm0 = a;
    if not(options_.occbin.likelihood.status && options_.occbin.likelihood.first_period_occbin_update==1)
        a_0_given_tm1 = T*a; %set state prediction for first Kalman step;
        if options_.Harvey_scale_factor==0
            Pstar = T*Pstar*T' + R*Q*R';
        end
    else
        a_0_given_tm1 = a; 
    end
    if options_.occbin.likelihood.status || (kalman_algo == 5)
        Z =zeros(length(bayestopt_.mf),size(T,1));
        for i = 1:length(bayestopt_.mf)
            Z(i,bayestopt_.mf(i))=1;
        end
        Zflag = 1;
    else
        Zflag = 0;
    end
  case 3% Diffuse Kalman filter (Durbin and Koopman)
    % Use standard Kalman filter except if the univariate filter is explicitly chosen.
    if kalman_algo == 0
        kalman_algo = 3;
    elseif ~( (kalman_algo==3) || (kalman_algo==4) )
        error(['The model requires Diffuse filter, but you specified a different Kalman filter. You must set options_.kalman_algo ' ...
               'to 0 (default), 3 or 4'])
    end
    [Pstar,Pinf] = compute_Pinf_Pstar(Z,T,R,Q,options_.qz_criterium);
    Z =zeros(length(bayestopt_.mf),size(T,1));
    for i = 1:length(bayestopt_.mf)
        Z(i,bayestopt_.mf(i))=1;
    end
    Zflag = 1;
    if options_.heteroskedastic_filter
        QQ=Qvec;
    else
        QQ=Q;
    end
    % Run diffuse Kalman filter on first periods.
    if (kalman_algo==3)
        % Multivariate Diffuse Kalman Filter
        a = zeros(mm,1);
        a = set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
        a_0_given_tm1 = T*a; %set state prediction for first Kalman step;
        Pstar0 = Pstar; % store Pstar
        if no_missing_data_flag
            [dLIK,dlik,a_0_given_tm1,Pstar] = kalman_filter_d(Y, 1, size(Y,2), ...
                                                  a_0_given_tm1, Pinf, Pstar, ...
                                                  kalman_tol, diffuse_kalman_tol, riccati_tol, options_.presample, ...
                                                  T,R,QQ,H,Z,mm,pp,rr);
        else
            [dLIK,dlik,a_0_given_tm1,Pstar] = missing_observations_kalman_filter_d(dataset_info.missing.aindex,dataset_info.missing.number_of_observations,dataset_info.missing.no_more_missing_observations, ...
                                                              Y, 1, size(Y,2), ...
                                                              a_0_given_tm1, Pinf, Pstar, ...
                                                              kalman_tol, diffuse_kalman_tol, riccati_tol, options_.presample, ...
                                                              T,R,QQ,H,Z,mm,pp,rr);
        end
        diffuse_periods = length(dlik);
        if isinf(dLIK)
            % Go to univariate diffuse filter if singularity problem.
            singular_diffuse_filter = 1;
            Pstar = Pstar0;
        end
    end
    if singular_diffuse_filter || (kalman_algo==4)
        % Univariate Diffuse Kalman Filter
        if isequal(H,0)
            H1 = zeros(pp,1);
            mmm = mm;
        else
            if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                H1 = diag(H);
                mmm = mm;
            else
                %Augment state vector (follows Section 6.4.3 of DK (2012))
                expanded_state_vector_for_univariate_filter=1;
                if Zflag
                    Z1=Z;
                else
                    Z1=zeros(pp,size(T,2));
                    for jz=1:length(Z)
                        Z1(jz,Z(jz))=1;
                    end
                end
                Z = [Z1, eye(pp)];
                Zflag=1;
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blkdiag(Pinf,zeros(pp));
                H1 = zeros(pp,1);
                mmm   = mm+pp;
                if options_.heteroskedastic_filter
                    clear QQ
                    for kv=1:size(Qvec,3)
                        QQ(:,:,kv) = blkdiag(Qvec(:,:,kv),H);
                    end
                    Qvec=QQ;
                else
                    QQ = Q;
                end
            end
        end

        a = zeros(mmm,1);
        a = set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
        a_0_given_tm1 = T*a;
        [dLIK,dlik,a_0_given_tm1,Pstar] = univariate_kalman_filter_d(dataset_info.missing.aindex,...
                                                         dataset_info.missing.number_of_observations,...
                                                         dataset_info.missing.no_more_missing_observations, ...
                                                         Y, 1, size(Y,2), ...
                                                         a_0_given_tm1, Pinf, Pstar, ...
                                                         kalman_tol, diffuse_kalman_tol, riccati_tol, options_.presample, ...
                                                         T,R,QQ,H1,Z,mmm,pp,rr);
        diffuse_periods = size(dlik,1);
    end
    if isnan(dLIK)
        fval = Inf;
        info(1) = 45;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end

  case 4% Start from the solution of the Riccati equation.
    if (kalman_algo~=2) && (kalman_algo~=5)
        % Use standard Kalman filter except if the univariate filter or pruned skewed filter is explicitly chosen.
        kalman_algo = 1;
    end
    try
        if isequal(H,0)
            Pstar = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))));
        else
            Pstar = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))),H);
        end
    catch ME
        disp(ME.message)
        disp('dsge_likelihood:: I am not able to solve the Riccati equation, so I switch to lik_init=1!');
        options_.lik_init = 1;
        Pstar=lyapunov_solver(T,R,Q,options_);
    end
    Pinf  = [];
    a     = zeros(mm,1);
    a = set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
    a_0_given_tm1 = T*a;
    if options_.occbin.likelihood.status || (kalman_algo == 5)
        Z =zeros(length(bayestopt_.mf),size(T,1));
        for i = 1:length(bayestopt_.mf)
            Z(i,bayestopt_.mf(i))=1;
        end
        Zflag = 1;
    else
        Zflag = 0;
    end
  case 5            % Old diffuse Kalman filter only for the non stationary variables
    [eigenvect, eigenv] = eig(T);
    eigenv = diag(eigenv);
    nstable = length(find(abs(abs(eigenv)-1) > 1e-7));
    V = eigenvect(:,abs(abs(eigenv)-1) < 1e-7);
    stable = find(sum(abs(V),2)<1e-5);
    nunit = length(eigenv) - nstable;
    Pstar = options_.Harvey_scale_factor*eye(nunit);
    if (kalman_algo~=2) && (kalman_algo~=5)
        % Use standard Kalman filter except if the univariate filter or pruned skewed filter is explicitly chosen.
        kalman_algo = 1;
    end
    R_tmp = R(stable, :);
    T_tmp = T(stable,stable);
    Pstar_tmp=lyapunov_solver(T_tmp,R_tmp,Q,options_);
    Pstar(stable, stable) = Pstar_tmp;
    Pinf  = [];
    a = zeros(mm,1);
    a = set_Kalman_starting_values(a,M_,dr,options_,bayestopt_);
    a_0_given_tm1 = T*a;
    if options_.occbin.likelihood.status || (kalman_algo == 5)
        Z =zeros(length(bayestopt_.mf),size(T,1));
        for i = 1:length(bayestopt_.mf)
            Z(i,bayestopt_.mf(i))=1;
        end
        Zflag = 1;
    else
        Zflag = 0;
    end
  otherwise
    error('dsge_likelihood:: Unknown initialization approach for the Kalman filter!')
end

if analytic_derivation
    offset = estim_params_.nvx;
    offset = offset+estim_params_.nvn;
    offset = offset+estim_params_.ncx;
    offset = offset+estim_params_.ncn;
    no_DLIK = 0;
    full_Hess = analytic_derivation==2;
    asy_Hess = analytic_derivation==-2;
    outer_product_gradient = analytic_derivation==-1;
    if asy_Hess
        analytic_derivation=1;
    end
    if outer_product_gradient
        analytic_derivation=1;
    end
    DLIK = [];
    iv = dr.restrict_var_list;
    if nargin<13 || isempty(derivatives_info)
        [~,~,~,~,dr, M_.params] = dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        if ~isempty(estim_params_.var_exo)
            indexo=estim_params_.var_exo(:,1);
        else
            indexo=[];
        end
        if ~isempty(estim_params_.param_vals)
            indparam=estim_params_.param_vals(:,1);
        else
            indparam=[];
        end
        old_order = options_.order;
        if options_.order > 1%not sure whether this check is necessary
            options_.order = 1; fprintf('Reset order to 1 for analytical parameter derivatives.\n');
        end
        old_analytic_derivation_mode = options_.analytic_derivation_mode;
        options_.analytic_derivation_mode = kron_flag;
        if full_Hess
            DERIVS = identification.get_perturbation_params_derivs(M_, options_, estim_params_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, indparam, indexo, [], true);
            indD2T = reshape(1:M_.endo_nbr^2, M_.endo_nbr, M_.endo_nbr);
            indD2Om = dyn_unvech(1:M_.endo_nbr*(M_.endo_nbr+1)/2);
            D2T = DERIVS.d2KalmanA(indD2T(iv,iv),:);
            D2Om = DERIVS.d2Om(dyn_vech(indD2Om(iv,iv)),:);
            D2Yss = DERIVS.d2Yss(iv,:,:);
        else
            DERIVS = identification.get_perturbation_params_derivs(M_, options_, estim_params_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, indparam, indexo, [], false);
        end
        DT = zeros(M_.endo_nbr, M_.endo_nbr, size(DERIVS.dghx,3));
        DT(:,M_.nstatic+(1:M_.nspred),:) = DERIVS.dghx;
        DT = DT(iv,iv,:);
        DOm = DERIVS.dOm(iv,iv,:);
        DYss = DERIVS.dYss(iv,:);
        options_.order = old_order; %make sure order is reset (not sure if necessary)
        options_.analytic_derivation_mode = old_analytic_derivation_mode;%make sure analytic_derivation_mode is reset (not sure if necessary)
    else
        DT = derivatives_info.DT(iv,iv,:);
        DOm = derivatives_info.DOm(iv,iv,:);
        DYss = derivatives_info.DYss(iv,:);
        if isfield(derivatives_info,'full_Hess')
            full_Hess = derivatives_info.full_Hess;
        end
        if full_Hess
            D2T = derivatives_info.D2T;
            D2Om = derivatives_info.D2Om;
            D2Yss = derivatives_info.D2Yss;
        end
        if isfield(derivatives_info,'no_DLIK')
            no_DLIK = derivatives_info.no_DLIK;
        end
        clear('derivatives_info');
    end
    DYss = [zeros(size(DYss,1),offset) DYss];
    DH=zeros([length(H),length(H),length(xparam1)]);
    DQ=zeros([size(Q),length(xparam1)]);
    DP=zeros([size(T),length(xparam1)]);
    if full_Hess
        for j=1:size(D2Yss,1)
            tmp(j,:,:) = blkdiag(zeros(offset,offset), squeeze(D2Yss(j,:,:)));
        end
        D2Yss = tmp;
        D2H=sparse(size(D2Om,1),size(D2Om,2)); %zeros([size(H),length(xparam1),length(xparam1)]);
        D2P=sparse(size(D2Om,1),size(D2Om,2)); %zeros([size(T),length(xparam1),length(xparam1)]);
        jcount=0;
    end
    if options_.lik_init==1
        for i=1:estim_params_.nvx
            k =estim_params_.var_exo(i,1);
            DQ(k,k,i) = 2*sqrt(Q(k,k));
            dum =  lyapunov_symm(T,DOm(:,:,i),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
            %         kk = find(abs(dum) < 1e-12);
            %         dum(kk) = 0;
            DP(:,:,i)=dum;
            if full_Hess
                for j=1:i
                    jcount=jcount+1;
                    dum =  lyapunov_symm(T,dyn_unvech(D2Om(:,jcount)),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
                    %             kk = (abs(dum) < 1e-12);
                    %             dum(kk) = 0;
                    D2P(:,jcount)=dyn_vech(dum);
                    %             D2P(:,:,j,i)=dum;
                end
            end
        end
    end
    offset = estim_params_.nvx;
    for i=1:estim_params_.nvn
        k = estim_params_.var_endo(i,1);
        DH(k,k,i+offset) = 2*sqrt(H(k,k));
        if full_Hess
            D2H(k,k,i+offset,i+offset) = 2;
        end
    end
    offset = offset + estim_params_.nvn;
    if options_.lik_init==1
        for j=1:estim_params_.np
            dum =  lyapunov_symm(T,DT(:,:,j+offset)*Pstar*T'+T*Pstar*DT(:,:,j+offset)'+DOm(:,:,j+offset),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
            %         kk = find(abs(dum) < 1e-12);
            %         dum(kk) = 0;
            DP(:,:,j+offset)=dum;
            if full_Hess
                DTj = DT(:,:,j+offset);
                DPj = dum;
                for i=1:j+offset
                    jcount=jcount+1;
                    DTi = DT(:,:,i);
                    DPi = DP(:,:,i);
                    D2Tij = reshape(D2T(:,jcount),size(T));
                    D2Omij = dyn_unvech(D2Om(:,jcount));
                    tmp = D2Tij*Pstar*T' + T*Pstar*D2Tij' + DTi*DPj*T' + DTj*DPi*T' + T*DPj*DTi' + T*DPi*DTj' + DTi*Pstar*DTj' + DTj*Pstar*DTi' + D2Omij;
                    dum = lyapunov_symm(T,tmp,options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
                    D2P(:,jcount) = dyn_vech(dum);
                end
            end
        end
    end
    if analytic_derivation==1
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP,asy_Hess};
    else
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P};
        clear DT DYss DOm DP D2T D2Yss D2Om D2H D2P
    end
else
    analytic_deriv_info={0};
end

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------
if options_.heteroskedastic_filter
    Q=Qvec;
end

singularity_has_been_detected = false;
% First test multivariate filter if specified; potentially abort and use univariate filter instead
if ((kalman_algo==1) || (kalman_algo==3)) || (kalman_algo == 5) % Multivariate Kalman Filter
    if no_missing_data_flag && ~options_.occbin.likelihood.status
        if options_.fast_kalman_filter
            if diffuse_periods
                %kalman_algo==3 requires no diffuse periods (stationary
                %observables) as otherwise FE matrix will not be positive
                %definite
                fval = Inf;
                info(1) = 55;
                info(4) = 0.1;
                exit_flag = 0;
                return
            end
            [LIK,lik] = kalman_filter_fast(Y,diffuse_periods+1,size(Y,2), ...
                                           a_0_given_tm1,Pstar, ...
                                           kalman_tol, ...
                                           options_.presample, ...
                                           T,H,Z,pp,Zflag,diffuse_periods);
        else
            if options_.kalman_algo == 5 % multivariate pruned skewed Kalman filter
                startTime = tic;
                [LIK, lik] = kalman_filter_pruned_skewed(Y, diffuse_periods+1, size(Y,2), options_.presample, ... % data, start, last, presample
                                                         a, Pstar, zeros(size(Pstar)), zeros(size(a)), eye(size(a,1)), ... % initialize CSN at Gaussian distribution
                                                         T, R, Z, ... % state space matrices
                                                         M_.csn.mu_e, M_.csn.Sigma_e, M_.csn.Gamma_e, M_.csn.nu_e, M_.csn.Delta_e, H, ... % shock CSN parameters, measurement error covariance
                                                         kalman_tol, options_.rescale_prediction_error_covariance,... % Gaussian Kalman filter options
                                                         options_.skewed_kalman.prune_tol, options_.skewed_kalman.mvnlogcdf, options_.skewed_kalman.rank_deficiency_transform, options_.debug); % skewed Kalman filter options
                fprintf('Time to evaluate likelihood with Pruned Skewed Kalman Filter with tol=%.4f: %s\n', options_.skewed_kalman.prune_tol, dynsec2hms(toc(startTime)));

            else % multivariate Gaussian Kalman filter
                if options_.kalman_filter_mex
                    [LIK,lik] = kalman_filter_mex(Y,a_0_given_tm1,Pstar, ...
                                                  kalman_tol, riccati_tol, ...
                                                  T,Q,R,Z,Zflag,H,diffuse_periods, ...
                                                  options_.presample);
                else
                    [LIK,lik] = kalman_filter(Y,diffuse_periods+1,size(Y,2), ...
                                              a_0_given_tm1,Pstar, ...
                                              kalman_tol, riccati_tol, ...
                                              options_.rescale_prediction_error_covariance, ...
                                              options_.presample, ...
                                              T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods, ...
                                              analytic_deriv_info{:});
                end
            end
        end
    else
        [LIK,lik] = missing_observations_kalman_filter(dataset_info.missing.aindex,dataset_info.missing.number_of_observations,dataset_info.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                                       a_0_given_tm1, Pstar, ...
                                                       kalman_tol, options_.riccati_tol, ...
                                                       options_.rescale_prediction_error_covariance, ...
                                                       options_.presample, ...
                                                       T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods, occbin_);
        if occbin_.status && isinf(LIK)
            fval = Inf;
            info(1) = 320;
            info(4) = 0.1;
            exit_flag = 0;
            return
        end
    end
    if analytic_derivation
        LIK1=LIK;
        LIK=LIK1{1};
        lik1=lik;
        lik=lik1{1};
    end
    if isinf(LIK)
        if options_.use_univariate_filters_if_singularity_is_detected
            singularity_has_been_detected = true;
            if kalman_algo == 1
                kalman_algo = 2;
            else
                kalman_algo = 4;
            end
        else
            fval = Inf;
            info(1) = 50;
            info(4) = 0.1;
            exit_flag = 0;
            return
        end
    else
        if options_.lik_init==3
            LIK = LIK + dLIK;
            if analytic_derivation==0 && nargout>3
                if ~singular_diffuse_filter
                    lik = [dlik; lik];
                else
                    lik = [sum(dlik,2); lik];
                end
            end
        end
    end
end

if (kalman_algo==2) || (kalman_algo==4)
    % Univariate Kalman Filter
    % resetting measurement error covariance matrix when necessary following DK (2012), Section 6.4.3                                                          %
    if isequal(H,0)
        H1 = zeros(pp,1);
        mmm = mm;
        if analytic_derivation
            DH = zeros(pp,length(xparam1));
        end
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H1 = diag(H);
            mmm = mm;
            clear('tmp')
            if analytic_derivation
                for j=1:pp
                    tmp(j,:)=DH(j,j,:);
                end
                DH=tmp;
            end
        else
            if ~expanded_state_vector_for_univariate_filter
                Z1=zeros(pp,size(T,2));
                for jz=1:length(Z)
                    Z1(jz,Z(jz))=1;
                end
                Z = [Z1, eye(pp)];
                Zflag=1;
                T = blkdiag(T,zeros(pp));
                if options_.heteroskedastic_filter
                    clear Q
                    for kv=1:size(Qvec,3)
                        Q(:,:,kv) = blkdiag(Qvec(:,:,kv),H);
                    end
                else
                    Q = blkdiag(Q,H);
                end
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blkdiag(Pinf,zeros(pp));
                H1 = zeros(pp,1);
                Zflag=1;
            end
            mmm   = mm+pp;
            if singularity_has_been_detected
                a_tmp = zeros(mmm,1);
                a_tmp(1:length(a_0_given_tm1)) = a_0_given_tm1;
                a_0_given_tm1 = a_tmp;
            elseif ~expanded_state_vector_for_univariate_filter
                a_0_given_tm1 = [a_0_given_tm1; zeros(pp,1)];
            end
        end
    end
    if analytic_derivation
        analytic_deriv_info{5}=DH;
    end
    [LIK, lik] = univariate_kalman_filter(dataset_info.missing.aindex,dataset_info.missing.number_of_observations,dataset_info.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                          a_0_given_tm1,Pstar, ...
                                          options_.kalman_tol, ...
                                          options_.riccati_tol, ...
                                          options_.presample, ...
                                          T,Q,R,H1,Z,mmm,pp,rr,Zflag,diffuse_periods,analytic_deriv_info{:});
    if analytic_derivation
        LIK1=LIK;
        LIK=LIK1{1};
        lik1=lik;
        lik=lik1{1};
    end
    if options_.lik_init==3
        LIK = LIK+dLIK;
        if analytic_derivation==0 && nargout>3
            lik = [dlik; lik];
        end
    end
end

if analytic_derivation
    if no_DLIK==0
        DLIK = LIK1{2};
    end
    if full_Hess
        Hess = -LIK1{3};
    end
    if asy_Hess
        Hess = LIK1{3};
    end
end

if isnan(LIK)
    fval = Inf; info(1) = 45; info(4) = 0.1; exit_flag = 0;
    return
end
if imag(LIK)~=0
    fval = Inf; info(1) = 46; info(4) = 0.1; exit_flag = 0;
    return
end
if isinf(LIK)
    fval = Inf; info(1) = 50; info(4) = 0.1; exit_flag = 0;
    return
end
likelihood = LIK;

% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------
if analytic_derivation
    if full_Hess
        [lnprior, dlnprior, d2lnprior] = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
        Hess = Hess - d2lnprior;
    else
        [lnprior, dlnprior] = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    end
    if no_DLIK==0
        DLIK = DLIK - dlnprior';
    end
    if outer_product_gradient
        dlik = lik1{2};
        dlik=[- dlnprior; dlik(start:end,:)];
        Hess = dlik'*dlik;
    end
else
    lnprior = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
end
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

if options_.endogenous_prior==1
    if options_.lik_init==2 || options_.lik_init==3
        error('Endogenous prior not supported with non-stationary models')
    else
        [lnpriormom]  = endogenous_prior(Y,dataset_info,Pstar,bayestopt_,H);
        fval    = (likelihood-lnprior-lnpriormom);
    end
elseif options_.init_state_endogenous_prior    
    if not(options_.lik_init==2 && options_.Harvey_scale_factor==0)
        error('Init state endogenous prior not supported without conditional likelihood')
    else
        [lnpriorendoinitstate, lnpriorinitstate]  = init_state_endogenous_prior(a_0_given_tm0,T,R,Q,xparam1,bayestopt_,options_);
    end
    fval    = likelihood-lnpriorendoinitstate - (lnprior-lnpriorinitstate);

else
    fval    = (likelihood-lnprior);
end

if options_.prior_restrictions.status
    tmp = feval(options_.prior_restrictions.routine, M_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, options_, dataset_, dataset_info);
    fval = fval - tmp;
end

if analytic_derivation==0 && nargout>3
    lik=lik(start:end,:);
    DLIK=[-lnprior; lik(:)];
end

function[occbin_options, occbin_filter_state_covariance] = set_occbin_options(options_)

% this builds the opts_simul options field needed by occbin.solver
occbin_options.opts_simul = options_.occbin.simul;
occbin_options.opts_simul.curb_retrench = options_.occbin.likelihood.curb_retrench;
occbin_options.opts_simul.maxit = options_.occbin.likelihood.maxit;
occbin_options.opts_simul.periods = options_.occbin.likelihood.periods;
occbin_options.opts_simul.check_ahead_periods = options_.occbin.likelihood.check_ahead_periods;
occbin_options.opts_simul.max_check_ahead_periods = options_.occbin.likelihood.max_check_ahead_periods;
occbin_options.opts_simul.periodic_solution = options_.occbin.likelihood.periodic_solution;
occbin_options.opts_simul.restrict_state_space = options_.occbin.likelihood.restrict_state_space;

occbin_options.opts_simul.full_output = options_.occbin.likelihood.full_output;
occbin_options.opts_simul.piecewise_only = options_.occbin.likelihood.piecewise_only;
occbin_options.opts_regime.init_binding_indicator = options_.occbin.likelihood.init_binding_indicator;
occbin_options.opts_regime.init_regime_history=options_.occbin.likelihood.init_regime_history;

% this checks particle filter options
if options_.occbin.filter.init_periods_using_particles || options_.occbin.filter.particle.status
    occbin_filter_state_covariance= true;
else
    occbin_filter_state_covariance= false;
end
