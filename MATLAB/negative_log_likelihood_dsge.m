function [negative_log_likelihood,exit_flag] = negative_log_likelihood_dsge(xparam1,DATAMAT,PARAM,MODEL,OPT)

% This function computes the negative log-likelihood of a DSGE model solved with
% perturbation at first order
% INPUTS
% - xparam1: vector of vales of estimated parameters (might be transformed to unbounded support for easier optimization)
% - DATAMAT: (nobs x nd) matrix with data
% - PARAM: Structure of Vector of all structural Parameters
% - MODEL: Structure with model information (also includes MODEL.Sigma and MODEL.Sigma_e)
% - OPT: Structure with options
% OUTPUTS
% - log_likelihood: value of log-likelihood function
% - exit_flag: equals to 1 if no errors occured

% initializations
exit_flag = 1;
large_number = Inf;
kf_variant = "gaussian"; % initialize

% Retransform parameters if needed

% Check bounds
if OPT.optimizer.bounds.penalize_objective
    if any(xparam1<OPT.optimizer.bounds.lb) || any(xparam1>OPT.optimizer.bounds.ub)
        negative_log_likelihood = large_number;
        exit_flag = 0;
        return
    end
end

% Update parameter structures
for jp = 1:MODEL.param_estim_nbr
    if any(ismember(MODEL.param_names(:,1),MODEL.param_estim_names{jp}))
        % update model parameters
        PARAM.(MODEL.param_estim_names{jp}) = xparam1(jp);
    elseif any(ismember("sqrt_Sigma_"+MODEL.varobs(:,1),MODEL.param_estim_names{jp}))
        idx = find(ismember("sqrt_Sigma_"+MODEL.endo_names_DR(MODEL.varobs_idx_DR),MODEL.param_estim_names{jp}));
        % update shock parameters
        MODEL.Sigma_eps(idx,idx) = xparam1(jp)^2;
    elseif any(ismember("sqrt_Sigma_"+MODEL.exo_names(:,1),MODEL.param_estim_names{jp}))
        idx = find(ismember("sqrt_Sigma_"+MODEL.exo_names(:,1),MODEL.param_estim_names{jp}));
        % update shock parameters
        MODEL.Sigma_eta(idx,idx) = xparam1(jp)^2;
    elseif any(ismember("Gamma_"+MODEL.exo_names(:,1),MODEL.param_estim_names{jp}))
        idx = find(ismember("Gamma_"+MODEL.exo_names(:,1),MODEL.param_estim_names{jp}));
        % update shock parameters
        MODEL.Gamma_eta(idx,idx) = xparam1(jp);        
    end
end
if any(diag(MODEL.Gamma_eta))
    kf_variant = "pruned_skewed";
end
MODEL.mu_eps = zeros(MODEL.varobs_nbr,1);
if kf_variant=="gaussian"
    MODEL.mu_eta  = zeros(MODEL.exo_nbr,1);
    MODEL.COV_eta = MODEL.Sigma_eta;
elseif kf_variant=="pruned_skewed"    
    MODEL.nu_eta    = zeros(MODEL.exo_nbr,1);
    MODEL.Delta_eta = eye(MODEL.exo_nbr);
    MODEL.mu_eta    = -csnMean(zeros(MODEL.exo_nbr,1),MODEL.Sigma_eta,MODEL.Gamma_eta,MODEL.nu_eta,MODEL.Delta_eta,OPT.cdfmvna_fct);
    MODEL.COV_eta   = csnVar(MODEL.Sigma_eta,MODEL.Gamma_eta,MODEL.nu_eta,MODEL.Delta_eta,OPT.cdfmvna_fct);
end

% Check if covariance matrices are positive definite
if isdiag(MODEL.COV_eta)
    if any(diag(MODEL.COV_eta)<0)
        negative_log_likelihood = large_number; exit_flag = 0; warning('COVeta not positive definite\n');
        return
    end
elseif ~ispd(MODEL.COV_eta)
    negative_log_likelihood = large_number; exit_flag = 0; warning('COVeta not positive definite\n');
    return
end
if isdiag(MODEL.Sigma_eps)
    if any(diag(MODEL.Sigma_eps)<0)
        negative_log_likelihood = large_number; exit_flag = 0; warning('Sigma_eps not positive definite\n');
        return
    end
elseif ~ispd(MODEL.Sigma_eps)
    negative_log_likelihood = large_number; exit_flag = 0; warning('Sigma_eps not positive definite\n');
    return
end

% check if theoretical bound on skewness coefficient is violated
for jp=1:MODEL.exo_nbr
    if abs(skewness_coef_theor(MODEL.Sigma_eta(jp,jp),MODEL.Gamma_eta(jp,jp))) > abs((sqrt(2)*(pi-4))/(pi-2)^(3/2))
        negative_log_likelihood = large_number; exit_flag = 0; warning('Skewness coefficient is out of bounds: %.4f\n',skewness_coef_theor(MODEL.Sigma_eta(jp,jp),MODEL.Gamma_eta(jp,jp)));
        return
    end
end

% compute linear approximation (first-order perturbation solution) around the deterministic steady state
[SOL,error_indicator] = get_first_order_perturbation_solution(MODEL,PARAM);

% Return, if no solution is found
if error_indicator
    negative_log_likelihood = large_number;
    exit_flag = 0;
    return
end

%% compute log likelihood

% initialize the state vector at the stationary distribution
mu_0 = zeros(MODEL.endo_nbr,1); % note that y are the model variables in deviation from steady-state, so the mean is zero by definition
%Sigma_0 = reshape( inv(eye(MODEL.endo_nbr*MODEL.endo_nbr) - kron(SOL.gx,SOL.gx))*reshape(SOL.gu*MODEL.COV_eta*SOL.gu',MODEL.endo_nbr*MODEL.endo_nbr,1) ,MODEL.endo_nbr,MODEL.endo_nbr); %analytical, but slow
Sigma_0 = dlyapdoubling(SOL.gx,SOL.gu*MODEL.COV_eta*SOL.gu'); % very fast and numerically accurate
if kf_variant == "gaussian"
    negative_log_likelihood = -1*kalman_gaussian(DATAMAT', mu_0,Sigma_0, SOL.gx,SOL.gu,MODEL.F, MODEL.mu_eta,MODEL.Sigma_eta, MODEL.mu_eps,MODEL.Sigma_eps, false, true);
elseif kf_variant == "pruned_skewed"
    Gamma_0 = zeros(MODEL.endo_nbr,MODEL.endo_nbr);
    nu_0    = zeros(MODEL.endo_nbr,1);
    Delta_0 = eye(MODEL.endo_nbr);
    %log_likelihood1 = kalman_csn(DATAMAT', mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, SOL.gx,SOL.gu,MODEL.F, MODEL.mu_eta,MODEL.Sigma_eta,MODEL.Gamma_eta,MODEL.nu_eta,MODEL.Delta_eta, MODEL.mu_eps,MODEL.Sigma_eps, OPT.cdfmvna_fct,OPT.prune_tol,false,true);
    negative_log_likelihood = -1*kalman_csn_dsge( ...
        SOL.gx, SOL.gu, MODEL.F, ...
        mu_0, Sigma_0, Gamma_0, nu_0, Delta_0, ...
        MODEL.mu_eps, MODEL.Sigma_eps, ...
        MODEL.mu_eta, MODEL.Sigma_eta, MODEL.Gamma_eta, MODEL.nu_eta, MODEL.Delta_eta, ...
        DATAMAT, true, OPT.prune_tol, OPT.cdfmvna_fct);
end

if isnan(negative_log_likelihood) || isinf(negative_log_likelihood) || ~isreal(negative_log_likelihood)
    negative_log_likelihood = large_number;
    exit_flag = 0;
end

end % main function end