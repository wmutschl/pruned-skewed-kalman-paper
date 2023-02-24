function [neg_log_likelihood,exit_flag] = dsge_negloglikelihood(xparam1,DATAMAT,PARAM,MODEL,OPT)

% This function computes the negative log-likelihood of a DSGE model solved with
% perturbation at first order
% INPUTS
% - xparam1: vector of vales of estimated parameters (might be transformed to unbounded support for easier optimization)
% - DATAMAT: (nobs x nd) matrix with data
% - PARAM: Structure of Vector of all structural Parameters
% - MODEL: Structure with model information (also includes MODEL.Sigma and MODEL.Sigma_e)
% - OPT: Structure with options
% OUTPUTS
% - neg_log_likelihood: value of negative log-likelihood function
% - exit_flag: equals to 1 if no errors occured

exit_flag = 1; % initialize

% Retransform parameters if needed
if OPT.optimizer.bounds.do_param_transform
    xparam1 = param_transform_bounded(xparam1,OPT.optimizer.bounds.lb,OPT.optimizer.bounds.ub);
end

% Update PARAM with xparam1
for j = 1:MODEL.param_estim_nbr
    PARAM.(MODEL.param_estim_names{j}) = xparam1(j);
end

% Check bounds
if OPT.optimizer.bounds.penalize_objective
    if any(xparam1<OPT.optimizer.bounds.lb) || any(xparam1>OPT.optimizer.bounds.ub)
        neg_log_likelihood = Inf;
        exit_flag = 0;
        return
    end
end
% Compute linear approximation around the deterministic steady state
[SOL,error_indicator] = get_first_order_perturbation_solution(MODEL,PARAM);

% Return, if no solution is found
if error_indicator
    neg_log_likelihood = Inf;
    exit_flag = 0;    
    return
end

[log_lik_t_tm1] = dsge_kalman_filter(DATAMAT,SOL.STEADY_STATE(MODEL.varobs_idx_DR),MODEL.varobs_selection_matrix,SOL.gx,SOL.gu,MODEL.Sigma_u,MODEL.Sigma_e);
neg_log_likelihood = -1*sum(log_lik_t_tm1);

if isnan(neg_log_likelihood) || isinf(neg_log_likelihood) || ~isreal(neg_log_likelihood)
    neg_log_likelihood = Inf;
    exit_flag = 0;
end

end % main function end