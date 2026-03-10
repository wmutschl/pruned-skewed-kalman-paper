function [exo_simul, diagnostics, best_draws, best_diag] = generate_exo_shocks_matched(periods, Sigma_e, Skew_e, csn)
% function [exo_simul, diagnostics] = generate_exo_shocks_matched(T, Sigma_e, Skew_e, varargin)
% -------------------------------------------------------------------------
% Generates independent exogenous shock series and redraws them until the
% sample moments are sufficiently close to the desired targets.
%
% The output exo_simul has dimension [n_exo by T], matching Dynare's
% convention for oo_.exo_simul.
%
% Supported distributions:
% - Gaussian: target skewness is ignored and set to zero
% - CSN: independent univariate CSN shocks with target variance/skewness
%
% -------------------------------------------------------------------------
% Name-value pairs:
% - 'distribution'   : 'auto' (default), 'gaussian', 'csn'
% - 'max_tries'      : redraw budget, default 200
% - 'mean_tol'       : max |sample mean| / target std, default 0.05
% - 'corr_tol'       : max absolute off-diagonal sample correlation,
%                      default 0.05
% - 'var_tol'        : max relative variance error, default 0.05
% - 'skew_tol'       : max absolute skewness error, default 0.10
% - 'seed'           : RNG seed, default []
% - 'antithetic'     : true/false, only for Gaussian, default false
% - 'verbose'        : true/false, default false
% - 'csn'            : optional struct with fields mu_e, Sigma_e, Gamma_e,
%                      nu_e, Delta_e. If omitted, a diagonal univariate CSN
%                      specification is constructed from Sigma_e and Skew_e.
% -------------------------------------------------------------------------

max_tries = 20000;
mean_tol = 0.05;
corr_tol = 0.05;
var_tol = 0.05;
skew_tol = 0.10;
seed = [];

target_var = diag(Sigma_e);
target_skew = zeros(size(Sigma_e,1),1);
target_skew(Skew_e(:,1),1) = Skew_e(:,4);

if ~isempty(seed)
    rng(seed);
end

best_score = inf;

for attempt = 1:max_tries
    if size(Skew_e, 1) > 0 % draw from skew normal distribution (special case of closed skew normal, see csn_update_specification.m for details)
        exo_candidate = rand_multivariate_csn(periods, csn.mu_e, csn.Sigma_e, csn.Gamma_e, csn.nu_e, csn.Delta_e);
    else % draw from Gaussian distribution
        exo_candidate = transpose(chol(Sigma_e))*randn(size(Sigma_e,1),periods);
    end 
    % recenter and rescale
    exo_candidate = exo_candidate - mean(exo_candidate,2);
    sample_std = std(exo_candidate,1,2);
    exo_candidate = exo_candidate ./ sample_std .* sqrt(target_var);
    % compute diagnostics
    sample_mean = mean(exo_candidate, 2);
    sample_var = var(exo_candidate, 1, 2);
    sample_skew = skewness(exo_candidate, 1, 2);
    if size(exo_candidate, 1) == 1
        corr_max_abs = 0;
    else
        sample_corr = corrcoef(exo_candidate');
        corr_offdiag = sample_corr - eye(size(exo_candidate, 1));
        corr_max_abs = max(abs(corr_offdiag(:)));
    end
    mean_max_rel = max(abs(sample_mean) ./ sqrt(target_var));
    var_rel_err = abs(sample_var - target_var) ./ target_var;
    skew_abs_err = abs(sample_skew - target_skew);

    diagnostics_candidate.mean_max_rel = max(mean_max_rel);
    diagnostics_candidate.var_max_rel = max(var_rel_err);
    diagnostics_candidate.skew_max_abs = max(skew_abs_err);
    diagnostics_candidate.corr_max_abs = corr_max_abs;
    diagnostics_candidate.accepted = diagnostics_candidate.mean_max_rel <= mean_tol ...
            && diagnostics_candidate.var_max_rel <= var_tol ...
            && diagnostics_candidate.skew_max_abs <= skew_tol ...
            && diagnostics_candidate.corr_max_abs <= corr_tol;

    score = diagnostics_candidate.mean_max_rel ...
          + diagnostics_candidate.corr_max_abs ...
          + diagnostics_candidate.var_max_rel ...
          + diagnostics_candidate.skew_max_abs;

    if score < best_score
        best_score = score;
        best_draws = exo_candidate;
        best_diag = diagnostics_candidate;
    end

    if diagnostics_candidate.accepted
        %fprintf('Found shocks....\n');
        exo_simul = exo_candidate;
        diagnostics = diagnostics_candidate;
        return
    end
    if attempt == max_tries
        error('Maximum number of tries reached....');
    end
end

end