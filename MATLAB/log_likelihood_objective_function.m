function [log_likelihood,exit_flag] = log_likelihood_objective_function(xparam1,y,OPT,PARAMS,kf_variant)

%% initializations
exit_flag = 1;
verysmall = -Inf;

%% update parameters with xparam1
idx = 0;
for jvar = OPT.param_names_estim(:,1)'
    if jvar == "G"
        G = reshape(xparam1(idx+(1:OPT.x_nbr^2)),OPT.x_nbr,OPT.x_nbr);
        idx = idx + OPT.x_nbr^2;
    end
    if jvar == "F"
        F = reshape(xparam1(idx+(1:OPT.y_nbr*OPT.x_nbr)),OPT.y_nbr,OPT.x_nbr);
        idx = idx + OPT.y_nbr*OPT.x_nbr;
    end
    if jvar == "R"
        R = reshape(xparam1(idx+(1:OPT.x_nbr*OPT.eta_nbr)),OPT.x_nbr,OPT.eta_nbr);
        idx = idx + OPT.x_nbr*OPT.eta_nbr;
    end
    if jvar == "mu_eps"
        mu_eps = xparam1(idx+(1:OPT.eps_nbr));
        idx = idx + OPT.eps_nbr;
    end
    if jvar == "diaglogSigma_eps"
        Sigma_eps = diag(exp(xparam1(idx+(1:OPT.eps_nbr)))); %undo log transform
        idx = idx + OPT.eps_nbr;
    end
    if jvar == "mu_eta"
        mu_eta = xparam1(idx+(1:OPT.eta_nbr));
        idx = idx + OPT.eta_nbr;
    end
    if jvar == "diaglogSigma_eta"
        Sigma_eta = diag(exp(xparam1(idx+(1:OPT.eta_nbr)))); %undo log transform
        idx = idx + OPT.eta_nbr;
    end    
    if strcmp(kf_variant{1},"pruned_skewed")
        if jvar == "diagGamma_eta"
            Gamma_eta = diag(xparam1(idx+(1:OPT.eta_nbr))); % make full matrix
            idx = idx + OPT.eta_nbr;
        end
    end
end

%% set calibrated parameters
for jvar = OPT.param_names_fixed(:,1)'
    if any(ismember(fieldnames(PARAMS),jvar))
        eval(sprintf('%s = PARAMS.%s;',jvar,jvar));
    end
end

%% Check stability and positive definitenes
if sum(abs(eig(G)) >= (1-1e-7))
    % disp("Some eigenvalues of G are outside the unit circle")
    log_likelihood = verysmall;
    exit_flag = 0;
    return
end

% Check if covariance matrices are positive definite
if strcmp(kf_variant{1},"pruned_skewed")
    COVeta = csnVar(Sigma_eta,Gamma_eta,nu_eta,Delta_eta,OPT.cdfmvna_fct);
elseif strcmp(kf_variant{1},"gaussian")
    COVeta = Sigma_eta;
end
if isdiag(COVeta)
    if any(diag(COVeta)<0)
        log_likelihood = verysmall;
        exit_flag = 0;
        return
    end
elseif ~ispd(COVeta)
    log_likelihood = verysmall;
    exit_flag = 0;
    return
end

if isdiag(Sigma_eps)
    if any(diag(Sigma_eps)<0)
        log_likelihood = verysmall;
        exit_flag = 0;
        return
    end
elseif ~ispd(Sigma_eps)
    log_likelihood = verysmall;
    exit_flag = 0;
    return
end

%% compute log likelihood
if strcmp(kf_variant{1},"gaussian")
    % initialize Kalman filter with a wide Normal prior: x_0 ~ CSN(0,Harvey_factor*eye(nx),0,0,eye(nx)) = N(0,10*eye(nx))
    mu_0     = zeros(OPT.x_nbr,1);
    Sigma_0  = OPT.Harvey_factor*eye(OPT.x_nbr);
    log_likelihood = kalman_gaussian(y, mu_0,Sigma_0, G,R,F, mu_eta,Sigma_eta, mu_eps,Sigma_eps, false, true);
elseif strcmp(kf_variant{1},"pruned_skewed")
    % initialize Kalman filter with a wide Normal prior: x_0 ~ CSN(0,Harvey_factor*eye(nx),0,0,eye(nx)) = N(0,10*eye(nx))
    mu_0     = zeros(OPT.x_nbr,1);
    Sigma_0  = OPT.Harvey_factor*eye(OPT.x_nbr);
    Gamma_0  = zeros(OPT.x_nbr,OPT.x_nbr);
    nu_0     = zeros(OPT.x_nbr,1);
    Delta_0  = eye(OPT.x_nbr);
    log_likelihood = kalman_csn(y, mu_0,Sigma_0,Gamma_0,nu_0,Delta_0, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, OPT.cdfmvna_fct,kf_variant{2},false,true);
end

if isnan(log_likelihood) || isinf(log_likelihood) || ~isreal(log_likelihood)
    log_likelihood = verysmall;
    exit_flag = 0;
end

end % main function end