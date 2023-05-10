function [log_lik, x_filt, x_smooth] = kalman_gaussian(Y, mu_tm1_tm1,Sigma_tm1_tm1, G,R,F, E_eta,V_eta, mu_eps,Sigma_eps, skip_lik,skip_loss,loss_fct)
% [log_lik, x_filt, x_smooth] = kalman_gaussian(Y, mu_tm1_tm1,Sigma_tm1_tm1, G,R,F, E_eta,V_eta, mu_eps,Sigma_eps, skip_lik,skip_loss,loss_fct)
% -------------------------------------------------------------------------
% Evaluate (1) log-likelihood value, (2) filtered states and (3) smoothed states
% of linear state space model with normally distributed innovations and normally distributed noise:
%   x(t) = G*x(t-1) + R*eta(t)    [state transition equation]
%   y(t) = F*x(t)   + eps(t)      [observation equation]
%   eta(t) ~ N(E_eta,V_eta)       [innovations, shocks]
%   eps(t) ~ N(mu_eps,Sigma_eps)  [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% -------------------------------------------------------------------------
% INPUTS
% - Y               [y_nbr by obs_nbr]             matrix with data
% - mu_tm1_tm1      [x_nbr by 1]                   initial value of location parameter of normally distributed states x (equals expectation vector)
% - Sigma_tm1_tm1   [x_nbr by x_nbr]               initial value of scale parameter of normally distributed states x (equals covariance matrix)
% - G               [x_nbr by x_nbr]               state transition matrix mapping previous states to current states
% - R               [x_nbr by eta_nbr]             state transition matrix mapping current innovations to current states
% - F               [y_nbr by x_nbr]               observation equation matrix mapping current states into current observables
% - E_eta           [eta_nbr by 1]                 location parameter of normally distributed innovations eta (equals expectation vector)
% - V_eta           [eta_nbr by eta_nbr]           scale parameter of CSN distributed innovations eta (equals covariance matrix
% - mu_eps          [y_nbr by 1]                   location parameter of normally distributed measurement errors eps (equals expectation vector)
% - Sigma_eps       [y_nbr by y_nbr]               scale parameter of normally distributed measurement errors eps (equals covariance matrix)
% - skip_lik        [boolean]                      1: skip log-likelihood computations (e.g. for doing filtering and smoothing only)
% - skip_loss       [boolean]                      1: skip filtering and smoothing computations (e.g. for doing log-likelihood only)
% - loss_fct        [structure]                    Underlying loss function for computing point estimate of filtered states x_t_t that minimizes the expected loss.
%                                                  loss_fct.type possible values:
%                                                  - "L1"   [boolean]   xtilde = median(x_t_t), i.e. absolute loss abs(xtilde-x)
%                                                  - "L2"   [boolean]   xtilde = E[x_t_t], i.e. squared loss (xtilde-x)^2
%                                                  - "La"   [boolean]   xtilde = quantile(x_t_t,a/(a+b)), i.e. asymmetric loss function a*abs(xtilde-x) for x>xtilde and b*abs(xtilde-x) for x<=xtilde, , a and b are taken from loss_fct.params.a and loss_fct.params.b
% -------------------------------------------------------------------------
% OUTPUTS
% - log_lik         [double]                       value of log likelihood
% - x_filt          [structure]                    filtered states according to different loss functions given by loss_fct.type
% - x_smooth        [structure]                    smoothed states according to different loss functions given by loss_fct.type
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler
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
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

% some settings (inspired by kalman_filter.m of Dynare)
kalman_tol     = 1e-10; % numerical tolerance for determining the singularity of the covariance matrix of the prediction errors during the Kalman filter (minimum allowed reciprocal of the matrix condition number)
Omega_singular = true;  % initialize

% get dimensions
[y_nbr,x_nbr] = size(F);
obs_nbr       = size(Y,2);

% initialize some matrices
E_eta    = R*E_eta;
V_eta    = R*V_eta*R';
const2pi = -0.5*y_nbr*log(2*pi);

if nargout > 1
    x_filt.L1 = nan(x_nbr,obs_nbr);
    x_filt.L2 = nan(x_nbr,obs_nbr);
    x_filt.La = nan(x_nbr,obs_nbr);
end
if nargout > 2
    % initialize "_pred" values to save parameters of predicted states
    mu_pred     = zeros(x_nbr,obs_nbr);
    Sigma_pred  = zeros(x_nbr, x_nbr, obs_nbr);
    % initialize "_filt" values to save parameters of predicted states
    mu_filt     = zeros(x_nbr,obs_nbr);
    Sigma_filt  = zeros(x_nbr,x_nbr,obs_nbr);
    % initialize structure for smoothed states
    x_smooth.L1 = nan(x_nbr,obs_nbr);
    x_smooth.L2 = nan(x_nbr,obs_nbr);
    x_smooth.La = nan(x_nbr,obs_nbr);    
end
log_lik_t = zeros(obs_nbr,1); % initialize vector of likelihood contributions
log_lik   = -Inf; % default value of log likelihood

for t=1:obs_nbr
    % prediction
    mu_t_tm1    = G*mu_tm1_tm1 + E_eta;
    Sigma_t_tm1 = G*Sigma_tm1_tm1*G' + V_eta;
    Sigma_t_tm1 = 0.5*(Sigma_t_tm1 + Sigma_t_tm1'); % ensure symmetry
    y_predicted = F*mu_t_tm1 + mu_eps;
    prediction_error = Y(:,t) - y_predicted;

    % Kalman gain
    Omega = F*Sigma_t_tm1*F' + Sigma_eps;
    Omega = 0.5*(Omega + Omega'); %ensure symmetry
    badly_conditioned_Omega = false;
    if rcond(Omega)<kalman_tol
        sig=sqrt(diag(Omega));
        if any(diag(Omega)<kalman_tol) || rcond(Omega./(sig*sig'))<kalman_tol
            badly_conditioned_Omega = true;
            warning('kalman_gaussian: badly_conditioned_Omega')
        end
    end
    if badly_conditioned_Omega
        if ~all(abs(Omega(:))<kalman_tol)
            % Use univariate filter (will remove observations with zero variance prediction error)
            error('kalman_gaussian: you should use an univariate filter, which is not in the replication codes')
        else
            % Pathological case, discard draw.
            warning('kalman_gaussian: discard draw due to badly_conditioned_Omega')
            return
        end
    else
        Omega_singular = false;
        log_detOmega = log(det(Omega));
        invOmega = inv(Omega);
        K_Gauss = Sigma_t_tm1*F'*invOmega;

        if ~skip_lik
            % log-likelihood contributions        
            % The conditional distribution of y(t) given y(t-1) is:
            % (y(t)|y(t-1)) ~ N(mu_y,Sigma_y) = mvnpdf(y(t),mu_y,Sigma_y)
            % where:
            %   mu_y    = F*mu_t_tm1 + mu_eps = y_predicted
            %   Sigma_y = F*Sigma_t_tm1*F' + Sigma_eps = Omega

            % evaluate Gaussian pdf
            log_lik_t(t) = const2pi - 0.5*log_detOmega - 0.5*transpose(prediction_error)*invOmega*prediction_error;
            % log_lik_t(t) = log(mvnpdf(Y(:,t), y_predicted, Omega));
            if isnan(log_lik_t(t))
                % penalize likelihood
                log_lik = -Inf;
                x_filt = nan;
                x_smooth = nan;
                return
            end
        end

        % filtering
        mu_t_t = mu_t_tm1 + K_Gauss*prediction_error;
        Sigma_t_t = Sigma_t_tm1 - K_Gauss*F*Sigma_t_tm1;
        if nargout > 1 && ~skip_loss % save the point estimate that minimizes loss of x_t_t
            if loss_fct.type.L1
                % absolue loss, i.e. compute median of Gaussian distributed x_t_t
                x_filt.L1(:,t) = mu_t_t;
            end
            if loss_fct.type.L2
                % squared loss, i.e. compute mean of Gaussian distributed x_t_t
                x_filt.L2(:,t) = mu_t_t;
            end
            if loss_fct.type.La
                % asymmetric loss, i.e. compute a/(a+b) and b/(a+b) quantiles of Gaussian distributed x_t_t
                x_filt.La(:,t) = gaussianQuantile(loss_fct.params.a/(loss_fct.params.a+loss_fct.params.b), mu_t_t, Sigma_t_t);
            end
        end

        % assign for next time step
        mu_tm1_tm1    = mu_t_t;
        Sigma_tm1_tm1 = Sigma_t_t;

        if nargout > 2
            % save the parameters of the predicted and filtered csn states for smoothing
            mu_pred(:,t)      = mu_t_tm1;
            Sigma_pred(:,:,t) = Sigma_t_tm1;
            mu_filt(:,t)      = mu_t_t;
            Sigma_filt(:,:,t) = Sigma_t_t;
        end
    end
end

if Omega_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% compute log-likelihood by summing individual contributions
log_lik = sum(log_lik_t);

if nargout > 2
    for t = obs_nbr:-1:1
        if t == obs_nbr
            % smoothing step for the last time point, i.e. 'T'
            mu_smooth(:,obs_nbr)      = mu_filt(:,obs_nbr);
            Sigma_smooth(:,:,obs_nbr) = Sigma_filt(:,:,obs_nbr);
        else
            Jt = Sigma_filt(:,:,t)*G'/Sigma_pred(:,:,t+1);
            % first two distribution parameters of smoothed value
            mu_smooth(:,t)      = mu_filt(:,t) + Jt*(mu_smooth(:,t+1)-mu_pred(:, t+1));
            Sigma_smooth(:,:,t) = Sigma_filt(:,:,t) + Jt*(Sigma_smooth(:,:,t+1) - Sigma_pred(:,:,t+1))*Jt';
        end
        if loss_fct.type.L1 && ~skip_loss
            % absolue loss, i.e. compute median of normally distributed x_t_t
            x_smooth.L1(:,t) = mu_smooth(:, t);
        end
        if loss_fct.type.L2 && ~skip_loss
            % squared loss, i.e. compute mean of normally distributed x_t_t
            x_smooth.L2(:,t) = mu_smooth(:, t);
        end
        if loss_fct.type.La && ~skip_loss
            % asymmetric loss, i.e. compute a/(a+b) quantile of normally distributed x_t_t
            x_smooth.La(:,t) = gaussianQuantile(loss_fct.params.a/(loss_fct.params.a+loss_fct.params.b), mu_smooth(:, t), Sigma_smooth(:, :, t));
        end
    end
end

end % kalman_gaussian