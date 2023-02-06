function [log_lik,x_filt,pred,filt] = kalman_csn(Y, mu_tm1_tm1,Sigma_tm1_tm1,Gamma_tm1_tm1,nu_tm1_tm1,Delta_tm1_tm1, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, cdfmvna_fct,prune_tol,skip_lik,skip_loss,loss_fct)
% [log_lik,x_filt,pred,filt] = kalman_csn(Y, mu_tm1_tm1,Sigma_tm1_tm1,Gamma_tm1_tm1,nu_tm1_tm1,Delta_tm1_tm1, G,R,F, mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta, mu_eps,Sigma_eps, cdfmvna_fct,prune_tol,skip_lik,skip_loss,loss_fct)
% -------------------------------------------------------------------------
% Evaluate (1) log-likelihood value, (2) filtered states and (3) smoothed states
% of linear state space model with csn distributed innovations and normally distributed noise:
%   x(t) = G*x(t-1) + R*eta(t)   [state transition equation]
%   y(t) = F*x(t)   + eps(t)     [observation equation]
%   eta(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta) [innovations, shocks]
%   eps(t) ~ N(mu_eps,Sigma_eps) [noise, measurement error]
% Dimensions:
%   x(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   eta(t) is (eta_nbr by 1) vector of innovations
%   eps(t) is (y_nbr by 1) vector of noise (measurement errors)
% -------------------------------------------------------------------------
% INPUTS
% - Y               [y_nbr by obs_nbr]             matrix with data
% - mu_tm1_tm1      [x_nbr by 1]                   initial value of 1st parameter of CSN distributed states x (location). Does not equal expectation vector unless Gamma_tm1_tm1=0.
% - Sigma_tm1_tm1   [x_nbr by x_nbr]               initial value of 2nd parameter of CSN distributed states x (scale). Does not equal covariance matrix unless Gamma_tm1_tm1=0.
% - Gamma_tm1_tm1   [skewx_dim by x_nbr]           initial value of 3rd parameter of CSN distributed states x (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu_tm1_tm1      [skewx_dim by x_nbr]           initial value of 4th parameter of CSN distributed states x (enables closure of CSN distribution under conditioning)
% - Delta_tm1_tm1   [skewx_dim by skewx_dim]       initial value of 5th parameter of CSN distributed states x (enables closure of CSN distribution und marginalization)
% - G               [x_nbr by x_nbr]               state transition matrix mapping previous states to current states
% - R               [x_nbr by eta_nbr]             state transition matrix mapping current innovations to current states
% - F               [y_nbr by x_nbr]               observation equation matrix mapping current states into current observables
% - mu_eta          [eta_nbr by 1]                 1st parameter of CSN distributed innovations eta (location). Does not equal expectation vector unless Gamma_eta=0.
% - Sigma_eta       [eta_nbr by eta_nbr]           2nd parameter of CSN distributed innovations eta (scale). Does not equal covariance matrix unless Gamma_eta=0.
% - Gamma_eta       [skeweta_dim by eta_nbr]       3rd parameter of CSN distributed innovations eta (regulates skewness from Gaussian (Gamma=0) to half normal)
% - nu_eta          [skeweta_dim by 1]             4th parameter of CSN distributed innovations eta (enables closure of CSN distribution under conditioning)
% - Delta_eta       [skeweta_dim by skeweta_dim]   5th parameter of CSN distributed innovations eta (enables closure of CSN distribution und marginalization)
% - mu_eps          [y_nbr by 1]                   location parameter of normally distributed measurement errors eps (equals expectation vector)
% - Sigma_eps       [y_nbr by y_nbr]               scale parameter of normally distributed measurement errors eps (equals covariance matrix)
% - cdfmvna_fct     [string]                       name of function to compute log Gaussian cdf, possible values: 'logmvncdf_ME', 'mvncdf', 'qsilatmvnv', 'qsimvnv'
% - prune_tol       [double]                       correlation threshold to prune redundant skewness dimensions, if set to 0 no pruning will be done
% - skip_lik        [boolean]                      1: skip log-likelihood computations (e.g. for doing filtering and smoothing only)
% - skip_loss       [boolean]                      1: skip filtering and smoothing computations (e.g. for doing log-likelihood only)
% - loss_fct        [structure]                    Underlying loss function for computing point estimate of filtered states x_t_t that minimizes the expected loss.
%                                                  loss_fct.type possible values:
%                                                  - "L1"   [boolean]   xtilde = median(x_t_t), i.e. absolute loss abs(xtilde-x)
%                                                  - "L2"   [boolean]   xtilde = E[x_t_t], i.e. squared loss (xtilde-x)^2
%                                                  - "La"   [boolean]   xtilde = quantile(x_t_t,a/(a+b)), i.e. asymmetric loss function a*abs(xtilde-x) for x>xtilde and b*abs(xtilde-x) for x<=xtilde, , a and b are taken from loss_fct.params.a and loss_fct.params.b
% -------------------------------------------------------------------------
% OUTPUTS
% - log_lik         [scalar]                       value of log likelihood
% - x_filt          [structure]                    filtered states according to different loss functions given by loss_fct.type
% - pred            [structure]                    CSN parameters of predicted states (mu_t_tm1,Sigma_t_tm1,Gamma_t_tm1,nu_t_tm1,Delta_t_tm1) needed e.g. for kalman_csn_smoother.m
% - filt            [structure]                    CSN parameters of filtered states (mu_t_t,Sigma_t_t,Gamma_t_t,nu_t_t,Delta_t_t) needed e.g. for kalman_csn_smoother.m
% =========================================================================
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================

% some settings (inspired by kalman_filter.m of Dynare)
kalman_tol = 1e-10;    % numerical tolerance for determining the singularity of the covariance matrix of the prediction errors during the Kalman filter (minimum allowed reciprocal of the matrix condition number)
Omega_singular = true; % initialize

% get dimensions
[y_nbr,x_nbr] = size(F);
obs_nbr       = size(Y,2);

% initialize some matrices
mu_eta    = R*mu_eta;
Sigma_eta = R*Sigma_eta*R';
Gamma_eta = Gamma_eta/(R'*R)*R';
Gamma_eta_X_Sigma_eta = Gamma_eta*Sigma_eta;
Delta22_common = Delta_eta + Gamma_eta_X_Sigma_eta*Gamma_eta';
const2pi = -0.5*y_nbr*log(2*pi);

if nargout > 1
    x_filt.L1 = nan(x_nbr,obs_nbr);
    x_filt.L2 = nan(x_nbr,obs_nbr);
    x_filt.La = nan(x_nbr,obs_nbr);
end
if nargout > 2
    % initialize "pred" structure to save parameters of predicted states
    pred.mu    = zeros(x_nbr,obs_nbr);
    pred.Sigma = zeros(x_nbr,x_nbr,obs_nbr);
    pred.Gamma = cell(obs_nbr,1); % use cell as skewness dimension is time-varying
    pred.nu    = cell(obs_nbr,1); % use cell as skewness dimension is time-varying
    pred.Delta = cell(obs_nbr,1); % use cell as skewness dimension is time-varying

    % initialize "filt" structure to save parameters of filtered states
    filt.mu    = zeros(x_nbr,obs_nbr);
    filt.Sigma = zeros(x_nbr,x_nbr,obs_nbr);
    filt.Gamma = cell(obs_nbr,1); % use cell as skewness dimension is time-varying
    filt.nu    = cell(obs_nbr,1); % use cell as skewness dimension is time-varying
    filt.Delta = cell(obs_nbr,1); % use cell as skewness dimension is time-varying
end
log_lik_t = zeros(obs_nbr,1); % initialize vector of likelihood contributions
log_lik   = -Inf; % default value of log likelihood

for t=1:obs_nbr
    % auxiliary matrices
    Gamma_tm1_tm1_X_Sigma_tm1_tm1      = Gamma_tm1_tm1*Sigma_tm1_tm1;
    Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT = Gamma_tm1_tm1_X_Sigma_tm1_tm1*G';

    % prediction
    mu_t_tm1  = G*mu_tm1_tm1 + mu_eta;
    Sigma_t_tm1 = G*Sigma_tm1_tm1*G' + Sigma_eta;
    Sigma_t_tm1 = 0.5*(Sigma_t_tm1 + Sigma_t_tm1'); % ensure symmetry
    invSigma_t_tm1 = inv(Sigma_t_tm1);
    Gamma_t_tm1 = [Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT; Gamma_eta_X_Sigma_eta]*invSigma_t_tm1;
    nu_t_tm1 = [nu_tm1_tm1; nu_eta];
    Delta11_t_tm1 = Delta_tm1_tm1 + Gamma_tm1_tm1_X_Sigma_tm1_tm1*Gamma_tm1_tm1' - Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT';
    Delta22_t_tm1 = Delta22_common - Gamma_eta_X_Sigma_eta*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
    Delta12_t_tm1 = -Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
    Delta_t_tm1 = [Delta11_t_tm1 , Delta12_t_tm1; Delta12_t_tm1' , Delta22_t_tm1];
    Delta_t_tm1 = 0.5*(Delta_t_tm1 + Delta_t_tm1'); % ensure symmetry
    y_predicted = F*mu_t_tm1 + mu_eps;
    prediction_error = Y(:,t) - y_predicted;

    % pruning redundant skewness dimensions to speed up filtering
    if prune_tol > 0
        [Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = csnPruneParams(Sigma_t_tm1,Gamma_t_tm1,nu_t_tm1,Delta_t_tm1,prune_tol);
    end

    % Kalman gains
    Omega = F*Sigma_t_tm1*F' + Sigma_eps;
    Omega = 0.5*(Omega + Omega'); %ensure symmetry
    badly_conditioned_Omega = false;
    if rcond(Omega)<kalman_tol
        sig=sqrt(diag(Omega));
        if any(diag(Omega)<kalman_tol) || rcond(Omega./(sig*sig'))<kalman_tol
            badly_conditioned_Omega = true;
            warning('kalman_csn: badly_conditioned_Omega')        
        end
    end
    if badly_conditioned_Omega
        if ~all(abs(Omega(:))<kalman_tol)
            % Use univariate filter (will remove observations with zero variance prediction error)
            error('kalman_csn: you should use an univariate filter, which is not in the replication codes')
        else
            % Pathological case, discard draw.
            warning('discard draw due to badly_conditioned_Omega')
            return
        end
    else
        Omega_singular = false;
        log_detOmega = log(det(Omega));
        invOmega = inv(Omega);
        K_Gauss = Sigma_t_tm1*F'*invOmega; % Gaussian Kalman Gain
        K_Skewed = Gamma_t_tm1*K_Gauss;    % Skewed Kalman Gain
    
        if ~skip_lik
            % log-likelihood contributions        
            % The conditional distribution of y(t) given y(t-1) is:
            % (y(t)|y(t-1)) ~ CSN(mu_y,Sigma_y,Gamma_y,nu_y,Delta_y)
            %               = mvncdf(Gamma_y*(y(t)-mu_y),nu_y,Delta_y) / mvncdf(0,nu_y,Delta_y+Gamma_y*Sigma_y*Gamma_y') * mvnpdf(y(t),mu_y,Sigma_y)
            % where:
            %   mu_y    = F*mu_t_tm1 + mu_eps = y_predicted
            %   Sigma_y = F*Sigma_t_tm1*F' + Sigma_eps = Omega
            %   Gamma_y = Gamma_t_tm1*Sigma_t_tm1*F'*inv(F*Sigma_t_tm1*F' + Sigma_eps) = K_Skewed
            %   nu_y    = nu_t_tm1
            %   Delta_y = Delta_t_tm1 + Gamma_t_tm1*Sigma_t_tm1*Gamma_t_tm1'...
            %             - Gamma_t_tm1*Sigma_t_tm1*F'*inv(F*Sigma_t_tm1*F')*F*Sigma_t_tm1*Gamma_t_tm1'...
            %             + (Gamma_t_tm1*Sigma_t_tm1*F'*inv(F*Sigma_t_tm1*F') - Gamma_t_tm1*Sigma_t_tm1*F'*inv(F*Sigma_t_tm1*F' + Sigma_eps))*F*Sigma_t_tm1*Gamma_t_tm1';
            %           = Delta_t_tm1 + (Gamma_t_tm1-K_Skewed*F)*Sigma_t_tm1*Gamma_t_tm1'
            Delta_y = Delta_t_tm1 + (Gamma_t_tm1-K_Skewed*F)*Sigma_t_tm1*Gamma_t_tm1';
            Delta_y = 0.5*(Delta_y + Delta_y'); %ensure symmetry

            % evaluate Gaussian cdfs, i.e.
            %  - bottom one: mvncdf(0,nu_y,Delta_y + Gamma_y*Sigma_y*Gamma_y')
            %  - top one: mvncdf(Gamma_y*(y(t)-mu_y),nu_y,Delta_y)
            cdf_bottom_cov = Delta_y + K_Skewed*Omega*K_Skewed';
            cdf_bottom_cov = 0.5*(cdf_bottom_cov + cdf_bottom_cov'); % ensure symmetry
            if strcmp(cdfmvna_fct,'logmvncdf_ME')
                % requires zero mean and correlation matrix as inputs
                normalization_bottom_cov = diag(1./sqrt(diag(cdf_bottom_cov)));
                cdf_bottom_cov = normalization_bottom_cov*cdf_bottom_cov*normalization_bottom_cov; % this is now a correlation matrix!
                cdf_bottom_cov = 0.5*(cdf_bottom_cov + cdf_bottom_cov'); % ensure symmetry
                if ~isempty(cdf_bottom_cov)
                    log_gaussian_cdf_bottom = logmvncdf_ME(-normalization_bottom_cov*nu_t_tm1, cdf_bottom_cov);
                else
                    log_gaussian_cdf_bottom = 0;
                end

                normalization_Delta_y = diag(1./sqrt(diag(Delta_y)));
                Delta_y = normalization_Delta_y*Delta_y*normalization_Delta_y; % this is now a correlation matrix!
                Delta_y = 0.5*(Delta_y + Delta_y');
                if ~isempty(Delta_y)
                    log_gaussian_cdf_top = logcdf_ME(normalization_Delta_y*(K_Skewed*prediction_error - nu_t_tm1), Delta_y);                    
                else
                    log_gaussian_cdf_top = 0;
                end
            elseif strcmp(cdfmvna_fct,'mvncdf')
                log_gaussian_cdf_bottom = log(mvncdf(zeros(size(nu_t_tm1,1),1), nu_t_tm1, cdf_bottom_cov));
                log_gaussian_cdf_top = log(mvncdf(K_Skewed*prediction_error, nu_t_tm1, Delta_y));
            elseif strcmp(cdfmvna_fct,'qsilatmvnv')
                k = size(nu_t_tm1,1);
                if k > 1
                    log_gaussian_cdf_bottom = log(qsilatmvnv( 1000*k, cdf_bottom_cov, repmat(-Inf,k,1)-nu_t_tm1, zeros(k,1)-nu_t_tm1 ));
                    log_gaussian_cdf_top = log(qsilatmvnv( 1000*k, Delta_y, repmat(-Inf,k,1)-nu_t_tm1, K_Skewed*prediction_error-nu_t_tm1 ));
                else
                    log_gaussian_cdf_bottom = log(normcdf(0, nu_t_tm1, cdf_bottom_cov));
                    log_gaussian_cdf_top = log(normcdf(K_Skewed*prediction_error, nu_t_tm1, Delta_y));
                end
            elseif strcmp(cdfmvna_fct,'qsimvnv')
                k = size(nu_t_tm1,1);
                if k > 1
                    log_gaussian_cdf_bottom = log(qsimvnv( 1000*k, cdf_bottom_cov, repmat(-Inf,k,1)-nu_t_tm1, zeros(k,1)-nu_t_tm1 ));
                    log_gaussian_cdf_top = log(qsimvnv( 1000*k, Delta_y, repmat(-Inf,k,1)-nu_t_tm1, K_Skewed*prediction_error-nu_t_tm1 ));
                else
                    log_gaussian_cdf_bottom = log(normcdf(0, nu_t_tm1, cdf_bottom_cov));
                    log_gaussian_cdf_top = log(normcdf(K_Skewed*prediction_error, nu_t_tm1, Delta_y));
                end
            end

            % evaluate Gaussian pdf
            % log_gaussian_pdf = log(mvnpdf(Y(:,t), y_predicted, Omega));
            log_gaussian_pdf = const2pi - 0.5*log_detOmega - 0.5*transpose(prediction_error)*invOmega*prediction_error;

            log_lik_t(t) = -log_gaussian_cdf_bottom + log_gaussian_pdf + log_gaussian_cdf_top;
            if isnan(log_lik_t(t))
                log_lik = -Inf;
                x_filt  = nan;
                return
            end
        end

        % filtering
        mu_t_t    = mu_t_tm1 + K_Gauss*prediction_error;
        Sigma_t_t = Sigma_t_tm1 - K_Gauss*F*Sigma_t_tm1;
        Gamma_t_t = Gamma_t_tm1;
        nu_t_t    = nu_t_tm1 - K_Skewed*prediction_error;
        Delta_t_t = Delta_t_tm1;
        Sigma_t_t = 0.5*(Sigma_t_t + Sigma_t_t'); % ensure symmetry
        Delta_t_t = 0.5*(Delta_t_t + Delta_t_t'); % ensure symmetry
        if nargout > 1 && ~skip_loss % save the point estimate that minimizes loss of x_t_t
            if loss_fct.type.L1
                % absolue loss, i.e. compute median of CSN distributed x_t_t
                x_filt.L1(:,t) = csnQuantile(0.5, mu_t_t, Sigma_t_t, Gamma_t_t, nu_t_t, Delta_t_t, cdfmvna_fct);
            end
            if loss_fct.type.L2
                % squared loss, i.e. compute mean of CSN distributed x_t_t
                x_filt.L2(:,t) = csnMean(mu_t_t, Sigma_t_t, Gamma_t_t, nu_t_t, Delta_t_t, cdfmvna_fct);
            end
            if loss_fct.type.La
                % asymmetric loss, i.e. compute a/(a+b) quantile of CSN distributed x_t_t
                x_filt.La(:,t) = csnQuantile(loss_fct.params.a/(loss_fct.params.a+loss_fct.params.b), mu_t_t, Sigma_t_t, Gamma_t_t, nu_t_t, Delta_t_t, cdfmvna_fct);
            end
        end

        % assign for next time step
        mu_tm1_tm1    = mu_t_t;
        Sigma_tm1_tm1 = Sigma_t_t;
        Gamma_tm1_tm1 = Gamma_t_t;
        nu_tm1_tm1    = nu_t_t;
        Delta_tm1_tm1 = Delta_t_t;

        if nargout > 2
            % save the parameters of the predicted and filtered csn states
            pred.mu(:,t)      = mu_t_tm1;
            pred.Sigma(:,:,t) = Sigma_t_tm1;
            pred.Gamma{t,1}   = Gamma_t_tm1;
            pred.nu{t,1}      = nu_t_tm1;
            pred.Delta{t,1}   = Delta_t_tm1;
            filt.mu(:,t)      = mu_t_t;
            filt.Sigma(:,:,t) = Sigma_t_t;
            filt.Gamma{t,1}   = Gamma_t_t;
            filt.nu{t,1}      = nu_t_t;
            filt.Delta{t,1}   = Delta_t_t;
        end
    end
end

if Omega_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% compute log-likelihood by summing individual contributions
log_lik = sum(log_lik_t);

end % main function end