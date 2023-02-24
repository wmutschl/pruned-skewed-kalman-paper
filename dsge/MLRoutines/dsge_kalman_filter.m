function [log_lik_t_tm1] = dsge_kalman_filter(d,dSS,H,gx,gu,SIGu,SIGe)

%==========================================================================
% This Function Implements the Kalman filter for the state space model:
%   d_{t} = dSS + H*y_{t} + e_{t} with e_{t} ~ N(0,SIGe) [Observation Equation]
%   y_{t} = gx*y_{t-1} + gu*u_{t}   with u_{t} ~ N(0,SIGu) [Transition Equation]
% 
% In the above system:
% - t=1,...,nobs measures the discrete time period,
% - d_{t} is a (nd x 1) vector of Gaussian observable variables (called controls)
% - y_{t} is a (ny x 1) vector of Gaussian latent variables (called states)
% - u_{t} is a (nu x 1) vector of Gaussian structural shocks (called innovations)
% - e_{t} is a (nd x 1) vector of Gaussian measurement errors (called noise)
%
% Note that in our DSGE model framework y_{t} corresponds to the vector of 
% all endogenous variables minus the corresponding steady-state. H is then
% simply a selection matrix which selects the variables that are
% observable; hence, we need to add the steady-state dSS
%
%==========================================================================
% INPUTS:
% - d:    (nobs x nd) matrix of observations for d(t)
% - dSS:  (nd x 1)    vector of steady-state values for the observables
% - H:    (nd x ny)   selection matrix that picks the observable variables from y
% - gx:   (ny x ny)   solution matrix with respect to states
% - gu:   (ny x nu)   solution matrix with respect to shocks
% - SIGu: (nu x nu)   covariance matrix of shocks
% - SIGe: (nd x nd)   covariance matrix of measurement errors
%  
% OUTPUTS:
% - log_lik_t_tm1 = is a (nobsx1) vector containing log(p(d_{t}|d_{t-1},...,d_{0})).
%   The first entry is based on the prediction of the state vector at its unconditional mean;
%
%
% Matrices in Kalman filter:
% - yhat_t_tm1:  forecast of y_{t} given d^{t-1}
% - yhat_tp1_t:  forecast of y_{t+1} given d^{t}
% - Sigma_t_tm1: mean-squared-error of y_t given d^{t-1}    
% - Sigma_tp1_t: mean-squared-error of y_{t+1} given d^{t}
% - K: Kalman gain
%==========================================================================

%% Get dimensions
[nobs,nd] = size(d);
ny        = size(gu,1);

%% Initialize the state vector at the stationary distribution
yhat_t_tm1 = zeros(ny,1); % note that y are the model variables in deviation from steady-state, so the mean is zero by definition
%Sigma_t_tm1 = reshape( inv(eye(ny*ny) - kron(gx,gx))*reshape(gu*SIGu*gu',ny*ny,1) ,ny,ny); %analytical, but slow
Sigma_t_tm1 = dlyapdoubling(gx,gu*SIGu*gu'); % very fast and numerically accurate
log_lik_t_tm1 = nan(nobs,1);
%% Kalman Filter Recursion
for t=1:nobs
    % Step 1: Compute Kalman Gain
    Omega = H*Sigma_t_tm1*H'+SIGe; % gain matrix
    det_Omega = det(Omega);
    if det_Omega<=0
        log_lik_t_tm1(t)=-10^8;
        return
    else
        K = gx*Sigma_t_tm1*H'/Omega; % Kalman gain
    end

    % Step 2: Compute the forecast error in the observations
    a = d(t,:)' - (dSS + H*yhat_t_tm1);
    
    % Step 3: Compute the state forecast for next period given today's information
    yhat_tp1_t = gx*yhat_t_tm1 + K*a;

    % Step 4: Update the covariance matrix
    Sigma_tp1_t = (gx-K*H)*Sigma_t_tm1*(gx'-H'*K') + gu*SIGu*gu' + K*SIGe*K';

    % contribution to log-likelihood
    log_lik_t_tm1(t) = -nd/2*log(2*pi) - 0.5*reallog(det_Omega) - 0.5*((a'/Omega*a)); % reallog produces errors for negative and complex numbers

    % reset values for next step
    yhat_t_tm1 = yhat_tp1_t;
    Sigma_t_tm1 = Sigma_tp1_t;    
end

end % main function end
