function [log_likelihood_val, filt_param_last, xfilt] = kalman_csn( ...
    G_mat, R_mat, F_mat, ...
    mu_t_t, Sigma_t_t, Gamma_t_t, nu_t_t, Delta_t_t, ...
    mu_eps, Sigma_eps, ...
    mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, ...
    data, eval_lik, cut_tol, cdf_eval_fnc, varargin ...
)

    % Briefly: 
    %       - More robust (also numerically) skewed Kalman filter recursions, ...
    %         log-likelihood evaluation
    %       - This skewed Kalman filter is tailored for DSGE models
    %
    % Structure:
    %       - State Space model is of the form:
    %               Yt = F_mat*X_t + eps_t        ,  eps_t N(mu_eps, Sigma_eps)
    %               Xt = G_mat*X_t-1 + R_mat*eta_t,  eta_t CSN(mu_eta, Sigma_eta, ...
    %                                                          Gamma_eta, ...
    %                                                          nu_eta, Delta_eta)
    % 
    % Inputs:
    %           G_mat       : Transition matrix in State equation
    %           F_mat       : Matrix in Observation equation
    %           ---
    %           mu_eps      : Mean of the Measurement error
    %           Sigma_eps   : Var-Cov matrix of the Measurement error
    %           ---
    %           R_mat       : Matrix multiplying structural shocks in State equation
    %           ---
    %           mu_eta      : Mean of the structural shock
    %           Sigma_eta   : Var-Cov matrix of the structural shock
    %           Gamma_eta   : Shape matrix which regulates skewness
    %           Delta_eta   : Last parameter of CSN distribution of structural shock, ...
    %                         usually identity matrix
    %           ---
    %           data        : Observations
    %           eval_lik    : Should the likelihood be evaluated, boolean value
    %           cdf_eval_fnc: Handle for the function with which multivariate normal cdf ...
    %                         is calculated or the string of the function_name ...
    %                         e.g. 'mvncdf' or @mvncdf
    %           varargin    : Additional inputs to "cdf_eval_fnc"
    %       
    % Outputs:
    %           log_likelihood_val: Sum of log-likelihood contributions
    %           filt_param_last   : Parameters of the last filtered values
    %           xfilt             : Mean of the filtered values
    %
    %
    % Literature:
    %       1. [HO.JO.2013]   : Horn, Johnson (2013) - ...
    %                           Matrix analysis -- Theorem 2.6.3 (page 150)
    %       2. [GO.DO.GU.2004]: Gonzales-Farias, Dominguez-Molina, Gupta (2004) - ...
    %                           Additive properties of skew normal random vectors -- ...
    %                           Theorem 2 (page 527)
    %
    % Moreover:
    %       - Codes by G.Guljanov 05.05.2022
    %       - All the vectors should be column vectors
    %       - Y should be as usual. Rows are observations(realizations) and ...
    %         columns are variables
    
    nn = length(G_mat);
    [TT, y_nbr] = size(data);

    if nargout > 1
        filt_param_last.mu = nan;
        filt_param_last.Sigma = nan;
        filt_param_last.Gamma = nan;
        filt_param_last.nu = nan;
        filt_param_last.Delta = nan;
    end

    if nargout > 2
        xfilt = nan(nn, TT);
    end
    
    likeli_contr = zeros(1, TT); %Likelihood contributions

    const2pi = -0.5 * y_nbr * log(2 * pi);

    
    try
        chol(Sigma_eta);
    catch
        warning("Sigma_eta is not positive definite")
        log_likelihood_val = -Inf;
        return
    end
    
    % Multiplying matrix of joint distribution
    G_mat_bar = [G_mat, R_mat];

    for tt = 1:TT
        
            %Prediction
                
        % Parameters of joint distribution, [X_t', eta_t']'
        mu_bar = [mu_t_t; mu_eta];
        Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
        Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
        nu_bar = [nu_t_t; nu_eta];
        Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
        
        % Linear transformation of the joint distribution
        [mu_t_tm1, Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = lintransf(...
            G_mat_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar ...
        );


        % Pruning / Cutting algorithm
        %[Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = dim_red6(Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1, cut_tol);
        [Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = csnPruneParams(Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1, cut_tol);


        prediction_error = data(tt, :)' - (F_mat * mu_t_tm1 + mu_eps);
        Sigma_data = F_mat * Sigma_t_tm1 * F_mat' + Sigma_eps;
        K_Gauss = Sigma_t_tm1 * F_mat' / Sigma_data;
        K_Skewed = Gamma_t_tm1 * K_Gauss;
        

            %Filtering
            
        Tmp_mat = K_Gauss * prediction_error;
        
        mu_t_t  = mu_t_tm1 + Tmp_mat;
        Sigma_t_t = Sigma_t_tm1 - K_Gauss * F_mat * Sigma_t_tm1;
        Gamma_t_t = Gamma_t_tm1;
        nu_t_t = nu_t_tm1 - Gamma_t_tm1 * Tmp_mat;
        Delta_t_t = Delta_t_tm1;
        
        % Evaluate the filtered values of latent states
        if nargout > 2
            xfilt(:, tt) = csnMean(mu_t_t, Sigma_t_t, Gamma_t_t, nu_t_t, Delta_t_t);
        end
        
        % Filtering parameters of the last step
        if tt == TT
            filt_param_last.mu = mu_t_t;
            filt_param_last.Sigma = Sigma_t_t;
            filt_param_last.Gamma = Gamma_t_t;
            filt_param_last.nu = nu_t_t;
            filt_param_last.Delta = Delta_t_t;
        end
        

            %Evaluate the likelihood

        if ~eval_lik
            continue
        end
        

        % Parameters of (Y_t|Y_t-1)
        Tmp_leval = Gamma_t_tm1 * Sigma_t_tm1;
        
        %mu_data is used in prediction error
        %Sigma_data is calculated above
        %Gamma_data is equal to K_Skewed
        %nu_data is equal to nu_t_tm1
        Delta_data = ( ...
            Delta_t_tm1 + Tmp_leval * Gamma_t_tm1' - K_Skewed * F_mat * Tmp_leval' ...
        );

        Sigma_data = 0.5 * (Sigma_data + Sigma_data');
        Delta_data = 0.5 * (Delta_data + Delta_data');

        
        %Evaluate CDF -- denominator
        Cov_mat_denom = Delta_data + K_Skewed * Sigma_data * K_Skewed';
        Cov_mat_denom = 0.5 * (Cov_mat_denom + Cov_mat_denom');
        
        covtocorr_mat = diag(1 ./ sqrt(diag(Cov_mat_denom)));

        Corr_mat_denom = covtocorr_mat * Cov_mat_denom * covtocorr_mat;
        Corr_mat_denom = 0.5 * (Corr_mat_denom + Corr_mat_denom');
        
        eval_point_denom = -covtocorr_mat * nu_t_tm1;
        try
            cdf_val_denom = feval( ...
                cdf_eval_fnc, eval_point_denom, Corr_mat_denom, varargin{:} ...
            );
        catch
            warning('kalman_can_dsge: cdf_val_denom something wrong')
            log_likelihood_val = -Inf;
            return
        end
        %Evaluate PDF
        %term9 = log(mvnpdf(data(tt, :), mu_data', Sigma_data));
        pdf_val = ( ...
            const2pi - 0.5 * log(det(Sigma_data)) ...
            - 0.5 * prediction_error' / Sigma_data * prediction_error ...
        );
        
        %Evaluate CDF -- numerator
        covtocorr_mat = diag(1 ./ sqrt(diag(Delta_data)));

        eval_point_num = covtocorr_mat * (K_Skewed * prediction_error - nu_t_tm1);
        
        Corr_mat_num = covtocorr_mat * Delta_data * covtocorr_mat;
        Corr_mat_num = 0.5 * (Corr_mat_num + Corr_mat_num');
        try
            cdf_val_num = feval(cdf_eval_fnc, eval_point_num, Corr_mat_num, varargin{:});
        catch
            warning('kalman_can_dsge: cdf_val_num something wrong')
            log_likelihood_val = -Inf;
            return
        end

        % Likelihood contribution
        likeli_contr(tt) = -cdf_val_denom + pdf_val + cdf_val_num;
        
        if isnan(likeli_contr(tt))
            warning('Likelihood value is NaN')
            log_likelihood_val = -Inf;
            return
        end
            
    end
    
    % Log-likelihood value is sum of likelihood contributions
    log_likelihood_val = sum(likeli_contr);

function res_mat = blkdiag_two(mat1, mat2)

    % Briefly: Makes a block diagonal matrix out of two matrices
    % 
    % Structure:
    % 
    % Inputs:
    %       mat1: First matrix
    %       mat2: Second matrix
    % 
    % Ouputs:
    %       res_mat: Resulting block diagonal matrix
    % 
    % Literature:
    % 
    % Moreover:
    %

    [nrow_mat1, ncol_mat1] = size(mat1);
    [nrow_mat2, ncol_mat2] = size(mat2);

    upper_mat = zeros(nrow_mat1, ncol_mat2);
    lower_mat = zeros(nrow_mat2, ncol_mat1);

    res_mat = [mat1, ...
               upper_mat; ...
               lower_mat, ...
               mat2];

end % blkdiag_two


end % kalman_csn_dsge




