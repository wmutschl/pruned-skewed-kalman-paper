function [Gamma, nu, Delta] = dim_red6(Sigma, Gamma, nu, Delta, cut_tol)
    
    % Briefly: 
    %       - Reduces the dimension of csn according to the correlations of the conditions
    %
    % Structure:
    %       - Get the lower left block of Covariance matrix of (W', Z')' random vector
    %       - Compute correlations out of the above computed matrix of covariances
    %       - Take absolute value of each element of above computed matrix of correlations
    %       - Look at the maximum values along the 2nd dimension
    %       - If max. values are below some threshold, delete corresponding dimension
    %       - Return newly cut matrices
    % 
    % Inputs:
    %       Sigma: 2nd parameter of the csn distributed random vector, a matrix of n x n
    %       Gamma: 3nd parameter of the csn distributed random vector, a matrix of m x n
    %       nu   : 4th parameter of the csn distributed random vector, a matrix of m x 1
    %       Delta: 5th parameter of the csn distributed random vector, a matrix of m x m
    %       
    % Outputs:
    %       Gamma: 3nd parameter of the csn distributed random vector, a matrix of m x n
    %       nu   : 4th parameter of the csn distributed random vector, a matrix of m x 1
    %       Delta: 5th parameter of the csn distributed random vector, a matrix of m x m
    %
    %
    % Literature:
    %       - Look at the pseudocode of our first paper about skewed Kalman filter
    %
    % Moreover:
    %       - Codes by G.Guljanov
    
    nn = length(Sigma);

    % Compute lower left block of the Correlation matrix
    Cov_mat_21 = Gamma * Sigma;

    diag_Sigma = sqrt(diag(Sigma));
    diag_aux_variance = sqrt(diag(Delta + Cov_mat_21 * Gamma'));

    Corr_mat_21 = abs(Cov_mat_21 ./ (diag_aux_variance * diag_Sigma'));
    
    % Hold-vector whose elements are "true" where the dimensions should be held
    max_vec = max(Corr_mat_21, [], 2);
    hold_vec = (max_vec > cut_tol);

    max_dim = 15;
    if sum(hold_vec) > max_dim
        hold_vec(:) = false;
        [~, max_ind] = maxk(max_vec, max_dim);
        hold_vec(max_ind) = true;
    end
    
    %Cut unnecessary dimensions
    Gamma = Gamma(hold_vec, :);
    nu = nu(hold_vec);
    Delta = Delta(hold_vec, hold_vec);

    Delta = 0.5 * (Delta + Delta'); %Ensure symmetry

    % All skewness dimensions are cut --> Gamma_t_tm1 as a zero row vector, i.e. gaussianity
    if isempty(Delta)
        Gamma = zeros(1, nn);
        nu = 0;
        Delta = 1;
    end

end