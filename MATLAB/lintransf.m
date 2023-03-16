function [ ...
    mu_A, Sigma_A, Gamma_A, nu_A, Delta_A ...
] = lintransf( ...
    A_mat, mu, Sigma, Gamma, nu, Delta ...
)
    
    % Briefly: 
    %       - Linearly transforms closed skew normally distributed random variables
    %
    % Structure:
    %           - 1) Make proper csn first, in case of singularity
    %           - 2) Either rank deficient or full rank linear transformation
    % 
    % Inputs:
    %           A_mat: Matrix multiplying the joint distribution of X_t and eta_t
    %           ---
    %           mu   : 1st parameter of csn distribution
    %           Sigma: 2nd parameter
    %           Gamma: 3rd parameter
    %           nu   : 4th parameter
    %           Delta: 5th parameter
    %           
    % Outputs:
    %           mu_A   : 1st parameter of linearly transformed csn variable
    %           Sigma_A: 2st parameter
    %           Gamma_A: 3st parameter
    %           nu_A   : 4st parameter
    %           Delta_A: 5st parameter
    %           
    % Literature:
    %
    % Moreover:
    %       - Codes by G.Guljanov 06.02.2023
    %       - All the vectors should be column vectors
    
    % First two parameters:
    mu_A = A_mat * mu;
    Sigma_A = A_mat * Sigma * A_mat';

    % The other parameters:
    % Make mu a zero vector
    mu = zeros(size(mu));

    %----------------------------------------------
    % Step 1 -- SSN to CSN
    %----------------------------------------------
   
    if rcond(Sigma) < 1e-10

        % Sigma_bar is rank deficient, use ssn2csn() fun
        [S_hat, ~, Sigma, Gamma] = ssn2csn(mu, Sigma, Gamma, 1e-8, 0);
        
        A_mat = A_mat * S_hat;

    end
    

    %--------------------------------
    % Step 2 -- Linear Transformation
    %--------------------------------

    % Check if A_mat is rank deficient        
    dimensions = size(A_mat);

    if dimensions(1) > dimensions(2)
        is_singular = rcond(A_mat' * A_mat) < 1e-8;
    else
        is_singular = rcond(Sigma_A) < 1e-8;
    end

    % If not full rank, use rank deficient linear transformation, i.e. rd_lintransf() fun
    if is_singular
        
        [Gamma_A, nu_A, Delta_A] = rd_lintransf(...
            A_mat, Sigma, Gamma, nu, Delta, 1e-8 ...
        );

        return
        
    end


    [Gamma_A, nu_A, Delta_A] = fr_lintransf(...
        A_mat, Sigma, Sigma_A, Gamma, nu, Delta ...
    );
        
end