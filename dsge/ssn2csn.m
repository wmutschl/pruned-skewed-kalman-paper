function [S_hat, mu_hat, Lambda_hat, Gamma_hat] = ssn2csn(mu, Sigma, Gamma, tol, hmdel)

    % Briefly: 
    %       - Write Singular Skew Normal (SSN) distribution as a linear transformation of
    %         proper Closed Skew Normal (CSN) distribution
    %
    % Structure:
    %       - The notation and the structure of this function closely ...
    %         follows [GU.MU.2023], esp. theorem about how to write ...
    %         SSN as linear transformation of CSN
    %       - For the somehow restrictive definition of SSN, ...
    %         see theorem 2 of [GO.DO.GU.2004]. ...
    %         But in this function and in [GU.MU.2023] we define SSN in ...
    %         more general way: whenever second parameter of the CSN is singular, ...
    %         it is SSN.
    % 
    % Inputs:
    %       mu   : First parameter of the SSN distribution
    %       Sigma: Second parameter of the SSN distribution
    %       Gamma: Third parameter of the SSN distribution
    %       tol  : Tolerance level showing which values should be assumed to be ...
    %              numerically zero
    %       hmdel: Minimum how many dimensions should be deleted, if known already
    %
    % Outputs:
    %       S_hat     : Multiplying matrix, which has more rows than columns
    %       mu_hat    : First parameter of the proper CSN
    %       Lambda_hat: Second parameter of the proper CSN
    %       Gamma_hat : Third parameter of the proper CSN
    %
    %
    % Literature:
    %       1. [GU.MU.2023]   : Our paper
    %       2. [GO.DO.GU.2004]: Gonzales-Farias, Dominguez-Molina, Gupta (2004) - ...
    %                           Additive properties of skew normal random vectors -- ...
    %                           Theorem 2 (page 527)
    %
    % Moreover:
    %       - Codes by G.Guljanov 06.02.2023
    %       - All the vectors should be column vectors

    if sum(mu ~= 0)
        error("mu should be a zero vector")
    end

    % Get Schur decomposition
    [S_mat, Lambda] = schur(Sigma);

    % Create a hold_vec, which tells which dimensions to hold
    if hmdel > 0
        
        % Reorder the diagonal elements of Lambda, in an increasing order
        len = length(Lambda);
        diag_el = diag(Lambda);
        [~, min_indices] = mink(diag_el, len);
    
        % Should we delete more than hmdel
        hmdel_tmp = hmdel;
    
        ii = 1;
        while diag_el(min_indices(hmdel_tmp + ii)) / diag_el(min_indices(len)) <= tol
            ii = ii + 1;
        end
        
        hmdel_tmp = hmdel_tmp + ii - 1;
        
        if hmdel_tmp >= len
            error("Sigma matrix should at least have rank one");
        end
            
        % Which vector to hold, i.e. to not cut
        hold_vec = true(len, 1);
        hold_vec(min_indices(1:hmdel_tmp)) = false;

    else

        hold_vec = diag(Lambda) > tol;

    end
    
    % Hold the relevant dimensions, i.e. delete the unnecessary ones
    S_hat = S_mat(:, hold_vec);
    
    mu_hat = zeros(sum(hold_vec), 1);
    
    Lambda_hat = Lambda(hold_vec, hold_vec);

    Gamma_hat = Gamma * S_hat;

end
