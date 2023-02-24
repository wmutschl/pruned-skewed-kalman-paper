function [Gamma_A, nu_A, Delta_A] = rd_lintransf(A_mat, Sigma, Gamma, nu, Delta, tol)

    % Briefly: 
    %       - Rank deficient linear transformation of CSN distributed random variable
    %       - Returns only skewness parameters, as other parameters are easy to obtain
    %
    % Structure:
    %       - The notation and the structure of this function closely ...
    %         follows [GU.MU.2023], esp. theorem about rank deficient linear transformation
    % 
    % Inputs:
    %       A_mat: Matrix to multiply CSN distributed random variable
    %
    %       Sigma: Second parameter of the CSN distribution
    %       Gamma: Third parameter of the CSN distribution
    %       nu   : Fourth parameter of the CSN distribution
    %       Delta: Fifth parameter of the CSN distribution
    %
    %       tol  : Tolerance level showing which values should be assumed to be ...
    %              numerically zero
    % Outputs:
    %       Gamma_A: Gamma parameter of the CSN after linear transformation
    %       nu_A   : nu parameter of the CSN after linear transformation
    %       Delta_A: Delta parameter of the CSN after linear transformation
    %
    %
    % Literature:
    %       1. [GU.MU.2023]   : Our paper
    %
    % Moreover:
    %       - Codes by G.Guljanov 06.02.2023
    %       - All the vectors should be column vectors

    [nrow, ncol] = size(A_mat);
    
    difference = abs(nrow - ncol);

    % Apply singular value decomposition
    [S_mat, Lambda_mat, T_mat] = svd(A_mat);

    if nrow >= ncol % more rows
        S_mat = S_mat(:, 1:end-difference);
        Lambda_mat = Lambda_mat(1:end-difference, :);
    else
        T_mat = T_mat(:, 1:end-difference);
        Lambda_mat = Lambda_mat(:, 1:end-difference);
    end

    % Which dimensions to delete
    hold_vec = abs(diag(Lambda_mat)) > tol;

    % Delete respective columns
    S_mat = S_mat(:, hold_vec);
    
    % Delete respective rows and columns
    Lambda_mat = Lambda_mat(hold_vec, hold_vec);

    % Delete respective columns
    T_mat = T_mat(:, hold_vec);

    % Apply the final formulas of our theorem
    Tmp_mat = T_mat' * Sigma * T_mat;
    Tmp_mat2 = Gamma * Sigma;
    Tmp_mat3 = Tmp_mat2 * T_mat / Tmp_mat;

    Gamma_A = (Tmp_mat3 / Lambda_mat) / (S_mat' * S_mat) * S_mat';
    nu_A = nu;
    Delta_A = ( ...
        Delta + Tmp_mat2 * Gamma' - Tmp_mat3 * T_mat' * Sigma * Gamma' ...
    );

end