function [Gamma_A, nu_A, Delta_A] = fr_lintransf(A_mat, Sigma, Sigma_A, Gamma, nu, Delta)

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
    
    %---------------
    %   More rows
    %---------------
    if nrow >= ncol

        nu_A = nu;
        Delta_A = Delta;
    
        if nrow == ncol
            Gamma_A = Gamma / A_mat;
            return
        end

        Gamma_A = Gamma / (A_mat' * A_mat) * A_mat';
        
        return
        
    end


    %-----------------
    %   More columns
    %-----------------
    Tmp_mat = Sigma * A_mat' / Sigma_A;

    Gamma_A = Gamma * Tmp_mat;
    nu_A = nu;
    Delta_A = Delta ...
            + Gamma * Sigma * Gamma' ...
            - Gamma * Tmp_mat * A_mat * Sigma * Gamma';

end