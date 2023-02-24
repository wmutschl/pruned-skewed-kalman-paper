function skewness = skewness_coef_theor(Sigma, Gamma)

    % Briefly: Evaluates the skewness coefficient of a univariate CSN distribution ...
    %          according to equation 3.6 of [Gr.Kl.Ko.2011] ...
    %          which assumes that delta2 = 1
    % 
    % Structure:
    % 
    % Inputs:
    %       sigma2: 2nd parameter of the CSN distribution
    %       gamma : 3rd parameter of the CSN distribution
    % 
    % Ouputs:
    %       skewness: Skewness coefficient
    %
    % Literature:
    %       1) [Gr.Kl.Ko.2011]: Grabek, Klos, Koloch (2011) ...
    %                           - Skew-normal shocks in the ...
    %                             linear state space form DSGE model
    % 
    % Moreover:
    %

    term1 = (4 - pi) / 2;
    term2 = (sqrt(2/pi) * Gamma * Sigma / sqrt(1 + Gamma^2 * Sigma))^3;
    term3 = (Sigma - (2/pi) * (Gamma^2 * Sigma^2) / (1 + Gamma^2 * Sigma))^(3/2);
    
    skewness = term1 * term2 / term3;

end