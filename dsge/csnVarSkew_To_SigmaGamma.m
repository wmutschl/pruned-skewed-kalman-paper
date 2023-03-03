function [Sigma,Gamma] = csnVarSkew_To_SigmaGamma(Var,Skew,check)
if nargin < 3
    check = false;
end
Sigma = (-2^(2/3) * (4-pi)^(1/3) * abs(Skew)^(2/3) + pi - 4) * Var / (pi - 4);
Gamma = sign(Skew) * ( 2^(1/3) * (4-pi)^(2/3)*sqrt(pi)*sqrt( abs(Skew)^(2/3) / ( ( -2 * 2^(1/3) * (4-pi)^(2/3) * (pi-2) * abs(Skew)^(4/3) + 2^(2/3) * (4-pi)^(7/3) * abs(Skew)^(2/3) + 2 * (pi - 4)^2 ) * Var ) ) );

if check 
    if (abs(csnVar(Sigma,Gamma,0,1) - Var)>1e-13) || (abs(skewness_coef_theor(Sigma,Gamma) - Skew) > 1e-13)
        error('something wrong parameter transformation')
    end
end