



Sigma_eta = 0.004^2;
gamma_eta_low = fsolve(@(x) skewness_coef_theor(Sigma_eta,x) + 0.99,-100,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_high = fsolve(@(x) skewness_coef_theor(Sigma_eta,x) - 0.99,-100,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
for Gamma_eta = gamma_eta_low:10:gamma_eta_high
    skew_eta = (4-pi)/2 * (sqrt(2/pi)*Gamma_eta*Sigma_eta/sqrt(1+Gamma_eta^2*Sigma_eta) )^3 / (Sigma_eta - 2/pi * Gamma_eta^2*Sigma_eta^2/(1+Gamma_eta^2*Sigma_eta))^(3/2);
    var_eta = Sigma_eta-2/pi*Gamma_eta^2*Sigma_eta^2/(1+Gamma_eta^2*Sigma_eta);
    stderr_eta = sqrt(var_eta);
    csnVar(Sigma_eta,Gamma_eta,0,1);
    skewness_coef_theor(Sigma_eta,Gamma_eta);
    Sigma_eta0 = (-2^(2/3) * (4-pi)^(1/3) * abs(skew_eta)^(2/3) + pi - 4) * (stderr_eta)^2 / (pi - 4);
    Gamma_eta0 = 2^(1/3) * (4-pi)^(2/3)*sqrt(pi)*sqrt( abs(skew_eta)^(2/3) / ( ( -2 * 2^(1/3) * (4-pi)^(2/3) * (pi-2) * abs(skew_eta)^(4/3) + 2^(2/3) * (4-pi)^(7/3) * abs(skew_eta)^(2/3) + 2 * (pi - 4)^2 ) * (stderr_eta)^2 ) );
    if skew_eta < 0 % -0.995 < skew_eta < 0, (sqrt(2)*(pi-4))/(pi-2)^(3/2)        
        Gamma_eta0 = -Gamma_eta0;
    end
    if norm([Sigma_eta-Sigma_eta0; Gamma_eta-Gamma_eta0],'Inf') > 1e-8
        error('stop')
    end
end    
