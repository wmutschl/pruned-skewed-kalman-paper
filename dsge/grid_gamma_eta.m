function gamma_eta_grid = grid_gamma_eta(sqrt_Sigma_eta_a,sqrt_Sigma_eta_e,sqrt_Sigma_eta_z,sqrt_Sigma_eta_r)
gamma_eta_a( 1) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_a^2,x) + 0.95,-200,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_a( 2) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_a^2,x) + 0.65,-70,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_a( 3) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_a^2,x) + 0.35,-40,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_a( 4) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_a^2,x) + 0.05,-15,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_a = [gamma_eta_a 0 -gamma_eta_a];

gamma_eta_e( 1) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_e^2,x) + 0.95,-7500,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_e( 2) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_e^2,x) + 0.65,-2500,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_e( 3) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_e^2,x) + 0.35,-1300,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_e( 4) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_e^2,x) + 0.05, -550,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_e = [gamma_eta_e 0 -gamma_eta_e];

gamma_eta_z( 1) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_z^2,x) + 0.95,-850,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_z( 2) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_z^2,x) + 0.65,-250,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_z( 3) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_z^2,x) + 0.35,-150,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_z( 4) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_z^2,x) + 0.05,-50,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_z = [gamma_eta_z 0 -gamma_eta_z];

gamma_eta_r( 1) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_r^2,x) + 0.95,-3000,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_r( 2) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_r^2,x) + 0.65,-1000,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_r( 3) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_r^2,x) + 0.35, -500,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_r( 4) = fsolve(@(x) skewness_coef_theor(sqrt_Sigma_eta_r^2,x) + 0.05, -200,optimset('display','off','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
gamma_eta_r = [gamma_eta_r 0 -gamma_eta_r];

gamma_eta_grid = combvec(gamma_eta_a,gamma_eta_e,gamma_eta_z,gamma_eta_r);