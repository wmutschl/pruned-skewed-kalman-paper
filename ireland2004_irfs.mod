% Simulates the model of Ireland (2004) with Gaussian and CSN distributed shocks
% Computes statistics on recessions based on simulated time series
% =========================================================================
% Copyright (C) 2025 Willi Mutschler
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
% Dynare always computes irfs with respect to a one standard deviation shock,
% so we abuse the shocks block to set the stderr to the quantile
% and adjust parameter SIGN_SHOCKS to get the correct sign on shocks (whether positive or negative)

@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GAUSSIAN                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use maximum likelihood estimates from ireland2004_ml_1_gaussian
OMEGA    = 0.058086387283355;
RHO_PI   = 0.386477141360770;
RHO_G    = 0.396013670264708;
RHO_X    = 0.165402611432957;
RHO_A    = 0.904795669581000;
RHO_E    = 0.990673584292426;
% shock parameters from ireland2004_ml_1_gaussian
gauss_std_eta_a = 3.016725412332537;
gauss_std_eta_e = 0.024764214899394;
gauss_std_eta_z = 0.886476699633342;
gauss_std_eta_r = 0.279030970775836;
% compute quantiles
gauss_q16 = [norminv(0.16,0,gauss_std_eta_a);
             norminv(0.16,0,gauss_std_eta_e);
             norminv(0.16,0,gauss_std_eta_z);
             norminv(0.16,0,gauss_std_eta_r);
            ];
gauss_q84 = [norminv(0.84,0,gauss_std_eta_a);
             norminv(0.84,0,gauss_std_eta_e);
             norminv(0.84,0,gauss_std_eta_z);
             norminv(0.84,0,gauss_std_eta_r);
            ];


% negative Gaussian shocks, i.e. 16th quantiles
SIGN_SHOCKS = -1; % flip sign of shocks for 16th quantiles as they are negative
shocks;
var eta_a; stderr (gauss_q16(1));
var eta_e; stderr (gauss_q16(2));
var eta_z; stderr (gauss_q16(3));
var eta_r; stderr (gauss_q16(4));
end;
stoch_simul(order=1,periods=0,irf=15,nodecomposition,nomoments,nocorr,nofunctions);
irfs_gaussian_neg = oo_.irfs;

% positive Gaussian shocks, i.e. 84th quantiles
SIGN_SHOCKS = 1; % don't flip sign of shocks for 84th quantiles as they are positive
shocks;
var eta_a; stderr (gauss_q84(1));
var eta_e; stderr (gauss_q84(2));
var eta_z; stderr (gauss_q84(3));
var eta_r; stderr (gauss_q84(4));
end;
stoch_simul(order=1,periods=0,irf=15,nodecomposition,nomoments,nocorr,nofunctions);
irfs_gaussian_pos = oo_.irfs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   CSN                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use maximum likelihood estimates from ireland2004_ml_3_csn
OMEGA    = 0.159645725481492;
RHO_PI   = 0.281005275686837;
RHO_G    = 0.338457069937301;
RHO_X    = 0.287099131411956;
RHO_A    = 0.913905604780658;
RHO_E    = 0.980518591389153;
% shock parameters from ireland2004_ml_3_csn
csn_std_eta_a = 2.433749588176093;
csn_std_eta_e = 0.020634966505142;
csn_std_eta_z = 0.791304331135172;
csn_std_eta_r = 0.285326737196001;
csn_skew_eta_a = -0.192359778594322;
csn_skew_eta_e = -0.217383469201676;
csn_skew_eta_z = -0.994999998739149;
csn_skew_eta_r =  0.817071416346745;
% get shock parameters from stderr and skew
[csn_Sigma_eta_a,csn_Gamma_eta_a] = csnVarSkew_To_SigmaGamma_univariate(csn_std_eta_a,csn_skew_eta_a,false);
[csn_Sigma_eta_e,csn_Gamma_eta_e] = csnVarSkew_To_SigmaGamma_univariate(csn_std_eta_e,csn_skew_eta_e,false);
[csn_Sigma_eta_z,csn_Gamma_eta_z] = csnVarSkew_To_SigmaGamma_univariate(csn_std_eta_z,csn_skew_eta_z,false);
[csn_Sigma_eta_r,csn_Gamma_eta_r] = csnVarSkew_To_SigmaGamma_univariate(csn_std_eta_r,csn_skew_eta_r,false);
csn_mu_eta_a = -csnMean(0,csn_Sigma_eta_a,csn_Gamma_eta_a,0,1);
csn_mu_eta_e = -csnMean(0,csn_Sigma_eta_e,csn_Gamma_eta_e,0,1);
csn_mu_eta_z = -csnMean(0,csn_Sigma_eta_z,csn_Gamma_eta_z,0,1);
csn_mu_eta_r = -csnMean(0,csn_Sigma_eta_r,csn_Gamma_eta_r,0,1);
% compute quantiles
csn_q16 = [csnQuantile(0.16,csn_mu_eta_a,csn_Sigma_eta_a,csn_Gamma_eta_a,0,1);
           csnQuantile(0.16,csn_mu_eta_e,csn_Sigma_eta_e,csn_Gamma_eta_e,0,1);
           csnQuantile(0.16,csn_mu_eta_z,csn_Sigma_eta_z,csn_Gamma_eta_z,0,1);
           csnQuantile(0.16,csn_mu_eta_r,csn_Sigma_eta_r,csn_Gamma_eta_r,0,1);
          ];
csn_q84 = [csnQuantile(0.84,csn_mu_eta_a,csn_Sigma_eta_a,csn_Gamma_eta_a,0,1);
           csnQuantile(0.84,csn_mu_eta_e,csn_Sigma_eta_e,csn_Gamma_eta_e,0,1);
           csnQuantile(0.84,csn_mu_eta_z,csn_Sigma_eta_z,csn_Gamma_eta_z,0,1);
           csnQuantile(0.84,csn_mu_eta_r,csn_Sigma_eta_r,csn_Gamma_eta_r,0,1);
          ];


% negative CSN shocks, i.e. 16th quantiles
SIGN_SHOCKS = -1; % flip sign of shocks for 16th quantiles as they are negative
shocks;
var eta_a; stderr (csn_q16(1));
var eta_e; stderr (csn_q16(2));
var eta_z; stderr (csn_q16(3));
var eta_r; stderr (csn_q16(4));
end;
stoch_simul(order=1,periods=0,irf=15,nodecomposition,nomoments,nocorr,nofunctions);
irfs_csn_neg = oo_.irfs;

% positive CSN shocks, i.e. 84th quantiles
SIGN_SHOCKS = 1; % don't flip sign of shocks for 84th quantiles as they are positive
shocks;
var eta_a; stderr (csn_q84(1));
var eta_e; stderr (csn_q84(2));
var eta_z; stderr (csn_q84(3));
var eta_r; stderr (csn_q84(4));
end;
stoch_simul(order=1,periods=0,irf=15,nodecomposition,nomoments,nocorr,nofunctions);
irfs_csn_pos = oo_.irfs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              STORE RESULTS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([M_.dname filesep 'Output' filesep M_.fname],'irfs_gaussian_neg','irfs_gaussian_pos','irfs_csn_neg','irfs_csn_pos','-v6');
