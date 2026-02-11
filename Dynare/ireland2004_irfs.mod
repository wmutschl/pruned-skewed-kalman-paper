% Computes impulse response functions of the Ireland (2004) with
% Gaussian and CSN distributed shocks using the 16th and 84th percentiles
% of the maximum likelihood estimates
% =========================================================================
% Copyright (C) 2025-2026 Willi Mutschler
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
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

% NOTE:
% Dynare always computes irfs with respect to a one standard deviation shock,
% which is provided by the user in the shocks block.
% Therefore, we use the shocks block to set the stderr to the quantiles
% and adjust the auxiliary parameter SIGN_SHOCKS to get the correct
% sign on shocks (whether positive or negative).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 MODEL                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#include "_ireland2004_common.inc"
var rAnnualized piAnnualized;
model;
rAnnualized = 4*rhat;
piAnnualized = 4*pihat;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GAUSSIAN                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use maximum likelihood estimates from ireland2004_ml_1_gaussian
set_param_value('OMEGA', 0.0581);
set_param_value('RHO_PI', 0.3865);
set_param_value('RHO_G', 0.3960);
set_param_value('RHO_X', 0.1654);
set_param_value('RHO_A', 0.9048);
set_param_value('RHO_E', 0.9907);
set_shock_stderr_value('eta_a',3.0167);
set_shock_stderr_value('eta_e',0.0248);
set_shock_stderr_value('eta_z',0.8865);
set_shock_stderr_value('eta_r',0.2790);

% compute quantiles of Gaussian distribution
gauss_q16 = [norminv(0.16, 0, sqrt(M_.Sigma_e(1,1)));
             norminv(0.16, 0, sqrt(M_.Sigma_e(2,2)));
             norminv(0.16, 0, sqrt(M_.Sigma_e(3,3)));
             norminv(0.16, 0, sqrt(M_.Sigma_e(4,4)));
            ];
gauss_q84 = [norminv(0.84, 0, sqrt(M_.Sigma_e(1,1)));
             norminv(0.84, 0, sqrt(M_.Sigma_e(2,2)));
             norminv(0.84, 0, sqrt(M_.Sigma_e(3,3)));
             norminv(0.84, 0, sqrt(M_.Sigma_e(4,4)));
            ];

% simulate negative Gaussian shocks, i.e. 16th quantiles
SIGN_SHOCKS = -1; % flip sign of shocks for 16th quantiles as they are negative
shocks;
var eta_a; stderr (gauss_q16(1));
var eta_e; stderr (gauss_q16(2));
var eta_z; stderr (gauss_q16(3));
var eta_r; stderr (gauss_q16(4));
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions);
irfs_gaussian_neg = oo_.irfs;

% simulate positive Gaussian shocks, i.e. 84th quantiles
SIGN_SHOCKS = 1; % don't flip sign of shocks for 84th quantiles as they are positive
shocks;
var eta_a; stderr (gauss_q84(1));
var eta_e; stderr (gauss_q84(2));
var eta_z; stderr (gauss_q84(3));
var eta_r; stderr (gauss_q84(4));
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions);
irfs_gaussian_pos = oo_.irfs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   CSN                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use maximum likelihood estimates from ireland2004_ml_3_csn
set_param_value('OMEGA', 0.15356);
set_param_value('RHO_PI', 0.28761);
set_param_value('RHO_G', 0.33795);
set_param_value('RHO_X', 0.28325);
set_param_value('RHO_A', 0.91704);
set_param_value('RHO_E', 0.98009);
set_shock_stderr_value('eta_a',2.5232);
set_shock_stderr_value('eta_e',0.021228);
set_shock_stderr_value('eta_z',0.79002);
set_shock_stderr_value('eta_r',0.28384);
set_shock_skew_value('eta_a',-0.1948);
set_shock_skew_value('eta_e',-0.21401);
set_shock_skew_value('eta_z',-0.99527);
set_shock_skew_value('eta_r',0.81275);
M_.csn = csn_update_specification(M_.Sigma_e, M_.Skew_e);

% compute quantiles of CSN distribution
csn_q16 = [csn_quantile(0.16, M_.csn.mu_e(1,1), M_.csn.Sigma_e(1,1), M_.csn.Gamma_e(1,1), 0, 1, 'mvncdf');
           csn_quantile(0.16, M_.csn.mu_e(2,1), M_.csn.Sigma_e(2,2), M_.csn.Gamma_e(2,2), 0, 1, 'mvncdf');
           csn_quantile(0.16, M_.csn.mu_e(3,1), M_.csn.Sigma_e(3,3), M_.csn.Gamma_e(3,3), 0, 1, 'mvncdf');
           csn_quantile(0.16, M_.csn.mu_e(4,1), M_.csn.Sigma_e(4,4), M_.csn.Gamma_e(4,4), 0, 1, 'mvncdf');
          ];
csn_q84 = [csn_quantile(0.84, M_.csn.mu_e(1,1), M_.csn.Sigma_e(1,1), M_.csn.Gamma_e(1,1), 0, 1, 'mvncdf');
           csn_quantile(0.84, M_.csn.mu_e(2,1), M_.csn.Sigma_e(2,2), M_.csn.Gamma_e(2,2), 0, 1, 'mvncdf');
           csn_quantile(0.84, M_.csn.mu_e(3,1), M_.csn.Sigma_e(3,3), M_.csn.Gamma_e(3,3), 0, 1, 'mvncdf');
           csn_quantile(0.84, M_.csn.mu_e(4,1), M_.csn.Sigma_e(4,4), M_.csn.Gamma_e(4,4), 0, 1, 'mvncdf');
          ];

% simulate negative CSN shocks, i.e. 16th quantiles
SIGN_SHOCKS = -1; % flip sign of shocks for 16th quantiles as they are negative
shocks;
var eta_a; stderr (csn_q16(1));
var eta_e; stderr (csn_q16(2));
var eta_z; stderr (csn_q16(3));
var eta_r; stderr (csn_q16(4));
% no need to specify for IRFs as the size of the shock is one standard deviation
skew eta_a = -0.1948;
skew eta_e = -0.21401;
skew eta_z = -0.99527;
skew eta_r =  0.81275;
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions);
irfs_csn_neg = oo_.irfs;

% simulate positive CSN shocks, i.e. 84th quantiles
SIGN_SHOCKS = 1; % don't flip sign of shocks for 84th quantiles as they are positive
shocks;
var eta_a; stderr (csn_q84(1));
var eta_e; stderr (csn_q84(2));
var eta_z; stderr (csn_q84(3));
var eta_r; stderr (csn_q84(4));
% no need to specify for IRFs as the size of the shock is one standard deviation
skew eta_a = -0.1948;
skew eta_e = -0.21401;
skew eta_z = -0.99527;
skew eta_r =  0.81275;
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions);
irfs_csn_pos = oo_.irfs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              STORE RESULTS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert IRF structures to a single tidy-format CSV table
irf_data   = {irfs_gaussian_neg, irfs_gaussian_pos, irfs_csn_neg, irfs_csn_pos};
irf_distrs = {'Gaussian',        'Gaussian',        'CSN',        'CSN'};
irf_signs  = {'neg',             'pos',             'neg',        'pos'};

fnames = fieldnames(irfs_gaussian_neg);
n_fields = length(fnames);
irf_length = length(irfs_gaussian_neg.(fnames{1}));
n_rows = 4 * n_fields * irf_length;

% preallocate
all_distribution = cell(n_rows, 1);
all_sign         = cell(n_rows, 1);
all_shock        = cell(n_rows, 1);
all_variable     = cell(n_rows, 1);
all_time         = zeros(n_rows, 1);
all_value        = zeros(n_rows, 1);

idx = 0;
for k = 1:4
    s = irf_data{k};
    for j = 1:n_fields
        fname = fnames{j};
        vals = s.(fname);
        parts = strsplit(fname, '_eta_');
        variable_name = parts{1};
        shock_name = ['eta_' parts{2}];
        for t = 1:irf_length
            idx = idx + 1;
            all_distribution{idx} = irf_distrs{k};
            all_sign{idx}         = irf_signs{k};
            all_shock{idx}        = shock_name;
            all_variable{idx}     = variable_name;
            all_time(idx)         = t;
            all_value(idx)        = vals(t);
        end
    end
end

tbl_irfs = table(all_distribution, all_sign, all_shock, all_variable, all_time, all_value, ...
                 'VariableNames', {'distribution', 'sign', 'shock', 'variable', 'time', 'value'});
writetable(tbl_irfs, sprintf('../results/ireland2004/irfs_%s_%s.csv', ARCH, MATLAB_VERSION));

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
