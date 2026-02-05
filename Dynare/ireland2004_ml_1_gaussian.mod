% =========================================================================
% Copyright (C) 2024-2026 Willi Mutschler
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
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMUM LIKELIHOOD ESTIMATION OF GAUSSIAN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD %
% NOTE THAT WE DON'T USE DYNARE'S KALMAN FILTER FUNCTION FOR A FAIR COMPARISON                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimated_params;
stderr eta_a, 3.0200,  0,    10;
stderr eta_e, 0.0200,  0,    10;
stderr eta_z, 0.8900,  0,    10;
stderr eta_r, 0.2800,  0,    10;
OMEGA,        0.0581,  0,     1;
RHO_PI,       0.3866,  0,     1;
RHO_G,        0.3960,  0,     1;
RHO_X,        0.1654,  0,     1;
RHO_A,        0.9048,  0,     1;
RHO_E,        0.9907,  0,     1;
end;

estimation(datafile = '../data/ireland2004_data.m'
          , mode_compute = 8 % 1,2,7,8,101 yield almost the same estimates
          , silent_optimizer % below we display optimization_info, so don't show intermediate optimization output
          , kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability;
                             % note that Dynare's implementation is much faster because it switches to the steady-state Kalman filter which we have not implemented yet for the PSKF
          , lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          );

fprintf('Optimization info:\n');
disp(oo_.posterior.optimization.optimization_info);
format long;
fprintf('MLE point estimates for structural parameters:\n');
disp(oo_.mle_mode.parameters);
fprintf('MLE point estimates for stderr parameters:\n');
disp(oo_.mle_mode.shocks_std);
format short;
fprintf('Final value of log-likelihood: %.15f\n', oo_.posterior.optimization.log_density);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE ESTIMATED SHOCK PARAMETERS AND SMOOTHED SHOCK VALUES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writetable(table(oo_.SmoothedShocks.eta_a, oo_.SmoothedShocks.eta_e, ...
                 oo_.SmoothedShocks.eta_z, oo_.SmoothedShocks.eta_r, ...
                 'VariableNames', {'eta_a', 'eta_e', 'eta_z', 'eta_r'}), ...
           ['../results/ireland2004/ml/' M_.fname '_smoothed_shocks.csv']);
writematrix(M_.csn.mu_e, ['../results/ireland2004/ml/' M_.fname '_csn_mu_e.csv']);
writematrix(M_.csn.Sigma_e, ['../results/ireland2004/ml/' M_.fname '_csn_Sigma_e.csv']);
writematrix(M_.csn.Gamma_e, ['../results/ireland2004/ml/' M_.fname '_csn_Gamma_e.csv']);
writematrix(M_.csn.nu_e, ['../results/ireland2004/ml/' M_.fname '_csn_nu_e.csv']);
writematrix(M_.csn.Delta_e, ['../results/ireland2004/ml/' M_.fname '_csn_Delta_e.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE LATEX TABLE 2 %
%%%%%%%%%%%%%%%%%%%%%%%%
addpath('_utils');
update_table_2('../results/ireland2004', 'ML_Gaussian', oo_, ARCH, MATLAB_VERSION);
rmpath('_utils');

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
