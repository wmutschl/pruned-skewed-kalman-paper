% Simulates the model of Ireland (2004) with Gaussian and CSN distributed shocks
% and computes statistics on recessions based on simulated time series.
% Recessions are defined by a period where output growth falls below -0.5% per quarter
% and stays negative for at least two quarters.
% Severe and mild recessions are determined by the top and bottom terciles of the implied
% distribution of peak-to-trough output losses.
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
% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
@#include "_ireland2004_common.inc"

% Calibration follows CSN-based maximum likelihood estimates from Table 2 in the paper
ALPHA_X  = 0;
ALPHA_PI = 0;
OMEGA    = 0.15356;
RHO_PI   = 0.28761;
RHO_G    = 0.33795;
RHO_X    = 0.28325;
RHO_A    = 0.91704;
RHO_E    = 0.98009;

shocks;
var eta_a; stderr 2.5232;
var eta_e; stderr 0.021228;
var eta_z; stderr 0.79002;
var eta_r; stderr 0.28384;
end;

CheckPath('graphs',M_.dname);

@#for j in ["GAUSSIAN","CSN"]
    @#if j == "CSN"
shocks;
skew eta_a = -0.1948;
skew eta_e = -0.21401;
skew eta_z = -0.99527;
skew eta_r = +0.81275;
end;
    @#endif
%set_dynare_seed('clock'); % uncomment for random draws
set_dynare_seed(75);
stoch_simul(order=1,periods=500000,irf=0,nodecomposition,nomoments,nocorr,nofunctions);
send_exogenous_variables_to_workspace;
send_endogenous_variables_to_workspace;

hh_fig_@{j} = dyn_figure(options_.nodisplay,'Name','Histogram of @{j} shocks');
subplot(2,2,1); histogram(eta_a,'normalization','pdf'); title('\eta_a');
subplot(2,2,2); histogram(eta_e,'normalization','pdf'); title('\eta_e');
subplot(2,2,3); histogram(eta_z,'normalization','pdf'); title('\eta_z');
subplot(2,2,4); histogram(eta_r,'normalization','pdf'); title('\eta_r');

addpath('_utils');
[meanDuration_@{j},meanOutputLoss_@{j},recessionNumber_@{j},effectiveSampleSize_@{j}] = compute_recessions_statistics(ghat, -0.005, 2);
rmpath('_utils');

@#endfor

%% print Latex table entries
% build data row strings using sprintf
row1 = sprintf('Number recessions (\\(\\hat{g}_{t} \\leq -0.5\\%% \\)) & %d & %d & %.2f & %d & %d & %.2f & %d & %d & %.2f \\\\', ...
    recessionNumber_GAUSSIAN(1), recessionNumber_CSN(1), 100*(recessionNumber_CSN(1)/recessionNumber_GAUSSIAN(1)-1), ...
    recessionNumber_GAUSSIAN(2), recessionNumber_CSN(2), 100*(recessionNumber_CSN(2)/recessionNumber_GAUSSIAN(2)-1), ...
    recessionNumber_GAUSSIAN(3), recessionNumber_CSN(3), 100*(recessionNumber_CSN(3)/recessionNumber_GAUSSIAN(3)-1));
row2 = sprintf('Frequency recessions (in \\%%)                     & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\', ...
    100*[recessionNumber_GAUSSIAN(1)/effectiveSampleSize_GAUSSIAN       recessionNumber_CSN(1)/effectiveSampleSize_CSN     ((recessionNumber_CSN(1)/effectiveSampleSize_CSN)/(recessionNumber_GAUSSIAN(1)/effectiveSampleSize_GAUSSIAN)-1) ...
         recessionNumber_GAUSSIAN(2)/effectiveSampleSize_GAUSSIAN       recessionNumber_CSN(2)/effectiveSampleSize_CSN     ((recessionNumber_CSN(2)/effectiveSampleSize_CSN)/(recessionNumber_GAUSSIAN(2)/effectiveSampleSize_GAUSSIAN)-1) ...
         recessionNumber_GAUSSIAN(3)/effectiveSampleSize_GAUSSIAN       recessionNumber_CSN(3)/effectiveSampleSize_CSN     ((recessionNumber_CSN(3)/effectiveSampleSize_CSN)/(recessionNumber_GAUSSIAN(3)/effectiveSampleSize_GAUSSIAN)-1)]);
row3 = sprintf('Mean duration (quarters)                         & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\', ...
    meanDuration_GAUSSIAN(1),   meanDuration_CSN(1),   100*(meanDuration_CSN(1)/meanDuration_GAUSSIAN(1)-1), ...
    meanDuration_GAUSSIAN(2),   meanDuration_CSN(2),   100*(meanDuration_CSN(2)/meanDuration_GAUSSIAN(2)-1), ...
    meanDuration_GAUSSIAN(3),   meanDuration_CSN(3),   100*(meanDuration_CSN(3)/meanDuration_GAUSSIAN(3)-1));
row4 = sprintf('Mean output loss (in \\%%)                         & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\', ...
    meanOutputLoss_GAUSSIAN(1), meanOutputLoss_CSN(1), 100*(meanOutputLoss_CSN(1)/meanOutputLoss_GAUSSIAN(1)-1), ...
    meanOutputLoss_GAUSSIAN(2), meanOutputLoss_CSN(2), 100*(meanOutputLoss_CSN(2)/meanOutputLoss_GAUSSIAN(2)-1), ...
    meanOutputLoss_GAUSSIAN(3), meanOutputLoss_CSN(3), 100*(meanOutputLoss_CSN(3)/meanOutputLoss_GAUSSIAN(3)-1));

% print to console and write data rows to tex file
fid_tex = fopen(sprintf('../results/ireland2004/table_3_%s_%s.tex', ARCH, MATLAB_VERSION), 'w');
for fid = [1, fid_tex] % 1 = stdout (console)
    fprintf(fid, '%s\n', row1);
    fprintf(fid, '%s\n', row2);
    fprintf(fid, '%s\n', row3);
    fprintf(fid, '%s\n', row4);
end
fclose(fid_tex);

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
