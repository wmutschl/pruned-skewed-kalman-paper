% Simulates the model of Ireland (2004) with Gaussian and CSN distributed shocks
% and computes statistics on recessions based on simulated time series
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

@#define TRANSFORM_PARAMETERS    = 0  // new feature to transform selected bounded parameters to unbounded domain using logit transform (only during optimization)
@#include "_ireland2004_common.inc"
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
dyn_saveas(hh_fig_@{j},[M_.dname, '/graphs/' M_.fname '_histogram_@{j}'],options_.nodisplay,options_.graph_format);

[meanDuration_@{j},meanOutputLoss_@{j},recessionNumber_@{j},effectiveSampleSize_@{j}] = ireland2004_recessions_statistics(ghat, -0.005, 2);

@#endfor

%% print Latex table entries
fprintf('Number recessions (\\(\\hat{g}_{t} \\leq -0.5\\%% \\)) & %d  & %d  & %d  & %d  & %d  & %d  \\\\ \n',[recessionNumber_GAUSSIAN recessionNumber_CSN round(recessionNumber_GAUSSIAN/3) round(recessionNumber_CSN/3) round(recessionNumber_GAUSSIAN/3) round(recessionNumber_CSN/3)]);
fprintf('Frequency recessions (in \\%%)                     & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n',100*[recessionNumber_GAUSSIAN/effectiveSampleSize_GAUSSIAN recessionNumber_CSN/effectiveSampleSize_CSN recessionNumber_GAUSSIAN/(3*effectiveSampleSize_GAUSSIAN) recessionNumber_CSN/(3*effectiveSampleSize_CSN) recessionNumber_GAUSSIAN/(3*effectiveSampleSize_GAUSSIAN) recessionNumber_CSN/(3*effectiveSampleSize_CSN)]);
fprintf('Mean duration (peak to trough, quarters)         & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n',meanDuration_GAUSSIAN(1), meanDuration_CSN(1), meanDuration_GAUSSIAN(2), meanDuration_CSN(2), meanDuration_GAUSSIAN(3), meanDuration_CSN(3));
fprintf('Mean output loss (peak to trough, in \\%%)         & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n', meanOutputLoss_GAUSSIAN(1), meanOutputLoss_CSN(1), meanOutputLoss_GAUSSIAN(2), meanOutputLoss_CSN(2), meanOutputLoss_GAUSSIAN(3), meanOutputLoss_CSN(3));