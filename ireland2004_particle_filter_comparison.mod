% This mod file compares the value of the log-likelihood computed by the
% Pruned Skewed Kalman filter with the value of the log-likelihood computed
% with the sequential importance particle filter (aka Bootstrap particle filter)
% The value is available in the log file under:
% "Initial value of the log posterior (or likelihood)"
% =========================================================================
% Copyright (C) 2026 Willi Mutschler
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
addpath('MATLAB/particle_overwrite_dynare');
@#include "_ireland2004_common.inc"

@#for DIST in ["gaussian", "csn"]

  @#for ME in [40, 20, 10, 5]

estimated_params(overwrite);
% CSN POSTERIOR MEAN
stderr eta_a, 3.2359;
stderr eta_e, 0.0603;
stderr eta_z, 0.7975;
stderr eta_r, 0.2920;
OMEGA,        0.1392;
RHO_PI,       0.4662;
RHO_G,        0.3586;
RHO_X,        0.2227;
RHO_A,        0.9216;
RHO_E,        0.9030;
    @#if DIST == "csn"
skew eta_a,  -0.0819;
skew eta_e,  -0.3584;
skew eta_z,  -0.3808;
skew eta_r,   0.5183;
    @#endif
% measurement errors set to @{ME}% of standard error of observables
stderr ghat,  @{ME}/100*0.007775644957988;
stderr rhat,  @{ME}/100*0.007740297075338;
stderr pihat, @{ME}/100*0.005142141163420;
end;

    @#for TOL in [0.1000, 0.0100, 0.0010, 0.0001]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOG-LIKELIHOOD WITH PRUNED SKEWED KALMAN FILTER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg = sprintf('# %s (ME=@{ME}%%): EVALUATE LOG-LIKELIHOOD VALUE WITH PRUNED SKEWED KALMAN FILTER (TOL=@{TOL}) #', upper('@{DIST}'));
fprintf('\n\n%s\n%s\n%s', repmat('#',1,length(msg)), msg, repmat('#',1,length(msg)));
options_.particle.status = false;
estimation(datafile = 'data/ireland2004_data.m', mode_compute = 0, cova_compute = 0, kalman_algo = 5, skewed_kalman_prune_tol = @{TOL}, skewed_kalman_smoother_skip);
    @#endfor // TOL

    @#for N in [2000, 20000, 200000, 2000000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOG-LIKELIHOOD WITH PARTICLE FILTER (N=@{N}) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg = sprintf('# %s (ME=@{ME}%%): EVALUATE LOG-LIKELIHOOD VALUE WITH PARTICLE FILTER (N=@{N}) #', upper('@{DIST}'));
fprintf('\n\n%s\n%s\n%s', repmat('#',1,length(msg)), msg, repmat('#',1,length(msg)));
options_.particle.status = true; options_.particle.use_reduced_rank_cholesky = true;
estimation(datafile = 'data/ireland2004_data.m', mode_compute = 0, cova_compute = 0, filter_algorithm = sis, number_of_particles = @{N});
    @#endfor // N

  @#endfor // ME

@#endfor // DIST

rmpath('MATLAB/particle_overwrite_dynare');