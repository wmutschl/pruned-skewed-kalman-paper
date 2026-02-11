% This mod file compares the value of the log-likelihood computed by the
% Pruned Skewed Kalman filter with the value of the log-likelihood computed
% with the sequential importance particle filter (aka Bootstrap particle filter)
% The value is available in the log file under:
% "Initial value of the log posterior (or likelihood)"
% Note that we overwrite the Dynare default particle filter implementation
% with our own in order to support the particle at first order approximation,
% using the CSN distribution for shock generation, and for simpler output to
% the log file; please refere to the functions inside the folder 
% "__particle_overwrite_dynare", which also contains the original files copied
% from the Dynare 7.0 source code (ending with *.m.orig).
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
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
addpath('_overwrite_dynare_particle_filter');
@#include "_ireland2004_common.inc"

@#define DISTRIBUTION = ["gaussian", "csn"]
@#define MEAS_ERROR = [40, 20, 10, 5]
@#define PSKF_TOL = [0.1000, 0.0100, 0.0010, 0.0001]
@#define PARTICLES = [2000, 20000, 200000, 2000000]

% initialize output structures
gaussian_pskf.fval     = nan(4, 4); gaussian_pskf.time     = nan(4, 4);
gaussian_particle.fval = nan(4, 4); gaussian_particle.time = nan(4, 4);
csn_pskf.fval          = nan(4, 4); csn_pskf.time          = nan(4, 4);
csn_particle.fval      = nan(4, 4); csn_particle.time      = nan(4, 4);

@#for DIST in 1:length(DISTRIBUTION)

  @#for ME in 1:length(MEAS_ERROR)

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
    @#if DISTRIBUTION[DIST] == "csn"
skew eta_a,  -0.0819;
skew eta_e,  -0.3584;
skew eta_z,  -0.3808;
skew eta_r,   0.5183;
    @#endif
% measurement errors set to @{MEAS_ERROR[ME]}% of standard error of observables
stderr ghat,  @{MEAS_ERROR[ME]}/100*0.007775644957988;
stderr rhat,  @{MEAS_ERROR[ME]}/100*0.007740297075338;
stderr pihat, @{MEAS_ERROR[ME]}/100*0.005142141163420;
end;

    @#for TOL in 1:length(PSKF_TOL)
msg = sprintf('# %s (ME=@{MEAS_ERROR[ME]}%%): EVALUATE LOG-LIKELIHOOD VALUE WITH PRUNED SKEWED KALMAN FILTER (TOL=@{PSKF_TOL[TOL]}) #', upper('@{DISTRIBUTION[DIST]}'));
fprintf('\n\n%s\n%s\n%s\n', repmat('#',1,length(msg)), msg, repmat('#',1,length(msg)));
options_.particle.status = false;
estimation(datafile = '../data/ireland2004_data.m', mode_compute = 0, cova_compute = 0, kalman_algo = 5, skewed_kalman_prune_tol = @{PSKF_TOL[TOL]}, skewed_kalman_smoother_skip) ghat rhat pihat;
@{DISTRIBUTION[DIST]}_pskf.fval(@{ME},@{TOL}) = oo_.particle_comparison.fval;
@{DISTRIBUTION[DIST]}_pskf.time(@{ME},@{TOL}) = oo_.particle_comparison.time;
    @#endfor // TOL

    @#for N in 1:length(PARTICLES)
msg = sprintf('# %s (ME=@{MEAS_ERROR[ME]}%%): EVALUATE LOG-LIKELIHOOD VALUE WITH PARTICLE FILTER (N=@{PARTICLES[N]}) #', upper('@{DISTRIBUTION[DIST]}'));
fprintf('\n\n%s\n%s\n%s\n', repmat('#',1,length(msg)), msg, repmat('#',1,length(msg)));
options_.particle.status = true; options_.particle.use_reduced_rank_cholesky = true;
estimation(datafile = '../data/ireland2004_data.m', mode_compute = 0, cova_compute = 0, filter_algorithm = sis, number_of_particles = @{PARTICLES[N]}) ghat rhat pihat;
@{DISTRIBUTION[DIST]}_particle.fval(@{ME},@{N}) = oo_.particle_comparison.fval;
@{DISTRIBUTION[DIST]}_particle.time(@{ME},@{N}) = oo_.particle_comparison.time;
    @#endfor // N

  @#endfor // ME

@#endfor // DIST

rmpath('_overwrite_dynare_particle_filter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate LaTeX table rows %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meas_error_labels = {'40\%', '20\%', '10\%', '5\%'};
panel_labels      = {'A', 'B'};
panel_pskf        = {gaussian_pskf, csn_pskf};
panel_particle    = {gaussian_particle, csn_particle};

msg = sprintf('# Log-Likelihood and CPU Time Comparison #');
fprintf('\n\n%s\n%s\n%s\n', repmat('#',1,length(msg)), msg, repmat('#',1,length(msg)));

for d = 1:2
    if d == 1
        fprintf('\nPanel A: Gaussian Shocks\n');
    else
        fprintf('\n\nPanel B: Skew-Normal (CSN) Shocks\n');
    end
    rows = cell(1, 5);
    % data rows: fval for each measurement error level
    for me = 1:4
        rows{me} = sprintf('%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\', ...
            meas_error_labels{me}, -1*panel_pskf{d}.fval(me,:), -1*panel_particle{d}.fval(me,:));
    end
    % average CPU time row
    time_strs = cell(1, 8);
    for j = 1:4
        avg_t = mean(panel_pskf{d}.time(:,j));
        if avg_t < 0.1,     time_strs{j} = '\(<\)0.1s';
        elseif avg_t < 1,   time_strs{j} = '\(<\)1s';
        elseif avg_t <= 60, time_strs{j} = sprintf('\\(\\approx \\)%ds', round(avg_t));
        else,               time_strs{j} = sprintf('\\(\\approx \\)%dm', round(avg_t/60));
        end
    end
    for j = 1:4
        avg_t = mean(panel_particle{d}.time(:,j));
        if avg_t < 0.1,     time_strs{4+j} = '\(<\)0.1s';
        elseif avg_t < 1,   time_strs{4+j} = '\(<\)1s';
        elseif avg_t <= 60, time_strs{4+j} = sprintf('\\(\\approx \\)%ds', round(avg_t));
        else,               time_strs{4+j} = sprintf('\\(\\approx \\)%dm', round(avg_t/60));
        end
    end
    rows{5} = sprintf('\\textit{Avg. CPU Time} & \\textit{%s} & \\textit{%s} & \\textit{%s} & \\textit{%s} & \\textit{%s} & \\textit{%s} & \\textit{%s} & \\textit{%s} \\\\', time_strs{:});
    % print to console and write to tex file
    fid_tex = fopen(sprintf('../results/ireland2004/online_appendix_table_5%s_%s_%s.tex', panel_labels{d}, ARCH, MATLAB_VERSION), 'w');
    for fid = [1, fid_tex]
        for r = 1:5
            if r < 5
                fprintf(fid, '%s\n', rows{r});
            else
                fprintf(fid, '%s', rows{r});
            end
        end
    end
    fclose(fid_tex);
end

%%%%%%%%%%%%%%%%
% Housekeeping %
%%%%%%%%%%%%%%%%
target_logfile = sprintf('../results/ireland2004/logs/%s_%s_%s.log', M_.fname, ARCH, MATLAB_VERSION);
