function update_montecarlo_table(results_folder, dist_label, true_values, param_names, sample_lengths, ARCH, MATLAB_VERSION, varargin)
% update_montecarlo_table - Create/update LaTeX table with Monte Carlo estimation results
% =========================================================================
% Creates a table similar to Table in online_appendix.tex for DGP (3),
% showing mean, [p5;p95] percentiles, and {NRMSE} for each parameter.
%
% INPUTS:
%   results_folder  - path to results folder (e.g., '../results/ireland2004/montecarlo')
%   dist_label      - 'csn' or 'gaussian'
%   true_values     - vector of true parameter values (nparam x 1)
%   param_names     - cell array of LaTeX parameter names (nparam x 1)
%   sample_lengths  - vector of sample lengths (e.g., [200, 500])
%   ARCH            - architecture string for filename
%   MATLAB_VERSION  - MATLAB version string for filename
%   varargin        - one matrix (nparam x ndatasets) per sample length
%
% The function creates a tex file with columns for CSN and Gaussian results.
% Each run (CSN or Gaussian) updates its own columns while preserving the others.
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

nT = length(sample_lengths);
nparam_csn = 8; % 4 stderr + 4 skew (full CSN model)
tex_file = sprintf('%s/montecarlo_table_%s_%s.tex', results_folder, ARCH, MATLAB_VERSION);

% Define all parameter rows (always 8 rows for both stderr and skew)
all_param_names = {
    '\(stderr(\eta_a)\)'
    '\(stderr(\eta_e)\)'
    '\(stderr(\eta_z)\)'
    '\(stderr(\eta_r)\)'
    '\(skew(\eta_a)\)'
    '\(skew(\eta_e)\)'
    '\(skew(\eta_z)\)'
    '\(skew(\eta_r)\)'
};
all_true_values = [2.5232; 0.0212; 0.7900; 0.2838; 0; -0.40; -0.95; +0.55];

% Parts layout: Param & Truth & CSN_T1 & CSN_T2 & Gaussian_T1 & Gaussian_T2
% CSN columns: parts{3,4}  |  Gaussian columns: parts{5,6}
nparts = 2 + 2*nT; % Param(1) + Truth(2) + CSN(nT) + Gaussian(nT)

% Try to load existing file, otherwise create template
if exist(tex_file, 'file')
    fid = fopen(tex_file, 'r');
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);
else
    % Create template
    lines = {};
    for jp = 1:nparam_csn
        if all_true_values(jp) == 0
            truth_str = sprintf('%.2f', all_true_values(jp));
        else
            truth_str = sprintf('%.4f', all_true_values(jp));
        end
        parts = cell(1, nparts);
        parts{1} = all_param_names{jp};
        parts{2} = [' ' truth_str];
        for jc = 3:nparts
            parts{jc} = ' ---';
        end
        lines{end+1} = [strjoin(parts, '&') ' \\']; %#ok<AGROW>
    end
end

% Determine which columns to update
if strcmp(dist_label, 'csn')
    col_offset = 2; % columns 3,4 (1-indexed in parts)
else
    col_offset = 2 + nT; % columns 5,6
end

% Update the table
nparam = length(true_values);
for i = 1:length(lines)
    line = lines{i};
    if isempty(line) || ~contains(line, '&')
        continue;
    end

    parts = strsplit(line, '&');
    if length(parts) < nparts
        continue;
    end

    % Preserve line ending (\\ and optional \midrule/\bottomrule) from last part
    line_ending = '';
    bs_pos = strfind(parts{end}, '\\');
    if ~isempty(bs_pos)
        line_ending = parts{end}(bs_pos(end):end);
        parts{end} = parts{end}(1:bs_pos(end)-1);
    end

    param_latex = strtrim(parts{1});

    % Find which parameter this row corresponds to
    row_idx = find(strcmp(all_param_names, param_latex), 1);
    if isempty(row_idx)
        continue;
    end

    % Check if this parameter is estimated in the current run
    % For Gaussian: only stderr params (first 4), skew rows get '---'
    if strcmp(dist_label, 'gaussian') && row_idx > 4
        for jT = 1:nT
            parts{col_offset + jT} = ' ---';
        end
    else
        % Map row_idx to the parameter index in xparam_results
        if strcmp(dist_label, 'gaussian')
            jparam = row_idx; % only stderr params, indices 1-4
        else
            jparam = row_idx; % all 8 params
        end

        if jparam <= nparam
            for jT = 1:nT
                vals = varargin{jT}(jparam, :);
                avg = mean(vals);
                p5 = prctile(vals, 5);
                p95 = prctile(vals, 95);
                tv = all_true_values(row_idx);
                if tv ~= 0
                    nrmse = sqrt(mean((vals - tv).^2)) / abs(tv);
                    nrmse_str = sprintf('%.3f', nrmse);
                else
                    nrmse_str = '~';
                end
                cell_str = sprintf(' \\(\\underset{\\{%s\\}}{\\underset{[%.2f;%.2f]}{%.4f}}\\)', ...
                    nrmse_str, p5, p95, avg);
                parts{col_offset + jT} = cell_str;
            end
        end
    end

    lines{i} = [strjoin(parts, '&') line_ending];
end

% Ensure the last line ends with \\ \bottomrule
lastline = lines{end};
lastline = regexprep(lastline, '\s*\\\\?\s*(\\bottomrule)?\s*%?\s*$', '');
lines{end} = [lastline ' \\ \bottomrule'];

% Write the updated tex file
fid = fopen(tex_file, 'w');
if fid == -1
    error('Cannot write to file: %s', tex_file);
end
for i = 1:length(lines)
    if i < length(lines)
        fprintf(fid, '%s\n', lines{i});
    else
        fprintf(fid, '%s', lines{i});
    end
end
fclose(fid);

fprintf('Monte Carlo table updated: %s (columns for %s)\n', tex_file, upper(dist_label));

end
