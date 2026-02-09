function update_table_2(results_folder, model_type, oo_, ARCH, MATLAB_VERSION)
% update_table_2 - Update columns in ireland2004 table2.tex with estimation results
% =========================================================================
% USAGE:
%   update_ireland2004_table2(tex_file, columns_to_update, oo_, log_likelihood)
%   update_ireland2004_table2(..., 'Precision', 4)
%
% INPUTS:
%   tex_file          - path to the tex file (e.g., '../results/ireland2004/table2.tex')
%   columns_to_update - cell array specifying which columns to update, e.g.:
%                       {'Mode', 'Std-dev'} for ML columns 2-3 (Gaussian ML)
%                       or column indices [2, 3] for Gaussian ML, [4, 5] for CSN ML
%                       For ML results: use [2,3] for Gaussian, [4,5] for CSN
%                       For Bayesian:   use [6,7,8] for Gaussian, [9,10,11] for CSN
%   oo_               - Dynare's oo_ structure containing estimation results
%   log_likelihood    - final log-likelihood value for Obj(mode) row
%
% OPTIONAL INPUTS (name-value pairs):
%   'Precision'       - number of decimal places (default: 4)
%   'HasSkewness'     - whether model has skewness parameters (default: auto-detect)
%
% EXAMPLE:
%   % For Gaussian ML (columns 2-3):
%   update_ireland2004_table2('../results/ireland2004/table2.tex', [2,3], oo_, loglik);
%
%   % For CSN ML (columns 4-5):
%   update_ireland2004_table2('../results/ireland2004/table2.tex', [4,5], oo_, loglik);
%
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

% Define parameter mapping from Dynare names to LaTeX symbols
param_map = {
    'OMEGA',  '\({\omega}\)'
    'RHO_PI', '\(\rho_{\pi}\)'
    'RHO_G',  '\(\rho_g\)'
    'RHO_X',  '\(\rho_x\)'
    'RHO_A',  '\(\rho_a\)'
    'RHO_E',  '\(\rho_e\)'
    'eta_a',  '\(stderr(\eta_a)\)'
    'eta_e',  '\(stderr(\eta_e)\)'
    'eta_z',  '\(stderr(\eta_z)\)'
    'eta_r',  '\(stderr(\eta_r)\)'
    'eta_a',  '\(skew(\eta_a)\)'
    'eta_e',  '\(skew(\eta_e)\)'
    'eta_z',  '\(skew(\eta_z)\)'
    'eta_r',  '\(skew(\eta_r)\)'
};

% Extract estimation results from oo_ structure
param_names = {};
bayes_mean = [];
bayes_mode = [];
bayes_hpdinf = [];
bayes_hpdsup = [];
ml_mode = [];
ml_std = [];

if strcmp(model_type, 'ML_GAUSSIAN') || strcmp(model_type, 'ML_CSN')
    % MLE results
    if isfield(oo_.mle_mode, 'parameters')
        pnames = fieldnames(oo_.mle_mode.parameters);
        for j = 1:length(pnames)
            param_names{end+1} = pnames{j};
            ml_mode(end+1) = oo_.mle_mode.parameters.(pnames{j});
            ml_std(end+1) = oo_.mle_std_at_mode.parameters.(pnames{j});
        end
    end
    if isfield(oo_.mle_mode, 'shocks_std')
        snames = fieldnames(oo_.mle_mode.shocks_std);
        for j = 1:length(snames)
            param_names{end+1} = ['stderr_' snames{j}];
            ml_mode(end+1) = oo_.mle_mode.shocks_std.(snames{j});
            ml_std(end+1) = oo_.mle_std_at_mode.shocks_std.(snames{j});
        end
    end
    if isfield(oo_.mle_mode, 'shocks_skew')
        sknames = fieldnames(oo_.mle_mode.shocks_skew);
        for j = 1:length(sknames)
            param_names{end+1} = ['skew_' sknames{j}];
            ml_mode(end+1) = oo_.mle_mode.shocks_skew.(sknames{j});
        end
    end
elseif strcmp(model_type, 'BAYESIAN_GAUSSIAN') || strcmp(model_type, 'BAYESIAN_CSN')
    % Bayesian results
    if isfield(oo_.posterior_mode, 'parameters')
        pnames = fieldnames(oo_.posterior_mode.parameters);
        for j = 1:length(pnames)
            param_names{end+1} = pnames{j};
            bayes_mean(end+1) = oo_.posterior_mean.parameters.(pnames{j});
            bayes_mode(end+1) = oo_.posterior_mode.parameters.(pnames{j});
            bayes_hpdinf(end+1) = oo_.posterior_hpdinf.parameters.(pnames{j});
            bayes_hpdsup(end+1) = oo_.posterior_hpdsup.parameters.(pnames{j});
        end
    end
    if isfield(oo_.posterior_mode, 'shocks_std')
        snames = fieldnames(oo_.posterior_mode.shocks_std);
        for j = 1:length(snames)
            param_names{end+1} = ['stderr_' snames{j}];
            bayes_mean(end+1) = oo_.posterior_mean.shocks_std.(snames{j});
            bayes_mode(end+1) = oo_.posterior_mode.shocks_std.(snames{j});
            bayes_hpdinf(end+1) = oo_.posterior_hpdinf.shocks_std.(snames{j});
            bayes_hpdsup(end+1) = oo_.posterior_hpdsup.shocks_std.(snames{j});
        end
    end
    if isfield(oo_.posterior_mode, 'shocks_skew')
        sknames = fieldnames(oo_.posterior_mode.shocks_skew);
        for j = 1:length(sknames)
            param_names{end+1} = ['skew_' sknames{j}];
            bayes_mean(end+1) = oo_.posterior_mean.shocks_skew.(sknames{j});
            bayes_mode(end+1) = oo_.posterior_mode.shocks_skew.(sknames{j});
            bayes_hpdinf(end+1) = oo_.posterior_hpdinf.shocks_skew.(sknames{j});
            bayes_hpdsup(end+1) = oo_.posterior_hpdsup.shocks_skew.(sknames{j});
        end
    end
end

% Read the tex file
tex_file = sprintf('%s/table_2_%s_%s.tex', results_folder, ARCH, MATLAB_VERSION);
fid = fopen(tex_file, 'r');
if fid == -1
    warning('Cannot open file %s, so I make a copy of latex2.tex first.', tex_file);
    copyfile([results_folder '/table_2.tex'], tex_file);
    fid = fopen(tex_file, 'r');
end
lines = {};
while ~feof(fid)
    lines{end+1} = fgetl(fid); %#ok<AGROW>
end
fclose(fid);

% Process each line
for i = 1:length(lines)
    line = lines{i};

    % Skip non-data lines (header, midrule, etc.)
    if isempty(line) || ~contains(line, '&') || contains(line, 'Parameter')
        continue;
    end

    parts = strsplit(line, '&');
    % Check if this is the Obj(mode) line
    if contains(line, 'Obj(mode)')
        if strcmp(model_type, 'ML_GAUSSIAN')
            multicolumn_idx = 2; span = 2;
        elseif strcmp(model_type, 'ML_CSN')
            multicolumn_idx = 3; span = 2;
        elseif strcmp(model_type, 'BAYESIAN_GAUSSIAN')
            multicolumn_idx = 4; span = 3;
        elseif strcmp(model_type, 'BAYESIAN_CSN')
            multicolumn_idx = 5; span = 3;
        end
        objval = sprintf('%.2f', oo_.posterior.optimization.log_density);
        parts{multicolumn_idx} = sprintf(' \\multicolumn{%d}{c}{%s} ', span, objval);
        lines{i} = strjoin(parts, '&');
        continue;
    end

    param_latex = strtrim(parts{1});
    param_idx = find(ismember(param_map(:,2), param_latex),1);
    param_dynare = param_map{param_idx,1};
    if strcmp(model_type, 'ML_GAUSSIAN')
        if contains(param_latex, 'skew')
            parts{2} = sprintf(' --- ');
            parts{3} = sprintf(' --- ');
        elseif contains(param_latex, 'stderr')
            parts{2} = sprintf(' %.4f ', oo_.mle_mode.shocks_std.(param_dynare));
            parts{3} = sprintf(' %.4f ', oo_.mle_std_at_mode.shocks_std.(param_dynare));
        else
            parts{2} = sprintf(' %.4f ', oo_.mle_mode.parameters.(param_dynare));
            parts{3} = sprintf(' %.4f ', oo_.mle_std_at_mode.parameters.(param_dynare));
        end
    elseif strcmp(model_type, 'ML_CSN')
        if contains(param_latex, 'skew')
            parts{4} = sprintf(' %.4f ', oo_.mle_mode.shocks_skew.(param_dynare));
            parts{5} = sprintf(' %.4f ', oo_.mle_std_at_mode.shocks_skew.(param_dynare));
        elseif contains(param_latex, 'stderr')
            parts{4} = sprintf(' %.4f ', oo_.mle_mode.shocks_std.(param_dynare));
            parts{5} = sprintf(' %.4f ', oo_.mle_std_at_mode.shocks_std.(param_dynare));
        else
            parts{4} = sprintf(' %.4f ', oo_.mle_mode.parameters.(param_dynare));
            parts{5} = sprintf(' %.4f ', oo_.mle_std_at_mode.parameters.(param_dynare));
        end
    elseif strcmp(model_type, 'BAYESIAN_GAUSSIAN')
        if contains(param_latex, 'skew')
            parts{6} = sprintf(' --- ');
            parts{7} = sprintf(' --- ');
            parts{8} = sprintf(' --- ');
        elseif contains(param_latex, 'stderr')
            parts{6} = sprintf(' %.4f ', oo_.posterior_mean.shocks_std.(param_dynare));
            parts{7} = sprintf(' %.4f ', oo_.posterior_mode.shocks_std.(param_dynare));
            parts{8} = sprintf(' [%.2f;%.2f] ', oo_.posterior_hpdinf.shocks_std.(param_dynare), oo_.posterior_hpdsup.shocks_std.(param_dynare));
        else
            parts{6} = sprintf(' %.4f ', oo_.posterior_mean.parameters.(param_dynare));
            parts{7} = sprintf(' %.4f ', oo_.posterior_mode.parameters.(param_dynare));
            parts{8} = sprintf(' [%.2f;%.2f] ', oo_.posterior_hpdinf.parameters.(param_dynare), oo_.posterior_hpdsup.parameters.(param_dynare));
        end
    elseif strcmp(model_type, 'BAYESIAN_CSN')
        if contains(param_latex, 'skew')
            parts{9} = sprintf(' %.4f ', oo_.posterior_mean.shocks_skew.(param_dynare));
            parts{10} = sprintf(' %.4f ', oo_.posterior_mode.shocks_skew.(param_dynare));
            parts{11} = sprintf(' [%.2f;%.2f] ', oo_.posterior_hpdinf.shocks_skew.(param_dynare), oo_.posterior_hpdsup.shocks_skew.(param_dynare));
        elseif contains(param_latex, 'stderr')
            parts{9} = sprintf(' %.4f ', oo_.posterior_mean.shocks_std.(param_dynare));
            parts{10} = sprintf(' %.4f ', oo_.posterior_mode.shocks_std.(param_dynare));
            parts{11} = sprintf(' [%.2f;%.2f] ', oo_.posterior_hpdinf.shocks_std.(param_dynare), oo_.posterior_hpdsup.shocks_std.(param_dynare));
        else
            parts{9} = sprintf(' %.4f ', oo_.posterior_mean.parameters.(param_dynare));
            parts{10} = sprintf(' %.4f ', oo_.posterior_mode.parameters.(param_dynare));
            parts{11} = sprintf(' [%.2f;%.2f] ', oo_.posterior_hpdinf.parameters.(param_dynare), oo_.posterior_hpdsup.parameters.(param_dynare));
        end
    end
    lines{i} = strjoin(parts, '&');

end

% Write the updated tex file
fid = fopen(tex_file, 'w');
if fid == -1
    error('Cannot write to file: %s', tex_file);
end
for i = 1:length(lines)
    if i < length(lines)
        fprintf(fid, '%s\n', lines{i});
    else
        fprintf(fid, '%s', lines{i}); % no newline at end
    end
end
fclose(fid);

end
