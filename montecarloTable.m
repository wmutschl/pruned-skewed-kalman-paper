%% Monte Carlo Table Generator
% Generates a LaTeX table for the Monte Carlo simulation results
% comparing T=100 and T=500 sample sizes

clear; clc;

% Load the estimation results
load('ireland2004_monte_carlo_xparam_opt.mat', 'xparam1_opt_T100', 'xparam1_opt_T500');
abs(xparam1_opt_T100(5:8,:) > 0.8)
% True parameter values (in the order they appear in the estimation)
param_names = {
    'stderr($\eta_a$)'
    'stderr($\eta_e$)'
    'stderr($\eta_z$)'
    'stderr($\eta_r$)'
    'skew($\eta_a$)'
    'skew($\eta_e$)'
    'skew($\eta_z$)'
    'skew($\eta_r$)'
};

true_values = [
    3.0167   % stderr eta_a
    0.0248   % stderr eta_e
    0.8865   % stderr eta_z
    0.2790   % stderr eta_r
   -0.3      % skew eta_a
    0        % skew eta_e
   -0.5      % skew eta_z
    0.8      % skew eta_r
];

n_params = length(true_values);
R = 1200;  % number of replications

% Compute statistics for T=100
mean_T100 = mean(xparam1_opt_T100, 2);
p05_T100 = prctile(xparam1_opt_T100, 5, 2);
p95_T100 = prctile(xparam1_opt_T100, 95, 2);
rmse_T100 = sqrt(mean((xparam1_opt_T100 - true_values).^2, 2));
nrmse_T100 = rmse_T100 ./ abs(true_values);  % normalized RMSE

% Compute statistics for T=500
mean_T500 = mean(xparam1_opt_T500, 2);
p05_T500 = prctile(xparam1_opt_T500, 5, 2);
p95_T500 = prctile(xparam1_opt_T500, 95, 2);
rmse_T500 = sqrt(mean((xparam1_opt_T500 - true_values).^2, 2));
nrmse_T500 = rmse_T500 ./ abs(true_values);  % normalized RMSE

%% Generate LaTeX table
fid = fopen('monte_carlo_table.tex', 'w');

% Table header
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\caption{Distribution of parameter estimates}\\label{tbl:MonteCarloEstimation}\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{l c c c}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '      &        & \\multicolumn{2}{c}{\\emph{Pruned Skewed KF ($10^{-2}$)}} \\\\\\cmidrule(lr){3-4}\n');
fprintf(fid, 'Param & Truth  & $T=100$ & $T=500$ \\\\\n');
fprintf(fid, '\\midrule\n');

% Generate rows for each parameter
for j = 1:n_params
    % Parameter name and true value
    fprintf(fid, '%s & %.4f\n', param_names{j}, true_values(j));
    
    % T=100 cell
    if true_values(j) == 0
        % NRMSE not defined for true value = 0
        fprintf(fid, '& $\\underset{\\{~\\}}{\\underset{[%.2f;%.2f]}{%.3f}}$\n', ...
            p05_T100(j), p95_T100(j), mean_T100(j));
    else
        fprintf(fid, '& $\\underset{\\{%.3f\\}}{\\underset{[%.2f;%.2f]}{%.3f}}$\n', ...
            nrmse_T100(j), p05_T100(j), p95_T100(j), mean_T100(j));
    end
    
    % T=500 cell
    if true_values(j) == 0
        % NRMSE not defined for true value = 0
        fprintf(fid, '& $\\underset{\\{~\\}}{\\underset{[%.2f;%.2f]}{%.3f}}$\n', ...
            p05_T500(j), p95_T500(j), mean_T500(j));
    else
        fprintf(fid, '& $\\underset{\\{%.3f\\}}{\\underset{[%.2f;%.2f]}{%.3f}}$\n', ...
            nrmse_T500(j), p05_T500(j), p95_T500(j), mean_T500(j));
    end
    
    fprintf(fid, '\\\\\n');
end

% Table footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\begin{tablenotes}\n');
fprintf(fid, '\\small\n');
fprintf(fid, 'Note: Cells contain average on top, [5,95] percentiles in square brackets, and \\{NRMSE\\} in curly brackets.\n');
fprintf(fid, 'For $\\gamma_{\\eta_e}$, NRMSE is not defined due to division by zero.\n');
fprintf(fid, '\\end{tablenotes}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

fprintf('LaTeX table saved to monte_carlo_table.tex\n');

%% Also display summary statistics in console
fprintf('\n=== Monte Carlo Results (R=%d replications) ===\n\n', R);
fprintf('%-20s %10s %12s %12s %12s %12s\n', 'Parameter', 'Truth', 'Mean T100', 'Mean T500', 'NRMSE T100', 'NRMSE T500');
fprintf('%s\n', repmat('-', 1, 80));
for j = 1:n_params
    if true_values(j) == 0
        fprintf('%-20s %10.4f %12.4f %12.4f %12s %12s\n', ...
            param_names{j}, true_values(j), mean_T100(j), mean_T500(j), 'N/A', 'N/A');
    else
        fprintf('%-20s %10.4f %12.4f %12.4f %12.4f %12.4f\n', ...
            param_names{j}, true_values(j), mean_T100(j), mean_T500(j), nrmse_T100(j), nrmse_T500(j));
    end
end