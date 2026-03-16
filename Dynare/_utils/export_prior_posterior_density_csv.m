function density_table = export_prior_posterior_density_csv(oo_, DISTRIB, ARCH, MATLAB_VERSION)

    density_table = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), double.empty(0,1), double.empty(0,1), ...
        'VariableNames', {'distribution', 'type', 'parameter', 'x', 'density'});

    sources = struct('posterior', oo_.posterior_density, 'prior', oo_.prior_density);
    source_names = fieldnames(sources);

    for sd = 1:length(source_names)
        dist_label = string(source_names{sd});
        src = sources.(source_names{sd});

        % Parameters
        if isfield(src, 'parameters')
            param_names = fieldnames(src.parameters);
            for jp = 1:length(param_names)
                pname = param_names{jp};
                density = src.parameters.(pname);
                n_rows = size(density, 1);
                density_table = [density_table;
                    table(repmat(dist_label, n_rows, 1), ...
                          repmat("deep", n_rows, 1), ...
                          repmat(string(pname), n_rows, 1), ...
                          density(:,1), density(:,2), ...
                        'VariableNames', {'distribution', 'type', 'parameter', 'x', 'density'})
                ]; %#ok<AGROW>
            end
        end

        % Shock standard deviations
        if isfield(src, 'shocks_std')
            shock_names = fieldnames(src.shocks_std);
            for js = 1:length(shock_names)
                sname = shock_names{js};
                density = src.shocks_std.(sname);
                n_rows = size(density, 1);
                density_table = [density_table;
                    table(repmat(dist_label, n_rows, 1), ...
                          repmat("stderr", n_rows, 1), ...
                          repmat(string(sname), n_rows, 1), ...
                          density(:,1), density(:,2), ...
                        'VariableNames', {'distribution', 'type', 'parameter', 'x', 'density'})
                ]; %#ok<AGROW>
            end
        end

        % Shock skewnesses
        if isfield(src, 'shocks_skew')
            shock_names = fieldnames(src.shocks_skew);
            for js = 1:length(shock_names)
                sname = shock_names{js};
                density = src.shocks_skew.(sname);
                n_rows = size(density, 1);
                density_table = [density_table;
                    table(repmat(dist_label, n_rows, 1), ...
                          repmat("skew", n_rows, 1), ...
                          repmat(string(sname), n_rows, 1), ...
                          density(:,1), density(:,2), ...
                        'VariableNames', {'distribution', 'type', 'parameter', 'x', 'density'})
                ]; %#ok<AGROW>
            end
        end
    end

    density_file = sprintf('../results/ireland2004/density_%s_%s_%s.csv', DISTRIB, ARCH, MATLAB_VERSION);
    writetable(density_table, density_file);
    fprintf('Density data saved to: %s\n', density_file);

end
