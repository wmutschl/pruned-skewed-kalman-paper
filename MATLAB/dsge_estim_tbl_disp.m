function dsge_estim_tbl_disp(xparams, estim_params_names, xstderr, log_likelihood, optim_names)


fprintf('POINT ESTIMATES\n')
disp(array2table([xparams; log_likelihood]...
                ,'RowNames',[estim_params_names;'Log-Lik']...
                ,'VariableNames', optim_names...
                ));
fprintf('STANDARD ERRORS\n')
disp(array2table([xstderr; log_likelihood]...
                ,'RowNames',[erase(estim_params_names,"transformed_");'Log-Lik']...
                ,'VariableNames', optim_names...
                ));
% 
% idx = 1;
% 
% tbl_skew_exo = nan(M_.exo_nbr,3);
% tbl_skew_exo(:,1) = M_.Skew_eta;
% for jexo = 1:estim_params_.nsx
%     tbl_skew_exo(estim_params_.skew_exo(jexo,1),2) = SE(idx);
%     tbl_skew_exo(estim_params_.skew_exo(jexo,1),3) = tbl_skew_exo(estim_params_.skew_exo(jexo,1),1)/SE(idx);
%     idx=idx+1;
% end
% 
% tbl_var_exo = nan(M_.exo_nbr,3);
% tbl_var_exo(:,1) = sqrt(diag(M_.Cov_eta));
% for jexo = 1:estim_params_.nvx
%     tbl_var_exo(estim_params_.var_exo(jexo,1),2) = SE(idx);
%     tbl_var_exo(estim_params_.var_exo(jexo,1),3) = tbl_var_exo(estim_params_.var_exo(jexo,1),1)/SE(idx);
%     idx=idx+1;
% end
% 
% tbl_var_endo = nan(length(M_.varobs),3);
% tbl_var_endo(:,1) = sqrt(diag(M_.Cov_eps));
% for jobs = 1:estim_params_.nvn
%     tbl_var_endo(estim_params_.var_endo(jendo,1),2) = SE(idx);
%     tbl_var_endo(estim_params_.var_endo(jendo,1),3) = tbl_var_endo(estim_params_.var_endo(jendo,1),1)/SE(idx);
%     idx=idx+1;
% end
% 
% tbl_param_vals = nan(M_.param_nbr,3);
% tbl_param_vals(:,1) = M_.params;
% for jp = 1:estim_params_.np
%     tbl_param_vals(estim_params_.param_vals(jp,1),2) = SE(idx);
%     tbl_param_vals(estim_params_.param_vals(jp,1),3) = tbl_param_vals(estim_params_.param_vals(jp,1),1)/SE(idx);
%     idx=idx+1;
% end
% 
% tbl = array2table([tbl_skew_exo;tbl_var_exo;tbl_var_endo;tbl_param_vals],...
%                   'RowNames',["skew("+M_.exo_names+")"; "stderr("+M_.exo_names+")"; "stderr("+M_.varobs+")"; M_.param_names],...
%                   'VariableNames',["estimate","stderr","t-stat"]);
