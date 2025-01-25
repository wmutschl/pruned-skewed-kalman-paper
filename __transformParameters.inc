% add columns for transform parameters
estim_params_.var_exo(:,11) = zeros(size(estim_params_.var_exo,1),1);
estim_params_.var_endo(:,11) = zeros(size(estim_params_.var_endo,1),1);
estim_params_.corrx(:,12) = zeros(size(estim_params_.corrx,1),1);
estim_params_.corrn(:,12) = zeros(size(estim_params_.corrn,1),1);
estim_params_.param_vals(:,11) = zeros(size(estim_params_.param_vals,1),1);

@#if TRANSFORM_PARAMETERS == 1
% optional new feature to use logit to transform bounded to unbounded parameters during optimization only
estim_params_.var_exo(estim_params_.var_exo(:,1)==find(ismember(M_.exo_names,'eta_a')),11) = 1;
estim_params_.var_exo(estim_params_.var_exo(:,1)==find(ismember(M_.exo_names,'eta_e')),11) = 1;
estim_params_.var_exo(estim_params_.var_exo(:,1)==find(ismember(M_.exo_names,'eta_z')),11) = 1;
estim_params_.var_exo(estim_params_.var_exo(:,1)==find(ismember(M_.exo_names,'eta_r')),11) = 1;
estim_params_.skew_exo(estim_params_.skew_exo(:,1)==find(ismember(M_.exo_names,'eta_a')),11) = 1;
estim_params_.skew_exo(estim_params_.skew_exo(:,1)==find(ismember(M_.exo_names,'eta_e')),11) = 1;
estim_params_.skew_exo(estim_params_.skew_exo(:,1)==find(ismember(M_.exo_names,'eta_z')),11) = 1;
estim_params_.skew_exo(estim_params_.skew_exo(:,1)==find(ismember(M_.exo_names,'eta_r')),11) = 1;
estim_params_.param_vals(estim_params_.param_vals(:,1)==find(ismember(M_.param_names,'OMEGA')),11) = 1;
estim_params_.param_vals(estim_params_.param_vals(:,1)==find(ismember(M_.param_names,'ALPHA_X')),11) = 1;
estim_params_.param_vals(estim_params_.param_vals(:,1)==find(ismember(M_.param_names,'ALPHA_PI')),11) = 1;
@#endif