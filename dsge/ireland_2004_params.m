function [xparam, PARAM, ESTIM_PARAM, MODEL, OPT] = ireland_2004_params(stage, MODEL, OPT, xparam_model, xparam_sqrt_diag_Sigma_eta, xparam_diag_Gamma_eta, xparam_stderr_eta, xparam_skew_eta)

%% initialization
PARAM = []; ESTIM_PARAM = [];
MODEL.idx_eta_a = find(ismember(MODEL.exo_names(:,1),"eta_a"),1);
MODEL.idx_eta_e = find(ismember(MODEL.exo_names(:,1),"eta_e"),1);
MODEL.idx_eta_z = find(ismember(MODEL.exo_names(:,1),"eta_z"),1);
MODEL.idx_eta_r = find(ismember(MODEL.exo_names(:,1),"eta_r"),1);

%% calibrated parameters
PARAM.BETA      = 0.99;   % discount factor
PARAM.PSI       = 0.1;    % output gap slope in Phillips curve
PARAM.ALPHA_PI  = 0.0001; % slope parameter in Phillips curve
%if OPT.datafile == "data_post_1980"
%    PARAM.ALPHA_X  = 0.0001; % slope parameter in IS curve
%end
if stage == 0
    MODEL.Sigma_eps = zeros(MODEL.varobs_nbr,MODEL.varobs_nbr);
    MODEL.Gamma_eta = zeros(MODEL.exo_nbr,MODEL.exo_nbr);
elseif stage == 1
    MODEL.Sigma_eps = zeros(MODEL.varobs_nbr,MODEL.varobs_nbr);
    PARAM.OMEGA     = xparam_model(ismember(MODEL.param_estim_names,'OMEGA'));   % scale parameter preference innovation in IS curve
    %if (OPT.datafile == "data_post_1980") == 0
    PARAM.ALPHA_X   = xparam_model(ismember(MODEL.param_estim_names,'ALPHA_X')); % slope parameter in IS curve
    %end
    PARAM.RHO_PI    = xparam_model(ismember(MODEL.param_estim_names,'RHO_PI'));  % feedback policy rule inflation
    PARAM.RHO_G     = xparam_model(ismember(MODEL.param_estim_names,'RHO_G'));   % feedback policy rule output growth
    PARAM.RHO_X     = xparam_model(ismember(MODEL.param_estim_names,'RHO_X'));   % feedback policy rule output gap
    PARAM.RHO_A     = xparam_model(ismember(MODEL.param_estim_names,'RHO_A'));   % persistence preference shock
    PARAM.RHO_E     = xparam_model(ismember(MODEL.param_estim_names,'RHO_E'));   % persistence cost-push shock
elseif stage == 2
    MODEL.Sigma_eps = zeros(MODEL.varobs_nbr,MODEL.varobs_nbr);
end

%% estimated parameters
% ESTIM_PARAM.PARAMNAME = {INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND}; % comment

if stage == 0
    if OPT.datafile == "data_full_sample"
        if OPT.use_stderr_skew_transform
            ESTIM_PARAM.stderr_eta_a = {0.0405, 0, 1};
            ESTIM_PARAM.stderr_eta_e = {0.0012, 0, 1};
            ESTIM_PARAM.stderr_eta_z = {0.0109, 0, 1};
            ESTIM_PARAM.stderr_eta_r = {0.0031, 0, 1};
        else
            ESTIM_PARAM.sqrt_diag_Sigma_eta_a = {0.0405, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_e = {0.0012, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_z = {0.0109, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_r = {0.0031, 0, 1};
        end
        ESTIM_PARAM.OMEGA   = {0.0617, 0, 1};
        ESTIM_PARAM.ALPHA_X = {0.0836, 0, 1};
        ESTIM_PARAM.RHO_PI  = {0.3597, 0, 1};
        ESTIM_PARAM.RHO_G   = {0.2536, 0, 1};
        ESTIM_PARAM.RHO_X   = {0.0347, 0, 1};
        ESTIM_PARAM.RHO_A   = {0.9470, 0, 1};
        ESTIM_PARAM.RHO_E   = {0.9625, 0, 1};
    elseif OPT.datafile == "data_post_1980"
        if OPT.use_stderr_skew_transform
            ESTIM_PARAM.stderr_eta_a = {0.0302, 0, 1};
            ESTIM_PARAM.stderr_eta_e = {0.0002, 0, 1};
            ESTIM_PARAM.stderr_eta_z = {0.0089, 0, 1};
            ESTIM_PARAM.stderr_eta_r = {0.0028, 0, 1};
        else
            ESTIM_PARAM.sqrt_diag_Sigma_eta_a = {0.0302, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_e = {0.0002, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_z = {0.0089, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_r = {0.0028, 0, 1};
        end
        ESTIM_PARAM.OMEGA   = {0.0581, 0, 1};
        ESTIM_PARAM.ALPHA_X = {0.0000, 0, 1};
        ESTIM_PARAM.RHO_PI  = {0.3866, 0, 1};
        ESTIM_PARAM.RHO_G   = {0.3960, 0, 1};
        ESTIM_PARAM.RHO_X   = {0.1654, 0, 1};
        ESTIM_PARAM.RHO_A   = {0.9048, 0, 1};
        ESTIM_PARAM.RHO_E   = {0.9907, 0, 1};
    elseif OPT.datafile == "data_pre_1980"
        if OPT.use_stderr_skew_transform
            ESTIM_PARAM.stderr_eta_a = {0.1538, 0, 1};
            ESTIM_PARAM.stderr_eta_e = {0.0035, 0, 1};
            ESTIM_PARAM.stderr_eta_z = {0.0104, 0, 1};
            ESTIM_PARAM.stderr_eta_r = {0.0033, 0, 1};
        else
            ESTIM_PARAM.sqrt_diag_Sigma_eta_a = {0.1538, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_e = {0.0035, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_z = {0.0104, 0, 1};
            ESTIM_PARAM.sqrt_diag_Sigma_eta_r = {0.0033, 0, 1};
        end
        ESTIM_PARAM.OMEGA   = {0.0000, 0, 1};
        ESTIM_PARAM.ALPHA_X = {0.2028, 0, 1};
        ESTIM_PARAM.RHO_PI  = {0.3053, 0, 1};
        ESTIM_PARAM.RHO_G   = {0.2365, 0, 1};
        ESTIM_PARAM.RHO_X   = {0.0000, 0, 1};
        ESTIM_PARAM.RHO_A   = {0.9910, 0, 1};
        ESTIM_PARAM.RHO_E   = {0.5439, 0, 1};
    end
elseif stage == 1
    if OPT.use_stderr_skew_transform
        ESTIM_PARAM.skew_eta_a   = {xparam_skew_eta(MODEL.idx_eta_a), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_e   = {xparam_skew_eta(MODEL.idx_eta_e), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_z   = {xparam_skew_eta(MODEL.idx_eta_z), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_r   = {xparam_skew_eta(MODEL.idx_eta_r), -0.95, 0.95};
        ESTIM_PARAM.stderr_eta_a = {xparam_stderr_eta(MODEL.idx_eta_a), 0, 1};
        ESTIM_PARAM.stderr_eta_e = {xparam_stderr_eta(MODEL.idx_eta_e), 0, 1};
        ESTIM_PARAM.stderr_eta_z = {xparam_stderr_eta(MODEL.idx_eta_z), 0, 1};
        ESTIM_PARAM.stderr_eta_r = {xparam_stderr_eta(MODEL.idx_eta_r), 0, 1};
    else
        ESTIM_PARAM.diag_Gamma_eta_a      = {xparam_diag_Gamma_eta(MODEL.idx_eta_a), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_e      = {xparam_diag_Gamma_eta(MODEL.idx_eta_e), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_z      = {xparam_diag_Gamma_eta(MODEL.idx_eta_z), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_r      = {xparam_diag_Gamma_eta(MODEL.idx_eta_r), -Inf, Inf};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_a = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_a), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_e = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_e), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_z = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_z), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_r = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_r), 0, 1};
    end    
elseif stage == 2
    if OPT.use_stderr_skew_transform
        ESTIM_PARAM.skew_eta_a   = {xparam_skew_eta(MODEL.idx_eta_a), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_e   = {xparam_skew_eta(MODEL.idx_eta_e), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_z   = {xparam_skew_eta(MODEL.idx_eta_z), -0.95, 0.95};
        ESTIM_PARAM.skew_eta_r   = {xparam_skew_eta(MODEL.idx_eta_r), -0.95, 0.95};
        ESTIM_PARAM.stderr_eta_a = {xparam_stderr_eta(MODEL.idx_eta_a), 0, 1};
        ESTIM_PARAM.stderr_eta_e = {xparam_stderr_eta(MODEL.idx_eta_e), 0, 1};
        ESTIM_PARAM.stderr_eta_z = {xparam_stderr_eta(MODEL.idx_eta_z), 0, 1};
        ESTIM_PARAM.stderr_eta_r = {xparam_stderr_eta(MODEL.idx_eta_r), 0, 1};
    else
        ESTIM_PARAM.diag_Gamma_eta_a      = {xparam_diag_Gamma_eta(MODEL.idx_eta_a), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_e      = {xparam_diag_Gamma_eta(MODEL.idx_eta_e), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_z      = {xparam_diag_Gamma_eta(MODEL.idx_eta_z), -Inf, Inf};
        ESTIM_PARAM.diag_Gamma_eta_r      = {xparam_diag_Gamma_eta(MODEL.idx_eta_r), -Inf, Inf};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_a = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_a), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_e = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_e), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_z = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_z), 0, 1};
        ESTIM_PARAM.sqrt_diag_Sigma_eta_r = {xparam_sqrt_diag_Sigma_eta(MODEL.idx_eta_r), 0, 1};
    end
    ESTIM_PARAM.OMEGA   = {xparam_model(ismember(MODEL.param_estim_names,'OMEGA')), 0, 1};
    %if (OPT.datafile == "data_post_1980") == 0
    ESTIM_PARAM.ALPHA_X = {xparam_model(ismember(MODEL.param_estim_names,'ALPHA_X')), 0, 1};
    %end
    ESTIM_PARAM.RHO_PI  = {xparam_model(ismember(MODEL.param_estim_names,'RHO_PI')), 0, 1};
    ESTIM_PARAM.RHO_G   = {xparam_model(ismember(MODEL.param_estim_names,'RHO_G')), 0, 1};
    ESTIM_PARAM.RHO_X   = {xparam_model(ismember(MODEL.param_estim_names,'RHO_X')), 0, 1};
    ESTIM_PARAM.RHO_A   = {xparam_model(ismember(MODEL.param_estim_names,'RHO_A')), 0, 1};
    ESTIM_PARAM.RHO_E   = {xparam_model(ismember(MODEL.param_estim_names,'RHO_E')), 0, 1};
end

% update fields
MODEL.param_estim_names = fieldnames(ESTIM_PARAM);
MODEL.param_estim_nbr   = size(MODEL.param_estim_names,1);
% create parameter vectors from structures
xparam = nan(MODEL.param_estim_nbr,1);
for j = 1:MODEL.param_estim_nbr
    OPT.optimizer.bounds.lb(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){2};
    OPT.optimizer.bounds.ub(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){3};
    xparam(j,1) = ESTIM_PARAM.(MODEL.param_estim_names{j}){1};
end