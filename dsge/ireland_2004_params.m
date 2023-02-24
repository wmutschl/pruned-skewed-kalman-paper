function [PARAM, ESTIM_PARAM, MODEL] = ireland_2004_params(MODEL,OPT)
ESTIM_PARAM = [];
%% Calibrate parameters
PARAM.BETA      = 0.99;   % discount factor
PARAM.PSI       = 0.1;    % output gap slope in Phillips curve
PARAM.ALPHA_PI  = 0.0001; % slope parameter in Phillips curve

MODEL.Sigma_eps = zeros(MODEL.varobs_nbr,MODEL.varobs_nbr); % initialize, calibrate or estimate below
MODEL.Sigma_eta = zeros(MODEL.varobs_nbr,MODEL.varobs_nbr); % initialize, calibrate or estimate below

idx_eta_a = find(ismember(MODEL.exo_names(:,1),"eta_a"),1);
idx_eta_e = find(ismember(MODEL.exo_names(:,1),"eta_e"),1);
idx_eta_z = find(ismember(MODEL.exo_names(:,1),"eta_z"),1);
idx_eta_r = find(ismember(MODEL.exo_names(:,1),"eta_r"),1);

%% Which parameters to estimate?
%%%%ESTIM_PARAM.PARAMNAME        = {INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND}; % comment
if OPT.datafile == "data_full_sample"
    
    % preference shock
    if OPT.estim.eta_a == 1
        ESTIM_PARAM.sqrt_Sigma_eta_a = {0.0405       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_a      = {0            , -Inf       , Inf        };
    elseif OPT.estim.eta_a == 2
        tmp_gamma = fsolve(@(x) skewness_coef_theor(0.0405^2,x) - 0.7,100,optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
        ESTIM_PARAM.sqrt_Sigma_eta_a = {0.0405       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_a      = {tmp_gamma    , -Inf       , Inf        };
    elseif OPT.estim.eta_a == 3
        tmp_gamma = fsolve(@(x) skewness_coef_theor(0.0405^2,x) + 0.7,100,optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
        ESTIM_PARAM.sqrt_Sigma_eta_a = {0.0405       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_a      = {tmp_gamma    , -Inf       , Inf        };
    elseif OPT.estim.eta_a == 4
        ESTIM_PARAM.sqrt_Sigma_eta_a = {0.0500884190745349     , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_a      = {-21.9389249481648       , -Inf      , Inf        };
    else
        if OPT.calib.eta_a == 1
            MODEL.Sigma_eta(idx_eta_a,idx_eta_a) = 0.0405^2;
            MODEL.Gamma_eta(idx_eta_a,idx_eta_a) = 0;
        elseif OPT.calib.eta_a == 2
            MODEL.Sigma_eta(idx_eta_a,idx_eta_a) = 0.0500884190745349^2;
            MODEL.Gamma_eta(idx_eta_a,idx_eta_a) = -21.9389249481648;
        end
    end
    
    % cost-push shock
    if OPT.estim.eta_e == 1
        ESTIM_PARAM.sqrt_Sigma_eta_e = {0.0012       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_e      = {0            , -Inf       , Inf        };
    elseif OPT.estim.eta_e == 2
        tmp_gamma = fsolve(@(x) skewness_coef_theor(0.0012^2,x) - 0.7,100,optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
        ESTIM_PARAM.sqrt_Sigma_eta_e = {0.0012       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_e      = {tmp_gamma    , -Inf       , Inf        };
    elseif OPT.estim.eta_e == 3
        tmp_gamma = fsolve(@(x) skewness_coef_theor(0.0012^2,x) + 0.7,100,optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6));
        ESTIM_PARAM.sqrt_Sigma_eta_e = {0.0012       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_e      = {tmp_gamma    , -Inf       , Inf        };
    elseif OPT.estim.eta_e == 4
        ESTIM_PARAM.sqrt_Sigma_eta_e = {0.00125796411613383     , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_e      = {31.95652980877       , -Inf       , Inf        };
    else
        if OPT.calib.eta_e == 1
            MODEL.Sigma_eta(idx_eta_e,idx_eta_e) = 0.0012^2;
            MODEL.Gamma_eta(idx_eta_e,idx_eta_e) = 0;
        elseif OPT.calib.eta_e == 2
            MODEL.Sigma_eta(idx_eta_e,idx_eta_e) = 0.00125796411613383^2;
            MODEL.Gamma_eta(idx_eta_e,idx_eta_e) = 31.95652980877;
        end
    end
    
    % tfp shock
    if OPT.estim.eta_z == 1
        ESTIM_PARAM.sqrt_Sigma_eta_z = {0.0109       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_z      = {0            , -Inf       , Inf        };
    elseif OPT.estim.eta_z == 2
        ESTIM_PARAM.sqrt_Sigma_eta_z = {0.0103800183811858     , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_z      = {-2.15180582728907       , -Inf       , Inf        };
    else
        if OPT.calib.eta_z == 1
            MODEL.Sigma_eta(idx_eta_z,idx_eta_z) = 0.0109^2;
            MODEL.Gamma_eta(idx_eta_z,idx_eta_z) = 0;
        elseif OPT.calib.eta_z == 2
            MODEL.Sigma_eta(idx_eta_z,idx_eta_z) = 0.0103800183811858^2;
            MODEL.Gamma_eta(idx_eta_z,idx_eta_z) = -2.15180582728907;
        end
    end

    % monetary policy shock
    if OPT.estim.eta_r == 1
        ESTIM_PARAM.sqrt_Sigma_eta_r = {0.0031       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_r      = {0            , -Inf       , Inf        };
    elseif OPT.estim.eta_r == 2
        ESTIM_PARAM.sqrt_Sigma_eta_r = {0.0031607       , 0          , 1          };
        ESTIM_PARAM.Gamma_eta_r      = {10000         , -Inf       , Inf        };
    else
        if OPT.calib.eta_r == 1
            MODEL.Sigma_eta(idx_eta_r,idx_eta_r) = 0.0031^2;
            MODEL.Gamma_eta(idx_eta_r,idx_eta_r) = 0;
        elseif OPT.calib.eta_r == 2
            MODEL.Sigma_eta(idx_eta_r,idx_eta_r) = 0.0031607^2;
            MODEL.Gamma_eta(idx_eta_r,idx_eta_r) = 70.828;
        end
    end

    % model parameters
    if OPT.estim.params
        ESTIM_PARAM.OMEGA            = {0.0617       , 0          , 1          }; % scale parameter preference innovation in IS curve
        ESTIM_PARAM.ALPHA_X          = {0.0836       , 0          , 1          }; % slope parameter in IS curve
        ESTIM_PARAM.RHO_PI           = {0.3597       , 0          , 1          }; % feedback policy rule inflation
        ESTIM_PARAM.RHO_G            = {0.2536       , 0          , 1          }; % feedback policy rule output growth
        ESTIM_PARAM.RHO_X            = {0.0347       , 0          , 1          }; % feedback policy rule output gap
        ESTIM_PARAM.RHO_A            = {0.9470       , 0          , 1          }; % persistence preference shock
        ESTIM_PARAM.RHO_E            = {0.9625       , 0          , 1          }; % persistence cost-push shock
    else
        PARAM.OMEGA     = 0.0617; % scale parameter preference innovation in IS curve
        PARAM.ALPHA_X   = 0.0836; % slope parameter in IS curve
        PARAM.RHO_PI    = 0.3597; % feedback policy rule inflation
        PARAM.RHO_G     = 0.2536; % feedback policy rule output growth
        PARAM.RHO_X     = 0.0347; % feedback policy rule output gap
        PARAM.RHO_A     = 0.9470; % persistence preference shock
        PARAM.RHO_E     = 0.9625; % persistence cost-push shock
    end
end
