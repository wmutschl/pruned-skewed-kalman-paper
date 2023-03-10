function DisplayResults_MakeLatexTables_ParameterEstimation(OPTIM_NAME,SAMPLE_SIZES,ARCH,ONLY_LATEX)
% DisplayResults_MakeLatexTables_ParameterEstimation(filename, ONLY_LATEX)
% -------------------------------------------------------------------------
% Collects and displays the results of the Monte-Carlo assessment of the 
% accuracy of parameter estimates of the CSN-distributed innovations in a 
% linear state-space system with both the Gaussian as well as Pruned Skewed
% Kalman filter in Parameter_Estimation.m
% Creates log files with results and also the Latex code for Table 4 of the
% paper "Pruned Skewed Kalman Filter and Smoother: With Application
% to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% -------------------------------------------------------------------------
% INPUTS
% - filename   [string]   load results created with Computing_Time_Log_Likliheood.m and saved in results/"filename".mat,
%                         for example: results_parameter_estimation_fminsearch_T50_R1000_maca64
% - ONLY_LATEX [boolean]  only print the latex rows in the command window
% -------------------------------------------------------------------------
% OUTPUTS
% - results are displayed in the command window
% - latex codes are displayed in the command window
% - a log file "filename.log" is created in the results folder
% =========================================================================
% Copyright (C) 2022-2023 Willi Mutschler
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================
if nargin < 4
    ONLY_LATEX = false;
end

load(['results/' filename],'OPT','PARAMS','RES');

%% start log file
if ~ONLY_LATEX
    diary off
    OPT.logfile = ['results/' filename,'.log'];
    if exist(OPT.logfile, 'file')
        delete(OPT.logfile)
    end
    diary(OPT.logfile)
end

if ~ONLY_LATEX
    format long
    disp(repmat('*',1,60));
    % print settings to logfile
    disp(OPT)
    disp(repmat('*',1,60));
    disp('Parameter Matrix G:');         disp(PARAMS.G);
    disp('Parameter Matrix F:');         disp(PARAMS.F);
    disp('Parameter Matrix R:');         disp(PARAMS.R);
    disp('Parameter Matrix mu_eps:');    disp(PARAMS.mu_eps);
    disp('Parameter Matrix Sigma_eps:'); disp(PARAMS.Sigma_eps);
    disp('Parameter Matrix mu_eta:');    disp(PARAMS.mu_eta);
    disp('Parameter Matrix Sigma_eta:'); disp(PARAMS.Sigma_eta);
    disp('Parameter Matrix nu_eta:');    disp(PARAMS.nu_eta);
    disp('Parameter Matrix Gamma_eta:'); disp(PARAMS.Gamma_eta);
    disp('Parameter Matrix Delta_eta:'); disp(PARAMS.Delta_eta);    
    disp(repmat('*',1,60));
    format short
end

%% get index for non-convergent runs
idx_failed = (RES.exitflags==0);
if ~ONLY_LATEX
    fprintf('Number of failed optimizations: %d\n',sum(idx_failed))
end

%% table with summary statistics
if ~ONLY_LATEX
    tblTitles.SampleStats = {};
    for j=1:size(PARAMS.F,1)
        tblTitles.SampleStats = [tblTitles.SampleStats sprintf('mean(y%d)',j) sprintf('sd(y%d)',j) sprintf('skew(y%d)',j)];
    end
    fprintf('%s\n%sSUMMARY STATISTICS\n',repmat('*',1,60),repmat(' ',1,10));
    disp(array2table(RES.SampleStats(~idx_failed,:),'VariableNames',tblTitles.SampleStats));
end

%% table with estimation results
eta_nbr = size(PARAMS.mu_eta,1);
Eeta_true = csnMean(PARAMS.mu_eta,PARAMS.Sigma_eta,PARAMS.Gamma_eta,PARAMS.nu_eta,PARAMS.Delta_eta,OPT.cdfmvna_fct);
COVeta_true = csnVar(PARAMS.Sigma_eta,PARAMS.Gamma_eta,PARAMS.nu_eta,PARAMS.Delta_eta,OPT.cdfmvna_fct);
theta_true = [Eeta_true; diag(COVeta_true); diag(PARAMS.Gamma_eta)];
estim_param_nbr = size(theta_true,1);
NRMSE = nan(estim_param_nbr,size(OPT.prune_tol,2)+1);
MEAN  = nan(estim_param_nbr,size(OPT.prune_tol,2)+1);
LOW   = nan(estim_param_nbr,size(OPT.prune_tol,2)+1);
HIGH  = nan(estim_param_nbr,size(OPT.prune_tol,2)+1);
theta_hat.gauss = [RES.EstimParams.gauss.mu_eta(:,~idx_failed);
                   RES.EstimParams.gauss.Sigma_eta(reshape(find(repmat(logical(eye(eta_nbr)),[1,1,sum(~idx_failed)])),eta_nbr,sum(~idx_failed)));
                   nan(eta_nbr,sum(~idx_failed));
                   ];

for jprune=1:size(OPT.prune_tol,2)
    mu_eta_j    = RES.EstimParams.csn.mu_eta(:,~idx_failed,jprune);
    Sigma_eta_j = RES.EstimParams.csn.Sigma_eta(:,:,~idx_failed,jprune);
    Gamma_eta_j = RES.EstimParams.csn.Gamma_eta(:,:,~idx_failed,jprune);
    E_eta_j = nan(size(mu_eta_j));
    V_eta_j = nan(size(Sigma_eta_j));
    for r = 1:sum(~idx_failed)
        E_eta_j(:,r)   = csnMean(mu_eta_j(:,r),Sigma_eta_j(:,:,r),Gamma_eta_j(:,:,r),PARAMS.nu_eta,PARAMS.Delta_eta,OPT.cdfmvna_fct);
        V_eta_j(:,:,r) = csnVar(Sigma_eta_j(:,:,r),Gamma_eta_j(:,:,r),PARAMS.nu_eta,PARAMS.Delta_eta,OPT.cdfmvna_fct);
    end    
    theta_hat.csn{jprune} = [E_eta_j;
                             V_eta_j(reshape(find(repmat(logical(eye(eta_nbr)),[1,1,sum(~idx_failed)])),eta_nbr,sum(~idx_failed)));
                             Gamma_eta_j(reshape(find(repmat(logical(eye(eta_nbr)),[1,1,sum(~idx_failed)])),eta_nbr,sum(~idx_failed)));
                            ];
    MEAN(:,jprune) = mean(theta_hat.csn{jprune},2);
    LOW(:,jprune) = quantile(theta_hat.csn{jprune},0.05,2);
    HIGH(:,jprune) = quantile(theta_hat.csn{jprune},0.95,2);
    NRMSE(:,jprune) = (1./theta_true).*sqrt(mean((theta_hat.csn{jprune}-theta_true).^2,2));
end
MEAN(:,3) = mean(theta_hat.gauss,2);
LOW(:,3) = quantile(theta_hat.gauss,0.05,2);
HIGH(:,3) = quantile(theta_hat.gauss,0.95,2);
NRMSE(:,3) = (1./theta_true).*sqrt(mean((theta_hat.gauss-theta_true).^2,2));
NRMSE(isinf(NRMSE)) = nan;


disp(['Total computing time : ' dynsec2hms(RES.time_Computing_Time) ]);




if ~ONLY_LATEX
    diary off
end