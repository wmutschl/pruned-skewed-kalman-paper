% Maximum Likelihood estimation of the DSGE model and data used in
% Ireland (2004) - "Technolog Shocks in The New Keynesian Model"
% The model is estimated first with the Gaussian Kalman filter to replicate
% the original paper; and then with the Pruned Skewed Kalman filter to allow
% for skewed innovations to preference, cost-push, productivity and monetary policy shocks.
% The model is estimated on the post-1980 sample.
% =========================================================================
% Copyright © 2023 Willi Mutschler
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
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

%% housekeeping
clearvars; clearvars -global; close all; clc;
addpath('MATLAB');
addpath('MATLAB/external');
addpath('MATLAB/ireland2004');

REDO = 1; % 0: load mat file with results, 1: redo estimation with initial values provided in estim_params file, 2: redo estimation with search for initial values

if REDO == 0
    load('results/results_ireland2004_stderrskew1_KalmanInit1_FixAx0_FixAp0_R2023a_maci64.mat');
else
    %% MODEL PREPROCESSING
    % create script files with analytic dynamic Jacobian, to compute the policy function numerically, inspired by Dynare's preprocessor, but written in MATLAB
    M_ = ireland2004_preprocessing;
    
    % common calibration, values taken from Ireland (2004)
    M_.params(ismember(M_.param_names,'BETA')) = 0.99;
    M_.params(ismember(M_.param_names,'PSI'))  = 0.1;
    
    %% DATA
    load data/ireland2004_gpr.dat; % original dataset
    ireland2004_gpr = array2timetable(ireland2004_gpr,...
           'RowTimes',datetime('1948-Q2','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'):calquarters(1):datetime('2003-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ'),...
           'VariableNames',{'ghat',...  % output growth (needs to be demeaned), based on quarterly changes in seasonally adjusted real GDP, converted to per capita terms by dividing by the civilian noninsitutional population aged 16 and over
                            'pihat',... % inflation (needs to be demeaned), based on quarterly changes in seasonally adjusted GDP deflator
                            'rhat'} ... % nominal interest rate (needs to be demeaned), based on quarterly averages of daily readings on the three-months US Treasury bill rate
           );
    ireland2004_post1980 = ireland2004_gpr(find(ireland2004_gpr.Time==datetime('1980-Q1','InputFormat','yyyy-QQQ','Format','yyyy-QQQ')):end,:);
    datamat = [];
    for jvarobs=1:length(M_.varobs)
        datamat = [datamat (ireland2004_post1980.(M_.varobs(jvarobs))-mean(ireland2004_post1980.(M_.varobs(jvarobs))))]; %demean data on actual sample
    end
    clear ireland2004_gpr ireland2004_post1980 jvarobs
    
    %% OPTIONS
    options_.computer_arch = computer('arch');
    options_.dsge = 1; % 1: we are interested in the structural model parameters and not the state-space parameters of the linearized solution
    options_.optim_opt = optimset('display','final','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-4,'TolX',1e-4); % optimization options
    options_.optim_opt.names = ["cmaes" "cmaes_dsge" "fminsearch" "fminsearchbnd" "fmincon" "fminunc" "simulannealbnd"]; % names of optimizer that will be used in parallel
    options_.kalman.csn.prune_tol = 0.01; % pruning threshold in Pruned Skewed Kalman filter
    options_.kalman.csn.cdfmvna_fct = "logmvncdf_ME"; % function to use to evaluate high-dimensional Gaussian log cdf, possible options: "logmvncdf_ME", "mvncdf", "qsilatmvnv", "qsimvnv"
    options_.kalman.lik_init = 1; % 1: stationary distribution, i.e. Gaussian distribution, where initial matrix of variance of the error of forecast is set equal to the unconditional variance of the state variables;
                                  % 2: wide prior, i.e. Gaussian with an initial matrix of variance of the error of forecast diagonal with 10 on the diagonal
                                  % 3: run Gaussian Kalman filter (without computing the likelihood) and take smoothed Sigma matrix
    options_.parameters.skewness_bounds = [-0.95 0.95]; % bounds for skewness coefficient used in estim_params file;
                                                        % note that the absolute value of the theoretical skewness coefficient of the CSN distribution is bounded by 0.995
                                                        % inference at the bound is problematic due to so-called pile ups; in our experience the likelihood becomes extremely flat with respect to Gamma_eta when close to the bound
    options_.parameters.use_stderr_skew = 1;            % 1: optimization algorithm optimizes over skewness coefficient and standard error of CSN shocks instead of Sigma_eta and Gamma_eta,
                                                        % 0: optimization algorithm optimizes over diag(Gamma_eta) and sqrt(diag(Sigma_eta)) instead of skewness coefficient and standard error of CSN shocks
                                                        % note that this has implications for the computation of standard errors
    % transform bounded parameters into unbounded ones for the optimization algorithm using the log-odds-transformation
    options_.parameters.fix.ALPHA_X            = false;
    options_.parameters.fix.ALPHA_PI           = false;
    options_.parameters.transform.OMEGA        = true;
    options_.parameters.transform.RHO_PI       = false;
    options_.parameters.transform.RHO_G        = false;
    options_.parameters.transform.RHO_X        = false;
    options_.parameters.transform.RHO_A        = false;
    options_.parameters.transform.RHO_E        = false;
    options_.parameters.transform.skew_eta_a   = true;
    options_.parameters.transform.skew_eta_e   = true;
    options_.parameters.transform.skew_eta_z   = true;
    options_.parameters.transform.skew_eta_r   = true;
    options_.parameters.transform.stderr_eta_a = true;
    options_.parameters.transform.stderr_eta_e = true;
    options_.parameters.transform.stderr_eta_z = true;
    options_.parameters.transform.stderr_eta_r = true;
    options_.parameters.transform.sqrt_diag_Sigma_eta_a = true;
    options_.parameters.transform.sqrt_diag_Sigma_eta_e = true;
    options_.parameters.transform.sqrt_diag_Sigma_eta_z = true;
    options_.parameters.transform.sqrt_diag_Sigma_eta_r = true;
    options_.parameters.transform.diag_Gamma_eta_a = true;
    options_.parameters.transform.diag_Gamma_eta_e = true;
    options_.parameters.transform.diag_Gamma_eta_z = true;
    options_.parameters.transform.diag_Gamma_eta_r = true;
    if isfield(options_.parameters.fix,'ALPHA_X') && (options_.parameters.fix.ALPHA_X==1)
        M_.params(ismember(M_.param_names,'ALPHA_X')) = 0;
    else
        options_.parameters.transform.ALPHA_X  = true;
    end
    if isfield(options_.parameters.fix,'ALPHA_PI') && (options_.parameters.fix.ALPHA_PI==1)
        M_.params(ismember(M_.param_names,'ALPHA_PI')) = 0;
    else    
        options_.parameters.transform.ALPHA_PI = true;
    end
    
    %% ESTIMATION ON POST 1980 SAMPLE
    options_.filename = sprintf('results_ireland2004_stderrskew%d_KalmanInit%d_FixAx%d_FixAp%d_R%s_%s',options_.parameters.use_stderr_skew,options_.kalman.lik_init,options_.parameters.fix.ALPHA_X,options_.parameters.fix.ALPHA_PI,version('-release'),options_.computer_arch);
    options_.kalman.csn.initval_search = 0; % 0: use initial values provided in estim_params file
    if REDO == 2
        options_.kalman.csn.initval_search = 1; % 1: try to find good initial values first from Gaussian estimation, then grid on csn shock parameters;
        options_.optim_opt.names = ["cmaes" "cmaes_dsge" "fminsearch" "fminsearchbnd" "fmincon" "fminunc"]; % names of optimizer that will be used in parallel ("simulannealbnd" is most time-consuming, so we skip this for initial value search; "cmaes" and "cmaes_dsge" are mildly time-consuming, fminsearch and fminsearchbnd can be mildly time-consuming, "fmincon" and "fminunc" are fast (but not as good) )
        options_.filename = [options_.filename '_initval'];
    end
    [oo_gauss, oo_csn, M_gauss, M_csn, estim_params_gauss, estim_params_csn] = dsge_maximum_likelihood_estimation_csn(M_,options_,datamat);
end % REDO

%% CREATE IMPULSE RESPONSE FUNCTIONS
options_.irfs = 16; % horizon for plots
% compute decision rule
[M_gauss.G, M_gauss.R] = dsge_perturbation_solution_order_1(M_gauss);
[M_csn.G  , M_csn.R  ] = dsge_perturbation_solution_order_1(M_csn  );
% which variables to plot
plot_vars_idx_DR = [find(M_csn.endo_names(M_csn.order_var)=="ghat")...  % output growth
                    find(M_csn.endo_names(M_csn.order_var)=="pihat")... % inflation
                    find(M_csn.endo_names(M_csn.order_var)=="rhat")...  % interest rate
                    find(M_csn.endo_names(M_csn.order_var)=="xhat")...  % output gap
                   ];
% which shocks to plot
plot_shocks_idx = 1:M_csn.exo_nbr; % plot all shocks
% initialize plotting variables
Yhat_gauss_lo = nan(length(plot_vars_idx_DR),options_.irfs,M_csn.exo_nbr);
Yhat_gauss_hi = nan(length(plot_vars_idx_DR),options_.irfs,M_csn.exo_nbr);
Yhat_csn_lo = nan(length(plot_vars_idx_DR),options_.irfs,M_csn.exo_nbr);
Yhat_csn_hi = nan(length(plot_vars_idx_DR),options_.irfs,M_csn.exo_nbr);
for jexo = plot_shocks_idx
    % set impulse to individual shock, all other shocks are zero
    eta_gauss_lo = zeros(M_gauss.exo_nbr,options_.irfs);
    eta_gauss_hi = zeros(M_gauss.exo_nbr,options_.irfs);
    eta_gauss_lo(jexo,1) = gaussianQuantile(0.16, 0, M_gauss.Cov_eta(jexo,jexo));
    eta_gauss_hi(jexo,1) = gaussianQuantile(0.84, 0, M_gauss.Cov_eta(jexo,jexo));
    eta_csn_lo = zeros(M_csn.exo_nbr,options_.irfs);
    eta_csn_hi = zeros(M_csn.exo_nbr,options_.irfs);
    eta_csn_lo(jexo,1) = csnQuantile(0.16, M_csn.mu_eta(jexo,1), M_csn.Sigma_eta(jexo,jexo), M_csn.Gamma_eta(jexo,jexo), M_csn.nu_eta(jexo,1), M_csn.Delta_eta(jexo,jexo), options_.kalman.csn.cdfmvna_fct);
    eta_csn_hi(jexo,1) = csnQuantile(0.84, M_csn.mu_eta(jexo,1), M_csn.Sigma_eta(jexo,jexo), M_csn.Gamma_eta(jexo,jexo), M_csn.nu_eta(jexo,1), M_csn.Delta_eta(jexo,jexo), options_.kalman.csn.cdfmvna_fct);
    
    % initialize at steady-state, note that we simulate for states in deviation from steady-state (so all zeros)
    Xhat_gauss_lo = zeros(M_gauss.endo_nbr,options_.irfs);
    Xhat_gauss_hi = zeros(M_gauss.endo_nbr,options_.irfs);
    Xhat_csn_lo = zeros(M_csn.endo_nbr,options_.irfs);
    Xhat_csn_hi = zeros(M_csn.endo_nbr,options_.irfs);
    for t = 2:options_.irfs
        % simulate using linear state-space system
        Xhat_gauss_lo(:,t) = M_gauss.G*Xhat_gauss_lo(:,t-1) + M_gauss.R*eta_gauss_lo(:,t-1);
        Xhat_gauss_hi(:,t) = M_gauss.G*Xhat_gauss_hi(:,t-1) + M_gauss.R*eta_gauss_hi(:,t-1);
        Xhat_csn_lo(:,t) = M_csn.G*Xhat_csn_lo(:,t-1) + M_csn.R*eta_csn_lo(:,t-1);
        Xhat_csn_hi(:,t) = M_csn.G*Xhat_csn_hi(:,t-1) + M_csn.R*eta_csn_hi(:,t-1);
    end
    % transform plotting variables
    Yhat_gauss_lo(:,:,jexo) = 100*Xhat_gauss_lo(plot_vars_idx_DR,:); % in percentage
    Yhat_gauss_hi(:,:,jexo) = 100*Xhat_gauss_hi(plot_vars_idx_DR,:); % in percentage
    Yhat_gauss_lo(2,:,jexo) = 4*Yhat_gauss_lo(2,:,jexo); %definition annualized inflation
    Yhat_gauss_hi(2,:,jexo) = 4*Yhat_gauss_hi(2,:,jexo); %definition annualized inflation
    Yhat_gauss_lo(3,:,jexo) = 4*Yhat_gauss_lo(3,:,jexo); %definition annualized interest rate
    Yhat_gauss_hi(3,:,jexo) = 4*Yhat_gauss_hi(3,:,jexo); %definition annualized interest rate
    Yhat_csn_lo(:,:,jexo) = 100*Xhat_csn_lo(plot_vars_idx_DR,:); % in percentage
    Yhat_csn_hi(:,:,jexo) = 100*Xhat_csn_hi(plot_vars_idx_DR,:); % in percentage
    Yhat_csn_lo(2,:,jexo) = 4*Yhat_csn_lo(2,:,jexo); %definition annualized inflation
    Yhat_csn_hi(2,:,jexo) = 4*Yhat_csn_hi(2,:,jexo); %definition annualized inflation
    Yhat_csn_lo(3,:,jexo) = 4*Yhat_csn_lo(3,:,jexo); %definition annualized interest rate
    Yhat_csn_hi(3,:,jexo) = 4*Yhat_csn_hi(3,:,jexo); %definition annualized interest rate
end

% options for plots
col_hi = "#800080";
col_lo = "#69b3a2";
linWidth = 3;
FontSize = 16;
strVars = ["output growth", "inflation", "interest rate", "output gap"];
strShocks = ["preference shock", "costpush shock", "productivity shock", "policy shock"];
if options_.computer_arch=="maca64" || options_.computer_arch=="maci64"
    setenv('PATH', [getenv('PATH') ':/usr/local/bin:$HOME/.local/bin:/Library/TeX/texbin']);
end

% create plots
for jexo = plot_shocks_idx
    for jvar=1:length(plot_vars_idx_DR)
        fig.(erase(strVars(jvar)," ")+"_"+erase(strShocks(jexo)," ")) = figure(name=strShocks(jexo) + ": effect on " + strVars(jvar));
        hold on;
        p_gauss_lo = plot(0:(options_.irfs-1),squeeze(Yhat_gauss_lo(jvar,:,jexo)),'--','LineWidth',linWidth,'Color',col_lo);
        p_gauss_hi = plot(0:(options_.irfs-1),squeeze(Yhat_gauss_hi(jvar,:,jexo)),'--','LineWidth',linWidth,'Color',col_hi);
        p_csn_lo = plot(0:(options_.irfs-1),squeeze(Yhat_csn_lo(jvar,:,jexo)),'-','LineWidth',linWidth,'Color',col_lo);
        p_csn_hi = plot(0:(options_.irfs-1),squeeze(Yhat_csn_hi(jvar,:,jexo)),'-','LineWidth',linWidth,'Color',col_hi);
        set(gca,'FontSize',FontSize);
        ylabel(sprintf('Percent'));
        xlabel('Quarters');
        if jexo==4
            ylim([-0.4, 0.4]);
            yticks(-0.4:0.1:0.4);
        end
        grid on;
        hold off;
        legend([p_csn_lo, p_gauss_lo, p_csn_hi, p_gauss_hi],["CSN (16th)" "Gaussian (16th)", "CSN (84th)", "Gaussian (84th)"],'NumColumns',2)
        filename="../Paper/plots/dsge_irfs_"+erase(strShocks(jexo)," ")+"_"+erase(strVars(jvar)," ")+".pdf";
        saveas(fig.(erase(strVars(jvar)," ")+"_"+erase(strShocks(jexo)," ")),filename,'pdf');
        system(sprintf('pdfcrop %s %s',filename,filename));
    end
end

%% SAVE PARAMETERS FOR PDF PLOTS DONE IN R
Sigma_gauss = M_gauss.Cov_eta;
mu_eta = M_csn.mu_eta;
Sigma_eta = M_csn.Sigma_eta;
Gamma_eta = M_csn.Gamma_eta;
nu_eta = M_csn.nu_eta;
Delta_eta = M_csn.Delta_eta;
save('results/results_ireland2004_params_for_pdfs','Sigma_gauss','mu_eta','Sigma_eta','Gamma_eta','nu_eta','Delta_eta','-V6');

%% CREATE ESTIMATION TABLE
% get best CSN estimates to print to table
[~,optim_col] = sort(oo_gauss.neg_log_likelihood); optim_col=optim_col(1);
xparams_gauss = oo_gauss.xparams(:,optim_col);
xparams_gauss_stderr = oo_gauss.xstderr(:,optim_col);
loglik_gauss = -1*oo_gauss.neg_log_likelihood(optim_col);
% get best Gaussian estimates to print to table
[~,optim_col] = sort(oo_csn.neg_log_likelihood); optim_col=optim_col(1);
xparams_csn = oo_csn.xparams(:,optim_col);
xparams_csn_stderr = oo_csn.xstderr(:,optim_col);
loglik_csn = -1*oo_csn.neg_log_likelihood(optim_col);

tbl_latex_model = string.empty; tbl_latex_shock = string.empty;
for jp=1:length(xparams_csn)
    switch erase(estim_params_csn.names(jp),'transformed_')
        case "skew_eta_a"
            latex_name = "$skew(\eta_a)$";
        case "skew_eta_e"
            latex_name = "$skew(\eta_e)$";
        case "skew_eta_z"
            latex_name = "$skew(\eta_z)$";
        case "skew_eta_r"
            latex_name = "$skew(\eta_r)$";
        case "stderr_eta_a"
            latex_name = "$stderr(\eta_a)$";
        case "stderr_eta_e"
            latex_name = "$stderr(\eta_e)$";
        case "stderr_eta_z"
            latex_name = "$stderr(\eta_z)$";
        case "stderr_eta_r"
            latex_name = "$stderr(\eta_r)$";
        case "OMEGA"
            latex_name = "$\omega$";
        case "ALPHA_X"
            latex_name = "$\alpha_x$";
        case "ALPHA_PI"
            latex_name = "$\alpha_\pi$";
        case "RHO_PI"
            latex_name = "$\rho_\pi$";
        case "RHO_G"
            latex_name = "$\rho_g$";
        case "RHO_X"
            latex_name = "$\rho_x$";
        case "RHO_A"
            latex_name = "$\rho_a$";
        case "RHO_E"
            latex_name = "$\rho_e$";
    end
    if contains(latex_name,"skew")
        tbl_latex_shock = [tbl_latex_shock;
                           sprintf('%s & $0$ & $\\underset{(%.4f)}{%.4f}$'...
                          ,latex_name ...
                          ,xparams_csn_stderr(jp),xparams_csn(jp) ...
                          )];
    elseif contains(latex_name,"stderr") || contains(latex_name,"skew")
        tbl_latex_shock = [tbl_latex_shock;
                           sprintf('%s & $\\underset{(%.4f)}{%.4f}$ & $\\underset{(%.4f)}{%.4f}$'...
                          ,latex_name ...
                          ,xparams_gauss_stderr(jp-4),xparams_gauss(jp-4) ...
                          ,xparams_csn_stderr(jp),xparams_csn(jp) ...
                          )];
    else
        tbl_latex_model = [tbl_latex_model;
                           sprintf('%s & $\\underset{(%.4f)}{%.4f}$ & $\\underset{(%.4f)}{%.4f}$'...
                           ,latex_name ...
                           ,xparams_gauss_stderr(jp-4),xparams_gauss(jp-4) ...
                           ,xparams_csn_stderr(jp),xparams_csn(jp) ...
                           )];
    end
end
tbl_latex = tbl_latex_model + "         &$\quad$&$\quad$&         " + tbl_latex_shock + " \\";
tbl_latex = [tbl_latex; "\midrule"; sprintf("\\multicolumn{6}{l}{Value of maximized Log-Likelihood function:} & %.2f & %.2f",loglik_gauss,loglik_csn)];
fprintf('\n\nLatex Table Entries:\n\n')
disp(tbl_latex);
fileID = fopen("../Paper/tables/dsge_ParamEstim.tex",'w');
for j=1:size(tbl_latex,1)
    fprintf(fileID,'%s\n',tbl_latex(j));
end
fclose(fileID);


%% LIKELIHOOD RATIO TEST
fprintf('\n\nLIKELIHOOD RATIO TEST: ')
lr_stat = 2 * (loglik_csn - loglik_gauss);
lr_pval = chi2cdf(lr_stat, 4, "upper");
lr_str = sprintf('test statistic of $%.2f$ and a p-value of $%.4f$.', lr_stat,lr_pval);
fprintf('%s\n',lr_str);
fileID = fopen("../Paper/values/dsge_LR_statistic.tex",'w');
fprintf(fileID,'%s', lr_str);
fclose(fileID);

%% HOUSEKEEPING
rmpath('MATLAB');
rmpath('MATLAB/external');
rmpath('MATLAB/ireland2004');