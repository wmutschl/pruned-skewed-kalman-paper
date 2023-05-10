function M_ = ireland2004_preprocessing
% function M_ = ireland2004_preprocessing
% -------------------------------------------------------------------------
% preprocesses the small scale New Keynesian model of Ireland (2004): 
% "Technology Shocks in The New Keynesian Model", The Review of Economics and Statistics
% based on Dynare replication codes kindly provided by Johannes Pfeifer at
% https://github.com/JohannesPfeifer/DSGE_mod/tree/master/Ireland_2004
% -------------------------------------------------------------------------
% Note that this is basically a stripped down MATLAB version of Dynare's
% preprocessor, using MATLAB's symbolic toolbox to set up the model and
% its dynamic Jacobian analytically and then writing the expressions out
% to script files. The same ordering as in Dynare is imposed, i.e. there
% is a distinction between declaration and DR order (which is used in the
% policy function), see the manual of Dynare for clarifications.
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% OUTPUTS
% - M_   [structure]   information on the model (inspired by Dynare's M_ structure)
% =========================================================================
% Copyright (C) 2023 Willi Mutschler
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

%% symbolic declarations
% Model name
M_.fname = "ireland2004";

% files are stored within the folder
if strcmp(computer('arch'),'win64')
    folder = erase(mfilename('fullpath'),"\"+M_.fname+"_preprocessing");
else
    folder = erase(mfilename('fullpath'),"/"+M_.fname+"_preprocessing");
end
oldfolder = cd(folder);

% declare names for endogenous variables as string array, second entry is assumption set for the sym command below
% similar to Dynare's var block
M_.endo_names = [
    "ahat",      "real"; % preference shock
    "ehat",      "real"; % cost push shock
    "zhat",      "real"; % TFP shock
    "xhat",      "real"; % output gap
    "pihat",     "real"; % inflation deviation from trend
    "yhat",      "real"; % output deviation from trend
    "ghat" ,     "real"; % output growth
    "rhat",      "real"; % interest deviations from trend
];

% Which variables are observables?
% similar to Dynare's varobs block
M_.varobs = ["ghat";"rhat";"pihat"];
varobs_nbr = length(M_.varobs);
M_.varobs_id = nan(1,varobs_nbr);
for jvarobs = 1:length(M_.varobs)
    M_.varobs_id(jvarobs) = find(contains(M_.endo_names(:,1),M_.varobs(jvarobs)));
end

% declare names for exogenous variables as string array, second entry is assumption set for the sym command below
% similar to Dynare's varexo block
M_.exo_names = [
    "eta_a" , "real"; % preference innovation
    "eta_e" , "real"; % cost push innovation
    "eta_z" , "real"; % TFP innovation
    "eta_r" , "real"; % monetary policy innovation
];

% declare names for model parameters as string array, second entry is assumption set for the sym command below
% similar to Dynare's parameters block
M_.param_names = [
    "BETA"     , "real"; % discount factor
    "PSI"      , "real"; % output gap slope in Phillips curve
    "ALPHA_PI" , "real"; % slope parameter in Phillips curve
    "ALPHA_X"  , "real"; % slope parameter in IS curve
    "RHO_A"    , "real"; % persistence preference shock
    "RHO_E"    , "real"; % persistence cost-push shock
    "OMEGA"    , "real"; % scale parameter preference innovation in IS curve
    "RHO_PI"   , "real"; % feedback policy rule inflation
    "RHO_G"    , "real"; % feedback policy rule output growth
    "RHO_X"    , "real"; % feedback policy rule output gap
];
aux_param_names = [];

% get sizes
M_.endo_nbr   = size(M_.endo_names,1);
M_.exo_nbr    = size(M_.exo_names,1);
M_.param_nbr  = size(M_.param_names,1);
aux_param_nbr = size(aux_param_names,1);

% create symbolic variables from declared names using MATLAB's symbolic toolbox
for j = 1:M_.endo_nbr
    eval(sprintf('%s_back = sym(''%s_back'',''%s'');',M_.endo_names(j,1),M_.endo_names(j,1),M_.endo_names(j,2))); %this evaluates expressions like x_back = sym('x_back','real');
    eval(sprintf('%s_curr = sym(''%s_curr'',''%s'');',M_.endo_names(j,1),M_.endo_names(j,1),M_.endo_names(j,2))); %this evaluates expressions like x_curr = sym('x_curr','real');
    eval(sprintf('%s_fwrd = sym(''%s_fwrd'',''%s'');',M_.endo_names(j,1),M_.endo_names(j,1),M_.endo_names(j,2))); %this evaluates expressions like x_fwrd = sym('x_fwrd','real');
end
for j = 1:M_.exo_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',M_.exo_names(j,1),M_.exo_names(j,1),M_.exo_names(j,2))); %this evaluates expressions like eta_x = sym('eta_x','real');
end
for j = 1:M_.param_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',M_.param_names(j,1),M_.param_names(j,1),M_.param_names(j,2))); %this evaluates expressions like param = sym('param','real');
end
for j = 1:aux_param_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',aux_param_names(j,1),aux_param_names(j,1),aux_param_names(j,2))); %this evaluates expressions like steady_state_x = sym('steady_state_x','real');
end

%% Model Equations
% similar to Dynare's model block, but with all equations written on the right-hand-side of the equations
% note the following convention with respect to Dynare: x_back for x(-1), x_curr for x(0), x_fwrd for x(+1) variables, steady_state_x for steady-state(x)

model_eqs(1) = -ahat_curr + RHO_A*ahat_back + eta_a; % temporary preference shock (15)
model_eqs(2) = -ehat_curr + RHO_E*ehat_back + eta_e; % temporary cost-push shock (16)
model_eqs(3) = -zhat_curr + eta_z; % technology shock (17)
model_eqs(4) = -xhat_curr + ALPHA_X*xhat_back + (1-ALPHA_X)*xhat_fwrd - (rhat_curr-pihat_fwrd) + (1-OMEGA)*(1-RHO_A)*ahat_curr; % New Keynesian IS curve (23)
model_eqs(5) = -pihat_curr + BETA*( ALPHA_PI*pihat_back + (1-ALPHA_PI)*pihat_fwrd ) + PSI*xhat_curr - ehat_curr; % New Keynesian PC curve (24)
model_eqs(6) = -xhat_curr + yhat_curr - OMEGA*ahat_curr; % output gap (20)
model_eqs(7) = -ghat_curr + yhat_curr - yhat_back + zhat_curr; % growth rate of output (21)
model_eqs(8) = -(rhat_curr-rhat_back) + RHO_PI*pihat_curr + RHO_G*ghat_curr + RHO_X*xhat_curr + eta_r; % policy rule (22)

%% DR ordering and different types of variables
% define and reorder variables
M_.lead_lag_incidence = zeros(3,M_.endo_nbr); % rows are time periods, columns are endogenous variables; this matrix encodes whether a varialbe appears at t-1, t or t+1
TIME = ["_back","_curr","_fwrd"];
idx = 1;
for jtime = 1:length(TIME)
    for jendo = 1:M_.endo_nbr
        if any(ismember(symvar(model_eqs),str2sym(M_.endo_names(jendo,1) + TIME(jtime))))
            M_.lead_lag_incidence(jtime,jendo) = idx; idx=idx+1;
        end
    end
end
% get group of variables as in Dynare's DR ordering
endo_static = M_.endo_names(ismember(transpose(M_.lead_lag_incidence)>0,[0 1 0],'rows'),1); % variable appears only at t
endo_pred   = M_.endo_names(ismember(transpose(M_.lead_lag_incidence([1 3],:))>0,[1 0],'rows'),1); % variable appears at t-1 but not at t+1 (don't care about t)
endo_both   = M_.endo_names(ismember(transpose(M_.lead_lag_incidence)>0,[1 1 1],'rows'),1); % variable appears at t-1, t, t+1
endo_fwrd   = M_.endo_names(ismember(transpose(M_.lead_lag_incidence([1 3],:))>0,[0 1],'rows'),1); % variable appears at t+1 but not at t-1 (don't care about t)
endo_names_DR = [endo_static;endo_pred;endo_both;endo_fwrd]; % corresponds to Dynare's M_.endo_names(oo_.dr.order_var)
if ~isequal(sort(endo_names_DR),sort(M_.endo_names(:,1)))
    error('Could not determine types of variables correctly');
end
M_.nstatic = size(endo_static,1);
M_.npred   = size(endo_pred,1);
M_.nboth   = size(endo_both,1);
M_.nfwrd   = size(endo_fwrd,1);
M_.nspred  = M_.npred + M_.nboth;
M_.nsfwrd  = M_.nfwrd + M_.nboth;
state_var  = [endo_pred;endo_both];
M_.state_var = nan(1,M_.nspred);
for jstate = 1:M_.nspred
    M_.state_var(jstate) = find(M_.endo_names(:,1)==state_var(jstate));
end
% get Dynare's useful indices to transform from declaration to DR order (and the other way around)
[~,M_.order_var] = ismember(endo_names_DR,M_.endo_names(:,1));
[~,M_.inv_order_var] = ismember(M_.endo_names(:,1),endo_names_DR);
% create selection matrix F such that varobs = F*y, where y is in DR order
M_.F = zeros(length(M_.varobs),M_.endo_nbr);
for jvarobs=1:length(M_.varobs)
    M_.F(jvarobs,find(ismember(M_.endo_names(M_.order_var),M_.varobs(jvarobs)),1)) = 1;
end

%% compute Jacobian of dynamic model
% take derivative with respect to [x_back;x_curr;x_fwrd;exo]
g1 = jacobian(model_eqs,[str2sym([endo_static;endo_pred;endo_both;endo_fwrd] + "_back");...
                         str2sym([endo_static;endo_pred;endo_both;endo_fwrd] + "_curr");...
                         str2sym([endo_static;endo_pred;endo_both;endo_fwrd] + "_fwrd");...
                         str2sym(M_.exo_names(:,1));...
                        ]);
g1 = simplify(g1); % simplify algebra

%% write out model equations to script files
nameOfFunction = M_.fname + "_f";
% Delete old version of file (if it exists)
if exist(nameOfFunction+".m",'file') > 0
    delete(nameOfFunction+".m")
end
fileID = fopen(nameOfFunction+".m",'w');
fprintf(fileID,'function f = %s(STEADY_STATE,PARAMS)\n',nameOfFunction);
fprintf(fileID,'%% function f = %s(STEADY_STATE,PARAMS)\n',nameOfFunction);
fprintf(fileID,'%% Model equations, automatically generated by %s_preprocessing.m\n',M_.fname);
fprintf(fileID,'\n%% Evaluate numerical values for parameters from PARAMS\n');
for j = 1:M_.param_nbr
    fprintf(fileID,'%s = PARAMS(%d);\n',M_.param_names{j,1},j);
end
fprintf(fileID,'\n%% Evaluate numerical values for steady-state from STEADY_STATE (this also updates endogenous parameters)\n');
for j = 1:M_.endo_nbr
    fprintf(fileID,'steady_state_%s = STEADY_STATE(%d);\n',M_.endo_names{j,1},j);
end
fprintf(fileID,'\n%% Evaluate exogenous variables\n');
for j = 1:M_.exo_nbr
    fprintf(fileID,'%s = 0;\n',M_.exo_names{j,1});
end
fprintf(fileID,'\n%% Set all endogenous variables to steady-state values\n');
for j = 1:M_.endo_nbr
    fprintf(fileID,'%s_back = steady_state_%s; ',M_.endo_names{j,1},M_.endo_names{j,1});
    fprintf(fileID,'%s_curr = steady_state_%s; ',M_.endo_names{j,1},M_.endo_names{j,1});
    fprintf(fileID,'%s_fwrd = steady_state_%s;\n',M_.endo_names{j,1},M_.endo_names{j,1});
end
fprintf(fileID,'\n%% Write out model equations\n');
for j = 1:M_.endo_nbr
    fprintf(fileID,'f(%d) = %s;\n',j,char(model_eqs(j)));
end
fprintf(fileID,'\nend %% function end \n');
fclose(fileID);

%% write out derivatives to script files
nameOfFunction = M_.fname + "_g1";
% Delete old version of file (if it exists)
if exist(nameOfFunction+".m",'file') > 0
    delete(nameOfFunction+".m")
end
fileID = fopen(nameOfFunction+".m",'w');
fprintf(fileID,'function g1 = %s(STEADY_STATE,PARAMS)\n',nameOfFunction);
fprintf(fileID,'%% function g1 = %s(STEADY_STATE,PARAMS)\n',nameOfFunction);
fprintf(fileID,'%% Dynamic model Jacobian, automatically generated by %s_preprocessing.m\n',M_.fname);
fprintf(fileID,'\n%% Evaluate numerical values for parameters from PARAMS\n');
for j = 1:M_.param_nbr
    fprintf(fileID,'%s = PARAMS(%d);\n',M_.param_names{j,1},j);
end
fprintf(fileID,'\n%% Evaluate numerical values for steady-state from STEADY_STATE (this also updates endogenous parameters)\n');
for j = 1:M_.endo_nbr
    fprintf(fileID,'steady_state_%s = STEADY_STATE(%d);\n',M_.endo_names{j,1},j);
end
fprintf(fileID,'\n%% Evaluate exogenous variables\n');
for j = 1:M_.exo_nbr
    fprintf(fileID,'%s = 0;\n',M_.exo_names{j,1});
end
fprintf(fileID,'\n%% Set all endogenous variables to steady-state values\n');
for j = 1:M_.endo_nbr
    fprintf(fileID,'%s_back = steady_state_%s; ',M_.endo_names{j,1},M_.endo_names{j,1});
    fprintf(fileID,'%s_curr = steady_state_%s; ',M_.endo_names{j,1},M_.endo_names{j,1});
    fprintf(fileID,'%s_fwrd = steady_state_%s;\n',M_.endo_names{j,1},M_.endo_names{j,1});
end
fprintf(fileID,'\n%% Initialize first dynamic Jacobian g1\n');
fprintf(fileID,'g1 = zeros(%d, %d);',size(g1,1),size(g1,2));
fprintf(fileID,'\n%% Evaluate non-zero entries in g1\n');
[nonzero_row,nonzero_col,nonzero_vals] = find(g1);
for j = 1:size(nonzero_vals,1)
    fprintf(fileID,'g1(%d,%d) = %s;\n',nonzero_row(j),nonzero_col(j),char(nonzero_vals(j)));
end
fprintf(fileID,'\nend %% function end \n');
fclose(fileID);

%% cleanup and initializations
% remove "real" from name structures
M_.endo_names = M_.endo_names(:,1);
M_.exo_names = M_.exo_names(:,1);
M_.param_names = M_.param_names(:,1);

% initialize parameters
M_.params    = zeros(M_.param_nbr,1);
M_.mu_eps    = zeros(length(M_.varobs),1);
M_.Sigma_eps = zeros(length(M_.varobs),length(M_.varobs));
M_.Cov_eps   = zeros(length(M_.varobs),length(M_.varobs));
M_.mu_eta    = zeros(M_.exo_nbr,1);
M_.nu_eta    = zeros(M_.exo_nbr,1);
M_.Sigma_eta = zeros(M_.exo_nbr,M_.exo_nbr);
M_.Gamma_eta = zeros(M_.exo_nbr,M_.exo_nbr);
M_.Delta_eta = eye(M_.exo_nbr);
M_.Cov_eta   = zeros(M_.exo_nbr,M_.exo_nbr);
M_.Skew_eta  = zeros(M_.exo_nbr,1);

% housekeeping
cd(oldfolder);