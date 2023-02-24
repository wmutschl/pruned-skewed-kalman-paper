function MODEL = ireland_2004_preprocessing
%% MINI DYNARE PREPROCESSING
MODEL.name = "ireland_2004";

%% Declare names of variables and parameters
% declare names for endogenous variables as string array, second entry is assumption set for the sym command below
MODEL.endo_names = [
    "a",         "real"; % preference shock
    "e",         "real"; % cost push shock
    "z",         "real"; % TFP shock
    "x",         "real"; % output gap
    "pihat",     "real"; % inflation deviation from trend
    "yhat",      "real"; % output deviation from trend
    "ghat" ,     "real"; % output growth
    "rhat",      "real"; % interest deviations from trend    
];

% Which variables are observables?
MODEL.varobs = ["ghat";"rhat";"pihat"];

% declare names for exogenous variables as string array, second entry is assumption set for the sym command below
MODEL.exo_names = [
    "eta_a" , "real"; % preference innovation
    "eta_e" , "real"; % cost push innovation
    "eta_z" , "real"; % TFP innovation
    "eta_r" , "real"; % monetary policy innovation
];

% declare names for model parameters as string array, second entry is assumption set for the sym command below
MODEL.param_names = [
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

% declare names for endogenous parameters (mostly steady-state) as string array, second entry is assumption set for the sym command below
% note that it is important to keep the same string as in endo_names and prepend a "steady_state_" to it, we also use this variable in the script files below
aux_param_names = [
];

% get sizes
MODEL.endo_nbr = size(MODEL.endo_names,1);
MODEL.exo_nbr = size(MODEL.exo_names,1);
MODEL.param_nbr = size(MODEL.param_names,1);
aux_param_nbr = size(aux_param_names,1);

%% create symbolic variables from declared names
for j = 1:MODEL.endo_nbr
    eval(sprintf('%s_back = sym(''%s_back'',''%s'');',MODEL.endo_names(j,1),MODEL.endo_names(j,1),MODEL.endo_names(j,2))); %this evaluates expressions like x_back = sym('x_back','real');
    eval(sprintf('%s_curr = sym(''%s_curr'',''%s'');',MODEL.endo_names(j,1),MODEL.endo_names(j,1),MODEL.endo_names(j,2))); %this evaluates expressions like x_curr = sym('x_curr','real');
    eval(sprintf('%s_fwrd = sym(''%s_fwrd'',''%s'');',MODEL.endo_names(j,1),MODEL.endo_names(j,1),MODEL.endo_names(j,2))); %this evaluates expressions like x_fwrd = sym('x_fwrd','real');
end
for j = 1:MODEL.exo_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',MODEL.exo_names(j,1),MODEL.exo_names(j,1),MODEL.exo_names(j,2))); %this evaluates expressions like eta_x = sym('eta_x','real');
end
for j = 1:MODEL.param_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',MODEL.param_names(j,1),MODEL.param_names(j,1),MODEL.param_names(j,2))); %this evaluates expressions like param = sym('param','real');
end
for j = 1:aux_param_nbr
    eval(sprintf('%s = sym(''%s'',''%s'');',aux_param_names(j,1),aux_param_names(j,1),aux_param_names(j,2))); %this evaluates expressions like steady_state_x = sym('steady_state_x','real');
end

%% Model Equations
% Note the following convention: x_back for t-1 variables, x_curr for t variables, x_fwrd for t+1 variables, steady_state_x for steady-state value

% temporary preference shock (15)
model_eqs(1) = -a_curr + RHO_A*a_back + eta_a;
% temporary cost-push shock (16)
model_eqs(2) = -e_curr + RHO_E*e_back + eta_e;
% technology shock (17)
model_eqs(3) = -z_curr + eta_z;
% New Keynesian IS curve (23)
model_eqs(4) = -x_curr + ALPHA_X*x_back + (1-ALPHA_X)*x_fwrd -(rhat_curr-pihat_fwrd)+(1-OMEGA)*(1-RHO_A)*a_curr;
% New Keynesian PC curve (24)
model_eqs(5) = -pihat_curr + BETA*(ALPHA_PI*pihat_back + (1-ALPHA_PI)*pihat_fwrd)+PSI*x_curr - e_curr;
% output gap (20)
model_eqs(6) = -x_curr + yhat_curr - OMEGA*a_curr;
% growth rate of output (21)
model_eqs(7) = -ghat_curr + yhat_curr - yhat_back + z_curr;
% policy rule (22)
model_eqs(8) = -(rhat_curr-rhat_back) + RHO_PI*pihat_curr + RHO_G*ghat_curr + RHO_X*x_curr + eta_r;

%% Define and reorder variables
MODEL.lead_lag_incidence = zeros(3,MODEL.endo_nbr); %rows are time periods, columns are endogenous variables
for j = 1:MODEL.endo_nbr
    MODEL.lead_lag_incidence(1,j) = any(ismember(symvar(model_eqs),str2sym([MODEL.endo_names{j,1}, '_back']))); %1 if variable appears at t-1
    MODEL.lead_lag_incidence(2,j) = any(ismember(symvar(model_eqs),str2sym([MODEL.endo_names{j,1}, '_curr']))); %1 if variable appears at t
    MODEL.lead_lag_incidence(3,j) = any(ismember(symvar(model_eqs),str2sym([MODEL.endo_names{j,1}, '_fwrd']))); %1 if variable appears at t+1  
end
MODEL.endo_static = MODEL.endo_names(ismember(transpose(MODEL.lead_lag_incidence),[0 1 0],'rows'),1); % variable appears only at t
MODEL.endo_pred   = MODEL.endo_names(ismember(transpose(MODEL.lead_lag_incidence([1 3],:)),[1 0],'rows'),1); % variable appears at t-1 but not at t+1 (don't care about t)
MODEL.endo_both   = MODEL.endo_names(ismember(transpose(MODEL.lead_lag_incidence),[1 1 1],'rows'),1); % variable appears at t-1, t, t+1
MODEL.endo_fwrd   = MODEL.endo_names(ismember(transpose(MODEL.lead_lag_incidence([1 3],:)),[0 1],'rows'),1); % variable appears at t+1 but not at t-1 (don't care about t)
MODEL.endo_names_DR = [MODEL.endo_static;MODEL.endo_pred;MODEL.endo_both;MODEL.endo_fwrd];
if ~isequal(sort(MODEL.endo_names_DR),sort(MODEL.endo_names(:,1)))
    error('Could not determine types of variables correctly');
end
MODEL.endo_static_nbr = size(MODEL.endo_static,1);
MODEL.endo_pred_nbr = size(MODEL.endo_pred,1);
MODEL.endo_both_nbr = size(MODEL.endo_both,1);
MODEL.endo_fwrd_nbr = size(MODEL.endo_fwrd,1);
MODEL.states = [MODEL.endo_pred;MODEL.endo_both];
[~,MODEL.DR_ORDER] = ismember(MODEL.endo_names_DR,MODEL.endo_names(:,1));
MODEL.varobs_nbr = size(MODEL.varobs,1);
for j = 1:MODEL.varobs_nbr
    MODEL.varobs_idx_DR(j) = find(ismember(MODEL.endo_names_DR,MODEL.varobs(j)));
end

%% Take derivative with respect to [x_back;x_curr;x_fwrd;exo]
g1 = jacobian(model_eqs,[str2sym([MODEL.endo_static;MODEL.endo_pred;MODEL.endo_both;MODEL.endo_fwrd] + "_back");...
                         str2sym([MODEL.endo_static;MODEL.endo_pred;MODEL.endo_both;MODEL.endo_fwrd] + "_curr");...
                         str2sym([MODEL.endo_static;MODEL.endo_pred;MODEL.endo_both;MODEL.endo_fwrd] + "_fwrd");...
                         str2sym(MODEL.exo_names(:,1));...
                        ]);
g1 = simplify(g1); % simplify algebra

%% Write out model equations to script files
nameOfFunction = MODEL.name + "_f";

% Delete old version of file (if it exists)
if exist(nameOfFunction+".m",'file') > 0
    delete(nameOfFunction+".m")
end
fileID = fopen(nameOfFunction+".m",'w');
fprintf(fileID,'function f = %s(STEADY_STATE,PARAM)\n',nameOfFunction);

fprintf(fileID,'\n%% Evaluate numerical values for parameters from PARAM\n');
for j = 1:MODEL.param_nbr
    fprintf(fileID,'%s = PARAM.%s;\n',MODEL.param_names{j,1},MODEL.param_names{j,1});
end

fprintf(fileID,'\n%% Evaluate numerical values for steady-state from STEADY_STATE (this also updates endogenous parameters)\n');
for j = 1:MODEL.endo_nbr
    fprintf(fileID,'steady_state_%s = STEADY_STATE.%s;\n',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
end

fprintf(fileID,'\n%% Evaluate exogenous variables from STEADY_STATE\n');
for j = 1:MODEL.exo_nbr
    fprintf(fileID,'%s = STEADY_STATE.%s;\n',MODEL.exo_names{j,1},MODEL.exo_names{j,1});
end

fprintf(fileID,'\n%% Set all endogenous variables to steady-state values\n');
for j = 1:MODEL.endo_nbr
    fprintf(fileID,'%s_back = steady_state_%s; ',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
    fprintf(fileID,'%s_curr = steady_state_%s; ',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
    fprintf(fileID,'%s_fwrd = steady_state_%s;\n',MODEL.endo_names{j,1},MODEL.endo_names{j,1});    
end

fprintf(fileID,'\n%% Write out model equations\n');
for j = 1:MODEL.endo_nbr
    fprintf(fileID,'f(%d) = %s;\n',j,char(model_eqs(j)));
end

fprintf(fileID,'\nend %% function end \n');
fclose(fileID);

%% Write out derivatives to script files
nameOfFunction = MODEL.name + "_g1";

% Delete old version of file (if it exists)
if exist(nameOfFunction+".m",'file') > 0
    delete(nameOfFunction+".m")
end
fileID = fopen(nameOfFunction+".m",'w');
fprintf(fileID,'function g1 = %s(STEADY_STATE,PARAM)\n',nameOfFunction);

fprintf(fileID,'\n%% Evaluate numerical values for parameters from PARAM\n');
for j = 1:MODEL.param_nbr
    fprintf(fileID,'%s = PARAM.%s;\n',MODEL.param_names{j,1},MODEL.param_names{j,1});
end

fprintf(fileID,'\n%% Evaluate numerical values for steady-state from STEADY_STATE (this also updates endogenous parameters)\n');
for j = 1:MODEL.endo_nbr
    fprintf(fileID,'steady_state_%s = STEADY_STATE.%s;\n',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
end

fprintf(fileID,'\n%% Evaluate exogenous variables from STEADY_STATE\n');
for j = 1:MODEL.exo_nbr
    fprintf(fileID,'%s = STEADY_STATE.%s;\n',MODEL.exo_names{j,1},MODEL.exo_names{j,1});
end

fprintf(fileID,'\n%% Set all endogenous variables to steady-state values\n');
for j = 1:MODEL.endo_nbr
    fprintf(fileID,'%s_back = steady_state_%s; ',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
    fprintf(fileID,'%s_curr = steady_state_%s; ',MODEL.endo_names{j,1},MODEL.endo_names{j,1});
    fprintf(fileID,'%s_fwrd = steady_state_%s;\n',MODEL.endo_names{j,1},MODEL.endo_names{j,1});    
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