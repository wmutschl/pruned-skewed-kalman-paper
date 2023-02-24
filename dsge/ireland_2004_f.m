function f = ireland_2004_f(STEADY_STATE,PARAM)

% Evaluate numerical values for parameters from PARAM
BETA = PARAM.BETA;
PSI = PARAM.PSI;
ALPHA_PI = PARAM.ALPHA_PI;
ALPHA_X = PARAM.ALPHA_X;
RHO_A = PARAM.RHO_A;
RHO_E = PARAM.RHO_E;
OMEGA = PARAM.OMEGA;
RHO_PI = PARAM.RHO_PI;
RHO_G = PARAM.RHO_G;
RHO_X = PARAM.RHO_X;

% Evaluate numerical values for steady-state from STEADY_STATE (this also updates endogenous parameters)
steady_state_a = STEADY_STATE.a;
steady_state_e = STEADY_STATE.e;
steady_state_z = STEADY_STATE.z;
steady_state_x = STEADY_STATE.x;
steady_state_pihat = STEADY_STATE.pihat;
steady_state_yhat = STEADY_STATE.yhat;
steady_state_ghat = STEADY_STATE.ghat;
steady_state_rhat = STEADY_STATE.rhat;

% Evaluate exogenous variables from STEADY_STATE
eta_a = STEADY_STATE.eta_a;
eta_e = STEADY_STATE.eta_e;
eta_z = STEADY_STATE.eta_z;
eta_r = STEADY_STATE.eta_r;

% Set all endogenous variables to steady-state values
a_back = steady_state_a; a_curr = steady_state_a; a_fwrd = steady_state_a;
e_back = steady_state_e; e_curr = steady_state_e; e_fwrd = steady_state_e;
z_back = steady_state_z; z_curr = steady_state_z; z_fwrd = steady_state_z;
x_back = steady_state_x; x_curr = steady_state_x; x_fwrd = steady_state_x;
pihat_back = steady_state_pihat; pihat_curr = steady_state_pihat; pihat_fwrd = steady_state_pihat;
yhat_back = steady_state_yhat; yhat_curr = steady_state_yhat; yhat_fwrd = steady_state_yhat;
ghat_back = steady_state_ghat; ghat_curr = steady_state_ghat; ghat_fwrd = steady_state_ghat;
rhat_back = steady_state_rhat; rhat_curr = steady_state_rhat; rhat_fwrd = steady_state_rhat;

% Write out model equations
f(1) = eta_a - a_curr + RHO_A*a_back;
f(2) = eta_e - e_curr + RHO_E*e_back;
f(3) = eta_z - z_curr;
f(4) = pihat_fwrd - rhat_curr - x_curr + ALPHA_X*x_back - x_fwrd*(ALPHA_X - 1) + a_curr*(OMEGA - 1)*(RHO_A - 1);
f(5) = PSI*x_curr - pihat_curr - e_curr + BETA*(ALPHA_PI*pihat_back - pihat_fwrd*(ALPHA_PI - 1));
f(6) = yhat_curr - x_curr - OMEGA*a_curr;
f(7) = yhat_curr - yhat_back - ghat_curr + z_curr;
f(8) = eta_r + rhat_back - rhat_curr + RHO_G*ghat_curr + RHO_PI*pihat_curr + RHO_X*x_curr;

end % function end 
