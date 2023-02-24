function g1 = ireland_2004_g1(STEADY_STATE,PARAM)

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

% Initialize first dynamic Jacobian g1
g1 = zeros(8, 28);
% Evaluate non-zero entries in g1
g1(1,3) = RHO_A;
g1(2,4) = RHO_E;
g1(7,5) = -1;
g1(8,6) = 1;
g1(4,7) = ALPHA_X;
g1(5,8) = ALPHA_PI*BETA;
g1(3,9) = -1;
g1(7,9) = 1;
g1(7,10) = -1;
g1(8,10) = RHO_G;
g1(1,11) = -1;
g1(4,11) = (OMEGA - 1)*(RHO_A - 1);
g1(6,11) = -OMEGA;
g1(2,12) = -1;
g1(5,12) = -1;
g1(6,13) = 1;
g1(7,13) = 1;
g1(4,14) = -1;
g1(8,14) = -1;
g1(4,15) = -1;
g1(5,15) = PSI;
g1(6,15) = -1;
g1(8,15) = RHO_X;
g1(5,16) = -1;
g1(8,16) = RHO_PI;
g1(4,23) = 1 - ALPHA_X;
g1(4,24) = 1;
g1(5,24) = -BETA*(ALPHA_PI - 1);
g1(1,25) = 1;
g1(2,26) = 1;
g1(3,27) = 1;
g1(8,28) = 1;

end % function end 
