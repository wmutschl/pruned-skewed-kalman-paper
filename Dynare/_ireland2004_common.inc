% Small scale New Keynesian model of Ireland (2004): "Technology Shocks in
% The New Keynesian Model", The Review of Economics and Statistics, based
% on Dynare replication codes kindly provided by Johannes Pfeifer at
% https://github.com/JohannesPfeifer/DSGE_mod/tree/master/Ireland_2004
% =========================================================================
% Copyright (C) 2024-2025 Willi Mutschler
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

var
ahat   ${\hat{a}}$    (long_name='preference shock')
ehat   ${\hat{e}}$    (long_name='cost push shock')
zhat   ${\hat{z}}$    (long_name='TFP shock')
xhat   ${\hat{x}}$    (long_name='output gap')
pihat  ${\hat{\pi}}$  (long_name='inflation deviation from trend')
yhat   ${\hat{y}}$    (long_name='output deviations from trend')
ghat   ${\hat{g}}$    (long_name='output growth')
rhat   ${\hat{r}}$    (long_name='interest deviations from trend')
;

varobs ghat rhat pihat;

varexo
eta_a  ${\varepsilon_a}$  (long_name='preference innovation')
eta_e  ${\varepsilon_e}$  (long_name='(negative) cost push innovation')
eta_z  ${\varepsilon_z}$  (long_name='TFP innovation')
eta_r  ${\varepsilon_r}$  (long_name='monetary policy innovation')
;

parameters
BETA      ${\beta}$         (long_name='discount factor')
PSI       ${\psi}$          (long_name='output gap coefficient in Phillips curve')
ALPHA_PI  ${\alpha_{\pi}}$  (long_name='slope parameter in Phillips curve')
ALPHA_X   ${\alpha_{x}}$    (long_name='slope parameter in IS curve')
RHO_A     ${\rho_{a}}$      (long_name='persistence preference shock')
RHO_E     ${\rho_{e}}$      (long_name='persistence (negative) cost-push shock')
OMEGA     ${\omega}$        (long_name='scale parameter preference innovation in IS curve')
RHO_PI    ${\rho_{\pi}}$    (long_name='feedback policy rule inflation')
RHO_G     ${\rho_{g}}$      (long_name='feedback policy rule output growth')
RHO_X     ${\rho_{x}}$      (long_name='feedback policy rule output gap')
SIGN_SHOCKS
;

% fixed parameters
BETA = 0.99;
PSI = 0.1;
ALPHA_X = 0;
ALPHA_PI = 0;
SIGN_SHOCKS = 1;

model;
[name='temporary preference shock (15)']
ahat = RHO_A*ahat(-1) + SIGN_SHOCKS*eta_a/100;
[name='temporary cost-push shock (16)']
ehat = RHO_E*ehat(-1) + SIGN_SHOCKS*eta_e/100;
[name='technology shock (17)']
zhat = SIGN_SHOCKS*eta_z/100;
[name='New Keynesian IS curve (23)']
xhat = ALPHA_X*xhat(-1) + (1-ALPHA_X)*xhat(+1) - (rhat-pihat(+1)) + (1-OMEGA)*(1-RHO_A)*ahat;
[name='New Keynesian PC curve (24)']
pihat = BETA * ( ALPHA_PI*pihat(-1) + (1-ALPHA_PI)*pihat(+1) ) + PSI*xhat - ehat;
[name='output gap (20)']
xhat = yhat - OMEGA*ahat;
[name='growth rate of output (21)']
ghat = yhat - yhat(-1) + zhat;
[name='policy rule (22)']
rhat - rhat(-1) = RHO_PI*pihat + RHO_G*ghat + RHO_X*xhat + SIGN_SHOCKS*eta_r/100;
end;

steady_state_model;
ahat = 0; ehat = 0; zhat = 0; xhat = 0; pihat = 0; yhat = 0; ghat = 0; rhat = 0;
end;

@#ifdef BAYESIAN
estimated_params;
stderr eta_a, ,  0, 10, inv_gamma_pdf, (sqrt(30)),   (sqrt(30)) ;
stderr eta_e, ,  0, 10, inv_gamma_pdf, (sqrt(0.08)), (sqrt(1))  ;
stderr eta_z, ,  0, 10, inv_gamma_pdf, (sqrt(5)),    (sqrt(15)) ;
stderr eta_r, ,  0, 10, inv_gamma_pdf, (sqrt(0.50)), (sqrt(2))  ;
OMEGA,        ,  0,  1, beta_pdf,      0.20,         0.10       ;
RHO_PI,       ,  0,  1, gamma_pdf,     0.30,         0.10       ;
RHO_G,        ,  0,  1, gamma_pdf,     0.30,         0.10       ;
RHO_X,        ,  0,  1, gamma_pdf,     0.25,         0.0625     ;
RHO_A,        ,  0,  1, beta_pdf,      0.85,         0.10       ;
RHO_E,        ,  0,  1, beta_pdf,      0.85,         0.10       ;
  @#if DISTRIB == "CSN"
skew eta_a,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_e,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_z,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
skew eta_r,   ,   ,   , beta_pdf,         0,         0.40       , -1, 1;
  @#endif
end;
@#endif