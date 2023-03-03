GAUSSIAN MAXIMUM LIKELIHOOD

BEST Optimizer fminsearch: value of maximimized log-likelihood function: 2648.4250160723.
Time required: 8.7791 seconds 
                        Estimate       s.d.       t-stat 
                        _________    _________    _______

    OMEGA                0.061727     0.064742    0.95343
    ALPHA_X              0.083579      0.10921    0.76529
    RHO_PI                0.35968      0.04687      7.674
    RHO_G                 0.25361     0.039082     6.4891
    RHO_X                0.034739     0.014961      2.322
    RHO_A                 0.94702     0.024666     38.394
    RHO_E                 0.96254     0.024575     39.168
    sqrt_Sigma_eta_a     0.040482      0.01539     2.6304
    sqrt_Sigma_eta_e     0.001238    0.0002607     4.7488
    sqrt_Sigma_eta_z      0.01086    0.0028167     3.8555
    sqrt_Sigma_eta_r    0.0031113    0.0003437     9.0524


Running grid search for CSN parameters:
100%[===================================================]

best_neg_log_likelihood_grid_1 =

   1.0e+03 *

   -4.9402   -4.1385   -4.1256   -3.9708


diag_Sigma_eta_grid_1 =

    0.0041    0.0024    0.0021    0.0021
    0.0000    0.0000    0.0000    0.0000
    0.0003    0.0002    0.0002    0.0002
    0.0000    0.0000    0.0000    0.0000


diag_Gamma_eta_grid_1 =

   1.0e+05 *

    0.0006   -0.0002   -0.0002   -0.0002
   -1.2021   -0.0055   -0.0055   -0.0070
    0.0017    0.0008    0.0008    0.0006
   -0.4653    0.4653    0.4653    0.4653


var_eta_grid_1 =

    0.0016    0.0016    0.0016    0.0016
    0.0000    0.0000    0.0000    0.0000
    0.0001    0.0001    0.0001    0.0001
    0.0000    0.0000    0.0000    0.0000


skew_eta_grid_1 =

    0.7819   -0.1422   -0.0711   -0.0711
   -0.9952   -0.0711   -0.0711   -0.1422
    0.6398    0.1422    0.1422    0.0711
   -0.9952    0.9952    0.9952    0.9952

 
OPTIMIZE OVER GAMMA_ETA AND SIGMA_ETA (KEEPING MODEL PARAMS FIXED)

Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -5363.576030 

 
Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -4653.088906 

 
Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -4607.463573 

 
Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -4432.413924 

IdleTimeout has been reached.
Parallel pool using the 'Processes' profile is shutting down.
neg_log_likelihood_1

neg_log_likelihood_1 =

   1.0e+03 *

   -5.3636   -4.6531   -4.6075   -4.4324

xparam_1

xparam_1 =

   1.0e+05 *

    0.0006   -0.0002   -0.0002   -0.0001
   -1.0587   -0.0062   -0.0061   -0.0103
    0.0018    0.0008    0.0008    0.0007
   -0.5114    0.5115    0.5970    0.4991
    0.0000    0.0000    0.0000    0.0000
    0.0000    0.0000    0.0000    0.0000
    0.0000    0.0000    0.0000    0.0000
    0.0000    0.0000    0.0000    0.0000

best_neg_log_likelihood_grid_1

best_neg_log_likelihood_grid_1 =

   1.0e+03 *

   -4.9402   -4.1385   -4.1256   -3.9708

MODEL_2 = MODEL; MODEL_2.param_estim_names = MODEL_0.param_estim_names; OPT_2 = OPT;
xparam_2 = nan(MODEL_0.param_estim_nbr+MODEL_2.exo_nbr,best_of) ; neg_log_likelihood_2 = nan(1,best_of);
parfor jbest=1:best_of    
    [xparam_2_j, PARAM_2_j, ESTIM_PARAM_2_j, MODEL_2_j, OPT_2_j] = feval(str2func(OPT_2.modelname + "_params"),  2, MODEL_2, OPT_2, xparam_0, sqrt(diag_Sigma_eta_grid_1(:,jbest)), diag_Gamma_eta_grid_1(:,jbest));    
    [xparam_2(:,jbest),neg_log_likelihood_2(jbest)] = fminsearch(@(x) negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2_j,MODEL_2_j,OPT_2_j),  xparam_2_j,OPT_2_j.optimizer.optim_options);
end
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to the parallel pool (number of workers: 24).
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
> In get_first_order_perturbation_solution (line 70)
In negative_log_likelihood_dsge (line 93)
In @(x)negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2_j,MODEL_2_j,OPT_2_j)
In fminsearch (line 267)
In parallel_function>make_general_channel/channel_general (line 837)
In remoteParallelFunction (line 67)
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
> In get_first_order_perturbation_solution (line 70)
In negative_log_likelihood_dsge (line 93)
In @(x)negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2_j,MODEL_2_j,OPT_2_j)
In fminsearch (line 325)
In parallel_function>make_general_channel/channel_general (line 837)
In remoteParallelFunction (line 67)
 
Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -5393.914963 

> In get_first_order_perturbation_solution (line 70)
In negative_log_likelihood_dsge (line 93)
In @(x)negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2_j,MODEL_2_j,OPT_2_j)
In fminsearch (line 267)
In parallel_function>make_general_channel/channel_general (line 837)
In remoteParallelFunction (line 67)
Warning: Blanchard-Khan order condition not fullfilled: No equilibrium exists.
> In get_first_order_perturbation_solution (line 70)
In negative_log_likelihood_dsge (line 93)
In @(x)negative_log_likelihood_dsge(x,DATA.MAT,PARAM_2_j,MODEL_2_j,OPT_2_j)
In fminsearch (line 325)
In parallel_function>make_general_channel/channel_general (line 837)
In remoteParallelFunction (line 67)
 
Exiting: Maximum number of iterations has been exceeded
         - increase MaxIter option.
         Current function value: -4543.210779 
