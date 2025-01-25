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

@#define TRANSFORM_PARAMETERS    = 0  // new feature to transform selected bounded parameters to unbounded domain using logit transform (only during optimization)
@#include "_ireland2004_common.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAYESIAN ESTIMATION WITH RWMH SAMPLER OF GAUSSIAN VERSION OF MODEL USING PSKF TO COMPUTE LIKELIHOOD %
% RWMH SAMPLER IS INITIALIZED AT MODE AND COVARIANCE MATRIX FROM OPTIMIZATION (MODE_COMPUTE=5)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end;
% structure estim_params_.skew_exo needs to be created by the preprocessor
estim_params_.skew_exo = zeros(0, 11);
@#include "__transformParameters.inc" // optional new feature to use logit transform of bounded parameters during optimization only

options_.kalman.pskf.use_in_gaussian_case = true; % use PSKF filter even for Gaussian shocks

estimation(datafile = 'data/ireland2004_data.m'
          , lik_init = 1
          , use_univariate_filters_if_singularity_is_detected = 0
          , keep_kalman_algo_if_singularity_is_detected
          , mh_replic = 0
          , mode_compute = 8
          , mode_file = 'ireland2004_bayes_3_mode_slice_gaussian/Output/ireland2004_bayes_3_mode_slice_gaussian_mh_mode.mat'
          , plot_priors = 0
          );
copyfile('ireland2004_bayes_7_mode_gaussian/Output/ireland2004_bayes_7_mode_gaussian_mode.mat','ireland2004_bayes_7_mode_gaussian/Output/ireland2004_bayes_7_mode_gaussian_mode8.mat');

estimation(datafile = 'data/ireland2004_data.m'
          , lik_init = 1
          , use_univariate_filters_if_singularity_is_detected = 0
          , keep_kalman_algo_if_singularity_is_detected
          , mh_replic = 0
          , mode_compute = 6
          , mode_file = 'ireland2004_bayes_7_mode_gaussian/Output/ireland2004_bayes_7_mode_gaussian_mode8.mat'
          , plot_priors = 0
          );
copyfile('ireland2004_bayes_7_mode_gaussian/Output/ireland2004_bayes_7_mode_gaussian_mode.mat','ireland2004_bayes_7_mode_gaussian/Output/ireland2004_bayes_7_mode_gaussian_mode6.mat');