function [LIK,lik] = auxiliary_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% [LIK,lik] = auxiliary_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% Evaluates the likelihood of a nonlinear model with the auxiliary particle filter
% allowing eventually resampling.
% INPUTS
%  - ReducedForm            [structure] decision rules
%  - Y                      [double]    dataset
%  - start                  [integer]   first observation for likelihood evaluation
%  - ParticleOptions        [structure] filter options
%  - ThreadsOptions         [structure] options for threading of mex files
%  - options_               [structure] describing the options
%  - M_                     [structure] describing the model
%
% OUTPUTS
% - LIK                [double]    scalar, likelihood
% - lik                [double]    (T-s+1)×1 vector, density of observations in each period.

% Copyright © 2011-2025 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% Set default
if isempty(start)
    start = 1;
end
if ParticleOptions.resampling.method.smooth
    error('auxiliary_particle_filter: resampling_method=smooth is not supported.')
end
% Get perturbation order
mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
sample_size = size(Y,2);
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);
number_of_particles = ParticleOptions.number_of_particles;

% Get initial condition for the state vector.
if options_.particle.use_reduced_rank_cholesky
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
else
    StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
end
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
Q_lower_triangular_cholesky = chol(ReducedForm.Q)';

% Set seed for randn().
set_dynare_seed_local_options([],false,'default');

%initialize output
lik  = NaN(sample_size,1);
LIK  = NaN;
if isempty(ReducedForm.H)
    ReducedForm.H = 0;
end
% filter out singular measurement error case
if rcond(ReducedForm.H) < 1e-12
    LIK = NaN;
    return
end

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables+log(det(ReducedForm.H));

% Initialization of the weights across particles.
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),ReducedForm.StateVectorMean);
if ParticleOptions.pruning && ~(options_.order==1)
    if options_.order == 2
        StateVectors_ = StateVectors;
        mf0_ = mf0;
    elseif options_.order == 3
        StateVectors_ = repmat(StateVectors,3,1);
        mf0_ = repmat(mf0,1,3); 
        mask2 = number_of_state_variables+1:2*number_of_state_variables;
        mask3 = 2*number_of_state_variables+1:3*number_of_state_variables;
        mf0_(mask2) = mf0_(mask2)+size(ghx,1);
        mf0_(mask3) = mf0_(mask3)+2*size(ghx,1);
    else
        error('Pruning is not available for orders > 3');
    end
else
    StateVectors_=[];
    mf0_ = mf0;
end

for t=1:sample_size
    tmp=iterate_law_of_motion(StateVectors,zeros(number_of_structural_innovations,number_of_particles),ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning,StateVectors_);
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    z = sum(PredictionError.*(ReducedForm.H\PredictionError),1) ;
    ddl = 3 ;
    tau_tilde = weights.*(exp(gammaln((ddl + 1) / 2) - gammaln(ddl/2))./(sqrt(ddl*pi).*(1 + (z.^2)./ddl).^((ddl + 1)/2))+1e-99) ;
    tau_tilde = tau_tilde/sum(tau_tilde) ;
    indx = resample(0,tau_tilde',ParticleOptions);
    weights_stage_1 = weights(indx)./tau_tilde(indx) ;
    if nnz(M_.Skew_e) > 0 % draw from skew normal distribution (special case of closed skew normal, see csn_update_specification.m for details)
        epsilon = rand_multivariate_csn(number_of_particles, M_.csn.mu_e, M_.csn.Sigma_e, M_.csn.Gamma_e, M_.csn.nu_e, M_.csn.Delta_e);
    else % draw from Gaussian distribution
        epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    end
    if ParticleOptions.pruning
        [tmp, tmp_]=iterate_law_of_motion(StateVectors(:,indx),epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning,StateVectors_(:,indx));
        StateVectors_ = tmp_(mf0_,:);
    else 
        [tmp]=iterate_law_of_motion(StateVectors(:,indx),epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning);
    end
    StateVectors = tmp(mf0,:);
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    weights_stage_2 = weights_stage_1.*(exp(-.5*(const_lik+sum(PredictionError.*(ReducedForm.H\PredictionError),1))) + 1e-99) ;
    lik(t) = log(mean(weights_stage_2)) ;
    weights = weights_stage_2/sum(weights_stage_2);
    if (ParticleOptions.resampling.status.generic && neff(weights)<ParticleOptions.resampling.threshold*sample_size) || ParticleOptions.resampling.status.systematic
        if ParticleOptions.pruning
            temp = resample([StateVectors' StateVectors_'],weights',ParticleOptions);
            StateVectors = temp(:,1:number_of_state_variables)';
            StateVectors_ = temp(:,number_of_state_variables+1:end)';
        else
            StateVectors = resample(StateVectors',weights',ParticleOptions)';
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    end
end

LIK = -sum(lik(start:end));
