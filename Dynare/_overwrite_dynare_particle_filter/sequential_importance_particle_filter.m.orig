function [LIK,lik] = sequential_importance_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% [LIK,lik] = sequential_importance_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% Evaluates the likelihood of a nonlinear model with a particle filter employing a sequential importance sampling approach 
%
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
%
% References:
% Implementation is e.g. described in Godsill/Doucet/West (2004): "Monte Carlo Smoothing for Nonlinear Time Series", 
% Journal of the American Statistical Association, March 2004, 99(465)


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

% Set default value for start
if isempty(start)
    start = 1;
end

sample_size = size(Y,2);
number_of_state_variables = length(ReducedForm.mf0);
number_of_observed_variables = length(ReducedForm.mf1);
number_of_structural_innovations = length(ReducedForm.Q);
number_of_particles = ParticleOptions.number_of_particles;

if isempty(ReducedForm.H)
    ReducedForm.H = 0;
end

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables +log(det(ReducedForm.H)) ;
lik  = NaN(sample_size,1);

% filter out singular measurement error case
if rcond(ReducedForm.H) < 1e-12
    LIK = NaN;
    return
end

StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';%reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2); % Get the rank of StateVectorVarianceSquareRoot

% Factorize the covariance matrix of the structural innovations
Q_lower_triangular_cholesky = chol(ReducedForm.Q)';

% Set seed for randn().
set_dynare_seed_local_options([],false,'default');

% Initialization of the weights across particles.
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),ReducedForm.StateVectorMean);
if ParticleOptions.pruning && ~(options_.order ==1)
    if options_.order == 2
        StateVectors_ = StateVectors;
        mf0_ = ReducedForm.mf0;
    elseif options_.order == 3
        StateVectors_ = repmat(StateVectors,3,1);
        mf0_ = repmat(ReducedForm.mf0,1,3); 
        mask2 = number_of_state_variables+1:2*number_of_state_variables;
        mask3 = 2*number_of_state_variables+1:3*number_of_state_variables;
        mf0_(mask2) = mf0_(mask2)+size(ReducedForm.ghx,1);
        mf0_(mask3) = mf0_(mask3)+2*size(ReducedForm.ghx,1);
    else
        error('Pruning is not available for orders > 3');
    end
else
    StateVectors_ = StateVectors;
    mf0_ = ReducedForm.mf0;
end

% Loop over observations
for t=1:sample_size
    epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    if ParticleOptions.pruning
        [tmp, tmp_]=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning,StateVectors_);
    else
        tmp=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning);
    end

    PredictionError = bsxfun(@minus,Y(:,t),tmp(ReducedForm.mf1,:));

    lnw = -.5*(const_lik+sum(PredictionError.*(ReducedForm.H\PredictionError),1));
    
    dfac = max(lnw);
    wtilde = weights.*exp(lnw-dfac);
    lik(t) = log(sum(wtilde))+dfac;
    weights = wtilde/sum(wtilde);
    if (ParticleOptions.resampling.status.generic && neff(weights)<ParticleOptions.resampling.threshold*sample_size) || ParticleOptions.resampling.status.systematic
        if ParticleOptions.pruning
            temp = resample([tmp(ReducedForm.mf0,:)' tmp_(mf0_,:)'],weights',ParticleOptions);
            StateVectors = temp(:,1:number_of_state_variables)';
            StateVectors_ = temp(:,number_of_state_variables+1:end)';
        else
            StateVectors = resample(tmp(ReducedForm.mf0,:)',weights',ParticleOptions)';
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif ParticleOptions.resampling.status.none
        StateVectors = tmp(ReducedForm.mf0,:);
        if ParticleOptions.pruning
            StateVectors_ = tmp_(mf0_,:);
        end
    end
end

LIK = -sum(lik(start:end));
