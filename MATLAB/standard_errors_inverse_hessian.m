function xstderr = standard_errors_inverse_hessian(objfct, xparams_untransformed, bounds_untransformed, varargin)
% function xstderr = standard_errors_inverse_hessian(objfct, xparams_untransformed, bounds_untransformed, varargin)
% -------------------------------------------------------------------------
% computes standard errors via inverse hessian method at different
% parameter vectors in parallel
% -------------------------------------------------------------------------
% INPUTS
% - objfct                 [function handle]  name of objective function
% - xparams_untransformed  [nparam x ncols]   different vectors of original (untransformed) parameters for which to compute standard errors
% - bounds_untransformed   [nparam x 2]       lower and upper bounds for original (untransformed) parameters
% - varargin               [cell]             additional arguments passed to objfct
% -------------------------------------------------------------------------
% OUTPUTS
% - xstderr                [nparam x ncols]   standard errors
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
xstderr = nan(size(xparams_untransformed));
parfor jx = 1:size(xparams_untransformed,2)
    xparams = xparams_untransformed(:,jx);
    gstep = max(abs(xparams), sqrt(1e-3)*ones(length(xparams), 1))*eps^(1/6);
    H = get_hessian_with_bounds(objfct, xparams, gstep, bounds_untransformed,    varargin{:});
    V = inv(reshape(H,length(xparams),length(xparams)));
    xstderr(:,jx) = sqrt(diag(V));
end