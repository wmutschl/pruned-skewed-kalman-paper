function dX_dY = transform_to_bounded_derivative(Y,a,b)
% function dX_dY = transform_to_bounded_derivative(Y,a,b)
% -------------------------------------------------------------------------
% Jacobian of inverse logit transform as described at
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
% -------------------------------------------------------------------------
% INPUTS
% - Y      [double vector]   unbounded parameter vector
% - a      [double vector]   lower bound of original parameters
% - b      [double vector]   upper bound of original parameters
% -------------------------------------------------------------------------
% OUTPUTS
% - dX_dY  [double vector]   computes Jacobian of original bounded parameter vector
%                            wrt to unbounded parameters
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
dX_dY = (b-a).*inv_logit(Y).*(1-inv_logit(Y));

function u = inv_logit(v)
    u = 1./(1+exp(-v));
end %inv_logit

end % transform_to_unbounded