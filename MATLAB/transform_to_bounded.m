function X = transform_to_bounded(Y,a,b)
% function X = transform_to_bounded(Y,a,b)
% -------------------------------------------------------------------------
% transforms an unbounded parameter vector back to original bounded parameters
% using an inverse log-odds transform as described at
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
% -------------------------------------------------------------------------
% INPUTS
% - Y      [double vector]   unbounded parameter vector
% - a      [double vector]   lower bound of original parameters
% - b      [double vector]   upper bound of original parameters
% -------------------------------------------------------------------------
% OUTPUTS
% - X      [double vector]   original bounded parameter vector
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

X = a + (b-a).*inv_logit(Y);

function u = inv_logit(v)
    u = 1./(1+exp(-v));
end %inv_logit

end % transform_to_unbounded