function Y = transform_to_unbounded(X,a,b)
% function Y = transform_to_unbounded(X,a,b)
% -------------------------------------------------------------------------
% transforms an original bounded parameter vector to unbounded parameters
% using a log-odds transform as described at
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
% -------------------------------------------------------------------------
% INPUTS
% - X      [double vector]   original bounded parameter vector
% - a      [double vector]   lower bound of original parameters
% - b      [double vector]   upper bound of original parameters
% -------------------------------------------------------------------------
% OUTPUTS
% - Y      [double vector]   unbounded parameter vector
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
Y = logit((X-a)./(b-a));

function logit_u = logit(u)
    logit_u = log(u./(1-u));
end % logit

end % transform_to_unbounded