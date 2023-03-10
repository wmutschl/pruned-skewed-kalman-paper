function SIGy = dlyapdoubling(A,SIGu)
% function SIGy = dlyapdoubling(A,SIGu)
% -------------------------------------------------------------------------
% Solves the Lyapunov equation SIGy = A*SIGy*A' + SIGu using the doubling algorithm
% -----------------------------------------------------------------------
% INPUTS
%   - A     [n by n]   square matrix, usually autoregressive or state space matrix
%   - SIGu  [n by n]   square matrix, usually usually covariance matrix
% -----------------------------------------------------------------------
% OUTPUT
%	- SIGy  [n by n]   square matrix that solves the Lyapunov equation
% =========================================================================
% Copyright (C) 2022-2023 Willi Mutschler
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
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================

max_iter   = 500;
A_old      = A;
SIGu_old   = SIGu;
SIGy_old   = eye(size(A));
difference = .1;
index1     = 1;
tol        = 1e-25;
while (difference > tol) && (index1 < max_iter)
    SIGy       = A_old*SIGy_old*transpose(A_old) + SIGu_old;
    difference = max(abs(SIGy(:)-SIGy_old(:)));
    SIGu_old   = A_old*SIGu_old*transpose(A_old) + SIGu_old;
    A_old      = A_old*A_old;
    SIGy_old   = SIGy;
    index1     = index1 + 1;
end

end %function end