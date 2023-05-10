function COMBOS = create_combos(n,m)
% function COMBOS = create_combos(n,m)
% -------------------------------------------------------------------------
% creates all n tupels (combinations) of numbers from 1 to m
% -------------------------------------------------------------------------
% INPUTS
% - n   [integer]        number of elements in tupel
% - m   [even integer]   highest number
% -------------------------------------------------------------------------
% OUTPUTS
% - COMBOS   [matrix]   all combinations
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
% - not efficient and only for n<=4, needs to be improved and made more general

COMBOS = zeros((m-1)^n,n,'int16');
idx = 1;
for j1=1:(m-1)
    if n==1
        COMBOS(idx,1) = j1;
        idx=idx+1;
    else
        for j2=1:(m-1)
            if n==2
                COMBOS(idx,1) = j1; COMBOS(idx,2) = j2;
                idx=idx+1;
            else
                for j3=1:(m-1)
                    if n==3
                        COMBOS(idx,1) = j1; COMBOS(idx,2) = j2; COMBOS(idx,3) = j3;
                        idx=idx+1;
                    else
                        for j4=1:(m-1)
                            if n==4
                                COMBOS(idx,1) = j1; COMBOS(idx,2) = j2; COMBOS(idx,3) = j3; COMBOS(idx,4) = j4;
                                idx=idx+1;
                            else
                                error('create_combos: please adjust the function for more shocks than 4; there is probably a better way than using for loops for this')
                            end
                        end
                    end
                end
            end
        end
    end
end
