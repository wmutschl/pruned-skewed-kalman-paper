function hms = dynsec2hms(secs)
% function hms = dynsec2hms(secs)
% -------------------------------------------------------------------------
% Converts a number of seconds (given e.g. by tic toc) into a hours-minutes-seconds string
% -------------------------------------------------------------------------
% INPUTS
% secs  [scalar]  seconds
% -------------------------------------------------------------------------
% OUTPUTS
% hms   [string]   human readable time in hours - minutes - seconds format
% =========================================================================
% Copyright (C) 2008-2009 Dynare Team
% Copyright (C) 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
secs = round(secs);
s = rem(secs, 60);
m = rem(floor(secs / 60), 60);
h = floor(secs / 3600);
hms = sprintf('%dh%02dm%02ds', h, m, s);