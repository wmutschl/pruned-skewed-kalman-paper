function log_cdf = logmvncdf_ME(Zj, R, cutoff)
% log_cdf = logmvncdf_ME(Zj, Corr_mat, cutoff)
% -------------------------------------------------------------------------
% Approximates Gaussian log(CDF) function according to Mendell and Elston (1974)
% -------------------------------------------------------------------------
% INPUTS
% - Zj       [n by 1]   column vector of points where Gaussian CDF is evaluated at
% - R        [n by n]   correlation matrix
% - cutoff   [2 by 1]   optional threshold points at which values in Zj are too low/high to be evaluated, defaults to [6 38]
% -------------------------------------------------------------------------
% OUTPUTS
% log_cdf    [double]   approximate value of Gaussian log(CDF)
% =========================================================================
% Copyright (C) 2015 Dietmar Bauer (original implementation)
% Copyright (C) 2022-2023 Gaygysyz Guljanov
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

if nargin < 3
    cutoff = [6 38]; % default values
end
% cutoff tails as these are too low to be evaluated
Zj(Zj>cutoff(1))  = cutoff(1);
Zj(Zj<-cutoff(1)) = -cutoff(1);
n = length(Zj); % remaining dimension of Z

% first element
cdf_val = phid(Zj(1));
pdf_val = phip(Zj(1));
log_cdf = log(cdf_val); % perform all calculations in logs

for jj = 1:(n-1)
    ajjm1 = pdf_val / cdf_val;

    % Update Zj and Rij
    tZ = Zj + ajjm1 * R(:, 1); % update Zj
    R_jj = R(:, 1) * R(1, :);
    tRij = R - R_jj * ( ajjm1 + Zj(1) ) * ajjm1; % update Rij

    % Convert Rij (i.e. Covariance matrix) to Correlation matrix
    COV_TO_CORR = sqrt( diag(tRij) );
    Zj = tZ ./ COV_TO_CORR;
    R = tRij ./ (COV_TO_CORR * COV_TO_CORR');

    % Cutoff those dimensions if they are too low to be evaluated
    Zj(Zj>cutoff(2))  = cutoff(2);  % use larger cutoff
    Zj(Zj<-cutoff(2)) = -cutoff(2); % use larger cutoff

    % Evaluate jj's probability
    cdf_val = phid(Zj(2));
    pdf_val = phip(Zj(2));

    % Delete unnecessary parts of updated Zj and Corr_mat
    last_el = n - jj + 1;
    Zj = Zj(2:last_el);
    R = R(2:last_el, 2:last_el);

    % Overall probability
    log_cdf = log_cdf + log(cdf_val);

end

end % logmvncdf_ME


%%%%%%%%%%%%%%
% Normal PDF %
%%%%%%%%%%%%%%
function p = phip(z)
    p = exp(-z.^2/2)/sqrt(2*pi);
end


%%%%%%%%%%%%%%
% Normal CDF %
%%%%%%%%%%%%%%
function p = phid(z)
    p = erfc( -z/sqrt(2) )/2;
end
