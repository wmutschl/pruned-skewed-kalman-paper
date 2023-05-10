% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede

function hessian_mat = get_hessian_with_bounds(func, x, gstep, bounds, varargin)
% Computes second order partial derivatives using two-sided finite
% difference method.
% based on get_hessian.m from Dynare with the following changes:
% - if a step up or down of a parameter passes the lower or upper bounds,
%   then the partial derivatives are computed using one-sided finite
%   differences.
%
% INPUTS
%    func        [string]   name of the function
%    x           [double]   vector, the Hessian of "func" is evaluated at x.
%    gstep       [double]   scalar, size of epsilon.
%    varargin    [void]     list of additional arguments for "func".
%
% OUTPUTS
%    hessian_mat [double]   Hessian matrix
%
% ALGORITHM
%    Uses Abramowitz and Stegun (1965) formulas 25.3.23
% \[
%     \frac{\partial^2 f_{0,0}}{\partial {x^2}} = \frac{1}{h^2}\left( f_{1,0} - 2f_{0,0} + f_{ - 1,0} \right)
% \]
% and 25.3.27 p. 884
%
% \[
%     \frac{\partial ^2f_{0,0}}{\partial x\partial y} = \frac{-1}{2h^2}\left(f_{1,0} + f_{-1,0} + f_{0,1} + f_{0,-1} - 2f_{0,0} - f_{1,1} - f_{-1,-1} \right)
% \]
%
% SPECIAL REQUIREMENTS
%    none
%

% Copyright (C) 2001-2017 Dynare Team
% Copyright (C) 2023 Willi Mutschler
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

if ~isa(func, 'function_handle')
    func = str2func(func);
end
ub = bounds(:,2);
lb = bounds(:,1);
n   = size(x,1);
%h1  = max(abs(x), sqrt(gstep(1))*ones(n, 1))*eps^(1/6)*gstep(2);
h1 = gstep;
h_1 = h1;
xh1 = x+h1;
h1  = xh1-x;
xh1 = x-h_1;
h_1 = x-xh1;
xh1 = x;
f0  = feval(func, x, varargin{:});
f1  = zeros(size(f0, 1), n);
f_1 = f1;

for i=1:n
    %do step up
    xh1(i)   = x(i)+h1(i);
    if xh1(i) <= ub(i)        
        f1(:,i)  = feval(func, xh1, varargin{:});
    else
        %fprintf('  - don''t do step up for parameter %d\n',i)
        f1(:,i)  = f0;
    end
    %do step down
    xh1(i)   = x(i)-h_1(i);
    if xh1(i) >= lb(i)        
        f_1(:,i) = feval(func, xh1, varargin{:});
    else
        %fprintf('  - don''t do step down for parameter %d\n',i)
        f_1(:,i) = f0;
    end
    %reset parameter
    xh1(i)   = x(i);
end

xh_1 = xh1;
temp = f1+f_1-f0*ones(1, n); %term f_(1,0)+f_(-1,0)-f_(0,0) used later

hessian_mat = zeros(size(f0,1), n*n);

for i=1:n
    if i > 1
        %fill symmetric part of Hessian based on previously computed results
        k = [i:n:n*(i-1)];
        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1) = hessian_mat(:,k);
    end
    hessian_mat(:,(i-1)*n+i) = (f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i)); %formula 25.3.23
    for j=i+1:n
        %step in up direction
        xh1(i) = x(i)+h1(i);  if xh1(i) > ub(i); xh1(i) = x(i); end %fprintf('  - don''t do cross step up for parameter %d\n',i); end
        xh1(j) = x(j)+h_1(j); if xh1(j) > ub(j); xh1(j) = x(j); end %fprintf('  - don''t do cross step up for parameter %d\n',j); end
        %step in down direction
        xh_1(i) = x(i)-h1(i);  if xh_1(i) < lb(i); xh_1(i) = x(i); end %fprintf('  - don''t do cross step down for parameter %d\n',i); end
        xh_1(j) = x(j)-h_1(j); if xh_1(j) < lb(j); xh_1(j) = x(j); end %fprintf('  - don''t do cross step down for parameter %d\n',j); end
        hessian_mat(:,(i-1)*n+j) =-(-feval(func, xh1, varargin{:})-feval(func, xh_1, varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j)); %formula 25.3.27
        %reset grid points
        xh1(i)  = x(i);
        xh1(j)  = x(j);
        xh_1(i) = x(i);
        xh_1(j) = x(j);
    end
end