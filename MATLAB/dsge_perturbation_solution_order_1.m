function [gx,gu,error_indicator] = dsge_perturbation_solution_order_1(M_)
% function [gx,gu,error_indicator] = dsge_perturbation_solution_order_1(M_)
% -------------------------------------------------------------------------
% computes the first-order perturbation solution of a DSGE model using
% the generalized Schur decomposition. Inspired by Dynare's first order
% solver, as described in Villemot (2011), but without taking special care
% of static variables and adding zero rows for non-state variables to have
% the following solution for all endogenous variables y(t):
%   y(t) = gx*y(t-1) + gu*eta(t)
% -------------------------------------------------------------------------
% INPUTS
% - M_         [structure]   model information
% - params_    [structure]   information on calibrated parameters
% - options_   [structure]   options
% -------------------------------------------------------------------------
% OUTPUTS
% - gx                [endo_nbr by nspred]    perturbation solution matrix with resepct to states
% - gu                [endo_nbr by exo_nbr]   perturbation solution matrix with resepct to exogenous
% - error_indicator   [boolean]     indicator if something was wrong
% =========================================================================
% Copyright © 2023 Willi Mutschler
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
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

%% Compute steady-state
[steadystate,params,error_indicator] = feval(str2func(M_.fname + '_steady_state'),[],[],M_);

if error_indicator
    gx = []; gu = [];
    return
end

%% Evaluate first dynamic Jacobian at steady-state
g1 = feval(str2func(M_.fname + '_g1'),steadystate,params);
if any(isnan(g1(:))) || any(isinf(g1(:))) || any(imag(g1(:)))
    error_indicator = 1;
    gx = []; gu = [];
    warning('Something wrong with the perturbation solution due to nan, inf or imaginary numbers');    
    return
end

%% Recovering gx
n = M_.endo_nbr;
D = [zeros(n,n)  g1(:,(2*n+1):(3*n));
    eye(n)       zeros(n,n)];
E = [-g1(:,1:n)  -g1(:,(n+1):(2*n));
    zeros(n,n)   eye(n)];

realsmall=1e-7;
try
    [S,T,Q,Z] = qz(D,E); % upper triangular factorization of the matrix pencil
catch
    error_indicator = 1;    
    warning('Error using qz');
    gx = []; gu = [];
    return
end
%norm(D-Q'*S*Z')
%norm(E-Q'*T*Z')
% Rule out that both s_ii and t_ii are zero
zxz = sum((abs(diag(S))<realsmall) & (abs(diag(T))<realsmall));
if ~(~zxz)
    error_indicator = 1;
	warning('Something wrong in the solution: Coincident zeros');
    gx = []; gu = [];
    return
end
%   If S is quasi-triangular, the diagonal elements of S and T,
%       sii = diag(S), tii = diag(T),
%   are the generalized eigenvalues that satisfy
%       D*V*diag(tii) = E*V*diag(sii)
%       diag(tii)*W'*D = diag(sii)*W'*E
%   where matrices V's and W' columns are the generalized eigenvectors
%   The eigenvalues produced by
%       lambda = eig(D,E)
%   are the ratios of the sii and tii.
%       [lambda  sii./tii]
% disp([eig(D,E) diag(S)./diag(T)]);
% reordering such that stable (smaller than one) generalized Eigenvalues 
% of E w.r.t. D is in the upper left corner of T and S, abs(tii)>abs(sii)
[S,T,Q,Z] = ordqz(S,T,Q,Z,'udo');
%disp(abs(diag(S))./abs(diag(T)));

% Blanchard-Khan order condition
if abs(T(n,n))>abs(S(n,n))
    error_indicator = 1;
	%warning('Blanchard-Khan order condition not fullfilled: No equilibrium exists.');
    gx = []; gu = [];
    return
elseif abs(T(n+1,n+1))<abs(S(n+1,n+1))
    error_indicator = 1;
	%warning('Blanchard-Khan order condition not fullfilled: Indeterminacy.');
    gx = []; gu = [];
    return
end

Z = Z'; % different to our notes
%z11 = Z(1:n,1:n);
%z12 = Z(1:n,n+1:end);
z22 = Z(n+1:end,n+1:end);
z21 = Z(n+1:end,1:n);

% Blanchard-Khan rank condition:
if rank(z22)<n
    error_indicator = 1;
    %warning('Blanchard-Khan rank condition not fullfilled: no stable solution');
    gx = []; gu = [];
    return
end
z22i = z22\eye(n);
gx = -z22i*z21;
gx = real(gx); % real function takes away very small imaginary parts of the solution

%% Recovering gu
nu = M_.exo_nbr;
gu = -(g1(:,(n+1):(2*n)) + g1(:,(2*n+1):(3*n)) * gx) \ g1(:,(3*n+1):(3*n+nu));