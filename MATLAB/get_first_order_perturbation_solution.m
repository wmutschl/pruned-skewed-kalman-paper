function [SOL,error_indicator] = get_first_order_perturbation_solution(MODEL,PARAM)
%% Compute steady-state
[STEADY_STATE,PARAM,error_indicator] = feval(str2func(MODEL.name + '_steady_state'),MODEL,PARAM);
if error_indicator
    SOL = [];
    return
else
    % put steady-state into DR order
    y_ss = struct2values(STEADY_STATE); % get steady-state values in declaration order
    y_ss = y_ss(1:MODEL.endo_nbr); % remove exogenous
    y_ss = y_ss(MODEL.DR_ORDER); % put into DR order
    SOL.STEADY_STATE = y_ss;
end

%% Evaluate first dynamic Jacobian at steady-state
g1 = feval(str2func(MODEL.name + '_g1'),STEADY_STATE,PARAM);
if any(isnan(g1(:))) || any(isinf(g1(:))) || any(imag(g1(:)))
    error_indicator = 1;
    SOL = [];
    warning('Something wrong with the perturbation solutiondue to nan, inf or imaginary numbers');
    return
end

%% Recovering gx
n = MODEL.endo_nbr;
D = [zeros(n,n)  g1(:,(2*n+1):(3*n));
    eye(n)       zeros(n,n)];
E = [-g1(:,1:n)  -g1(:,(n+1):(2*n));
    zeros(n,n)   eye(n)];

realsmall=1e-7;
try
    [S,T,Q,Z] = qz(D,E); % upper triangular factorization of the matrix pencil
catch
    error_indicator = 1;
    SOL = [];
    warning('Error using qz');
    return
end
%norm(D-Q'*S*Z')
%norm(E-Q'*T*Z')
% Rule out that both s_ii and t_ii are zero
zxz = sum((abs(diag(S))<realsmall) & (abs(diag(T))<realsmall));
if ~(~zxz)
    error_indicator = 1;
    SOL = [];
	warning('Something wrong in the solution: Coincident zeros');
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
    SOL = [];
	warning('Blanchard-Khan order condition not fullfilled: No equilibrium exists.');
    return
elseif abs(T(n+1,n+1))<abs(S(n+1,n+1))
    error_indicator = 1;
    SOL = [];
	warning('Blanchard-Khan order condition not fullfilled: Indeterminacy.');
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
    SOL = [];
    warning('Blanchard-Khan rank condition not fullfilled: no stable solution');
    return
end
z22i = z22\eye(n);
gx = -z22i*z21;
gx = real(gx); % real function takes away very small imaginary parts of the solution

%% Recovering gu
nu = MODEL.exo_nbr;
gu = -(g1(:,(n+1):(2*n)) + g1(:,(2*n+1):(3*n)) * gx) \ g1(:,(3*n+1):(3*n+nu));

%% Write out to structure
SOL.gx = gx;
SOL.gu = gu;
