function [ p, e ] = qsimvnv( m, r, a, b )
%
%  [ P E ] = QSIMVNV( M, R, A, B )
%    uses a randomized quasi-random rule with m points to estimate an
%    MVN probability for positive definite covariance matrix r,
%     with lower integration limit column vector a and upper
%     integration limit column vector b. 
%   Probability p is output with error estimate e.
%    Example:
%     r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
%     a = -inf*[1 1 1 1 ]'; b = [ 1 2 3 4 ]';
%     [ p e ] = qsimvnv( 5000, r, a, b ); disp([ p e ])
%
%   This function uses an algorithm given in the paper
%      "Numerical Computation of Multivariate Normal Probabilities", in
%      J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : alangenz@wsu.edu
%  The primary references for the numerical integration are 
%   "On a Number-Theoretical Integration Method"
%   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11, and
%   "Randomization of Number Theoretic Methods for Multiple Integration"
%    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%
%   Alan Genz is the author of this function and following Matlab functions.
%
%
% Copyright (C) 2013, Alan Genz,  All rights reserved.               
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided the following conditions are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the 
%      distribution.
%   3. The contributor name(s) may not be used to endorse or promote 
%      products derived from this software without specific prior 
%      written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Initialization
%
[ch as bs] = chlrdr(r,a,b); ct = ch(1,1); ai = as(1); bi = bs(1); 
if ai > -9*ct, if ai < 9*ct, c = Phi(ai/ct); else, c=1; end, else c=0; end  
if bi > -9*ct, if bi < 9*ct, d = Phi(bi/ct); else, d=1; end, else d=0; end
[n, n] = size(r); ci = c; dci = d - ci; p = 0; e = 0;
ps = sqrt(primes(5*n*log(n+1)/4)); q = ps(1:n-1)'; % Richtmyer generators
ns = 12; nv = fix( max( [ m/ns 1 ] ) ); on = ones(1,nv); y = zeros(n-1,nv);
%
% Randomization loop for ns samples
%
for j = 1 : ns, c = ci*on; dc = dci*on; pv = dc;
  for i = 2 : n, x = abs( 2*mod( q(i-1)*[1:nv] + rand, 1 ) - 1 );
    y(i-1,:) = Phinv( c + x.*dc ); s = ch(i,1:i-1)*y(1:i-1,:); 
    ct = ch(i,i); ai = as(i) - s; bi = bs(i) - s; c = on; d = c; 
    c(find( ai < -9*ct )) = 0; d(find( bi < -9*ct )) = 0; 
    tstl = find( abs(ai) < 9*ct ); c(tstl) = Phi( ai(tstl)/ct ); 
    tstl = find( abs(bi) < 9*ct ); d(tstl) = Phi( bi(tstl)/ct ); 
    dc = d - c; pv = pv.*dc; 
  end, d = ( mean(pv) - p )/j; p = p + d; e = ( j - 2 )*e/j + d^2; 
end, e = 3*sqrt(e); % error estimate is 3 x standard error with ns samples.
%
% end qsimvnv
%
%
%  Standard statistical normal distribution functions
%
function p =   Phi(z), p =  erfc( -z/sqrt(2) )/2;
%function z = Phinv(p), z = norminv( p );
function z = Phinv(p), z = -sqrt(2)*erfcinv( 2*p ); % use if no norminv
%
function [ c, ap, bp ] = chlrdr( R, a, b )
%
%  Computes permuted lower Cholesky factor c for R which may be singular, 
%   also permuting integration limit vectors a and b.
%
ep = 1e-10; % singularity tolerance;
%
[n,n] = size(R); c = R; ap = a; bp = b; d = sqrt(max(diag(c),0));
for i = 1 :  n
  if d(i) > 0, c(:,i) = c(:,i)/d(i); c(i,:) = c(i,:)/d(i); 
    ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i);
  end
end
y = zeros(n,1); sqtp = sqrt(2*pi);
for k = 1 : n, im = k; ckk = 0; dem = 1; s = 0; 
  for i = k : n 
    if c(i,i) > eps, cii = sqrt( max( [c(i,i) 0] ) ); 
      if i > 1, s = c(i,1:k-1)*y(1:k-1); end
      ai = ( ap(i)-s )/cii; bi = ( bp(i)-s )/cii; de = Phi(bi) - Phi(ai);
      if de <= dem, ckk = cii; dem = de; am = ai; bm = bi; im = i; end
    end
  end
  if im > k, c(im,im) = c(k,k); 
    ap([im k]) = ap([k im]); bp([im k]) = bp([k im]);
    c([im;k],1:k-1) = c([k;im],1:k-1); c(im+1:n,[im k]) = c(im+1:n,[k im]); 
    t = c(k+1:im-1,k); c(k+1:im-1,k) = c(im,k+1:im-1)'; c(im,k+1:im-1) = t'; 
  end, c(k,k+1:n) = 0;
  if ckk > ep*k, c(k,k) = ckk; 
    for i = k+1 : n
      c(i,k) = c(i,k)/ckk; c(i,k+1:i) = c(i,k+1:i) - c(i,k)*c(k+1:i,k)';
    end
    if abs(dem) > ep, y(k) = ( exp(-am^2/2) - exp(-bm^2/2) )/(sqtp*dem); 
    else
      if am < -10, y(k) = bm;
      elseif bm > 10, y(k) = am;
      else, y(k) = ( am + bm )/2;
      end
    end
  else, c(k:n,k) = 0; y(k) = 0;
  end
end
%
% end chlrdr
%
